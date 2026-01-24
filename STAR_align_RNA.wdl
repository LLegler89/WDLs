version 1.0

workflow STAR_align_paired_rnaseq {
  input {
    String sample_id

    File read1_fastq
    File read2_fastq

    # Uncompressed tar containing STAR_index_dir/ at its root (can also be .tar.gz)
    File star_index_tar

    Int threads = 8
    String memory = "64G"

    # Optional: to match whatever your STAR-Fusion WDL uses
    String docker_image = "trinityctat/starfusion:latest"

    # Optional extra STAR args (ex: "--outSAMattributes NH HI AS nM --chimOutType Junctions")
    String extra_star_args = ""

    # --- Cutadapt options ---
    # If empty, runs cutadapt in "minimal" mode (quality trimming only if you specify it via extra_cutadapt_args)
    # You can also pass adapters through extra_cutadapt_args, e.g.:
    #   "--adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    String extra_cutadapt_args = ""

    # If true, delete large intermediates where safe
    Boolean cleanup_intermediates = true
  }

  call Cutadapt_Trim_Paired {
    input:
      sample_id = sample_id,
      read1_fastq = read1_fastq,
      read2_fastq = read2_fastq,
      threads = threads,
      memory = memory,
      docker_image = docker_image,
      extra_cutadapt_args = extra_cutadapt_args
  }

  call STAR_Align_SortedBam {
    input:
      sample_id = sample_id,
      # feed trimmed reads into STAR
      read1_fastq = Cutadapt_Trim_Paired.trimmed_read1_fastq,
      read2_fastq = Cutadapt_Trim_Paired.trimmed_read2_fastq,
      star_index_tar = star_index_tar,
      threads = threads,
      memory = memory,
      docker_image = docker_image,
      extra_star_args = extra_star_args,
      cleanup_intermediates = cleanup_intermediates
  }

  output {
    File trimmed_read1 = Cutadapt_Trim_Paired.trimmed_read1_fastq
    File trimmed_read2 = Cutadapt_Trim_Paired.trimmed_read2_fastq
    File cutadapt_report = Cutadapt_Trim_Paired.cutadapt_report

    File bam = STAR_Align_SortedBam.sorted_bam
    File bam_bai = STAR_Align_SortedBam.sorted_bam_bai
    File star_log_final = STAR_Align_SortedBam.log_final
    File star_log_out = STAR_Align_SortedBam.log_out
    File star_log_progress = STAR_Align_SortedBam.log_progress

    # STAR gene counts (requires index built with a GTF)
    File gene_counts = STAR_Align_SortedBam.gene_counts
  }
}

task Cutadapt_Trim_Paired {
  input {
    String sample_id
    File read1_fastq
    File read2_fastq
    Int threads
    String memory
    String docker_image
    String extra_cutadapt_args = ""
  }

  Int disk_gb = ceil(
     size(read1_fastq, "GB") +
     size(read2_fastq, "GB")
  ) + 50

  command <<<
    set -euo pipefail

    echo "Checking for cutadapt..."
    if ! command -v cutadapt >/dev/null 2>&1; then
      echo "cutadapt not found; attempting install via python/pip..."
      if command -v python3 >/dev/null 2>&1; then
        python3 -m pip install --user --no-cache-dir cutadapt
        export PATH="$HOME/.local/bin:$PATH"
      elif command -v python >/dev/null 2>&1; then
        python -m pip install --user --no-cache-dir cutadapt
        export PATH="$HOME/.local/bin:$PATH"
      else
        echo "ERROR: Neither cutadapt nor python is available in the container." >&2
        exit 1
      fi
    fi

    echo "cutadapt version:"
    cutadapt --version

    # Build input commands (support gz or plain)
    R1="~{read1_fastq}"
    R2="~{read2_fastq}"

    OUT1="~{sample_id}.cutadapt.R1.fastq.gz"
    OUT2="~{sample_id}.cutadapt.R2.fastq.gz"

    # Note: cutadapt auto-detects gzip by filename, so we always write .gz outputs.
    # You can pass adapters/quality/trimming parameters through extra_cutadapt_args.
    echo "Running cutadapt..."
    cutadapt \
      -j ~{threads} \
      -o "${OUT1}" \
      -p "${OUT2}" \
      ~{extra_cutadapt_args} \
      "${R1}" "${R2}" \
      > "~{sample_id}.cutadapt.report.txt"

    # Sanity checks
    if [[ ! -s "${OUT1}" ]] || [[ ! -s "${OUT2}" ]]; then
      echo "ERROR: cutadapt did not produce trimmed FASTQs." >&2
      ls -lah >&2
      exit 1
    fi
  >>>

  output {
    File trimmed_read1_fastq = "~{sample_id}.cutadapt.R1.fastq.gz"
    File trimmed_read2_fastq = "~{sample_id}.cutadapt.R2.fastq.gz"
    File cutadapt_report = "~{sample_id}.cutadapt.report.txt"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{memory}"
    disks: "local-disk ~{disk_gb} SSD"
  }
}

task STAR_Align_SortedBam {
  input {
    String sample_id
    File read1_fastq
    File read2_fastq
    File star_index_tar
    Int threads
    String memory
    String docker_image
    String extra_star_args = ""
    Boolean cleanup_intermediates = true
  }

  Int disk_gb = ceil(
     size(read1_fastq, "GB") +
     size(read2_fastq, "GB") +
     size(star_index_tar, "GB")
  ) + 200

  command <<<
    set -euo pipefail

    # Support .tar or .tar.gz; tar will usually autodetect with -xf, but -xzf will fail for plain .tar.
    # So: try gzip mode first, then fallback.
    if tar -tzf "~{star_index_tar}" >/dev/null 2>&1; then
      tar -xzf "~{star_index_tar}" --no-overwrite-dir --no-same-owner --no-same-permissions -C .
    else
      tar -xf "~{star_index_tar}" --no-overwrite-dir --no-same-owner --no-same-permissions -C .
    fi

    GENOME_DIR=$(find . -maxdepth 3 -type f -name SA -printf '%h\n' | head -n 1)
    if [[ -z "${GENOME_DIR}" ]]; then
      echo "ERROR: Could not find STAR index (SA file) after untarring." >&2
      ls -lah >&2
      exit 1
    fi

    if [[ "~{cleanup_intermediates}" == "true" ]]; then
      rm -f "~{star_index_tar}" || true
    fi

    READ_CMD=""
    if [[ "~{read1_fastq}" == *.gz ]]; then
      READ_CMD="--readFilesCommand zcat"
    fi

    STAR \
      --runThreadN ~{threads} \
      --genomeDir "${GENOME_DIR}" \
      --readFilesIn "~{read1_fastq}" "~{read2_fastq}" \
      ${READ_CMD} \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "~{sample_id}.star." \
      --quantMode GeneCounts \
      ~{extra_star_args}

    if [[ ! -f ~{sample_id}.star.Aligned.sortedByCoord.out.bam ]]; then
      echo "ERROR: STAR did not produce ~{sample_id}.star.Aligned.sortedByCoord.out.bam" >&2
      ls -lah >&2
      exit 1
    fi

    if [[ ! -f ~{sample_id}.star.ReadsPerGene.out.tab ]]; then
      echo "ERROR: STAR did not produce ~{sample_id}.star.ReadsPerGene.out.tab (gene counts)." >&2
      echo "This typically happens if the STAR genome index was not built with a GTF (sjdbGTFfile)." >&2
      ls -lah >&2
      exit 1
    fi

    samtools index -@ ~{threads} ~{sample_id}.star.Aligned.sortedByCoord.out.bam
  >>>

  output {
    File sorted_bam = "~{sample_id}.star.Aligned.sortedByCoord.out.bam"
    File sorted_bam_bai = "~{sample_id}.star.Aligned.sortedByCoord.out.bam.bai"
    File log_final = "~{sample_id}.star.Log.final.out"
    File log_out = "~{sample_id}.star.Log.out"
    File log_progress = "~{sample_id}.star.Log.progress.out"
    File gene_counts = "~{sample_id}.star.ReadsPerGene.out.tab"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{memory}"
    disks: "local-disk ~{disk_gb} SSD"
  }
}
