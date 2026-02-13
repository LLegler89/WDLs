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

    # Container used for both fastp and STAR
    String docker_image = "trinityctat/starfusion:latest"

    # Optional extra STAR args
    String extra_star_args = ""

    # --- fastp options ---
    # Examples:
    #   "--detect_adapter_for_pe --cut_front --cut_tail --cut_mean_quality 20 --length_required 20"
    #   "--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    String extra_fastp_args = "--detect_adapter_for_pe"

    Boolean cleanup_intermediates = true
  }

  call Fastp_Trim_Paired {
    input:
      sample_id = sample_id,
      read1_fastq = read1_fastq,
      read2_fastq = read2_fastq,
      threads = threads,
      memory = memory,
      docker_image = docker_image,
      extra_fastp_args = extra_fastp_args
  }

  call STAR_Align_SortedBam {
    input:
      sample_id = sample_id,
      read1_fastq = Fastp_Trim_Paired.trimmed_read1_fastq,
      read2_fastq = Fastp_Trim_Paired.trimmed_read2_fastq,
      star_index_tar = star_index_tar,
      threads = threads,
      memory = memory,
      docker_image = docker_image,
      extra_star_args = extra_star_args,
      cleanup_intermediates = cleanup_intermediates
  }

  output {
    File trimmed_read1 = Fastp_Trim_Paired.trimmed_read1_fastq
    File trimmed_read2 = Fastp_Trim_Paired.trimmed_read2_fastq
    File fastp_json = Fastp_Trim_Paired.fastp_json
    File fastp_html = Fastp_Trim_Paired.fastp_html

    File bam = STAR_Align_SortedBam.sorted_bam
    File bam_bai = STAR_Align_SortedBam.sorted_bam_bai
    File star_log_final = STAR_Align_SortedBam.log_final
    File star_log_out = STAR_Align_SortedBam.log_out
    File star_log_progress = STAR_Align_SortedBam.log_progress
    File gene_counts = STAR_Align_SortedBam.gene_counts
  }
}

task Fastp_Trim_Paired {
  input {
    String sample_id
    File read1_fastq
    File read2_fastq
    Int threads
    String memory
    String docker_image
    String extra_fastp_args = "--detect_adapter_for_pe"
  }

  Int disk_gb = ceil(
     size(read1_fastq, "GB") +
     size(read2_fastq, "GB")
  ) + 500

  command <<<
    set -euo pipefail

    echo "Checking for fastp..."
    if ! command -v fastp >/dev/null 2>&1; then
      echo "fastp not found; attempting install..."

      # Try apt-get if present (Debian/Ubuntu-based images)
      if command -v apt-get >/dev/null 2>&1; then
        apt-get update -y
        apt-get install -y fastp || true
      fi

      # If still missing, fetch static binary (requires network)
      if ! command -v fastp >/dev/null 2>&1; then
        echo "fastp still not found; downloading static binary..."
        curl -fsSL -o fastp \
          https://github.com/OpenGene/fastp/releases/latest/download/fastp
        chmod +x fastp
        export PATH="$PWD:$PATH"
      fi
    fi

    fastp --version

    OUT1="~{sample_id}.fastp.R1.fastq.gz"
    OUT2="~{sample_id}.fastp.R2.fastq.gz"
    JSON="~{sample_id}.fastp.json"
    HTML="~{sample_id}.fastp.html"

    # fastp can read gz and write gz; no special readFilesCommand stuff needed here.
    fastp \
      -i "~{read1_fastq}" \
      -I "~{read2_fastq}" \
      -o "${OUT1}" \
      -O "${OUT2}" \
      -w ~{threads} \
      -j "${JSON}" \
      -h "${HTML}" \
      ~{extra_fastp_args}

    if [[ ! -s "${OUT1}" ]] || [[ ! -s "${OUT2}" ]]; then
      echo "ERROR: fastp did not produce trimmed FASTQs." >&2
      ls -lah >&2
      exit 1
    fi

    if [[ ! -s "${JSON}" ]] || [[ ! -s "${HTML}" ]]; then
      echo "ERROR: fastp did not produce JSON/HTML reports." >&2
      ls -lah >&2
      exit 1
    fi
  >>>

  output {
    File trimmed_read1_fastq = "~{sample_id}.fastp.R1.fastq.gz"
    File trimmed_read2_fastq = "~{sample_id}.fastp.R2.fastq.gz"
    File fastp_json = "~{sample_id}.fastp.json"
    File fastp_html = "~{sample_id}.fastp.html"
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

    # Support .tar or .tar.gz
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
