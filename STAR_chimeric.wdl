version 1.0

workflow wf_star {
    meta {
        version: 'v0.2'
        description: 'STAR RNA-seq alignment workflow with gene counting. Compatible with Terra Bio.'
    }

    input {
        # Original mandatory inputs
        Array[File] read1
        Array[File]? read2
        File idx_tar
        String prefix = "rna-project"
        String genome_name

        # Optional inputs
        File? gtf_annotation
        Int cpus = 16
        Int disk_gb = 100
        Int mem_gb = 128
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_star"
        Int preemptible_tries = 2
        Int max_retries = 1
        String extra_star_args = ""
    }

    call rna_align {
        input:
            fastq_R1            = read1,
            fastq_R2            = read2,
            genome_index_tar    = idx_tar,
            gtf_annotation      = gtf_annotation,
            genome_name         = genome_name,
            prefix              = prefix,
            cpus                = cpus,
            disk_gb             = disk_gb,
            mem_gb              = mem_gb,
            docker_image        = docker_image,
            preemptible_tries   = preemptible_tries,
            max_retries         = max_retries,
            extra_star_args     = extra_star_args
    }

    output {
        File rna_alignment_raw   = rna_align.rna_alignment
        File rna_alignment_index = rna_align.rna_alignment_index
        File rna_alignment_log   = rna_align.rna_alignment_log
        File rna_gene_counts     = rna_align.rna_gene_counts
    }
}

task rna_align {
    meta {
        version: 'v0.2'
        description: 'Align RNA-seq reads with STAR and produce gene counts.'
    }

    input {
        Array[File] fastq_R1
        Array[File]? fastq_R2
        File genome_index_tar
        File? gtf_annotation
        String genome_name
        String prefix = "rna"
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_star"
        Int cpus = 16
        Int disk_gb = 100
        Int mem_gb = 128
        Int preemptible_tries = 2
        Int max_retries = 1
        String extra_star_args = ""
    }

    Int samtools_cpus = 6
    Int samtools_mem_gb = 8

    String sorted_bam = "${default="rna" prefix}.rna.align.${genome_name}.sorted.bam"
    String sorted_bai = "${default="rna" prefix}.rna.align.${genome_name}.sorted.bam.bai"
    String star_prefix = "${default="rna" prefix}.rna.align.${genome_name}."

    command {
        set -euo pipefail

        # Extract genome index
        tar xvzf ${genome_index_tar} --no-overwrite-dir --no-same-owner --no-same-permissions -C ./

        mkdir -p out

        $(which STAR) \
            --runThreadN ${cpus} \
            --chimOutType WithinBAM \
            --genomeDir ./ \
            --readFilesIn ${sep=',' fastq_R1} ${sep=',' fastq_R2} \
            --outFileNamePrefix out/${star_prefix} \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --outSAMattributes NH HI AS nM MD \
            --limitOutSJcollapsed 2000000 \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
            ${"--sjdbGTFfile " + gtf_annotation} \
            ${extra_star_args}

        $(which samtools) sort \
            -@ ${samtools_cpus} \
            -m ${samtools_mem_gb}G \
            -o out/${sorted_bam} \
            out/${star_prefix}Aligned.out.bam

        $(which samtools) index \
            -@ ${cpus} \
            out/${sorted_bam}
    }

    output {
        File rna_alignment       = "out/${sorted_bam}"
        File rna_alignment_index = "out/${sorted_bai}"
        File rna_alignment_log   = glob('out/*.Log.final.out')[0]
        File rna_gene_counts     = glob('out/*ReadsPerGene.out.tab')[0]
    }

    runtime {
        cpu: cpus
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_gb + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    parameter_meta {
        fastq_R1: {
            description: 'Read1 fastq',
            help: 'Processed fastq for read1. Supports gs:// paths for Terra.',
            example: 'gs://my-bucket/sample.R1.fq.gz'
        }
        fastq_R2: {
            description: 'Read2 fastq (optional)',
            help: 'Processed fastq for read2 (paired-end). Supports gs:// paths.',
            example: 'gs://my-bucket/sample.R2.fq.gz'
        }
        genome_index_tar: {
            description: 'STAR genome index',
            help: 'Pre-built STAR index files in tar.gz. Supports gs:// paths.',
            example: 'gs://my-bucket/star_index.tar.gz'
        }
        gtf_annotation: {
            description: 'GTF annotation (optional)',
            help: 'GTF file for gene counting. Only needed if not already built into the genome index.',
            example: 'gs://my-bucket/annotation.gtf'
        }
        genome_name: {
            description: 'Reference name',
            help: 'The name of the reference genome used for output file naming.',
            example: ['hg38', 'mm10', 'both']
        }
        prefix: {
            description: 'Prefix for output files',
            help: 'Prefix that will be used to name the output files.',
            example: 'MyExperiment'
        }
        cpus: {
            description: 'Number of cpus',
            help: 'Set the number of cpus used by STAR.',
            example: '16'
        }
        docker_image: {
            description: 'Docker image',
            help: 'Docker image for alignment. Dependencies: STAR, samtools.',
            example: ['us.gcr.io/buenrostro-share-seq/share_task_star']
        }
        extra_star_args: {
            description: 'Additional STAR arguments',
            help: 'Any additional STAR options as a single string, appended to the command.',
            example: '--outFilterType BySJout --alignIntronMax 1000000 --twopassMode Basic'
        }
    }
}
