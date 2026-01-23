version 1.0

workflow wf_star {
    meta {
        version: 'v0.2'
    }

    input {
        # Original mandatory/default inputs preserved
        Array[File] read1
        Array[File]? read2
        File idx_tar
        String prefix = "rna-project"
        String genome_name
        Int cpus = 16
        Int disk_gb = 100
        Int mem_gb = 128
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_star"

        # New optional inputs
        File? gtf_annotation        # Only needed if GTF was not baked into the index
        String extra_star_args = "" # Any additional STAR flags as a single string
        Int preemptible = 2
    }

    call rna_align {
        input:
            fastq_R1         = read1,
            fastq_R2         = read2,
            genome_index_tar = idx_tar,
            genome_name      = genome_name,
            prefix           = prefix,
            cpus             = cpus,
            disk_gb          = disk_gb,
            mem_gb           = mem_gb,
            docker_image     = docker_image,
            gtf_annotation   = gtf_annotation,
            extra_star_args  = extra_star_args,
            preemptible      = preemptible
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
    }

    input {
        Array[File] fastq_R1
        Array[File]? fastq_R2
        File genome_index_tar
        String genome_name
        String prefix = "rna"
        String docker_image = "us.gcr.io/buenrostro-share-seq/share_task_star"
        Int cpus = 16
        Int disk_gb = 100
        Int mem_gb = 128
        File? gtf_annotation
        String extra_star_args = ""
        Int preemptible = 2
    }

    Int samtools_cpus = 6
    Int samtools_mem_gb = 8

    String sorted_bam  = "~{prefix}.rna.align.~{genome_name}.sorted.bam"
    String sorted_bai  = "~{prefix}.rna.align.~{genome_name}.sorted.bam.bai"
    String star_prefix = "~{prefix}.rna.align.~{genome_name}."

    command <<<
        set -e

        # Untar the genome
        tar xvzf ~{genome_index_tar} --no-overwrite-dir --no-same-owner --no-same-permissions -C ./

        mkdir -p out

        STAR \
            --runThreadN ~{cpus} \
            --chimOutType WithinBAM \
            --genomeDir ./ \
            --readFilesIn ~{sep="," fastq_R1} ~{if defined(fastq_R2) then sep(",", select_first([fastq_R2])) else ""} \
            --outFileNamePrefix out/~{star_prefix} \
            --outFilterMultimapNmax 20 \
            --outFilterScoreMinOverLread 0.3 \
            --outFilterMatchNminOverLread 0.3 \
            --outSAMattributes NH HI AS nM MD \
            --limitOutSJcollapsed 2000000 \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
            ~{"--sjdbGTFfile " + gtf_annotation} \
            ~{extra_star_args}

        samtools sort \
            -@ ~{samtools_cpus} \
            -m ~{samtools_mem_gb}G \
            -o out/~{sorted_bam} \
            out/~{star_prefix}Aligned.out.bam

        samtools index \
            -@ ~{cpus} \
            out/~{sorted_bam}
    >>>

    output {
        File rna_alignment       = "out/~{sorted_bam}"
        File rna_alignment_index = "out/~{sorted_bai}"
        File rna_alignment_log   = glob("out/*.Log.final.out")[0]
        File rna_gene_counts     = glob("out/*ReadsPerGene.out.tab")[0]
    }

    runtime {
        cpu: cpus
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        docker: docker_image
        preemptible: preemptible
        maxRetries: 0
    }

    parameter_meta {
        fastq_R1: {
            description: 'Read1 fastq',
            help: 'Processed fastq for read1.',
            example: 'gs://my-bucket/processed.rna.R1.fq.gz'
        }
        genome_index_tar: {
            description: 'STAR indexes',
            help: 'Index files for STAR to use during alignment in tar.gz.',
            example: 'gs://my-bucket/star_index.tar.gz'
        }
        genome_name: {
            description: 'Reference name',
            help: 'The name of the reference genome used by the aligner.',
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
            example: 'us.gcr.io/buenrostro-share-seq/share_task_star'
        }
        gtf_annotation: {
            description: 'GTF annotation file',
            help: 'GTF file required for gene counting if not already included in the genome index.',
            example: 'gs://my-bucket/gencode.v38.annotation.gtf'
        }
        extra_star_args: {
            description: 'Additional STAR arguments',
            help: 'Any additional STAR options as a single string, appended to the command.',
            example: '--alignIntronMax 1000000 --outFilterType BySJout --twopassMode Basic'
        }
    }
}
