
include: "../rules/Fastq/Merge_fastq_WP3.smk"

#PU:X_210122_TWIST / PU:{sample}
#-c 250 / 10000
#-M	Mark shorter split hits as secondary (for Picard compatibility).
bwa_mem_input = ["fastq/{sample}_R1.fastq.gz", "fastq/{sample}_R2.fastq.gz"]
bwa_params_extra = r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:illumina\tPU:{sample}' -v 1 -c 250 -M "
include: "../rules/Alignment/bwa-mem.smk"
include: "../rules/Alignment/index_bam.smk"


#GATK4 Haplotype caller
#The Genome Analysis Toolkit (GATK) v4.1.3.0 / v4.1.7.0
#HTSJDK Version: 2.20.1
#Picard Version: 2.20.5
include: "../rules/SNV/Haplotypecaller_WP3.smk"
include: "../rules/VCF_fix/VEP_WP3.smk"
include: "../rules/VCF_fix/Filter_vcf_WP3.smk"
include: "../rules/QC/sex_and_pedegree_WP3.smk"


#Uses vcfanno / VEP med dbsnp
#vcfanno -p 16 /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/dbsnp.conf /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D20-07453.vcf.gz |  bgzip -c > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpi5mocnf_/D20-07453-annotated.vcf.gz
#cat /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D21-00083-annotated-nomissingalt.vcf  | /sw/pipelines/bcbio-nextgen/1.1.5/anaconda/bin/bgzip --threads 16 -c > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmp_j3lp91j/D21-00083-annotated-nomissingalt.vcf.gz

#Tabix and bcftools filtering
#tabix -f -p vcf /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpm2_4q4ss/D21-00083-annotated-nomissingalt.vcf.gz
#bcftools filter -O v --soft-filter 'GATKCutoffSNP' -e 'TYPE="snp" && (MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QD < 2.0 || FS > 60.0 || (QD < 10.0 && AD[0:1] / (AD[0:1] + AD[0:0]) < 0.25 && ReadPosRankSum < 0.0) || MQ < 30.0)' -m '+' /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D21-00083-annotated-nomissingalt.vcf.gz | sed 's/\\"//g' | bgzip -c > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpk8rszybt/D21-00083-annotated-nomissingalt-filterSNP.vcf.gz
#bcftools filter -O v --soft-filter 'GATKCutoffIndel' -e 'TYPE="indel" && (ReadPosRankSum < -20.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || (QD < 10.0 && AD[0:1] / (AD[0:1] + AD[0:0]) < 0.25 && ReadPosRankSum < 0.0))' -m '+' /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D20-08020-annotated-nomissingalt-filterSNP.vcf.gz | sed 's/\\"//g' | bgzip -c > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpzx_ym6m3/D20-08020-annotated-nomissingalt-filterSNP-filterINDEL.vcf.gz

#Peddy - kör ej denna version
#/sw/pipelines/bcbio-nextgen/1.1.5/usr/local/bin/peddy -p 16  --plot --prefix /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpeozl22n6/D21-00083 /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D21-00083-annotated-nomissingalt-filterSNP-filterINDEL.vcf.gz /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D21-00083-annotated-nomissingalt-filterSNP-filterINDEL.ped 2> /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpeozl22n6/run-stderr.log

#hts_nim_tools count-reads -t 16 -F 1804 /beegfs-scratch/wp3/TE41_210122/bcbio/coverage/D21-00083/counts/fullgenome.bed /beegfs-scratch/wp3/TE41_210122/bcbio/align/D21-00083/D21-00083-sort.bam > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmp3bb0c8g0/fullgenome-1804-counts.txt

#FastQC
#fastqc -d /beegfs-scratch/wp3/TE41_210122/bcbio/qc/D21-00083/bcbiotx/tmpzgfu_uyq -t 16 --extract -o /beegfs-scratch/wp3/TE41_210122/bcbio/qc/D21-00083/bcbiotx/tmpzgfu_uyq -f bam /beegfs-scratch/wp3/TE41_210122/bcbio/qc/D21-00083/D21-00083-sort-downsample.bam

#Callable
#bedtools slop -i /beegfs-scratch/wp3/TE41_210122/bcbio/align/D21-00083/D21-00083-sort-callable_sample.bed -g /data/ref_genomes/bcbio-nextgen/sam/hg19.with.mt.fasta.fai -b 200 | bedtools merge -i - > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmp6xee6wmb/D21-00083-sort-callable_sample-padded.bed

#Bcftools stats
#bcftools stats -s D20-07453 -f PASS,. /beegfs-scratch/wp3/TE41_210122/bcbio/gatk-haplotype/D20-07453-annotated-nomissingalt-filterSNP-filterINDEL.vcf.gz > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpj0bar4fk/orig_D20-07453_bcftools_stats.txt

#Uses mosdepth
#export MOSDEPTH_Q0=NO_COVERAGE && export MOSDEPTH_Q1=LOW_COVERAGE && export MOSDEPTH_Q2=CALLABLE && mosdepth -t 16 -F 1804 -Q 1 --no-per-base --by /beegfs-scratch/wp3/TE41_210122/bcbio/coverage/D20-07453/target-genome.bed --quantize 0:1:4: /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmphhl_0f2h/D20-07453-variant_regions /beegfs-scratch/wp3/TE41_210122/bcbio/align/D20-07453/D20-07453-sort.bam

#Uses samtools stats
#/sw/pipelines/bcbio-nextgen/1.1.5/galaxy/../anaconda/bin/samtools stats -@ 16 /beegfs-scratch/wp3/TE41_210122/bcbio/align/D20-07453/D20-07453-sort.bam > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpt7yayxf6/D20-07453.txt

#Uses samtools idxstats
#/sw/pipelines/bcbio-nextgen/1.1.5/galaxy/../anaconda/bin/samtools idxstats /beegfs-scratch/wp3/TE41_210122/bcbio/align/D20-07453/D20-07453-sort.bam > /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpdxq_eq6d/D20-07453-idxstats.txt


#QC enligt ariells lina

#MultiQC
#multiqc -f -l /beegfs-scratch/wp3/TE41_210122/bcbio/qc/multiqc/list_files.txt  -o /beegfs-scratch/wp3/TE41_210122/bcbio/bcbiotx/tmpo3qtk1u7
