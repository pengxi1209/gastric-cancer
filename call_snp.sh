#!/bin/bash

fastp -w 10 -z 6 -i ${sample}_1.fq.gz -o cleanData/${sample}_1.clean.fq.gz -I ${sample}_2.fq.gz -O cleanData/${sample}_2.clean.fq.gz

fastqc -o qcReport -t 10 cleanData/${sample}_1.clean.fq.gz cleanData/${sample}_2.clean.fq.gz

bwa mem -t 10 -K 100000000 -R "@RG\tID:${sample}\tSM:${sample}\tPL:DNBseq" -Y hg38full.fa cleanData/${sample}_1.clean.fq.gz cleanData/${sample}_2.clean.fq.gz | samtools view -@ 10 -1 - > alignData2/${sample}.bam

samtools sort -@ 10 -o alignData2/${sample}.sort.bam alignData2/${sample}.bam
samtools index -@ 10 alignData2/${sample}.sort.bam

mkdir ${sample}_tamp
gatk --java-options "-Xmx10G  -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${sample}_tamp" MarkDuplicates -I alignData2/${sample}.sort.bam -O alignData2/${sample}.sort.rmdup.bam -M alignData2/${sample}.redup.metrics.txt
samtools index alignData2/${sample}.sort.rmdup.bam


## BaseRecalibrator
gatk --java-options "-Xmx2G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={sample}_tamp" BaseRecalibrator \
    -R g38full.fa \
    -I ${sample}.sort.rmdup.bam \
    --known-sites 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --known-sites dbsnp_146.hg38.vcf.gz \
    --known-sites Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O {sample}.recal_data.table

## ApplyBQSR
gatk --java-options "-Xmx2G -XX:ParallelGCThreads=2 -Djava.io.tmpdir={sample}_tamp" ApplyBQSR \
    -R g38full.fa \
    -I {sample}.sort.rmdup.bam \
    -bqsr {sample}.recal_data.table \
    -O {sample}.recal.bam


gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" Mutect2 \
    -R g38full.fa \
    -I {sample}/{sample}.recal.bam \
    -tumor {sample} \
    --panel-of-normals 1000g_pon.hg38.vcf.gz \
    -germline-resource gastric_somaticall_pon.vcf.gz \
    --f1r2-tar-gz f1r2.tar.gz \
    -O {sample}.mutect2.unfiltered.vcf
#
gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O read-orientation-model.tar.gz

gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" GetPileupSummaries \
    -I {sample}/{sample}.recal.bam \
    -V small_exac_common_3.hg38.vcf.gz \
    -L small_exac_common_3.hg38.vcf.gz \
    -O getpileupsummaries.table

gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" CalculateContamination \
    -I getpileupsummaries.table \
    -O calculatecontamination.table

gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" FilterMutectCalls \
    -V {sample}.mutect2.unfiltered.vcf \
    --ob-priors read-orientation-model.tar.gz \
    --contamination-table calculatecontamination.table \
    -R hg38full.fa \
    -O {sample}.mutect2.filtered.vcf

gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" SelectVariants \
    -V {sample}.mutect2.filtered.vcf \
    -O {sample}.mutect2.selected.vcf \
    --select-type-to-include SNP \
    --select-type-to-include INDEL

gatk --java-options "-Xmx10G -XX:ParallelGCThreads=10 -Djava.io.tmpdir={sample}_tamp" SelectVariants     -V {sample}.mutect2.selected.vcf     -O {sample}.mutect2.selected.PASS.vcf     --exclude-filtered


rm -rf {sample}_tamp



 perl table_annovar.pl ${sample}.mutect2.selected.PASS.vcf hunmadb_annovar --outfile ${sample} --buildver hg38 --protocol refGene,avsnp150,exac03,clinvar_20210501,icgc28,cosmic70,dbnsfp42a,gnomad41_genome --operation g,f,f,f,f,f,f,f --vcfinput --dot2underline --nastring . --remove

