#!/bin/bash


#################### WGBSseq pipeline ####################
workdir="your_workdir"
cd $workdir
mkdir 1.rawdata 2.cleandata 3.bismark_index 4.bismark 5.lambda 6.Methylation
## vi Sample ##sequence data list


###### 1.rawdata ######
cd ${workdir}/1.rawdata
### download sequence 


###### 2.cleandata ######
cd ${workdir}/2.cleandata
for id in `cat ../Sample`; do
~/biosoft/fastp/fastp -i ${workdir}/1.rawdata/${id}_1.fastq.gz -o ./${id}_1.fq.gz -I ${workdir}/1.rawdata/${id}_2.fastq.gz -O ./${id}_2.fq.gz
done


###### 3.bismark_index ######
cd ${workdir}/3.index
## copy ref-genome 'zs11.genome.fa' to this folder
~/biosoft/bismark_v0.22.1/bismark_genome_preparation --bowtie2 --verbose ${workdir}/3.bismark_index


###### 4.bismark ######
cd ${workdir}/4.zs11_bismark
for id in `cat ../Sample`;do
bismark --genome ${workdir}/3.bismark_index -1 ${workdir}/2.cleandata/${id}_1.fq.gz -2 ${workdir}/2.cleandata/${id}_2.fq.gz -p 20 -o ./ 2> ${id}_bowtie2.log;
deduplicate_bismark --bam -p ./${id}_1_bismark_bt2_pe.bam --samtools_path ~/biosoft/samtools-1.9/samtools;
bismark_methylation_extractor --multicore 20 --gzip --bedGraph --buffer_size 20G --CX --cytosine_report --genome_folder ${workdir}/3.bismark_index ./${id}_1_bismark_bt2_pe.deduplicated.bam 2>${id}_extracor.log;
done


###### 5.lambda ######  
cd ${workdir}/5.lambda
## To calculate the methylation conversion rate
## download ref-genome 'lambda.fasta' to this folder  ## https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1?report=fasta
~/biosoft/bismark_v0.22.1/bismark_genome_preparation --bowtie2 --verbose ${workdir}/5.lambda
for id in `cat ../Sample`;do
bismark --genome ./ -1 ${workdir}/2.cleandata/${id}_1.fq.gz -2 ${workdir}/2.cleandata/${id}_2.fq.gz -p 20 -o ./ 2> ${id}_bowtie2.log;
done


###### 6.Methylation ######
cd ${workdir}/6.Methylation
# cp calculate_methylation.sh ./
# cp Get_methy.py ./
# cp Merge_methy_site.py ./
for id in `cat ../Sample`;do
sh calculate_methylation.sh ${id}
done