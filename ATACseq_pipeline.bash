#!/bin/bash

#################### ATACseq pipeline ####################
workdir="your_workdir"
cd $workdir
mkdir 1.rawdata 2.cleandata 3.index 4.mapping 5.macs
## vi Sample ##sequence data list


###### 1.rawdata ######
cd ${workdir}/1.rawdata
### download sequence 


###### 2.cleandata ######
cd ${workdir}/2.cleandata
for id in `cat ../Sample`; do
~/biosoft/fastp/fastp -i ${workdir}/1.rawdata/${id}_1.fastq.gz -o ./${id}_1.fq.gz -I ${workdir}/1.rawdata/${id}_2.fastq.gz -O ./${id}_2.fq.gz
done


###### 3.index Bowtie2 ######
cd ${workdir}/3.index
module load Bowtie2/2.4.4
bowtie2-build -f zs11.genome.fa --threads 10 ZS11_bowtie2_index


###### 4.mapping ######
cd ${workdir}/4.mapping
for id in `cat ../Sample`; do
bowtie2 -q -p 10 --very-sensitive -X 1000 --fr -x ${workdir}/3.index/ZS11_bowtie2_index -1 ${workdir}/2.cleandata/${id}_1.fq.gz -2 ${workdir}/2.cleandata/${id}_2.fq.gz -S ${id}.sam 2> ${id}.bowtie2.log
samtools view -@ 10 -b -S -h -F 12  ${id}.sam > ${id}.bam
samtools sort -o ${id}.sorted.bam -@ 10 ${id}.bam
samtools index -b -@ 10 ${id}.sorted.bam
java -jar ~/biosoft/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${id}.sorted.bam O=${id}.sorted.rmdup.bam M=${id}.sorted.rmdup.log
samtools index -b -@ 10 ${id}.sorted.rmdup.bam
done


###### 5.macs ######
cd ${workdir}/5.macs
for id in `cat ../Sample`; do
macs2 callpeak -t ${workdir}/4.mapping/${id}.sorted.rmdup.bam -f BAM --nomodel -q 0.01 --extsize 200 --shift -100 -g 1010887457 --keep-dup all -n ${id} -B --call-summit --outdir ./ 2>> ./${id}.macs2.log
done



