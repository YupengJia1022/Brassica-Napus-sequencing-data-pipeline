#!/bin/bash


#################### RNAseq pipeline ####################
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
java -jar ~/biosoft/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 5 ./${id}.1.fastq.gz ./${id}.2.fastq.gz ./${id}_1.fq.gz ./${id}_unpaired_1.fq.gz ./${id}_2.fq.gz ./${id}_unpaired_2.fq.gz ILLUMINACLIP:~/Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
done


###### 3.index Bowtie2 ######
cd ${workdir}/3.index
module load Bowtie2/2.4.4
hisat2-build –p 10 zs11.genome.fa ZS11_hisat2_index


###### 4.mapping ######
module load hisat2/2.2.0
module load SAMtools/1.9
cd ${workdir}/4.mapping
for id in `cat ../Sample`; do
hisat2 -p 10 --dta -x ${workdir}/3.index/ZS11_hisat2_index -1 ${workdir}/2.cleandata/${id}_1.fq.gz -2 ${workdir}/2.cleandata/${id}_2.fq.gz -S ${id}.sam --summary-file ./${id}.hisat2.log
samtools view -@ 10 -bS ${id}.sam -o ${id}.bam
samtools sort -@ 10 ${id}.bam -o ./${id}.sorted.bam
samtools index ${id}.sorted.bam
bamCoverage --bam ${id}.sorted.bam -o ${id}.sorted.RPKM.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 1010887457 --extendReads
~/biosoft/stringtie-2.0.Linux_x86_64/stringtie ${id}.sorted.bam -p 10 -B -e -G ${workdir}/3.index/zs11.gtf -o ${id}/${id}.gtf -A ${id}/${id}.tab -l ${id}/${id} 
done



###### 5.re_mapping ######
cd ${workdir}/5.macs
for id in `cat ../Sample`; do
echo "${workdir}/4.mapping/${i}/${i}.gtf" >> ./mergelist.txt
done

~/biosoft/stringtie-2.0.Linux_x86_64/stringtie --merge -p 10 -G ${workdir}/3.index/zs11.gtf -o ./stringtie_merged.gtf ./mergelist.txt

for id in `cat ../Sample`; do
~/biosoft/stringtie-2.0.Linux_x86_64/stringti ${workdir}/4.mapping/${id}.sorted.bam -p 24 -B -e -G ../stringtie_merged.gtf -o ${id}/${i}.gtf -A ${id}/${i}.tab -l ${id}/${i}
done


