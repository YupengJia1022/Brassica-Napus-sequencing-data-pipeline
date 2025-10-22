#!/bin/bash


#################### Hi-C pipeline ####################
workdir="your_workdir"
cd $workdir
ln -s ~/biosoft/juicer/CPU scripts
mkdir fastq references matrix compartment

###### fastq ######
cd ${workdir}/fastq
### download sequence ##
### rename ${id}_R1.fastq.gz ${id}_R2.fastq.gz


###### references ######
cd ${workdir}/references
## https://github.com/aidenlab/juicer/edit/master/misc/generate_site_positions.py
module load BWA/0.7.15
bwa index ./zs11.genome.fa 
python2.7 generate_site_positions.py DpnII ZS11 zs11.genome.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ZS11_DpnII.txt | head -n 19 > ./ZS11.sizes


###### juicer ######
cd $workdir
./scripts/juicer.sh -d $workdir -s DpnII -p $workdir/references/ZS11.sizes -y $workdir/references/ZS11_DpnII.txt -z $workdir/references/zs11.genome.fa -D $workdir


###### matrix ######
CHR=(scaffoldA01 scaffoldA02 scaffoldA03 scaffoldA04 scaffoldA05 scaffoldA06 scaffoldA07 scaffoldA08 scaffoldA09 scaffoldA10 scaffoldC01 scaffoldC02 scaffoldC03 scaffoldC04 scaffoldC05 scaffoldC06 scaffoldC07 scaffoldC08 scaffoldC09)
for i in $(seq 0 18);do
for j in $(seq 0 18);do
java -jar ~/biosoft/juicer/CPU/juicer_tools.jar dump observed NONE $workdir/aligned/inter_30.hic ${CHR[$i]} ${CHR[$j]} BP 500000 $workdir/matrix/${CHR[$i]}_${CHR[$j]}.txt;
done
done


###### compartment ######
for i in scaffoldA{01..10} scaffoldC{01..09} ;do
java -jar ~/biosoft/juicer/CPU/juicer_tools.jar eigenvector NONE -p $workdir/aligned/inter_30.hic ${i} BP 500000 > $workdir/compartment/${i}.cigenvetor
done


###### downstream analysis ######
python matrix2h5.py -i $workdir/matrix/ -b zs11.500k.bed -r 500000 -o zs11.500K
python h5_to_triple.py