
#!/bin/bash

ATACseq_workdir="./"
WGBSseq_workdir="./"

###################### ATACseq OCRs ######################
echo "Processing ATACseq OCRs..."


cp ${ATACseq_workdir}/5.macs/*_peaks.narrowPeak ./
if [ $? -ne 0 ]; then
    echo "Error: Failed to copy ATACseq files" >&2
    exit 1
fi

for sample in 9Bao22 B262 NY10 SWU99; do
    echo "Processing sample: $sample"

    if [[ ! -f "${sample}_r1_peaks.narrowPeak" || ! -f "${sample}_r2_peaks.narrowPeak" ]]; then
        echo "Warning: Missing files for sample $sample, skipping..." >&2
        continue
    fi
    
    bedtools intersect \
        -a <(cut -f1-3 "${sample}_r1_peaks.narrowPeak") \
        -b <(cut -f1-3 "${sample}_r2_peaks.narrowPeak") \
        -wa -wb | \
    awk '{
        start = ($2 < $5) ? $2 : $5
        end = ($3 > $6) ? $3 : $6
        print $1 "\t" start "\t" end "\t1"
    }' > "1.${sample}.intersection.OCR.bdg"
done


for sample in 9Bao22 B262 NY10 SWU99; do
    if [[ ! -f "1.${sample}.intersection.OCR.bdg" ]]; then
        continue
    fi
    
    echo "Merging regions for sample: $sample"
    > "2.${sample}.intersection.OCR.bdg" 
    
    for chrom in scaffoldA{01..10} scaffoldC{01..09}; do
        grep -E "^${chrom}[[:space:]]" "1.${sample}.intersection.OCR.bdg" > "1.${chrom}" 2>/dev/null
        
        if [[ ! -s "1.${chrom}" ]]; then
            continue  
        fi

        if python Merge.py --n "${chrom}" 2>/dev/null && [[ -s "2.${chrom}" ]]; then
            awk '{print $0 "\t1"}' "2.${chrom}" >> "2.${sample}.intersection.OCR.bdg"
        fi
        
        rm -f "1.${chrom}" "2.${chrom}"
    done
done

ocr_files=()
for sample in 9Bao22 B262 NY10 SWU99; do
    if [[ -s "2.${sample}.intersection.OCR.bdg" ]]; then
        ocr_files+=("2.${sample}.intersection.OCR.bdg")
    fi
done

if [[ ${#ocr_files[@]} -eq 0 ]]; then
    echo "Error: No OCR files to merge" >&2
    exit 1
fi

cat "${ocr_files[@]}" | sort -k1,1 -k2,2n -k3,3n | \
    bedtools merge | \
    awk '{print $0 "\t1"}' > ATACseq_mergeOCR.bdg

###################### WGBSseq UMRs ######################
echo "Processing WGBSseq UMRs..."
cp ${WGBSseq_workdir}/6.Methylation/7.*.UMR ./
if [ $? -ne 0 ]; then
    echo "Error: Failed to copy WGBSseq files" >&2
    exit 1
fi

for sample in 9Bao22 B262 NY10 SWU99; do
    echo "Processing sample: $sample"

    if [[ ! -f "7.${sample}_rep1.UMR" || ! -f "7.${sample}_rep2.UMR" || ! -f "7.${sample}_rep3.UMR" ]]; then
        echo "Warning: Missing UMR files for sample $sample, skipping..." >&2
        continue
    fi

    bedtools intersect \
        -a "7.${sample}_rep1.UMR" \
        -b "7.${sample}_rep2.UMR" | \
    bedtools intersect \
        -a - \
        -b "7.${sample}_rep3.UMR" | \
    awk '{print $0 "\t1"}' > "1.${sample}.intersection.UMR.bdg"
done

for sample in 9Bao22 B262 NY10 SWU99; do
    if [[ ! -f "1.${sample}.intersection.UMR.bdg" ]]; then
        continue
    fi
    
    echo "Merging UMR regions for sample: $sample"
    > "2.${sample}.intersection.UMR.bdg"
    
    for chrom in scaffoldA{01..10} scaffoldC{01..09}; do
        grep -E "^${chrom}[[:space:]]" "1.${sample}.intersection.UMR.bdg" > "1.${chrom}" 2>/dev/null
        
        if [[ ! -s "1.${chrom}" ]]; then
            continue
        fi
        
        if python Merge.py --n "${chrom}" 2>/dev/null && [[ -s "2.${chrom}" ]]; then
            awk '{print $0 "\t1"}' "2.${chrom}" >> "2.${sample}.intersection.UMR.bdg"
        fi
        
        rm -f "1.${chrom}" "2.${chrom}"
    done
done

umr_files=()
for sample in 9Bao22 B262 NY10 SWU99; do
    if [[ -s "2.${sample}.intersection.UMR.bdg" ]]; then
        umr_files+=("2.${sample}.intersection.UMR.bdg")
    fi
done

if [[ ${#umr_files[@]} -eq 0 ]]; then
    echo "Error: No UMR files to merge" >&2
    exit 1
fi

cat "${umr_files[@]}" | sort -k1,1 -k2,2n -k3,3n | \
    bedtools merge | \
    awk '{print $0 "\t1"}' > WGBSseq_mergeUMR.bdg

###################### ACR+UMRs REs ######################
echo "Identifying ACR+UMR regulatory elements..."

if [[ ! -s "ATACseq_mergeOCR.bdg" || ! -s "WGBSseq_mergeUMR.bdg" ]]; then
    echo "Error: Missing input files for RE analysis" >&2
    exit 1
fi


{
    bedtools intersect -a ATACseq_mergeOCR.bdg -b WGBSseq_mergeUMR.bdg -v | \
    awk '{len=$3-$2+1; print $1 "\t" $2 "\t" $3 "\t" len "\t0"}'
    
    bedtools intersect -a ATACseq_mergeOCR.bdg -b WGBSseq_mergeUMR.bdg -wo | \
    awk '{
        key = $1 "\t" $2 "\t" $3
        len = $3 - $2 + 1
        overlap[key] += $7
        region[key] = $1 "\t" $2 "\t" $3 "\t" len
    }
    END {
        for (key in overlap) {
            split(key, parts, "\t")
            len = parts[3] - parts[2] + 1
            if (len > 1) {  # 避免除零错误
                coverage = overlap[key] / (len - 1)
                print region[key] "\t" overlap[key] "\t" coverage
            }
        }
    }'
} | sort -k1,1 -k2,2n -k3,3n | awk '$6 > 0' > RE.txt

echo "Analysis completed. Results saved to RE.txt"
echo "Done!"
