
#!/bin/bash

ATACseq_workdir="./"
WGBSseq_workdir="./"

###################### ATACseq OCRs ######################
echo "Processing ATACseq OCRs..."

# 复制文件并检查是否存在
cp ${ATACseq_workdir}/5.macs/*_peaks.narrowPeak ./
if [ $? -ne 0 ]; then
    echo "Error: Failed to copy ATACseq files" >&2
    exit 1
fi

# 处理每个样本的OCR
for sample in 9Bao22 B262 NY10 SWU99; do
    echo "Processing sample: $sample"
    
    # 检查输入文件是否存在
    if [[ ! -f "${sample}_r1_peaks.narrowPeak" || ! -f "${sample}_r2_peaks.narrowPeak" ]]; then
        echo "Warning: Missing files for sample $sample, skipping..." >&2
        continue
    fi
    
    # 一步完成交集和区域合并，减少临时文件
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

# 按染色体合并区域
for sample in 9Bao22 B262 NY10 SWU99; do
    if [[ ! -f "1.${sample}.intersection.OCR.bdg" ]]; then
        continue
    fi
    
    echo "Merging regions for sample: $sample"
    > "2.${sample}.intersection.OCR.bdg"  # 清空输出文件
    
    for chrom in scaffoldA{01..10} scaffoldC{01..09}; do
        # 提取当前染色体的区域
        grep -E "^${chrom}[[:space:]]" "1.${sample}.intersection.OCR.bdg" > "1.${chrom}" 2>/dev/null
        
        if [[ ! -s "1.${chrom}" ]]; then
            continue  # 跳过空文件
        fi
        
        # 合并重叠区域
        if python Merge.py --n "${chrom}" 2>/dev/null && [[ -s "2.${chrom}" ]]; then
            awk '{print $0 "\t1"}' "2.${chrom}" >> "2.${sample}.intersection.OCR.bdg"
        fi
        
        # 清理临时文件
        rm -f "1.${chrom}" "2.${chrom}"
    done
done

# 合并所有样本的OCR
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

# 复制文件
cp ${WGBSseq_workdir}/6.Methylation/7.*.UMR ./
if [ $? -ne 0 ]; then
    echo "Error: Failed to copy WGBSseq files" >&2
    exit 1
fi

# 处理每个样本的UMR
for sample in 9Bao22 B262 NY10 SWU99; do
    echo "Processing sample: $sample"
    
    # 检查三个重复文件是否存在
    if [[ ! -f "7.${sample}_rep1.UMR" || ! -f "7.${sample}_rep2.UMR" || ! -f "7.${sample}_rep3.UMR" ]]; then
        echo "Warning: Missing UMR files for sample $sample, skipping..." >&2
        continue
    fi
    
    # 一步完成三个重复的交集
    bedtools intersect \
        -a "7.${sample}_rep1.UMR" \
        -b "7.${sample}_rep2.UMR" | \
    bedtools intersect \
        -a - \
        -b "7.${sample}_rep3.UMR" | \
    awk '{print $0 "\t1"}' > "1.${sample}.intersection.UMR.bdg"
done

# 按染色体合并UMR区域
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

# 合并所有样本的UMR
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

# 检查输入文件
if [[ ! -s "ATACseq_mergeOCR.bdg" || ! -s "WGBSseq_mergeUMR.bdg" ]]; then
    echo "Error: Missing input files for RE analysis" >&2
    exit 1
fi

# 使用awk一次性处理，减少临时文件
{
    # 处理非重叠区域
    bedtools intersect -a ATACseq_mergeOCR.bdg -b WGBSseq_mergeUMR.bdg -v | \
    awk '{len=$3-$2+1; print $1 "\t" $2 "\t" $3 "\t" len "\t0"}'
    
    # 处理重叠区域并计算覆盖率
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

# 清理临时文件（可选）
# rm -f 1.*.bdg 2.*.bdg

echo "Done!"
