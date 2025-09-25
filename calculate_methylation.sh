#!/bin/bash

# 配置参数
#workdir="${workdir}"

sample="${sample}"
input_file="${workdir}/4.bismark/${sample}/${sample}_1_bismark_bt2_pe.deduplicated.CX_report.txt.gz"
output_dir="${workdir}/6.Methylation"
genome_fai="${workdir}/3.index/zs11.genome.fa.fai"

# 创建输出目录
mkdir -p "$output_dir"
cd "$output_dir"

echo "Starting methylation analysis pipeline..."

# Step 1: 生成基因组bins
echo "Step 1: Generating 100bp genomic bins..."
cat "$genome_fai" | cut -f 1-2 > zs11.size
bedtools makewindows -g zs11.size -w 100 > zs11.100bp.Bin

# 为每个scaffold创建单独的bin文件
chromosomes=($(seq -f "scaffoldA%02g" 1 10) $(seq -f "scaffoldC%02g" 1 9))
for chrom in "${chromosomes[@]}"; do
    grep -E "^${chrom}" zs11.100bp.Bin | awk -v chrom="$chrom" '{print $0 "\t" int($2/100)}' > "${chrom}_100bp.bin"
done

# Step 2: 并行处理每个染色体
process_chromosome() {
    local chrom=$1
    local short_name=$(echo "$chrom" | sed 's/scaffold//')
  i  
    echo "Processing $chrom..."
    
    # 提取染色体数据并计算统计量
    zcat "$input_file" | awk -v chr="$chrom" '
        $1 == chr && ($4 + $5 != 0) {
            total = $4 + $5
            rate = ($4 + $5 > 0) ? $4 / total : 0
            print $1 "\t" $2 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" total "\t" rate "\t" $6 "\t" $7
        }' > "2.${short_name}"
    
    # 区间交集分析
    bedtools intersect -a "${chrom}_100bp.bin" -b "2.${short_name}" -wa -wb -loj > "3.${short_name}"
    
    # 计算甲基化统计
    python3 Get_methy.py --name "$short_name"
    
    # 筛选低甲基化区域
    awk '($4>=2 && $5>=2 && $7<0.1 && $8<0.1 && $9<0.1) || ($4>=2 && $5<2 && $7<0.1 && $9<0.1)' "4.${short_name}" > "5.${short_name}"
    
    # 合并相邻区域
    python3 Merge_methy_site.py --name "$short_name"
   
}

export -f process_chromosome
export input_file output_dir

# 使用parallel并行处理（如果系统支持）
if command -v parallel &> /dev/null; then
    printf "%s\n" "${chromosomes[@]}" | parallel -j 4 process_chromosome
else
    # 顺序处理
    for chrom in "${chromosomes[@]}"; do
        process_chromosome "$chrom"
    done
fi

# Step 3: 合并所有结果
echo "Merging results from all chromosomes..."
rm -f 7.UMR
for id in A{01..10} C{01..09}; do
    if [[ -f "6.${id}" ]]; then
        cat "6.${id}" >> 7.${sample}.UMR
    fi
done

echo "Pipeline completed! Final results in 7.UMR"
