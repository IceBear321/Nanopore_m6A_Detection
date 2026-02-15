#!/bin/bash
# =============================================================================
# 02_alignment.sh
# 序列比对 - 将reads比对到参考基因组
# =============================================================================
# 输入: .fastq 序列文件
# 输出: .bam 比对文件
# =============================================================================

set -euo pipefail

# 配置参数
INPUT_FASTQ="${1:-./fastq/reads.fastq}"
REF_GENOME="${2:-./reference/genome.fa}"
OUTPUT_BAM="${3:-./aligned.bam}"
THREADS="${4:-16}"

echo "============================================"
echo "Nanopore Alignment Pipeline"
echo "============================================"
echo "Input FASTQ: $INPUT_FASTQ"
echo "Reference: $REF_GENOME"
echo "Output: $OUTPUT_BAM"
echo "============================================"

# 检查输入文件
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input FASTQ not found: $INPUT_FASTQ"
    exit 1
fi

if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome not found: $REF_GENOME"
    exit 1
fi

# 检查minimap2索引
MMIDX="${REF_GENOME}.mmi"
if [ ! -f "$MMIDX" ]; then
    echo ""
    echo ">>> Creating minimap2 index..."
    minimap2 -d "$MMIDX" "$REF_GENOME"
fi

# 创建samtools索引
FAI="${REF_GENOME}.fai"
if [ ! -f "$FAI" ]; then
    echo ""
    echo ">>> Creating samtools index..."
    samtools faidx "$REF_GENOME"
fi

# 临时文件
TMP_BAM="${OUTPUT_BAM}.tmp"

# 使用minimap2进行比对
echo ""
echo ">>> Running Minimap2 alignment..."
minimap2 \
    -ax map-ont \
    -t "$THREADS" \
    -K 1G \
    --MD \
    -r 0 \
    "$REF_GENOME" \
    "$INPUT_FASTQ" | \
    samtools sort -@ "$THREADS" -o "$TMP_BAM"

# 移动到最终位置
mv "$TMP_BAM" "$OUTPUT_BAM"

# 创建索引
echo ""
echo ">>> Creating BAM index..."
samtools index "$OUTPUT_BAM"

# 统计比对结果
TOTAL_READS=$(samtools view -c "$OUTPUT_BAM")
MAPPED_READS=$(samtools view -c -F 4 "$OUTPUT_BAM")
UNMAPPED_READS=$(samtools view -c -f 4 "$OUTPUT_BAM")

echo ""
echo ">>> Alignment Complete!"
echo ">>> Total reads: $TOTAL_READS"
echo ">>> Mapped reads: $MAPPED_READS"
echo ">>> Unmapped reads: $UNMAPPED_READS"

# 检查修饰标签
MOD_COUNT=$(samtools view "$OUTPUT_BAM" | grep -c "MM:Z:" || echo "0")
echo ">>> Reads with modification tags: $MOD_COUNT"

echo ""
echo "============================================"
echo "Output: $OUTPUT_BAM"
echo "============================================"
