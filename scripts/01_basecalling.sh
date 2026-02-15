#!/bin/bash
# =============================================================================
# 01_basecalling.sh
# 纳米孔测序碱基识别 - 将原始信号转为序列
# =============================================================================
# 输入: .pod5 原始信号文件
# 输出: .fastq 序列文件
# =============================================================================

set -euo pipefail

# 配置参数
INPUT_POD5="${1:-./reads.pod5}"
OUTPUT_DIR="${2:-./fastq}"
MODEL="${3:-dna_r10.4.1_e8.2_400bps_hac}"
THREADS="${4:-16}"

echo "============================================"
echo "Nanopore Basecalling Pipeline"
echo "============================================"
echo "Input: $INPUT_POD5"
echo "Model: $MODEL"
echo "Output: $OUTPUT_DIR"
echo "============================================"

mkdir -p "$OUTPUT_DIR"

# 检查Dorado是否安装
if ! command -v dorado &> /dev/null; then
    echo "Error: Dorado not found. Please install Dorado:"
    echo "  https://github.com/nanoporetech/dorado"
    exit 1
fi

# 检查pod5文件
if [ ! -d "$INPUT_POD5" ] && [ ! -f "$INPUT_POD5" ]; then
    echo "Error: Input pod5 not found: $INPUT_POD5"
    exit 1
fi

# Dorado碱基识别
# 注意: 使用--emit-moves保留修饰碱基信息
echo ""
echo ">>> Running Dorado Basecalling..."

dorado basecaller \
    --emit-moves \
    --mod-base-models "$MODEL" \
    "$MODEL" \
    "$INPUT_POD5" \
    --threads "$THREADS" \
    > "$OUTPUT_DIR/basecalls.bam"

# 转换为fastq
echo ""
echo ">>> Converting to FASTQ..."
samtools fastq "$OUTPUT_DIR/basecalls.bam" > "$OUTPUT_DIR/reads.fastq"

# 统计reads数量
READ_COUNT=$(grep -c "^@" "$OUTPUT_DIR/reads.fastq")
echo ""
echo ">>> Basecalling Complete!"
echo ">>> Total reads: $READ_COUNT"

# 清理中间文件
rm -f "$OUTPUT_DIR/basecalls.bam"

echo ""
echo "============================================"
echo "Output: $OUTPUT_DIR/reads.fastq"
echo "============================================"
