#!/bin/bash
# =============================================================================
# run_full_pipeline.sh
# 完整m6A检测流程 - 一键运行所有步骤
# =============================================================================
# 配置参数 (修改这里以适应你的数据)
# =============================================================================

set -euo pipefail

# 路径配置
WORK_DIR="$(cd "$(dirname "$0")" && pwd)"
POD5_INPUT="${1:-$WORK_DIR/../raw/reads.pod5}"
REF_GENOME="${2:-/path/to/genome.fa}"
OUTPUT_DIR="${3:-$WORK_DIR/results}"

# 参数配置
THREADS="${4:-16}"
MODEL="${5:-dna_r10.4.1_e8.2_400bps_hac}"
MOD_THRESHOLD="${6:-0.1}"
MIN_COVERAGE="${7:-10}"

echo "============================================"
echo "m6A Detection Full Pipeline"
echo "============================================"
echo "Work Directory: $WORK_DIR"
echo "Input POD5: $POD5_INPUT"
echo "Reference: $REF_GENOME"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Model: $MODEL"
echo "============================================"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/fastq"
mkdir -p "$OUTPUT_DIR/m6a"

# 检查参考基因组
if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome not found: $REF_GENOME"
    exit 1
fi

cd "$WORK_DIR"

# Step 1: 碱基识别
echo ""
echo "============================================"
echo "Step 1: Basecalling"
echo "============================================"
bash scripts/01_basecalling.sh \
    "$POD5_INPUT" \
    "$OUTPUT_DIR/fastq" \
    "$MODEL" \
    "$THREADS"

# Step 2: 序列比对
echo ""
echo "============================================"
echo "Step 2: Alignment"
echo "============================================"
bash scripts/02_alignment.sh \
    "$OUTPUT_DIR/fastq/reads.fastq" \
    "$REF_GENOME" \
    "$OUTPUT_DIR/aligned.bam" \
    "$THREADS"

# Step 3: m6A修饰分析
echo ""
echo "============================================"
echo "Step 3: m6A Modification Analysis"
echo "============================================"
bash scripts/03_modkit_pileup.sh \
    "$OUTPUT_DIR/aligned.bam" \
    "$REF_GENOME" \
    "$OUTPUT_DIR/m6a" \
    "$THREADS" \
    "$MOD_THRESHOLD" \
    "$MIN_COVERAGE"

echo ""
echo "============================================"
echo "Pipeline Complete!"
echo "============================================"
echo "Results:"
echo "  - FASTQ: $OUTPUT_DIR/fastq/reads.fastq"
echo "  - BAM: $OUTPUT_DIR/aligned.bam"
echo "  - m6A sites: $OUTPUT_DIR/m6a/m6A_pileup.csv"
echo ""
ls -lh "$OUTPUT_DIR/m6a/"
