#!/bin/bash
# =============================================================================
# 03_modkit_pileup.sh
# m6A修饰位点分析 - 使用modkit进行pileup分析
# =============================================================================
# 输入: .bam 比对文件
# 输出: .csv 修饰位点文件
# =============================================================================

set -euo pipefail

# 配置参数
INPUT_BAM="${1:-./aligned.bam}"
REF_GENOME="${2:-./reference/genome.fa}"
OUTPUT_DIR="${3:-./m6a_results}"
THREADS="${4:-16}"
MOD_THRESHOLD="${5:-0.1}"
MIN_COVERAGE="${6:-10}"

echo "============================================"
echo "m6A Modification Detection Pipeline"
echo "============================================"
echo "Input BAM: $INPUT_BAM"
echo "Reference: $REF_GENOME"
echo "Output: $OUTPUT_DIR"
echo "Mod Threshold: $MOD_THRESHOLD"
echo "Min Coverage: $MIN_COVERAGE"
echo "============================================"

mkdir -p "$OUTPUT_DIR"

# 检查输入文件
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome not found: $REF_GENOME"
    exit 1
fi

# 检查modkit是否安装
if ! command -v modkit &> /dev/null; then
    echo "Error: modkit not found. Please install:"
    echo "  pip install modkit"
    exit 1
fi

# 运行modkit pileup
echo ""
echo ">>> Running modkit pileup analysis..."

modkit pileup \
    "$INPUT_BAM" \
    "$OUTPUT_DIR" \
    --ref "$REF_GENOME" \
    --threads "$THREADS" \
    --mod-threshold "$MOD_THRESHOLD" \
    --min-coverage "$MIN_COVERAGE" \
    --prefix "m6A" \
    --log-level "INFO"

# 检查输出文件
OUTPUT_CSV="$OUTPUT_DIR/m6A_pileup.csv"
if [ -f "$OUTPUT_CSV" ]; then
    # 统计结果
    echo ""
    echo ">>> Analysis Complete!"
    
    # 总修饰位点数
    TOTAL_SITES=$(tail -n +2 "$OUTPUT_CSV" | wc -l)
    echo ">>> Total m6A sites: $TOTAL_SITES"
    
    # 高置信度位点 (freq > 0.5)
    HIGH_CONF=$(awk -F',' 'NR>1 && $6>0.5' "$OUTPUT_CSV" | wc -l)
    echo ">>> High confidence sites (freq>0.5): $HIGH_CONF"
    
    # 平均修饰频率
    AVG_FREQ=$(awk -F',' 'NR>1 {sum+=$6; count++} END {print sum/count}' "$OUTPUT_CSV")
    echo ">>> Average modification frequency: $AVG_FREQ"
    
    # 平均覆盖度
    AVG_COV=$(awk -F',' 'NR>1 {sum+=$7; count++} END {print sum/count}' "$OUTPUT_CSV")
    echo ">>> Average coverage: $AVG_COV"
else
    echo "Warning: Output CSV not found. Checking for alternative formats..."
    ls -la "$OUTPUT_DIR/"
fi

echo ""
echo "============================================"
echo "Output: $OUTPUT_DIR/"
echo "============================================"
