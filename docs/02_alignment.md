# Step 2: 序列比对 (Alignment)

## 目的

将碱基识别产生的序列reads比对到参考基因组，生成包含位置信息的BAM文件。此步骤需要保留修饰碱基的MM/MI标签信息，同时进行质量控制。

## 完整工作流程

```
输入FASTQ → pychopper修剪 → NanoFilt过滤 → minimap2比对 → SAM→BAM → 排序 → 去重 → 索引
```

## 工具

### Minimap2 (推荐)

专为长读长序列设计的快速比对工具。

**安装:**
```bash
conda install -c bioconda minimap2
# 或编译安装
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```

## 命令详解

### 完整比对流程

```bash
# 1. pychopper修剪 (去除低质量序列和接头)
pychopper -t 8 input.fastq output.fq

# 2. NanoFilt过滤 (质量过滤和长度过滤)
NanoFilt -l 50 -q 7 output.fq > output.clear.fq

# 3. minimap2比对
minimap2 -ax map-ont -t 8 -uf -k 14 \
    reference.fa \
    output.clear.fq > output.sam

# 4. SAM转BAM
samtools view -@ 8 -Sb output.sam > output.unsorted.bam

# 5. 排序
samtools sort -@ 8 -o output.sorted.bam -T tmp output.unsorted.bam

# 6. 去重
samtools markdup -@ 8 -r output.sorted.bam output.rmdup.bam

# 7. 索引
samtools index output.rmdup.bam
```

### minimap2 关键参数

```bash
minimap2 -ax map-ont -t 8 -uf -k 14 \
    reference.fa \
    reads.fastq > output.sam
```

| 参数 | 说明 | 推荐值 | 解释 |
|------|------|--------|------|
| `-ax map-ont` | 比对模式 | 必须 | 专为ONT长读长优化 |
| `-t` | 线程数 | 8-16 | 根据CPU核心数调整 |
| `-uf` | 只比对正向链 | 可选 | 用于RNA-seq |
| `-k` | k-mer大小 | 14 | R10.4用14,R9.4用15 |
| `-K` | 批处理大小 | 1G | 影响内存和速度 |
| `--MD` | 保留MD标签 | 推荐 | 用于修饰检测 |
| `-r` | 误差率 | 0 | 使用实际误差率 |

#### 参数详解

**-ax map-ont**
```
可用模式:
- map-ont:    ONT长读长 (推荐)
- map-hifi:   PacBio HiFi读长
- map-pb:     PacBio CLR读长
- asm20:      基因组组装 (低一致性)
- asm5:       基因组组装 (中等一致性)
- splice:     有剪接感知 (RNA-seq)
```

**-k 参数选择**
```bash
# R10.4 测序仪 (更高精度)
minimap2 -k 14 ...

# R9.4 测序仪
minimap2 -k 15 ...

# 更短的k-mer可提高灵敏度，但可能增加错误
minimap2 -k 12 ...
```

**-uf 参数 (Forward-only)**
```
-uf: 只比对到参考基因组的正向链
用途:
- RNA-seq:  只保留转录本比对
- 避免比对到反义链的干扰
```

### pychopper 修剪

**作用:** 去除接头序列和低质量区域

```bash
pychopper -t 8 input.fastq output.fq
```

| 参数 | 说明 |
|------|------|
| `-t` | 线程数 |
| `-r` | 报告文件 |
| `-S` | 统计文件 |

### NanoFilt 过滤

**作用:** 质量过滤和长度过滤

```bash
NanoFilt -l 50 -q 7 input.fq > output.fq
```

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `-l` | 最小长度 | 50-100 |
| `-q` | 最小质量分数 | 7-10 |
| `--head` | 保留前N条reads | - |
| `--tail` | 保留后N条reads | - |

### SAM处理

```bash
# SAM转BAM
samtools view -@ 8 -Sb input.sam > input.bam

# 排序
samtools sort -@ 8 -o sorted.bam -T tmp input.bam

# 去重
samtools markdup -@ 8 -r sorted.bam rmdup.bam

# 索引
samtools index rmdup.bam
```

## 索引准备

### 创建参考基因组索引

```bash
# minimap2索引 (推荐预建)
minimap2 -d genome.fa.mmi genome.fa

# samtools索引
samtools faidx genome.fa
```

### 检查索引

```bash
# 检查索引文件
ls -lh genome.fa.mmi genome.fa.fai

# 验证索引可用
samtools faidx genome.fa chr1:1-1000
```

## 输出格式

### BAM文件结构

| 字段 | 说明 | 示例 |
|------|------|------|
| QNAME | Read名称 | read_001 |
| FLAG | 比对状态 | 0 (正向) / 16 (反向) |
| RNAME | 染色体 | chr1 |
| POS | 起始位置 | 1234567 |
| MAPQ | 比对质量 | 60 |
| CIGAR | 比对详情 | 150M |
| SEQ | 序列 | ATCG... |
| QUAL | 质量分数 | IIII... |
| TAGS | 修饰标签 | MM:Z:m+6A |

### FLAG说明

| FLAG | 含义 |
|------|------|
| 0 | 正向比对 |
| 16 | 反向比对 |
| 4 | 未比对 |
| 2048 | PCR重复 |

### 修饰标签

```
MM:Z:m+6A      # m6A修饰类型
ML:B:C,200     # 置信度 (0-255)
MP:A:+,10,15   # 修饰位置
```

## 统计信息

### 比对率统计

```bash
# 总reads数
samtools view -c aligned.bam

# 比对上的reads
samtools view -c -F 4 aligned.bam

# 未比对的reads
samtools view -c -f 4 aligned.bam

# 比对率计算
echo "scale=2; $(samtools view -c -F 4 aligned.bam) / $(samtools view -c aligned.bam) * 100" | bc
```

### 覆盖度统计

```bash
# 平均覆盖度
samtools depth aligned.bam | awk '{sum+=$3} END {print sum/NR}'

# 覆盖度分布
samtools depth aligned.bam | awk '{print $3}' | sort -n | uniq -c

# 使用bamCoverage
bamCoverage -b aligned.bam -of bigwig -o coverage.bw
```

### 读长分布

```bash
# 平均读长
samtools view aligned.bam | awk '{print length($10)}' | awk '{sum+=$1} END {print sum/NR}'

# 读长直方图
samtools view aligned.bam | awk '{print length($10)}' | \
    sort | uniq -c | awk '{print $2"\t"$1}' | head -20
```

### 去重统计

```bash
# 使用flagstat
samtools flagstat aligned.rmdup.bam

# 输出示例:
# 1000000 + 0 in total (QC-passed reads + QC-failed reads)
# 5000 + 0 secondary
# 0 + 0 supplementary
# 950000 + 0 mapped (95.00% : N/A)
# 10000 + 0 duplicates
```

## 常见问题

### 问题1: 比对率低

**可能原因:**
1. 测序质量差
2. 参考基因组不匹配
3. 序列太短

**解决方案:**
```bash
# 使用更宽松参数
minimap2 -ax map-ont -k 12 -w 5 ...

# 降低最小长度过滤
NanoFilt -l 30 -q 5 ...
```

### 问题2: 内存不足

**解决方案:**
```bash
# 减小线程数
minimap2 -t 4 ...

# 减小batch-size
minimap2 -K 100M ...

# 分批处理
split -l 10000 reads.fastq reads_
for f in reads_*; do
    minimap2 ... $f >> output.sam
done
```

### 问题3: 修饰标签丢失

**检查:**
```bash
# 检查MM标签
samtools view aligned.bam | grep "MM:Z:" | head

# 检查MD标签
samtools view aligned.bam | grep "MD:Z:" | head
```

**原因与解决:**
- minimap2版本过旧 → 升级到最新版本
- 未使用--MD参数 → 添加--MD参数
- 比对后格式转换丢失 → 检查转换步骤

### 问题4: 去重后覆盖度降低太多

**说明:** 纳米孔数据本身PCR重复率较低，过多去重可能是因为:
- 同一reads多次比对
- 建库PCR扩增导致

**调整:**
```bash
# 只标记不删除
samtools markdup -@ 8 aligned.bam output.bam

# 或跳过去重步骤
```

### 问题5: CIGAR错误

**说明:** minimap2使用简化CIGAR，部分工具不兼容

**验证:**
```bash
# 检查CIGAR格式
samtools view aligned.bam | head | awk '{print $6}'
```

## 完整参数模板

### for RNA-seq (m6A检测)

```bash
#!/bin/bash
# ONT RNA-seq m6A检测比对流程

REF="genome.fa"
THREADS=16

# 1. 预处理
pychopper -t "$THREADS" input.fastq tmp.fq
NanoFilt -l 50 -q 7 tmp.fq > clean.fq

# 2. 比对
minimap2 -ax map-ont -t "$THREADS" -uf -k 14 \
    "$REF" \
    clean.fq > aligned.sam

# 3. 处理BAM
samtools view -@ "$THREADS" -Sb aligned.sam | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" -r - \
    aligned.rmdup.bam

# 4. 索引
samtools index aligned.rmdup.bam

# 清理
rm -f tmp.fq clean.fq aligned.sam
```

### for DNA-seq (直接测序)

```bash
#!/bin/bash
# ONT DNA-seq 比对流程

REF="genome.fa"
THREADS=16

# 1. 过滤
NanoFilt -l 100 -q 10 input.fastq > clean.fq

# 2. 比对 (不使用-uf)
minimap2 -ax map-ont -t "$THREADS" -k 14 \
    "$REF" \
    clean.fq > aligned.sam

# 3. 处理BAM
samtools view -@ "$THREADS" -Sb aligned.sam | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" -r - \
    aligned.rmdup.bam

samtools index aligned.rmdup.bam
```

## 验证检查清单

```bash
# 1. 检查BAM文件
samtools flagstat aligned.rmdup.bam

# 2. 检查修饰标签
echo "MM tags:"
samtools view aligned.rmdup.bam | grep -c "MM:Z:"

# 3. 检查覆盖度
samtools depth aligned.rmdup.bam | awk '{sum+=$3} END {print "Avg coverage:", sum/NR}'

# 4. 检查bigWig
ls -lh processed/*.bw

# 5. 查看reads分布
samtools view -c aligned.rmdup.bam
```

## 与其他工具比较

| 工具 | 优点 | 缺点 |
|------|------|------|
| Minimap2 | 快速，专为ONT优化 | RNA比对需额外参数 |
| BWA-MEM | 兼容性好 | 对长读长不够优化 |
| GraphMap | 高精度 | 速度慢 |
| NGMLR | 适合变异检测 | 速度慢 |
