# Step 1: 碱基识别 (Basecalling)

## 目的

将纳米孔测序产生的原始电信号数据(.pod5)转换为核苷酸序列(.fastq)。现代碱基识别器可以直接识别并标记修饰碱基(m6A等)。

## 工具选择

### Dorado (推荐)

ONT官方推荐的碱基识别工具，支持GPU加速。

**安装:**
```bash
# 下载预编译版本
wget https://github.com/nanoporetech/dorado/releases/download/v0.8.1/dorado-0.8.1-linux-x64.tar.gz
tar -xzf dorado-0.8.1-linux-x64.tar.gz
export PATH=$PATH:/path/to/dorado-0.8.1-linux-x64/bin

# 验证安装
dorado --version
```

### Guppy

ONT官方工具，GPU和CPU版本都可用。

```bash
# conda安装
conda install -c conda-forge ont-guppy

# 或使用docker
docker pull ontporetech/guppy:7.0.0
```

## 模型选择

### R10.4.1 (最新)

| 模型 | 精度 | 速度 | 修饰检测 |
|------|------|------|----------|
| dna_r10.4.1_e8.2_400bps_hac | 高 | 中 | 支持 |
| dna_r10.4.1_e8.2_400bps_sup | 最高 | 慢 | 支持 |
| dna_r10.4.1_e8.2_400bps_fast | 中 | 快 | 支持 |

### R9.4.1 (经典)

| 模型 | 精度 | 速度 |
|------|------|------|
| dna_r9.4.1_e8_hac | 高 | 中 |
| dna_r9.4.1_e8_sup | 最高 | 慢 |
| dna_r9.4.1_e8_fast | 中 | 快 |

## 命令详解

### 基础命令

```bash
dorado basecaller \
    --emit-moves \
    --mod-base-models dna_r10.4.1_e8.2_400bps_hac \
    dna_r10.4.1_e8.2_400bps_hac \
    reads.pod5 > basecalls.bam
```

### 关键参数

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `--emit-moves` | 保留修饰碱基信息 | 必须加 |
| `--mod-base-models` | 启用修饰检测模型 | dna_r10.4.1_e8.2_400bps_hac |
| `--threads` | 线程数 | 16-32 |
| `--batch-size` | 批处理大小 | 1000 |
| `--min-qscore` | 最小质量分数 | 10 |

### 输出格式转换

```bash
# BAM转FASTQ
samtools fastq basecalls.bam > reads.fastq

# 直接输出FASTQ (需较新版本)
dorado basecaller ... --fastq > reads.fastq
```

## 输入文件格式

### .pod5 文件

**文件结构:**
- 二进制格式
- 包含原始电流信号
- 每个read有对应的时间戳和通道信息

**查看信息:**
```bash
pod5 view --stats reads.pod5
pod5 count reads.pod5
```

**文件大小:**
- 每个read约2-5KB
- 100万reads约2-5GB

## 修饰检测原理

### Mm/MI标签

Dorado在碱基识别过程中会生成修饰信息:

```
MM:Z:m+6A    # m6A修饰
ML:B:C,255   # 置信度分数
```

**标签说明:**
- `MM`: 修饰碱基类型和位置
- `ML`: 修饰置信度 (0-255)

## 常见问题

### 问题1: GPU内存不足

**症状:**
```
CUDA out of memory
```

**解决方案:**
```bash
# 减少batch-size
dorado basecaller --batch-size 500 ...

# 或使用CPU模式
dorado basecaller --device cpu ...
```

### 问题2: 速度太慢

**解决方案:**
```bash
# 使用fast模型
dorado basecaller dna_r10.4.1_e8.2_400bps_fast ...

# 增加batch-size
dorado basecaller --batch-size 2000 ...
```

### 问题3: 修饰检测失败

**检查项:**
1. 确认使用--emit-moves参数
2. 确认使用支持modification的模型
3. 检查Dorado版本(需>=0.7.0)

```bash
# 检查版本
dorado --version

# 验证修饰标签
samtools view basecalls.bam | head | grep "MM:Z:"
```

### 问题4: pod5文件损坏

**检测:**
```bash
pod5 view --stats reads.pod5

# 修复(如有必要)
pod5 validate reads.pod5
```

## 性能基准

| 测序仪 | GPU | 模型 | 速度(reads/s) |
|--------|-----|------|---------------|
| PromethION | A100 | hac | ~3000 |
| PromethION | A100 | sup | ~1500 |
| MinION | RTX3090 | hac | ~450 |
| MinION | RTX3090 | fast | ~800 |

## 输出验证

```bash
# 统计reads数量
samtools view -c basecalls.bam

# 查看序列质量分布
samtools fastq basecalls.bam | \
    awk 'NR%4==0' | \
    tr -d '\n' | \
    fold -w 1 | \
    sort | \
    uniq -c | \
    sort -rn

# 检查修饰标签
samtools view basecalls.bam | \
    grep -o 'MM:Z:[^ ]*' | \
    sort | \
    uniq -c
```
