# Step 3: m6A修饰检测 (modkit)

## 目的

使用modkit对纳米孔测序数据进行m6A等修饰碱基的检测和分析。通过pileup分析统计每个潜在修饰位点的覆盖度和修饰频率。

## 工具

### modkit

ONT官方推荐的修饰检测工具。

**安装:**
```bash
# 方法1: pip安装 (推荐)
pip install modkit

# 方法2: conda安装
conda install -c bioconda modkit

# 方法3: 源码安装
git clone https://github.com/nanoporetech/modkit.git
cd modkit && pip install .
```

## 命令详解

### 基础命令

```bash
modkit pileup \
    aligned.bam \
    output_dir \
    --ref genome.fa \
    --threads 16 \
    --mod-threshold 0.1 \
    --min-coverage 10
```

### 关键参数

| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `--ref` | 参考基因组FASTA | (必需) |
| `--threads` | 线程数 | 16 |
| `--mod-threshold` | 最小修饰频率 | 0.1 |
| `--min-coverage` | 最小覆盖度 | 10 |
| `--prefix` | 输出文件前缀 | "modkit" |
| `--cpg` | 专门分析CpG位点 | 可选 |
| `--chunks` | 分块数 | 自动 |

### 参数详解

#### --mod-threshold

修饰频率阈值，低于此值的位点不报告。

```bash
# 严格模式 (高置信度)
modkit pileup ... --mod-threshold 0.5

# 平衡模式 (推荐)
modkit pileup ... --mod-threshold 0.1

# 宽松模式 (捕获更多)
modkit pileup ... --mod-threshold 0.05
```

#### --min-coverage

最小覆盖度要求，过滤低覆盖度位点。

```bash
# 高覆盖度要求
modkit pileup ... --min-coverage 50

# 默认值
modkit pileup ... --min-coverage 10

# 低覆盖度要求
modkit pileup ... --min-coverage 5
```

### pileup 子命令

```bash
# 标准pileup
modkit pileup input.bam output_dir --ref genome.fa

# 批量pileup
modkit pileup batch input_dir output_dir --ref genome.fa

# 生成统计摘要
modkit summary aligned.bam --ref genome.fa
```

### probs 子命令

提取修饰概率信息:

```bash
modkit probs aligned.bam output.tsv --ref genome.fa
```

## 输出文件

### m6A_pileup.csv

**列说明:**

| 列名 | 描述 | 示例 |
|------|------|------|
| chrom | 染色体 | chr1 |
| start | 起始位置 | 1234567 |
| end | 终止位置 | 1234568 |
| strand | 链方向 | + |
| mod_base | 修饰碱基 | m6A |
| count_modified | 修饰reads数 | 85 |
| count_unmodified | 未修饰reads数 | 15 |
| coverage | 总覆盖度 | 100 |
| freq | 修饰频率 | 0.85 |
| score | 置信度分数 | 25.5 |
| motif | 序列motif | DRACH |

### m6A_pileup.stats

**内容:**

```
Total sites: 10000
High confidence (freq > 0.5): 5000
Medium (0.1 < freq <= 0.5): 3000
Low (freq <= 0.1): 2000
```

## 修饰类型

### 支持的修饰碱基

| 修饰 | 标签 | 说明 |
|------|------|------|
| m6A | m6A | N6-甲基腺嘌呤 |
| m5C | m5C | 5-甲基胞嘧啶 |
| m4C | m4C | N4-甲基胞嘧啶 |
| m3C | m3C | N3-甲基胞嘧啶 |
| hm5C | hm5C | 5-羟甲基胞嘧啶 |

### 指定修饰类型

```bash
# 只分析m6A
modkit pileup ... --mod-threshold 0.1

# 只分析m5C
modkit pileup ... --only-mods m5C

# 分析所有修饰
modkit pileup ... --only-mods m6A,m5C
```

## 结果解读

### 修饰频率分布

```bash
# 统计修饰频率
awk -F',' 'NR>1 {print $9}' m6A_pileup.csv | \
    sort | uniq -c | sort -k2 -n

# 高置信度位点
awk -F',' 'NR>1 && $9>0.5' m6A_pileup.csv | wc -l
```

### 覆盖度分析

```bash
# 平均覆盖度
awk -F',' 'NR>1 {sum+=$11} END {print sum/NR}' m6A_pileup.csv

# 高覆盖度位点 (>100x)
awk -F',' 'NR>1 && $11>100' m6A_pileup.csv | wc -l
```

### Motif分析

```bash
# 提取motif分布
awk -F',' 'NR>1 {print $12}' m6A_pileup.csv | \
    sort | uniq -c | sort -rn | head -10

# 典型m6A motif: DRACH (D=A/G/T, R=A/G, A=A, C=C, H=A/C/T)
```

## 可视化

### 曼哈顿图

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('m6A_pileup.csv')

# 简化曼哈顿图
plt.figure(figsize=(12, 4))
plt.scatter(df['start'], df['freq'], c='blue', alpha=0.5, s=1)
plt.xlabel('Position')
plt.ylabel('Modification Frequency')
plt.title('m6A Distribution')
plt.savefig('manhattan.png', dpi=300)
```

### 覆盖度分布

```python
plt.figure(figsize=(12, 4))
plt.hist(df['coverage'], bins=50, alpha=0.7)
plt.xlabel('Coverage')
plt.ylabel('Count')
plt.title('Coverage Distribution')
plt.savefig('coverage.png', dpi=300)
```

## 常见问题

### 问题1: 无修饰位点输出

**检查项:**
1. 确认BAM中有MM标签
2. 检查覆盖度是否足够
3. 确认阈值设置

```bash
# 检查MM标签
samtools view aligned.bam | grep -c "MM:Z:"

# 降低阈值
modkit pileup ... --mod-threshold 0.05 --min-coverage 5
```

### 问题2: 修饰频率过低

**原因:**
- 覆盖度不够
- 样本本身修饰较少
- 碱基识别模型问题

**解决方案:**
```bash
# 增加覆盖度要求
modkit pileup ... --mod-threshold 0.05

# 检查碱基识别是否正确
samtools view aligned.bam | head | grep "MM:Z:"
```

### 问题3: 运行时间过长

**解决方案:**
```bash
# 增加线程数
modkit pileup ... --threads 32

# 减少分块
modkit pileup ... --chunks 4
```

### 问题4: 内存不足

**解决方案:**
```bash
# 减少线程数
modkit pileup ... --threads 8

# 增加分块数(减少内存)
modkit pileup ... --chunks 16
```

### 问题5: 与已知m6A位点不匹配

**可能原因:**
1. 参考基因组版本不同
2. 阈值设置过于严格
3. 覆盖度不够

**解决:**
- 使用相同基因组版本
- 降低阈值重新分析
- 增加测序深度

## 参数优化

### 高精度模式

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.5 \
    --min-coverage 30 \
    --threads 32
```

### 灵敏度模式

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.05 \
    --min-coverage 10 \
    --threads 16
```

### 平衡模式

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.1 \
    --min-coverage 20 \
    --threads 16
```

## 验证

### 与其他方法比较

- **MeRIP-seq**: 传统免疫沉淀方法
- **m6A-seq**: m6A特异性测序
- **SCARLET**: 单分子验证金标准

### 统计检验

```python
# 计算修饰率置信区间
import scipy.stats as stats

coverage = df['coverage'].values
freq = df['freq'].values

# 二项分布置信区间
for i in range(len(coverage)):
    n = coverage[i]
    k = int(n * freq[i])
    ci = stats.binom.interval(0.95, n, freq[i])
    print(f"Site {i}: {freq[i]:.3f} (95% CI: {ci[0]/n:.3f}-{ci[1]/n:.3f})")
```

## 引用

- **modkit**: Oxford Nanopore Technologies. https://github.com/nanoporetech/modkit
- **原始文献**: Liu, H. et al. (2024). Accurate detection of m6A RNA modifications in native nanopore sequencing. Nature Methods.
