# Nanopore m6A Detection Pipeline

针对纳米孔测序数据的m6A修饰检测完整流程，从原始信号数据到修饰位点鉴定的全流程分析方案。

## 概述

本流程实现从Nanopore原始信号数据(.pod5)到m6A修饰位点鉴定的完整分析:

```
原始信号(.pod5) → 碱基识别 → 序列比对 → Mm/MI标记 → modkit分析 → m6A位点
```

### 流程步骤

| 步骤 | 软件 | 输入 | 输出 | 描述 |
|------|------|------|------|------|
| 1 | Dorado/Guppy | .pod5 | .fastq | 碱基识别 |
| 2 | Minimap2 | .fastq + genome | .bam | 序列比对 |
| 3 | modkit | .bam | .csv | 修饰位点鉴定 |

## 目录结构

```
m6a_detection_pipeline/
├── scripts/
│   ├── 01_basecalling.sh          # 碱基识别
│   ├── 02_alignment.sh           # 序列比对
│   ├── 03_modkit_pileup.sh      # m6A修饰分析
│   └── run_full_pipeline.sh     # 一键运行
├── docs/
│   ├── 01_basecalling.md         # 碱基识别详解
│   ├── 02_alignment.md           # 比对详解
│   ├── 03_modkit.md             # modkit详解
│   └── troubleshooting.md       # 故障排除
├── reference/                     # 参考基因组
├── README.md
└── .gitignore
```

## 安装依赖

### 方式一: Conda环境 (推荐)

```bash
# 创建环境
conda create -n m6a python=3.10
conda activate m6a

# 安装核心软件
conda install -c conda-forge numpy pandas
conda install -c bioconda minimap2 samtools bedtools

# 安装Dorado (ONT官方)
# 下载地址: https://github.com/nanoporetech/dorado
# 或使用容器版本
```

### 方式二: Docker/Singularity

```bash
# 使用ONT官方容器
docker pull nanotechpore/dorado:0.8.1

# 或Singularity
singularity pull docker://nanotechpore/dorado:0.8.1
```

### 软件版本要求

| 软件 | 最低版本 | 推荐版本 |
|------|----------|----------|
| Dorado | 0.7.0 | 0.8.1 |
| Guppy | 6.0.0 | 7.0.0 |
| minimap2 | 2.24 | 2.26 |
| samtools | 1.15 | 1.19 |
| modkit | 0.2.0 | 0.5.0 |

## 快速开始

### 完整流程

```bash
# 修改配置参数
vim scripts/run_full_pipeline.sh

# 一键运行
bash scripts/run_full_pipeline.sh

# 或分步运行
bash scripts/01_basecalling.sh
bash scripts/02_alignment.sh
bash scripts/03_modkit_pileup.sh
```

### 分步运行

```bash
# Step 1: 碱基识别
bash scripts/01_basecalling.sh \
    /path/to/reads.pod5 \
    /path/to/output/reads \
    dna_r10.4.1_e8.2_400bps_hac

# Step 2: 序列比对
bash scripts/02_alignment.sh \
    /path/to/reads.fastq \
    /path/to/reference/genome.fa \
    /path/to/output/aligned.bam

# Step 3: m6A修饰分析
bash scripts/03_modkit_pileup.sh \
    /path/to/aligned.bam \
    /path/to/reference/genome.fa \
    /path/to/output/m6A_sites.csv
```

## 输入文件格式

### 原始信号文件 (.pod5)

Nanopore测序仪产生的电信号数据文件。

**格式说明:**
- 扩展名: `.pod5`
- 内容: 纳米孔电信号强度数据
- 生成: MinKNOW软件或MinKNOW API

**获取方式:**
- 直接从MinKNOW输出目录获取
- 从ONT Rerio数据库下载

```bash
# 查看文件信息
pod5 view --stats reads.pod5

# 统计reads数量
pod5 count reads.pod5
```

### 参考基因组 (.fa/.fasta)

标准FASTA格式的参考基因组序列。

**格式要求:**
```
>chr1
ATGCGCTAGCTAGCTAGCTAGCTAGCTAGCTA...
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA...
```

**准备步骤:**
```bash
# 创建索引
samtools faidx genome.fa

# 创建mmi索引 (用于minimap2)
minimap2 -d genome.fa.mmi genome.fa
```

## 输出文件格式

### 碱基识别结果 (.fastq)

**格式:**
```
@read_id
ATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIII
```

**质量分数:** Phred score格式

### 比对结果 (.bam)

**格式:** 二进制SAM格式

**关键标签:**
- `MM`: 修饰碱基信息 (如 `MM:Z:m6A`)
- `ML`: 修饰置信度分数
- `MP`: 修饰位置信息

```bash
# 查看修饰标签
samtools view -e 'MM:Z:*' aligned.bam | head

# 查看m6A标记
samtools view aligned.bam | grep -o 'MM:Z:m6A[^,]*'
```

### m6A修饰位点 (.csv)

**输出列说明:**

| 列名 | 描述 | 示例 |
|------|------|------|
| chrom | 染色体 | chr1 |
| start | 起始位置 | 1234567 |
| end | 终止位置 | 1234568 |
| strand | 链方向 | + |
| mod_base | 修饰碱基 | m6A |
| freq | 修饰频率 | 0.85 |
| coverage | 覆盖度 | 150 |
| score | 置信度分数 | 25.5 |
| motif | 序列motif | DRACH |

**示例:**
```
chrom,start,end,strand,mod_base,freq,coverage,score,motif
chr1,1234567,1234568,+,m6A,0.85,150,25.5,GGA
chr1,2345678,2345679,+,m6A,0.72,98,18.3,GGAC
```

## 脚本详解

### 01_basecalling.sh

碱基识别脚本，将原始信号转为序列。

**参数:**
| 参数 | 说明 | 示例 |
|------|------|------|
| `$1` | 输入pod5文件/目录 | `./reads/` |
| `$2` | 输出目录 | `./fastq/` |
| `$3` | 模型名称 | `dna_r10.4.1_e8.2_400bps_hac` |

**模型选择:**
- `dna_r10.4.1_e8.2_400bps_hac` - 高精度(推荐)
- `dna_r10.4.1_e8.2_400bps_sup` - 超高精度
- `dna_r9.4.1_e8_hac` - R9.4.1通用

**启用修饰检测:**
```bash
# Dorado (>=0.7.0)
dorado basecaller \
    --emit-moves \
    --mod-base-models dna_r10.4.1_e8.2_400bps_hac \
    $MODEL \
    reads.pod5 > basecalls.bam

# 转换为fastq
samtools fastq basecalls.bam > reads.fastq
```

**注意事项:**
- 确保GPU内存充足(推荐24GB+)
- 使用--emit-moves保留修饰信息
- 高modification model会显著增加运行时间

---

### 02_alignment.sh

序列比对脚本，将reads比对到参考基因组。

**参数:**
| 参数 | 说明 |
|------|------|
| `$1` | 输入fastq文件 |
| `$2` | 参考基因组 |
| `$3` | 输出bam文件 |

**命令示例:**
```bash
minimap2 -ax map-ont \
    -t 16 \
    genome.fa \
    reads.fastq | \
    samtools sort -@ 16 -o aligned.bam

# 创建索引
samtools index aligned.bam
```

**关键参数:**
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `-ax map-ont` | 模式:ONT长读长 | - |
| `-t` | 线程数 | 16-32 |
| `-K` | 批处理大小 | 1G |
| `--MD` | 记录MD标签 | - |

**注意事项:**
- 使用`-ax map-ONT`而非`map-ont`
- 必须对bam文件建索引
- 建议保留MD标签用于后续分析

---

### 03_modkit_pileup.sh

m6A修饰位点分析脚本。

**参数:**
| 参数 | 说明 |
|------|------|
| `$1` | 输入bam文件 |
| `$2` | 参考基因组 |
| `$3` | 输出目录 |

**命令:**
```bash
modkit pileup \
    aligned.bam \
    output_dir \
    --ref genome.fa \
    --threads 16 \
    --cpg \
    --mod-threshold 0.1
```

**关键参数:**
| 参数 | 说明 | 推荐值 |
|------|------|--------|
| `--ref` | 参考基因组 | (必需) |
| `--threads` | 线程数 | 16 |
| `--cpg` | 分析CpG位点 | - |
| `--mod-threshold` | 最小修饰频率 | 0.1 |
| `--prefix` | 输出文件前缀 | - |

**输出文件:**
- `modkit_pileup.csv` - 修饰位点列表
- `modkit_pileup.stats` - 统计信息

## 常见问题

### 1. 碱基识别速度慢

**原因:**
- GPU性能不足
- 模型太大
- pod5文件太多小文件

**解决方案:**
```bash
# 使用更快模式
dorado basecaller --fast $MODEL reads.pod5 > fastq

# 批量处理
ls *.pod5 | xargs -I{} dorado basecaller $MODEL {} > {}.fastq
```

### 2. 修饰检测失败

**检查项:**
```bash
# 确认MM标签存在
samtools view aligned.bam | grep -c "MM:Z:"

# 检查modkit版本
modkit --version
```

**解决方案:**
- 使用支持modification的模型
- 确保使用--emit-moves参数

### 3. 修饰频率过低

**原因分析:**
- 覆盖度太低
- 阈值设置过高
- 样本本身m6A较少

**调整方法:**
```bash
# 降低阈值
modkit pileup input.bam output --mod-threshold 0.05

# 增加覆盖度要求
modkit pileup input.bam output --min-coverage 20
```

### 4. 内存不足

**解决方案:**
```bash
# 减少线程数
minimap2 -t 8 ...

# 分批处理
split -l 10000 reads.fastq reads_chunk_
for f in reads_chunk_*; do minimap2 ... $f; done
```

### 5. BAM文件过大

**解决方案:**
```bash
# 使用CRAM格式
samtools view -C aligned.bam -o aligned.cram

# 压缩
samtools sort -@ 16 -o aligned.sorted.bam aligned.bam
samtools index aligned.sorted.bam
```

## 结果解读

### 修饰频率分布

```bash
# 统计修饰频率
awk -F',' 'NR>1 {print $6}' m6A_sites.csv | \
    sort | uniq -c | sort -k2 -n
```

### 覆盖度分布

```bash
# 统计不同覆盖度区间
awk -F',' 'NR>1 {if($7>100) print "high"; else if($7>50) print "med"; else print "low"}' \
    m6A_sites.csv | sort | uniq -c
```

### 保守motif分析

```bash
# 提取motif分布
awk -F',' 'NR>1 {print $9}' m6A_sites.csv | \
    sort | uniq -c | sort -rn | head -10
```

## 验证与优化

### 与其他方法比较

可与以下方法结果比较验证:
- **MeRIP-seq**: 传统m6A检测方法
- **m6A-seq**: 免疫沉淀方法
- **SCARLET**: 单分子验证

### 参数优化建议

```bash
# 提高精度
modkit pileup input.bam output \
    --mod-threshold 0.2 \
    --min-coverage 30 \
    --refine-boundaries

# 平衡精度与覆盖
modkit pileup input.bam output \
    --mod-threshold 0.1 \
    --min-coverage 20
```

## 引用

- **Dorado**: Oxford Nanopore Technologies. https://github.com/nanoporetech/dorado
- **modkit**: Foundation for Applied Bioinformatics. https://github.com/nanoporetech/modkit
- **minimap2**: Li H. (2018). Minimap2: pairwise alignment for genomic sequences. Journal of Computational Biology.

## 许可证

MIT License

## 联系方式

问题反馈: https://github.com/IceBear321/Nanopore_m6A_Detection/issues
