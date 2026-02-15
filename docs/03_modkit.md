# Step 3: m6A Modification Detection (modkit)

## Purpose

Use modkit to detect and analyze m6A and other modified bases in Nanopore sequencing data. Pileup analysis calculates coverage and modification frequency at each potential modification site.

## Tool

**modkit**: OFFicially recommended modification detection tool by ONT.

**Installation:**
```bash
# Method 1: pip installation (recommended)
pip install modkit

# Method 2: conda installation
conda install -c bioconda modkit

# Method 3: Source installation
git clone https://github.com/nanoporetech/modkit.git
cd modkit && pip install .
```

## Command Details

### Basic Command

```bash
modkit pileup \
    aligned.bam \
    output_dir \
    --ref genome.fa \
    --threads 16 \
    --mod-threshold 0.1 \
    --min-coverage 10
```

### Key Parameters

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `--ref` | Reference genome FASTA | (required) |
| `--threads` | Thread count | 16 |
| `--mod-threshold` | Minimum modification frequency | 0.1 |
| `--min-coverage` | Minimum coverage | 10 |
| `--prefix` | Output file prefix | "modkit" |
| `--cpg` | Analyze CpG sites specifically | optional |
| `--chunks` | Number of chunks | auto |

### Parameter Details

#### --mod-threshold

Modification frequency threshold, sites below this value are not reported.

```bash
# Strict mode (high confidence)
modkit pileup ... --mod-threshold 0.5

# Balanced mode (recommended)
modkit pileup ... --mod-threshold 0.1

# Loose mode (capture more)
modkit pileup ... --mod-threshold 0.05
```

#### --min-coverage

Minimum coverage requirement, filters low-coverage sites.

```bash
# High coverage requirement
modkit pileup ... --min-coverage 50

# Default value
modkit pileup ... --min-coverage 10

# Low coverage requirement
modkit pileup ... --min-coverage 5
```

### pileup Subcommand

```bash
# Standard pileup
modkit pileup input.bam output_dir --ref genome.fa

# Batch pileup
modkit pileup batch input_dir output_dir --ref genome.fa

# Generate summary statistics
modkit summary aligned.bam --ref genome.fa
```

### probs Subcommand

Extract modification probability information:

```bash
modkit probs aligned.bam output.tsv --ref genome.fa
```

## Output Files

### m6A_pileup.csv

**Column Description:**

| Column | Description | Example |
|--------|-------------|---------|
| chrom | Chromosome | chr1 |
| start | Start position | 1234567 |
| end | End position | 1234568 |
| strand | Strand direction | + |
| mod_base | Modified base | m6A |
| count_modified | Number of modified reads | 85 |
| count_unmodified | Number of unmodified reads | 15 |
| coverage | Total coverage | 100 |
| freq | Modification frequency | 0.85 |
| score | Confidence score | 25.5 |
| motif | Sequence motif | DRACH |

### m6A_pileup.stats

**Content:**

```
Total sites: 10000
High confidence (freq > 0.5): 5000
Medium (0.1 < freq <= 0.5): 3000
Low (freq <= 0.1): 2000
```

## Modification Types

### Supported Modified Bases

| Modification | Label | Description |
|-------------|-------|-------------|
| m6A | m6A | N6-methyladenosine |
| m5C | m5C | 5-methylcytosine |
| m4C | m4C | N4-methylcytosine |
| m3C | m3C | N3-methylcytosine |
| hm5C | hm5C | 5-hydroxymethylcytosine |

### Specify Modification Type

```bash
# Analyze only m6A
modkit pileup ... --mod-threshold 0.1

# Analyze only m5C
modkit pileup ... --only-mods m5C

# Analyze all modifications
modkit pileup ... --only-mods m6A,m5C
```

## Result Interpretation

### Modification Frequency Distribution

```bash
# Calculate modification frequency
awk -F',' 'NR>1 {print $9}' m6A_pileup.csv | \
    sort | uniq -c | sort -k2 -n

# High confidence sites
awk -F',' 'NR>1 && $9>0.5' m6A_pileup.csv | wc -l
```

### Coverage Analysis

```bash
# Average coverage
awk -F',' 'NR>1 {sum+=$11} END {print sum/NR}' m6A_pileup.csv

# High coverage sites (>100x)
awk -F',' 'NR>1 && $11>100' m6A_pileup.csv | wc -l
```

### Motif Analysis

```bash
# Extract motif distribution
awk -F',' 'NR>1 {print $12}' m6A_pileup.csv | \
    sort | uniq -c | sort -rn | head -10

# Typical m6A motif: DRACH (D=A/G/T, R=A/G, A=A, C=C, H=A/C/T)
```

## Visualization

### Manhattan Plot

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('m6A_pileup.csv')

# Simplified Manhattan plot
plt.figure(figsize=(12, 4))
plt.scatter(df['start'], df['freq'], c='blue', alpha=0.5, s=1)
plt.xlabel('Position')
plt.ylabel('Modification Frequency')
plt.title('m6A Distribution')
plt.savefig('manhattan.png', dpi=300)
```

### Coverage Distribution

```python
plt.figure(figsize=(12, 4))
plt.hist(df['coverage'], bins=50, alpha=0.7)
plt.xlabel('Coverage')
plt.ylabel('Count')
plt.title('Coverage Distribution')
plt.savefig('coverage.png', dpi=300)
```

## Common Issues

### Issue 1: No Modification Sites Output

**Check:**
1. Confirm MM tags in BAM
2. Check coverage sufficiency
3. Confirm threshold settings

```bash
# Check MM tags
samtools view aligned.bam | grep -c "MM:Z:"

# Lower threshold
modkit pileup ... --mod-threshold 0.05 --min-coverage 5
```

### Issue 2: Modification Frequency Too Low

**Cause:**
- Insufficient coverage
- Few modifications in sample
- Basecalling model issue

**Solution:**
```bash
# Increase coverage requirement
modkit pileup ... --mod-threshold 0.05

# Check basecalling
samtools view aligned.bam | head | grep "MM:Z:"
```

### Issue 3: Runtime Too Long

**Solution:**
```bash
# Increase thread count
modkit pileup ... --threads 32

# Reduce chunks
modkit pileup ... --chunks 4
```

### Issue 4: Insufficient Memory

**Solution:**
```bash
# Reduce thread count
modkit pileup ... --threads 8

# Increase chunks (reduce memory)
modkit pileup ... --chunks 16
```

## Parameter Optimization

### High Precision Mode

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.5 \
    --min-coverage 30 \
    --threads 32
```

### Sensitivity Mode

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.05 \
    --min-coverage 10 \
    --threads 16
```

### Balanced Mode

```bash
modkit pileup input.bam output \
    --ref genome.fa \
    --mod-threshold 0.1 \
    --min-coverage 20 \
    --threads 16
```

## Verification

### Compare with Other Methods

- **MeRIP-seq**: Traditional immunoprecipitation method
- **m6A-seq**: m6A-specific sequencing
- **SCARLET**: Single-molecule validation gold standard

### Statistical Testing

```python
# Calculate confidence interval of modification rate
import scipy.stats as stats

coverage = df['coverage'].values
freq = df['freq'].values

# Binomial distribution confidence interval
for i in range(len(coverage)):
    n = coverage[i]
    k = int(n * freq[i])
    ci = stats.binom.interval(0.95, n, freq[i])
    print(f"Site {i}: {freq[i]:.3f} (95% CI: {ci[0]/n:.3f}-{ci[1]/n:.3f})")
```

## Citation

- **modkit**: Oxford Nanopore Technologies. https://github.com/nanoporetech/modkit
- **Original Paper**: Liu, H. et al. (2024). Accurate detection of m6A RNA modifications in native nanopore sequencing. Nature Methods.
