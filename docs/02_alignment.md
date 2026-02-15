# Step 2: Sequence Alignment

## Purpose

Align sequence reads from basecalling to the reference genome to generate BAM files containing position information. This step needs to preserve MM/MI tag information for modified bases while performing quality control.

## Complete Workflow

```
Input FASTQ → pychopper → NanoFilt → minimap2 → SAM→BAM → Sort → Deduplicate → Index
```

## Tools

### Minimap2 (Recommended)

Fast alignment tool designed for long-read sequences.

**Installation:**
```bash
conda install -c bioconda minimap2
# Or compile from source
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```

## Command Details

### Complete Alignment Workflow

```bash
# 1. pychopper (remove low-quality sequences and adapters)
pychopper -t 8 input.fastq output.fq

# 2. NanoFilt (quality and length filtering)
NanoFilt -l 50 -q 7 output.fq > output.clear.fq

# 3. minimap2 alignment
minimap2 -ax map-ont -t 8 -uf -k 14 \
    reference.fa \
    output.clear.fq > output.sam

# 4. SAM to BAM
samtools view -@ 8 -Sb output.sam > output.unsorted.bam

# 5. Sort
samtools sort -@ 8 -o output.sorted.bam -T tmp output.unsorted.bam

# 6. Deduplicate
samtools markdup -@ 8 -r output.sorted.bam output.rmdup.bam

# 7. Index
samtools index output.rmdup.bam
```

### minimap2 Key Parameters

```bash
minimap2 -ax map-ont -t 8 -uf -k 14 \
    reference.fa \
    reads.fastq > output.sam
```

| Parameter | Description | Recommended | Explanation |
|-----------|-------------|------------|-------------|
| `-ax map-ont` | Alignment mode | Required | Optimized for ONT long reads |
| `-t` | Thread count | 8-16 | Adjust based on CPU cores |
| `-uf` | Forward strand only | Optional | For RNA-seq |
| `-k` | k-mer size | 14 | Use 14 for R10.4, 15 for R9.4 |
| `-K` | Batch size | 1G | Affects memory and speed |
| `--MD` | Preserve MD tags | Recommended | For modification detection |
| `-r` | Error rate | 0 | Use actual error rate |

#### Parameter Details

**-ax map-ont**
```
Available modes:
- map-ont:    ONT long reads (Recommended)
- map-hifi:   PacBio HiFi reads
- map-pb:     PacBio CLR reads
- asm20:      Genome assembly (low consistency)
- asm5:       Genome assembly (medium consistency)
- splice:    Splice-aware (RNA-seq)
```

**-k Parameter Selection**
```bash
# R10.4 sequencer (higher accuracy)
minimap2 -k 14 ...

# R9.4 sequencer
minimap2 -k 15 ...

# Shorter k-mer increases sensitivity but may add errors
minimap2 -k 12 ...
```

**-uf Parameter (Forward-only)**
```
-uf: Only align to forward strand of reference genome
Use cases:
- RNA-seq:  Only keep transcript alignments
- Avoid antisense strand interference
```

### pychopper

**Purpose:** Remove adapter sequences and low-quality regions

```bash
pychopper -t 8 input.fastq output.fq
```

| Parameter | Description |
|-----------|-------------|
| `-t` | Thread count |
| `-r` | Report file |
| `-S` | Statistics file |

### NanoFilt

**Purpose:** Quality filtering and length filtering

```bash
NanoFilt -l 50 -q 7 input.fq > output.fq
```

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `-l` | Minimum length | 50-100 |
| `-q` | Minimum quality score | 7-10 |
| `--head` | Keep first N reads | - |
| `--tail` | Keep last N reads | - |

### SAM Processing

```bash
# SAM to BAM
samtools view -@ 8 -Sb input.sam > input.bam

# Sort
samtools sort -@ 8 -o sorted.bam -T tmp input.bam

# Deduplicate
samtools markdup -@ 8 -r sorted.bam rmdup.bam

# Index
samtools index rmdup.bam
```

## Index Preparation

### Create Reference Genome Index

```bash
# minimap2 index (recommended to pre-build)
minimap2 -d genome.fa.mmi genome.fa

# samtools index
samtools faidx genome.fa
```

### Check Index

```bash
# Check index files
ls -lh genome.fa.mmi genome.fa.fai

# Verify index is usable
samtools faidx genome.fa chr1:1-1000
```

## Output Format

### BAM File Structure

| Field | Description | Example |
|-------|-------------|---------|
| QNAME | Read name | read_001 |
| FLAG | Alignment status | 0 (forward) / 16 (reverse) |
| RNAME | Chromosome | chr1 |
| POS | Start position | 1234567 |
| MAPQ | Alignment quality | 60 |
| CIGAR | Alignment details | 150M |
| SEQ | Sequence | ATCG... |
| QUAL | Quality score | IIII... |
| TAGS | Modification tags | MM:Z:m+6A |

### FLAG Description

| FLAG | Meaning |
|------|---------|
| 0 | Forward alignment |
| 16 | Reverse alignment |
| 4 | Not aligned |
| 2048 | PCR duplicate |

### Modification Tags

```
MM:Z:m+6A      # m6A modification type
ML:B:C,200     # Confidence (0-255)
MP:A:+,10,15   # Modification position
```

## Statistics

### Alignment Rate

```bash
# Total reads
samtools view -c aligned.bam

# Aligned reads
samtools view -c -F 4 aligned.bam

# Unaligned reads
samtools view -c -f 4 aligned.bam

# Calculate alignment rate
echo "scale=2; $(samtools view -c -F 4 aligned.bam) / $(samtools view -c aligned.bam) * 100" | bc
```

### Coverage Statistics

```bash
# Average coverage
samtools depth aligned.bam | awk '{sum+=$3} END {print sum/NR}'

# Coverage distribution
samtools depth aligned.bam | awk '{print $3}' | sort -n | uniq -c

# Using bamCoverage
bamCoverage -b aligned.bam -of bigwig -o coverage.bw
```

### Read Length Distribution

```bash
# Average read length
samtools view aligned.bam | awk '{print length($10)}' | awk '{sum+=$1} END {print sum/NR}'

# Read length histogram
samtools view aligned.bam | awk '{print length($10)}' | \
    sort | uniq -c | awk '{print $2"\t"$1}' | head -20
```

### Deduplication Statistics

```bash
# Using flagstat
samtools flagstat aligned.rmdup.bam

# Output example:
# 1000000 + 0 in total (QC-passed reads + QC-failed reads)
# 5000 + 0 secondary
# 0 + 0 supplementary
# 950000 + 0 mapped (95.00% : N/A)
# 10000 + 0 duplicates
```

## Common Issues

### Issue 1: Low Alignment Rate

**Possible causes:**
1. Poor sequencing quality
2. Reference genome mismatch
3. Sequences too short

**Solutions:**
```bash
# Use more lenient parameters
minimap2 -ax map-ont -k 12 -w 5 ...

# Lower minimum length filtering
NanoFilt -l 30 -q 5 ...
```

### Issue 2: Out of Memory

**Solutions:**
```bash
# Reduce thread count
minimap2 -t 4 ...

# Reduce batch-size
minimap2 -K 100M ...

# Process in batches
split -l 10000 reads.fastq reads_
for f in reads_*; do
    minimap2 ... $f >> output.sam
done
```

### Issue 3: Modification Tags Lost

**Check:**
```bash
# Check MM tags
samtools view aligned.bam | grep "MM:Z:" | head

# Check MD tags
samtools view aligned.bam | grep "MD:Z:" | head
```

**Solutions:**
- minimap2 version too old → Upgrade to latest version
- Didn't use --MD parameter → Add --MD parameter
- Format conversion lost tags → Check conversion steps

### Issue 4: Too Much Coverage Reduction After Deduplication

**Note:** ONT data has inherently low PCR duplication rate, excessive deduplication may indicate:
- Same read aligned multiple times
- Library PCR amplification

**Adjustment:**
```bash
# Only mark, don't delete
samtools markdup -@ 8 aligned.bam output.bam

# Or skip deduplication step
```

### Issue 5: CIGAR Errors

**Note:** minimap2 uses simplified CIGAR, some tools may be incompatible

**Verify:**
```bash
# Check CIGAR format
samtools view aligned.bam | head | awk '{print $6}'
```

## Complete Parameter Templates

### for RNA-seq (m6A Detection)

```bash
#!/bin/bash
# ONT RNA-seq m6A detection alignment workflow

REF="genome.fa"
THREADS=16

# 1. Preprocessing
pychopper -t "$THREADS" input.fastq tmp.fq
NanoFilt -l 50 -q 7 tmp.fq > clean.fq

# 2. Alignment
minimap2 -ax map-ont -t "$THREADS" -uf -k 14 \
    "$REF" \
    clean.fq > aligned.sam

# 3. BAM processing
samtools view -@ "$THREADS" -Sb aligned.sam | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" -r - \
    aligned.rmdup.bam

# 4. Index
samtools index aligned.rmdup.bam

# Cleanup
rm -f tmp.fq clean.fq aligned.sam
```

### for DNA-seq (Direct Sequencing)

```bash
#!/bin/bash
# ONT DNA-seq alignment workflow

REF="genome.fa"
THREADS=16

# 1. Filter
NanoFilt -l 100 -q 10 input.fastq > clean.fq

# 2. Alignment (don't use -uf)
minimap2 -ax map-ont -t "$THREADS" -k 14 \
    "$REF" \
    clean.fq > aligned.sam

# 3. BAM processing
samtools view -@ "$THREADS" -Sb aligned.sam | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" -r - \
    aligned.rmdup.bam

samtools index aligned.rmdup.bam
```

## Verification Checklist

```bash
# 1. Check BAM file
samtools flagstat aligned.rmdup.bam

# 2. Check modification tags
echo "MM tags:"
samtools view aligned.rmdup.bam | grep -c "MM:Z:"

# 3. Check coverage
samtools depth aligned.rmdup.bam | awk '{sum+=$3} END {print "Avg coverage:", sum/NR}'

# 4. Check bigWig
ls -lh processed/*.bw

# 5. Check reads distribution
samtools view -c aligned.rmdup.bam
```

## Comparison with Other Tools

| Tool | Advantages | Disadvantages |
|------|------------|---------------|
| Minimap2 | Fast, optimized for ONT | RNA alignment needs extra parameters |
| BWA-MEM | Better compatibility | Not optimized for long reads |
| GraphMap | High precision | Slow |
| NGMLR | Good for variant detection | Slow |
