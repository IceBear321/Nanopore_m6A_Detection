# Step 1: Basecalling

## Purpose

Convert raw Nanopore sequencing signal data (.pod5) to nucleotide sequences (.fastq). Modern basecallers can directly identify and label modified bases (m6A, etc.).

## Tool Selection

### Dorado (Recommended)

ONT-official basecalling tool with GPU acceleration support.

**Installation:**
```bash
# Download precompiled version
wget https://github.com/nanoporetech/dorado/releases/download/v0.8.1/dorado-0.8.1-linux-x64.tar.gz
tar -xzf dorado-0.8.1-linux-x64.tar.gz
export PATH=$PATH:/path/to/dorado-0.8.1-linux-x64/bin

# Verify installation
dorado --version
```

### Guppy

ONT official tool, both GPU and CPU versions available.

```bash
# conda installation
conda install -c conda-forge ont-guppy

# Or use docker
docker pull ontporetech/guppy:7.0.0
```

## Model Selection

### R10.4.1 (Latest)

| Model | Accuracy | Speed | Modification Detection |
|-------|----------|--------|----------------------|
| dna_r10.4.1_e8.2_400bps_hac | High | Medium | Supported |
| dna_r10.4.1_e8.2_400bps_sup | Highest | Slow | Supported |
| dna_r10.4.1_e8.2_400bps_fast | Medium | Fast | Supported |

### R9.4.1 (Classic)

| Model | Accuracy | Speed |
|-------|----------|--------|
| dna_r9.4.1_e8_hac | High | Medium |
| dna_r9.4.1_e8_sup | Highest | Slow |
| dna_r9.4.1_e8_fast | Medium | Fast |

## Command Details

### Basic Command

```bash
dorado basecaller \
    --emit-moves \
    --mod-base-models dna_r10.4.1_e8.2_400bps_hac \
    dna_r10.4.1_e8.2_400bps_hac \
    reads.pod5 > basecalls.bam
```

### Key Parameters

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `--emit-moves` | Preserve modification information | Must include |
| `--mod-base-models` | Enable modification detection model | dna_r10.4.1_e8.2_400bps_hac |
| `--threads` | Thread count | 16-32 |
| `--batch-size` | Batch size | 1000 |
| `--min-qscore` | Minimum quality score | 10 |

### Output Format Conversion

```bash
# BAM to FASTQ
samtools fastq basecalls.bam > reads.fastq

# Direct FASTQ output (requires newer version)
dorado basecaller ... --fastq > reads.fastq
```

## Input File Format

### .pod5 File

**File Structure:**
- Binary format
- Contains raw current signals
- Each read has corresponding timestamp and channel information

**View Information:**
```bash
pod5 view --stats reads.pod5
pod5 count reads.pod5
```

**File Size:**
- ~2-5KB per read
- ~2-5GB for 1 million reads

## Modification Detection Principle

### Mm/MI Tags

Dorado generates modification information during basecalling:

```
MM:Z:m+6A    # m6A modification
ML:B:C,255   # Confidence score
```

**Tag Description:**
- `MM`: Modification type and position
- `ML`: Modification confidence (0-255)

## Common Issues

### Issue 1: GPU Out of Memory

**Symptoms:**
```
CUDA out of memory
```

**Solutions:**
```bash
# Reduce batch-size
dorado basecaller --batch-size 500 ...

# Or use CPU mode
dorado basecaller --device cpu ...
```

### Issue 2: Too Slow

**Solutions:**
```bash
# Use fast model
dorado basecaller dna_r10.4.1_e8.2_400bps_fast ...

# Increase batch-size
dorado basecaller --batch-size 2000 ...
```

### Issue 3: Modification Detection Failed

**Check Items:**
1. Confirm using --emit-moves parameter
2. Confirm using modification-supported model
3. Check Dorado version (must be >=0.7.0)

```bash
# Check version
dorado --version

# Verify modification tags
samtools view basecalls.bam | head | grep "MM:Z:"
```

### Issue 4: pod5 File Corrupted

**Detection:**
```bash
pod5 view --stats reads.pod5

# Repair if necessary
pod5 validate reads.pod5
```

## Performance Benchmarks

| Sequencer | GPU | Model | Speed (reads/s) |
|-----------|-----|-------|-----------------|
| PromethION | A100 | hac | ~3000 |
| PromethION | A100 | sup | ~1500 |
| MinION | RTX3090 | hac | ~450 |
| MinION | RTX3090 | fast | ~800 |

## Output Verification

```bash
# Count reads
samtools view -c basecalls.bam

# View quality distribution
samtools fastq basecalls.bam | \
    awk 'NR%4==0' | \
    tr -d '\n' | \
    fold -w 1 | \
    sort | \
    uniq -c | \
    sort -rn

# Check modification tags
samtools view basecalls.bam | \
    grep -o 'MM:Z:[^ ]*' | \
    sort | \
    uniq -c
```
