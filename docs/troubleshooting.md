# Troubleshooting Guide

This document lists common issues and solutions when using the Nanopore m6A detection pipeline.

## Table of Contents

1. [Basecalling Issues](#basecalling-issues)
2. [Alignment Issues](#alignment-issues)
3. [Modification Detection Issues](#modification-detection-issues)
4. [Performance Issues](#performance-issues)

---

## Basecalling Issues

### Issue 1: Dorado Installation Failed

**Symptoms:**
```
Command 'dorado' not found
```

**Solutions:**
```bash
# Method 1: Download precompiled version
wget https://github.com/nanoporetech/dorado/releases/download/v0.8.1/dorado-0.8.1-linux-x64.tar.gz
tar -xzf dorado-0.8.1-linux-x64.tar.gz
export PATH=$PATH:/path/to/dorado-0.8.1-linux-x64/bin

# Method 2: Use Docker
docker pull nanotechpore/dorado:0.8.1

# Method 3: Use Singularity
singularity pull docker://nanotechpore/dorado:0.8.1
```

### Issue 2: GPU Out of Memory

**Symptoms:**
```
CUDA out of memory
```

**Solutions:**
```bash
# Reduce batch-size
dorado basecaller --batch-size 500 ...

# Use CPU mode (slower but works)
dorado basecaller --device cpu ...

# Or use Guppy CPU version
guppy_basecaller -c dna_r10.4.1_e8.2_400bps_hac.cfg -i reads.pod5 -o output
```

### Issue 3: Modification Detection Failed

**Symptoms:**
```
No modification detected
```

**Check:**
1. Confirm using modification-supported model
2. Confirm using --emit-moves parameter
3. Check Dorado version

```bash
# Check version
dorado --version

# Use modification model
dorado basecaller \
    --emit-moves \
    --mod-base-models dna_r10.4.1_e8.2_400bps_hac \
    dna_r10.4.1_e8.2_400bps_hac \
    reads.pod5 > basecalls.bam
```

### Issue 4: Basecalling Too Slow

**Symptoms:**
```
Running time too long
```

**Solutions:**
```bash
# Use fast model
dorado basecaller dna_r10.4.1_e8.2_400bps_fast ...

# Increase batch-size
dorado basecaller --batch-size 2000 ...

# Process multiple pod5 files in parallel
ls *.pod5 | xargs -I{} dorado basecaller $MODEL {} > {}.fastq
```

---

## Alignment Issues

### Issue 5: minimap2 Installation Failed

**Symptoms:**
```
Command 'minimap2' not found
```

**Solutions:**
```bash
# conda installation
conda install -c bioconda minimap2

# Or compile from source
git clone https://github.com/lh3/minimap2.git
cd minimap2 && make -j8
```

### Issue 6: Low Alignment Rate

**Symptoms:**
```
Alignment rate below 50%
```

**Possible causes:**
1. Poor sequencing quality
2. Reference genome mismatch
3. Sequences too short

**Solutions:**
```bash
# Lower filtering criteria
NanoFilt -l 30 -q 5 input.fq > clean.fq

# Use more lenient alignment parameters
minimap2 -ax map-ont -k 12 -w 5 -r 0.5 ...

# Check reference genome
minimap2 -d genome.fa.mmi genome.fa  # Rebuild index
```

### Issue 7: Modification Tags Lost

**Symptoms:**
```
No MM tags in BAM file
```

**Check:**
```bash
# Check for MM tags
samtools view aligned.bam | grep "MM:Z:" | head

# Check if alignment succeeded
samtools flagstat aligned.bam
```

**Solutions:**
1. Check if Dorado used --emit-moves
2. Check if minimap2 used --MD parameter
3. Confirm no improper format conversion was performed

### Issue 8: Out of Memory

**Symptoms:**
```
Killed (out of memory)
```

**Solutions:**
```bash
# Reduce thread count
minimap2 -t 4 ...

# Reduce batch-size
minimap2 -K 100M ...

# Process in batches
split -l 5000 reads.fastq reads_chunk_
for f in reads_chunk_*; do
    minimap2 ... $f >> aligned.sam
done
```

### Issue 9: pychopper Not Found

**Symptoms:**
```
Command 'pychopper' not found
```

**Solutions:**
```bash
# conda installation
conda install -c bioconda pychopper

# Or pip installation
pip install pychopper
```

### Issue 10: NanoFilt Error

**Symptoms:**
```
NanoFilt: error: too few arguments
```

**Solutions:**
```bash
# Check parameters
NanoFilt -l 50 -q 7 input.fq > output.fq

# Parameter description:
# -l: minimum length
# -q: minimum quality score
# Input file at the end
```

---

## Modification Detection Issues

### Issue 11: modkit Installation Failed

**Symptoms:**
```
Command 'modkit' not found
```

**Solutions:**
```bash
# pip installation
pip install modkit

# conda installation
conda install -c bioconda modkit

# Verify
modkit --version
```

### Issue 12: No Modification Sites Output

**Symptoms:**
```
Output directory is empty
```

**Check:**
```bash
# 1. Confirm MM tags in BAM
samtools view aligned.bam | grep -c "MM:Z:"

# 2. Check modification types
samtools view aligned.bam | grep -o "MM:Z:[^ ]*" | head

# 3. Lower threshold and re-analyze
modkit pileup aligned.bam output \
    --ref genome.fa \
    --mod-threshold 0.05 \
    --min-coverage 5
```

### Issue 13: Modification Frequency Too Low

**Symptoms:**
```
Modification frequency generally below 10%
```

**Possible causes:**
1. Sample itself has few m6A
2. Coverage insufficient
3. Basecalling model issue

**Solutions:**
```bash
# Lower threshold
modkit pileup ... --mod-threshold 0.05

# Increase sequencing depth
# Or check raw data quality
```

### Issue 14: Doesn't Match Known Sites

**Symptoms:**
```
Detected sites differ significantly from literature reports
```

**Solutions:**
1. Confirm reference genome version is consistent
2. Use more lenient thresholds
3. Increase coverage
4. Cross-validate with other methods

### Issue 15: Runtime Too Long

**Symptoms:**
```
modkit runs for several hours
```

**Solutions:**
```bash
# Increase thread count
modkit pileup ... --threads 32

# Reduce chunks
modkit pileup ... --chunks 4
```

---

## Performance Issues

### Issue 16: Insufficient Disk Space

**Symptoms:**
```
No space left on device
```

**Solutions:**
```bash
# Clean temporary files
rm -f *.tmp *.sam

# Remove intermediate files
rm -f *.unsorted.bam

# Use CRAM format
samtools view -C aligned.bam -o aligned.cram

# Clean conda cache
conda clean -a
```

### Issue 17: Large File Processing Crash

**Symptoms:**
```
Program not responding or crashes
```

**Solutions:**
```bash
# Batch processing
# Split large pod5 files into smaller files

# Use GNU parallel
ls *.pod5 | parallel -j 2 'dorado basecaller {} > {.}.fastq'

# Batch alignment
split -l 10000 large.fastq chunk_
for f in chunk_*; do minimap2 ... $f; done
```

---

## Debugging Tips

### Enable Debug Mode

```bash
# View detailed output
dorado basecaller -v ...

# bash debug mode
bash -x run_pipeline.sh
```

### Check Intermediate Results

```bash
# Check FASTQ
head -20 reads.fastq
wc -l reads.fastq

# Check SAM/BAM
samtools view aligned.bam | head
samtools flagstat aligned.bam

# Check modification tags
samtools view aligned.bam | awk '{print $NF}' | grep "MM" | head
```

### Test with Small Dataset

```bash
# Extract small number of reads for testing
head -400 input.fastq > test.fastq

# Run full pipeline
bash run_pipeline.sh test.fastq

# Verify before processing full data
```

---

## Getting Help

If the above solutions don't resolve your issue:

1. Check all software versions
2. Check detailed error messages in log files
3. Search GitHub issues: https://github.com/IceBear321/Nanopore_m6A_Detection/issues
4. When opening an issue, include:
   - Complete error message
   - Command used
   - System environment (`uname -a`)
   - Software versions (`conda list` or `dorado --version`)
