# 故障排除指南

本文档列出使用Nanopore m6A检测流程时可能遇到的常见问题及其解决方案。

## 目录

1. [碱基识别问题](#碱基识别问题)
2. [比对问题](#比对问题)
3. [修饰检测问题](#修饰检测问题)
4. [性能问题](#性能问题)

---

## 碱基识别问题

### 问题1: Dorado安装失败

**症状:**
```
Command 'dorado' not found
```

**解决方案:**
```bash
# 方法1: 下载预编译版本
wget https://github.com/nanoporetech/dorado/releases/download/v0.8.1/dorado-0.8.1-linux-x64.tar.gz
tar -xzf dorado-0.8.1-linux-x64.tar.gz
export PATH=$PATH:/path/to/dorado-0.8.1-linux-x64/bin

# 方法2: 使用Docker
docker pull nanotechpore/dorado:0.8.1

# 方法3: 使用Singularity
singularity pull docker://nanotechpore/dorado:0.8.1
```

### 问题2: GPU内存不足

**症状:**
```
CUDA out of memory
```

**解决方案:**
```bash
# 减少batch-size
dorado basecaller --batch-size 500 ...

# 使用CPU模式 (慢但可用)
dorado basecaller --device cpu ...

# 或使用Guppy CPU版本
guppy_basecaller -c dna_r10.4.1_e8.2_400bps_hac.cfg -i reads.pod5 -o output
```

### 问题3: 修饰检测失败

**症状:**
```
No modification detected
```

**检查项:**
1. 确认使用支持modification的模型
2. 确认使用--emit-moves参数
3. 检查Dorado版本

```bash
# 检查版本
dorado --version

# 使用modification模型
dorado basecaller \
    --emit-moves \
    --mod-base-models dna_r10.4.1_e8.2_400bps_hac \
    dna_r10.4.1_e8.2_400bps_hac \
    reads.pod5 > basecalls.bam
```

### 问题4: 碱基识别速度慢

**症状:**
```
运行时间过长
```

**解决方案:**
```bash
# 使用fast模型
dorado basecaller dna_r10.4.1_e8.2_400bps_fast ...

# 增加batch-size
dorado basecaller --batch-size 2000 ...

# 并行处理多个pod5文件
ls *.pod5 | xargs -I{} dorado basecaller $MODEL {} > {}.fastq
```

---

## 比对问题

### 问题5: minimap2安装失败

**症状:**
```
Command 'minimap2' not found
```

**解决方案:**
```bash
# conda安装
conda install -c bioconda minimap2

# 或编译安装
git clone https://github.com/lh3/minimap2.git
cd minimap2 && make -j8
```

### 问题6: 比对率低

**症状:**
```
比对率低于50%
```

**可能原因:**
1. 测序质量差
2. 参考基因组不匹配
3. 序列太短

**解决方案:**
```bash
# 降低过滤标准
NanoFilt -l 30 -q 5 input.fq > clean.fq

# 使用更宽松的比对参数
minimap2 -ax map-ont -k 12 -w 5 -r 0.5 ...

# 检查参考基因组
minimap2 -d genome.fa.mmi genome.fa  # 重建索引
```

### 问题7: 修饰标签丢失

**症状:**
```
BAM文件中无MM标签
```

**检查:**
```bash
# 查看是否有MM标签
samtools view aligned.bam | grep "MM:Z:" | head

# 检查比对是否成功
samtools flagstat aligned.bam
```

:**
1. 检查Dorado是否**解决方案使用--emit-moves
2. 检查minimap2是否使用--MD参数
3. 确认没有进行不当的格式转换

### 问题8: 内存不足

**症状:**
```
Killed (out of memory)
```

**解决方案:**
```bash
# 减少线程数
minimap2 -t 4 ...

# 减小batch-size
minimap2 -K 100M ...

# 分批处理
split -l 5000 reads.fastq reads_chunk_
for f in reads_chunk_*; do
    minimap2 ... $f >> aligned.sam
done
```

### 问题9: pychopper找不到

**症状:**
```
Command 'pychopper' not found
```

**解决方案:**
```bash
# conda安装
conda install -c bioconda pychopper

# 或pip安装
pip install pychopper
```

### 问题10: NanoFilt报错

**症状:**
```
NanoFilt: error: too few arguments
```

**解决方案:**
```bash
# 检查参数
NanoFilt -l 50 -q 7 input.fq > output.fq

# 参数说明:
# -l: 最小长度
# -q: 最小质量分数
# 输入文件在最后
```

---

## 修饰检测问题

### 问题11: modkit安装失败

**症状:**
```
Command 'modkit' not found
```

**解决方案:**
```bash
# pip安装
pip install modkit

# conda安装
conda install -c bioconda modkit

# 验证
modkit --version
```

### 问题12: 无修饰位点输出

**症状:**
```
Output directory is empty
```

**检查项:**
```bash
# 1. 确认BAM中有MM标签
samtools view aligned.bam | grep -c "MM:Z:"

# 2. 检查修饰类型
samtools view aligned.bam | grep -o "MM:Z:[^ ]*" | head

# 3. 降低阈值重新分析
modkit pileup aligned.bam output \
    --ref genome.fa \
    --mod-threshold 0.05 \
    --min-coverage 5
```

### 问题13: 修饰频率过低

**症状:**
```
修饰频率普遍低于10%
```

**可能原因:**
1. 样本本身m6A较少
2. 覆盖度不够
3. 碱基识别模型问题

**解决方案:**
```bash
# 降低阈值
modkit pileup ... --mod-threshold 0.05

# 增加测序深度
# 或检查原始数据质量
```

### 问题14: 与已知位点不匹配

**症状:**
```
检测到的位点与文献报道差异大
```

**解决方案:**
1. 确认参考基因组版本一致
2. 使用更宽松的阈值
3. 增加覆盖度
4. 与其他方法结果交叉验证

### 问题15: 运行时间过长

**症状:**
```
modkit运行超过数小时
```

**解决方案:**
```bash
# 增加线程数
modkit pileup ... --threads 32

# 减少chunks
modkit pileup ... --chunks 4
```

---

## 性能问题

### 问题16: 磁盘空间不足

**症状:**
```
No space left on device
```

**解决方案:**
```bash
# 清理临时文件
rm -f *.tmp *.sam

# 删除中间文件
rm -f *.unsorted.bam

# 使用CRAM格式
samtools view -C aligned.bam -o aligned.cram

# 清理conda缓存
conda clean -a
```

### 问题17: 处理大文件崩溃

**症状:**
```
程序无响应或崩溃
```

**解决方案:**
```bash
# 分批处理
# 将大pod5文件分成多个小文件

# 使用GNU parallel
ls *.pod5 | parallel -j 2 'dorado basecaller {} > {.}.fastq'

# 分批比对
split -l 10000 large.fastq chunk_
for f in chunk_*; do minimap2 ... $f; done
```

---

## 调试技巧

### 启用调试模式

```bash
# 查看详细输出
dorado basecaller -v ...

# bash调试模式
bash -x run_pipeline.sh
```

### 检查中间结果

```bash
# 检查FASTQ
head -20 reads.fastq
wc -l reads.fastq

# 检查SAM/BAM
samtools view aligned.bam | head
samtools flagstat aligned.bam

# 检查修饰标签
samtools view aligned.bam | awk '{print $NF}' | grep "MM" | head
```

### 测试小数据集

```bash
# 提取少量reads测试
head -400 input.fastq > test.fastq

# 运行完整流程
bash run_pipeline.sh test.fastq

# 验证后再处理完整数据
```

---

## 获取帮助

如果以上方案都不能解决问题:

1. 检查所有软件的版本
2. 查看日志文件中的详细错误信息
3. 搜索GitHub issues: https://github.com/IceBear321/Nanopore_m6A_Detection/issues
4. 提issue时附上:
   - 完整的错误信息
   - 使用的命令
   - 系统环境 (`uname -a`)
   - 软件版本 (`conda list` 或 `dorado --version`)
