# wgs全基因组序列比对流程
### 用到的软件
>  [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/"%3Ehttps://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/%3C/a%3E)
>  [bgzip](
http://www.htslib.org/doc/bgzip.html)
>  [tabix](
http://www.htslib.org/doc/tabix.html)
>  [bwa](
http://bio-bwa.sourceforge.net/)
>  [samtools](
http://samtools.sourceforge.net/)
>  [GATK4.0](
https://software.broadinstitute.org/gatk/)
 [gatk 官方流程](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
)
***
### 过程步骤
##### 一. 下载准备需要的文件
1. 下载参考序列基因组文件
    > 1.建立索引
     `bwa  index ref.fasta`
     >> 完成之后 会看到几个ref.fasta为前缀的文件
    > 2. 为参考序列生成dict文件
     ` gatk CreateSequenceDictionary -R ref.fasta -O ref.dict` 
    > 3. samtools 建索引
    `samtools faidx ref.fasta`
2. 下载测序文件
   ` fastaq-dump --split-files SRR*****`
>下载的文件是双末端测序从两端读的read1和read2
     >> 用bgzip压缩
        ` bgzip seq1_.fasta`
       ` bgzip seq2_.fasta`
****
##### 二.处理文件
###### 将read比对到参考基因组
 `bwa mem -t 4 -R '@RG\tID:foo\tPL:illumina\tSM:seq' ref.fasta seq_1.fastaq.gz seq_2.fastaq.gz | samtools view -Sb -> seq.bam
`
1. 首先利用bwa mem比对模块将E.coli K12质控后的测序数据定位到其参考基因组上（我们这里设置了4个线程来完成比对，根据电脑性能可以适当调大），同时通过管道（'|' 操作符）将比对数据流引到samtools转换为BAM格式（SAM的二进制压缩格式），然后重定向('>'操作符)输出到文件中保存下来。

2. -R 设置Read Group信息，它是read数据的组别标识，并且其中的ID，PL和SM信息在正式的项目中是不能缺少的(如果样本包含多个测序文库的话，LB信息也不要省略)，另外由于考虑到与GATK的兼容关系，PL（测序平台）信息不能随意指定，必须是：ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，PACBIO，IONTORRENT，CAPILLARY，HELICOS或UNKNOWN这12个中的一个

 ###### 排序
``samtools sort -@ 4 -O bam -o seq.sorted.bam  seq.bam ``
___
##### 标记PCR重复
```
gatk MarkDuplicates -I seqsorted.bam \
 -O seq.sorted.markdup.bam \
-M E_seq.sorted.markdup_metrics.txt
```

> 在制备文库的过程中，由于PCR扩增过程中会存在一些偏差，也就是说有的序列会被过量扩增。这样，在比对的时候，这些过量扩增出来的完全相同的序列就会比对到基因组的相同位置。而这些过量扩增的reads并不是基因组自身固有序列，不能作为变异检测的证据，因此，要尽量去除这些由PCR扩增所形成的duplicates。去重复的过程是给这些序列设置一个flag以标志它们，方便GATK的识别。还可以设置 REMOVE_DUPLICATES=true 来丢弃duplicated序列。对于是否选择标记或者删除，对结果应该没有什么影响，GATK官方流程里面给出的例子是仅做标记不删除。这里定义的重复序列是这样的：如果两条reads具有相同的长度而且比对到了基因组的同一位置，那么就认为这样的reads是由PCR扩增而来，就会被GATK标记

###### 创建比对索引文件 
`samtools index seq.sorted.markdup.bam`
___
#### 三. 变异检测
如果是多样本，为每一个样本都生成一个中间文件GVCF
```
gatk HaplotypeCaller \
-R ref.fa \
--emit-ref-confidence GVCF \
-I  seq.sorted.markdup.bam \
-O  seq.g.vcf
```
> 这里非常费时间,每个样本需要很长时间 而且gatk的HaplotypeCaller的spark功能目前还没有开发完成，可以python脚本多进程运行并行计算，如果是一个样本高深度测序，也可以指定 -L 参数分每条染色体并行跑会快很多
> 
```
#此处为写的python脚本代码 开多个进程
from multiprocessing import Pool
import subprocess

s = ['C001','C002','C003']
def run(num):
    cmd = """
             gatk HaplotypeCaller \
            -R ref.fa \
            --emit-ref-confidence GVCF \
            -I  seq_{}.sorted.markdup.bam \
            -O  seq{}.g.vcf.gz """.format(num,num)
    subprocess.getoutput(cmd)
if __name__ == '__main__':
    p = Pool(5)
    for i in s:
        p.apply_async(run,(i,))
    p.close()
    p.join()
```
> 这里的每一个样本对比排序完成之后的文件名是seq_C00.sorted.markdup.bam

####  合并GVCF文件
(在大量产生的gvcf文件有的文件索引可能会出问题报错，可删除索引文件重新生成索引
> `bgzip sample.g.vcf`
> `tabix -p vcf sample.g.vcf.gz`)
也有可能是系统内存不够，这里就需要升级配置


这里有两种方式  一种是传统的combineGVCFs
```
gatk  CombineGVCFs \
 -V seq1.g.vcf \
-V seq2.g.vcf  \
-V seq3.g.vcf  \
 -O combine.g.vcf
```
第二种方式是为每条染色体建立一个数据库  ，速度快
```
gatk  GenomicsDBImport  \
    -V C001.g.vcf  \
    -V C002.g.vcf  \ 
    -V C003.g.vcf  \
    --genomicsdb-workspace-path database_chr1 \
    -L chr1
```



通过gvcf筛选变异
```
gatk GenotypeGVCFs \
-R ref.fa \
-V combine.g.vcf \(-V gendb://database_chr01)
-O seq.vcf
```


> 用数据库筛选变异，只能每一条染色体分开筛选，这样效率高
***
#### 四. 对筛选到的变异进行过滤分离
```
1.筛选SNP
gatk SelectVariants \
-select-type SNP \
-V seq.vcf \
--restrict-alleles-to BIALLELIC  \ #限制只有两种基因型
-O SNP.vcf

2.过滤SNP
gatk VariantFiltration \
-V     SNP.vcf
--filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "Filter" \
-O     SNP.filter.vcf
```
```
1.筛选INDEL
gatk SelectVariants \
-select-type INDEL\
--restrict-alleles-to BIALLELIC  \ #限制只有两种基因型
-V seq.vcf \
-O INDEL.vcf

2.过滤INDEL
gatk VariantFiltration \
-V     INDEL.vcf.gz \
--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "Filter" \
-O         INDEL.filter.vcf


3. 提取通过filter的结果
gatk SelectVariants  -V SNP.filter.gz --exclude-filtered true -O filter.vcf.gz
```



> QualByDepth(QD): 变异位点可信度除以未过滤的非参考read数 
> FisherStrand (FS): Fisher精确检验评估当前变异是strand bias的可能性，这个值在0-60间 
> RMSMappingQuality (MQ): 所有样本中比对质量的平方根
> MappingQualityRankSumTest (MQRankSum): 根据REF和ALT的read的比对质量来评估可信度
> ReadPosRankSumTest (ReadPosRankSum) : 通过变异在read的位置来评估变异可信度，通常在read的两端的错误率比较高
> StrandOddsRatio (SOR) : 综合评估strand bias的可能性




>> 对结果进行压缩
`bgzip SNP.filter.vcf`
`bgzip INDEL.filter.vcf`
____
#### 五. 用snpeff变异注释
[官方教程点这里](
http://snpeff.sourceforge.net/SnpEff_manual.html)
1. 下载注释文件（gff3）
2. 修改配置文件snp.Eff.config 
 >添加一行 Rice.genome : Rice
3. 在snpeff的安装目录下data目录创建Rice目录和genomes
> 把基因组序列文件放在genomes下面  Rice下存放注释文件（gff3）
> gff3文件改名为genes.gff
4. 执行建数据库命令
`java -jar snpEff.jar build -gff3 -v Rice`
5. 进行序列注释
`java -jar snpEff.jar Rice SNP.filter.vcf.gz > rice.ann.vcf`
#### 六. 通过分割任务来提高速度
> 在HaplotypeCaller和GenomicsDBImport这两步中，通过指定 -L 参数来分割每次运行的区间，从而多进程多节点运行提高速度

1. 策略一：按照染色体进行分割
2. 策略二：按照染色体的NNN区域来进行分割
> 这里用到了seqkit locate 定位染色体上NNN， 然后用bedtools complement 得到后续输入的bed文件
```
seqkit locate -j 20 -Pp '"N{10,}"' ref.fasta | tail -n+2 | awk '{print $1"\t"$5"\t"$6}' > ath_n_pos.bed
bioawk -c fastx '{print $name"\t"length($seq)}' ref.fasta > tair10.txt
bedtools complement -i ath_n_pos.bed -g tair10.txt  > tair10_interval.bed
```
最后产生了大量vcf文件通过gatk提供的MergeVcfs工具来合并结果
