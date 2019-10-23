# DNA
dna 分析流程

# 使用方法
    nextflow run main.nf -profile {x86|arm} --reads1=sample_R1.fq.gz  --reads2=sample.R2.fq.gz --sample={SAMPLE_NAME} --output={OUTPUT DIR} 
    参数：
    profile|string          必填    指定计算平台（arm|x86)
    reads1|string           必填    测序下机数据中reads1文件
    reads2|string           必填    测序下机数据中reads2文件
    sample|string           必填    样本名称
    output|string           必填    输出结果所在的目录
    bwa_cpu|int             可选    bwa软件进程数（默认：2）
    gatk_cpu|int            可选    gatk软件进程数（默认：2）
    bwa_db_prefix|string    可选    bwa 用到的参考基因组索引文件（默认：/data/wgs_db/hg19/bwa/chr1-x）
    ref_sequence|string     可选    参考基因组文件（默认：/data/wgs_db/hg19/bwa/chr1-x.fasta）

# 依赖环境

    java:1.8
    nextflow:lastest
    docker:lastest
