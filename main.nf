#!/usr/bin/env nextflow



params.output='output'
params.tmp='/tmp'
params.database="/DG/database/genomes/Homo_sapiens/hg19"
params.target="/DG/programs/beta/scripts/DNA/genetest/adult/dgadulttarget.bed"
params.input="/DG/project0/DGgene/rawdata/20180621L3-SR16008/*_S74_L003_R{1,2}_001.fastq.gz"

Channel.fromFilePairs(params.input).into{samples0;samples1;samples2}
samples2.println()
log.info "DNA - NF ~ version 1.0"
log.info "=============================="
log.info "fastq path    :  ${params.input}"
log.info "database home :  ${params.database}"
log.info "target file   :  ${params.target}"

log.info "\n"

process soapnuke{
    tag { sample_name }
    container 'registry.cn-hangzhou.aliyuncs.com/bio/soapnuke'
    input:
        set sample_name,files from samples0
    
    output:
        set sample_name,file("*.fastq.gz") into clean_samples

    script:
    """
    SOAPnuke filter -1 ${files[0]} -2 ${files[1]} -l 15 -q 0.5 -Q 2 -G -o . \
        -C sample.clean1.fastq.gz -D sample.clean2.fastq.gz
    """
}

process fastqc{
    tag { sample_name }
    container 'biocontainers/fastqc:v0.11.5'
    //conda "fastqc"
    publishDir { "${params.output}/fastqc/"+ sample_name }
    input:
        set sample_name,files from samples1
    output:
        file "*_fastqc/Images/*.png"
    script:
        """
        echo \$PWD
        fastqc --extract -o . ${files[0]} ${files[1]}
        """
}

process bwa{
    container 'biocontainers/bwa:0.7.15'
    input:
        set sample_name,files from clean_samples
    
    output:
        set sample_name,file('sample.clean.sam') into sam_res

    script:
    """
    bwa mem -t 10 -M -T 30 ${params.database}/bwa/chr1-x \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:DG\tPU:illumina" \
        ${files[0]} ${files[1]} > sample.clean.sam
    """     
}

process sort{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //picard=1.1150
    input:
        set sample_name,file('sample.clean.sam') from sam_res
    
    output:
        set sample_name,file('sample.clean.bam') into bam_res

    script:
    """
    java -Xmx2g -Djava.io.tmpdir=${params.tmp} -jar /usr/gitc/picard.jar  SortSam \
        INPUT=sample.clean.sam \
        OUTPUT=sample.clean.bam \
        MAX_RECORDS_IN_RAM=500000 \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=true
    """
}

bam_res.into{bam_res0;bam_res1;bam_res2;bam_res3}

process align_stat{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //picard=1.1150
    input:
        set sample_name,file('sample.clean.bam') from bam_res0

    script:
        """
        java -Xmx2g -Djava.io.tmpdir=${params.tmp} -jar /usr/gitc/picard.jar CollectAlignmentSummaryMetrics \
            INPUT=sample.clean.bam \
            OUTPUT=sample.mapped.stat \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT \
            REFERENCE_SEQUENCE=${params.database}/bwa/chr1-x.fasta
        """
}

process depth{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.clean.bam') from bam_res1
    output:
        file('sample.target.basedepth.sample_interval_summary') into depth_res

    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T DepthOfCoverage \
            -R ${params.database}/bwa/chr1-x.fasta \
            -L ${params.target} \
            -I sample.clean.bam \
            -ct 1 -ct 10 -ct 20 -ct 30 -ct 50 -ct 100 -ct 200 -ct 1000 \
            -o sample.target.basedepth
        """
}

process relign{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.clean.bam') from bam_res2
    output:
        file('sample.realigner.dedupped.clean.intervals') into target_intervals
    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T RealignerTargetCreator \
            -R ${params.database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -o sample.realigner.dedupped.clean.intervals \
            -I sample.clean.bam \
            -known ${params.database}/GATK/1000G_phase1.indels.hg19.vcf \
            -known ${params.database}/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
        """
}

process IndelRealigner{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.clean.bam') from bam_res3
        file('sample.realigner.dedupped.clean.intervals') from target_intervals
    
    output:
        set sample_name,file('sample.realigned.clean.bam') into realigned_bam_res
    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T IndelRealigner \
            -filterNoBases \
            -R ${params.database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -I sample.clean.bam \
            -targetIntervals sample.realigner.dedupped.clean.intervals \
            -o sample.realigned.clean.bam \
            -known ${params.database}/GATK/1000G_phase1.indels.hg19.vcf \
            -known ${params.database}/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
        """
}


realigned_bam_res.into{realigned_bam_res0;realigned_bam_res1}

process BQSR{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.realigned.clean.bam') from realigned_bam_res0
    
    output:
        file('sample.recal.table') into bsqr_res

    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T BaseRecalibrator \
            -R ${params.database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -I sample.realigned.clean.bam \
            -knownSites ${params.database}/GATK/dbsnp_137.hg19.vcf \
            -knownSites ${params.database}/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf \
            -knownSites ${params.database}/GATK/1000G_phase1.indels.hg19.vcf \
            -o sample.recal.table
        """
}

process print_reads{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.realigned.clean.bam') from realigned_bam_res1
        file('sample.recal.table') from bsqr_res
    output:
        set sample_name,file('sample.recal.final.clean.bam') into recal_bam_res
    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T PrintReads \
            -R ${params.database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -I sample.realigned.clean.bam \
            -BQSR sample.recal.table \
            -o sample.recal.final.clean.bam
        """
}

recal_bam_res.into{recal_bam_res0;recal_bam_res1}

process UnifiedGenotyper_snp{

    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    publishDir { "${params.output}/snp/"}
    input:
        set sample_name,file('sample.recal.final.clean.bam') from recal_bam_res0
    output:
        set sample_name,file('*.snp.vcf') into snp_vcf_res
    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T UnifiedGenotyper \
            -R ${params.database}/bwa/chr1-x.fasta \
            -I sample.recal.final.clean.bam \
            -glm SNP \
            -D ${params.database}/GATK/dbsnp_137.hg19.vcf \
            -o ${sample_name}.snp.vcf \
            -stand_call_conf 30 \
            -baqGOP 30 \
            -L /DG/programs/beta/scripts/DNA/genetest/adult/dgadultsnp.bed \
            -nct 12 \
            -dcov  10000 \
            -U ALLOW_SEQ_DICT_INCOMPATIBILITY -A VariantType -A QualByDepth \
            -A HaplotypeScore -A BaseQualityRankSumTest \
            -A MappingQualityRankSumTest -A ReadPosRankSumTest \
            -A FisherStrand -A DepthPerAlleleBySample \
            -A ClippingRankSumTest \
            --output_mode EMIT_ALL_SITES
        """
}

process UnifiedGenotyper_indel{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    publishDir { "${params.output}/indel/"}
    input:
        set sample_name,file('sample.recal.final.clean.bam') from recal_bam_res1
    output:
        set sample_name,file('*.indel.vcf') into indel_vcf_res
    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T UnifiedGenotyper \
            -R ${params.database}/bwa/chr1-x.fasta \
            -I sample.recal.final.clean.bam \
            -glm INDEL \
            -D ${params.database}/GATK/dbsnp_137.hg19.vcf \
            -o ${sample_name}.indel.vcf \
            -stand_call_conf 30 \
            -baqGOP 30 \
            -L /DG/programs/beta/scripts/DNA/genetest/adult/dgadultindel.bed \
            -nct 12 \
            -U ALLOW_SEQ_DICT_INCOMPATIBILITY -A VariantType -A QualByDepth \
            -A HaplotypeScore -A BaseQualityRankSumTest \
            -A MappingQualityRankSumTest -A ReadPosRankSumTest \
            -A FisherStrand -A DepthPerAlleleBySample \
            -A ClippingRankSumTest
        """
}

process genotype{
    container 'registry.cn-hangzhou.aliyuncs.com/bio/script-tool'
    publishDir {"${params.output}/genotype/"}

    input:
        set sample_name,file('snp.vcf') from snp_vcf_res
        set sample_name_1,file('index.vcf') from indel_vcf_res
        file('sample.target.basedepth.sample_interval_summary') from depth_res

    output:
        file("${sample_name}.geno") into genotype_res

    script:
        """
        Rscript $baseDir/bin/dgadultgenotype.R --args -o snp.vcf,index.vcf,sample.target.basedepth.sample_interval_summary,./${sample_name}
        """
}

workflow.onComplete {
    def msg="""
Pipeline execution summary
---------------------------
ScriptId        :   ${workflow.scriptId}
ScriptName      :   ${workflow.scriptName}
scriptFile      :   ${workflow.scriptFile}
Repository      :   ${workflow.repository?:'-'}
Revision        :   ${workflow.revision?:'-'}
ProjectDir      :   ${workflow.projectDir}
LaunchDir       :   ${workflow.launchDir}
ConfigFiles     :   ${workflow.configFiles}
Container       :   ${workflow.container}
CommandLine     :   ${workflow.commandLine}
Profile         :   ${workflow.profile}
RunName         :   ${workflow.runName}
SessionId       :   ${workflow.sessionId}
Resume          :   ${workflow.resume}
Start           :   ${workflow.start}

Completed at    :   ${workflow.complete}
Duration        :   ${workflow.duration}
Success         :   ${workflow.success}
Exit status     :   ${workflow.exitStatus}
ErrorMessage    :   -
Error report    :   -
"""
    log.info(msg)
    sendMail(
        to: 'panyunlai@126.com',
        subject: 'dna workflow run complete！',
        body:msg,
        attach:"${workflow.launchDir}/report.html"
    )
}