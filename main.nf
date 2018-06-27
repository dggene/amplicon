#!/usr/bin/env nextflow



params.output='output'
params.tmp='/tmp'
database="/DG/database/genomes/Homo_sapiens/hg19"
target="/DG/programs/beta/scripts/DNA/genetest/adult/dgadulttarget.bed"

Channel.fromFilePairs('/DG/project0/DGgene/rawdata/20180621L3-SR16008/*_S74_L003_R{1,2}_001.fastq.gz').into{samples0;samples1;samples2}
samples2.println()

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
    bwa mem -t 10 -M -T 30 $database/bwa/chr1-x \
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
            REFERENCE_SEQUENCE=${database}/bwa/chr1-x.fasta
        """
}

process depth{
    container 'broadinstitute/genomes-in-the-cloud:2.3.1-1512499786'
    //gatk=3.6
    input:
        set sample_name,file('sample.clean.bam') from bam_res1
    

    script:
        """
        java -Xmx15g -jar /usr/gitc/GATK36.jar \
            -T DepthOfCoverage \
            -R ${database}/bwa/chr1-x.fasta \
            -L $target \
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
            -R ${database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -o sample.realigner.dedupped.clean.intervals \
            -I sample.clean.bam \
            -known $database/GATK/1000G_phase1.indels.hg19.vcf \
            -known $database/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
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
            -R ${database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -I sample.clean.bam \
            -targetIntervals sample.realigner.dedupped.clean.intervals \
            -o sample.realigned.clean.bam \
            -known $database/GATK/1000G_phase1.indels.hg19.vcf \
            -known $database/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
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
            -R ${database}/bwa/chr1-x.fasta \
            -L /DG/project0/DGgene/adult/20180621/dgadult.bed \
            -I sample.realigned.clean.bam \
            -knownSites ${database}/GATK/dbsnp_137.hg19.vcf \
            -knownSites ${database}/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf \
            -knownSites ${database}/GATK/1000G_phase1.indels.hg19.vcf \
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
            -R ${database}/bwa/chr1-x.fasta \
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
            -R ${database}/bwa/chr1-x.fasta \
            -I sample.recal.final.clean.bam \
            -glm SNP \
            -D ${database}/GATK/dbsnp_137.hg19.vcf \
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
            -R ${database}/bwa/chr1-x.fasta \
            -I sample.recal.final.clean.bam \
            -glm INDEL \
            -D ${database}/GATK/dbsnp_137.hg19.vcf \
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

}

workflow.onComplete {
    log.info """
Pipeline execution summary
---------------------------
ScriptId    :   ${workflow.scriptId}
ScriptName  :   ${workflow.scriptName}
scriptFile  :   ${workflow.scriptFile}
Repository  :   ${workflow.repository?:'-'}
Revision    :   ${workflow.revision?:'-'}
ProjectDir  :   ${workflow.projectDir}
LaunchDir   :   ${workflow.launchDir}
ConfigFiles :   ${workflow.configFiles}
Container   :   ${workflow.container}
CommandLine :   ${workflow.commandLine}
Profile     :   ${workflow.profile}
RunName     :   ${workflow.runName}
SessionId   :   ${workflow.sessionId}
Resume      :   ${workflow.resume}
Start       :   ${workflow.start}

Completed at:   ${workflow.complete}
Duration    :   ${workflow.duration}
Success     :   ${workflow.success}
Exit status :   ${workflow.exitStatus}
ErrorMessage:   -
Error report:   -
"""
}