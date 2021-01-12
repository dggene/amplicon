#!/usr/bin/env nextflow



params.output='output'
params.tmp='/tmp'
params.reads1=null
params.reads2=null
params.sample=null

params.bwa_db_prefix=null
params.picard_jar_path=null
params.gatk_jar_path=null
params.ref_sequence=null
params.depth_target=null
params.gatk_default_target=null
params.gatk_snp_target=null
params.gatk_indel_target=null
params.v1000G_phase1_indels_hg19_vcf=null
params.Mills_and_1000G_gold_standard_indels_hg19_vcf=null
params.dbsnp_137_hg19_vcf=null
params.genotype_bed=null
params.bwa_cpu=2
params.gatk_cpu=2


params.help=null

version="v1.0.6"

log.info ""
log.info "------------------------------------------"
log.info "        DNA pipeline "
log.info "        version  ${version}"
log.info "run 'nextflow run main.nf --help' for get help "
log.info "------------------------------------------"
log.info ""

if(params.help){
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow main.nf --reads1=sample.R1.fq.gz --reads2=sample.R2.fq.gz --sample=sample_name --output=output_dir"
    log.info "Mandatory arguments:"
    log.info "--reads1                 String                reads1"
    log.info "--reads2                 String                reads2"
    log.info "--sample                 String                sample name"
    log.info "--output                 String                output dir"
    log.info "----------------------------------------------------------------"
    log.info "----------------------------------------------------------------"
   	log.info "Optional arguments (General):"
   	log.info "--tmp                    String                Name of folder that will contain the tmp file"
    log.info ""
    log.info "--------------------------------------------------------"
    exit 0
}else{
    log.info "fastq reads1          :  ${params.reads1}"
    log.info "fastq reads2          :  ${params.reads2}"
    log.info "sample name           :  ${params.sample}"          
    log.info "bwa_db_prefix         :  ${params.bwa_db_prefix}"
    log.info "ref_sequence          :  ${params.ref_sequence}"
    log.info "depth_target          :  ${params.depth_target}"
    log.info "gatk_default_target   :  ${params.gatk_default_target}"
    log.info "gatk_snp_target       :  ${params.gatk_snp_target}"
    log.info "gatk_indel_target     :  ${params.gatk_indel_target}"
    log.info "v1000G_phase1_indels_hg19_vcf   :  ${params.v1000G_phase1_indels_hg19_vcf}"
    log.info "Mills_and_1000G_gold_standard_indels_hg19_vcf   :  ${params.Mills_and_1000G_gold_standard_indels_hg19_vcf}"
    log.info "dbsnp_137_hg19_vcf    :  ${params.dbsnp_137_hg19_vcf}"
}
log.info "\n"

if(!params.reads1){
    exit 1,"Need params reads1"
}
if(!params.reads2){
    exit 1,"Need params reads2"
}
if(!params.sample){
    exit 1,"Need params sample"
}

reads1_file=file(params.reads1)
reads2_file=file(params.reads2)

if( !reads1_file.exists() ) {
  exit 1, "The specified input file does not exist: ${params.reads1}"
}

if( !reads2_file.exists() ) {
  exit 1, "The specified input file does not exist: ${params.reads2}"
}

sample=Channel.from(params.sample)
reads1=Channel.fromPath(params.reads1)
reads2=Channel.fromPath(params.reads2)

process soapnuke{
    tag { sample_name }
    input:
        val sample_name from sample
        file 'sample_R1.fq.gz' from reads1
        file 'sample_R2.fq.gz' from reads2

    output:
        set sample_name,file("*.fastq.gz") into clean_samples

    script:
    """
    SOAPnuke filter -1 sample_R1.fq.gz -2 sample_R2.fq.gz -l 15 -q 0.5 -Q 2 -o . \
        -C sample.clean1.fastq.gz -D sample.clean2.fastq.gz
    """
}

process bwa{
    input:
        set sample_name,files from clean_samples
    
    output:
        set sample_name,file('sample.clean.sam') into sam_res

    script:
    """
    bwa mem -t ${params.bwa_cpu} -M -T 30 ${params.bwa_db_prefix} \
        -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:DG\tPU:illumina" \
        ${files[0]} ${files[1]} > sample.clean.sam
    """     
}

process sort{
    //picard=1.1150
    input:
        set sample_name,file('sample.clean.sam') from sam_res
    
    output:
        set sample_name,file('sample.clean.bam') into bam_res

    script:
    """
    java -Xmx2g -Djava.io.tmpdir=${params.tmp} -jar ${params.picard_jar_path}  SortSam \
        INPUT=sample.clean.sam \
        OUTPUT=sample.clean.bam \
        MAX_RECORDS_IN_RAM=500000 \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=true
    """
}

bam_res.into{bam_res0;bam_res1;bam_res2;bam_res3}



process depth{
    label 'gatk'
    //gatk=3.6
    input:
        set sample_name,file('sample.clean.bam') from bam_res1
    output:
        file('sample.target.basedepth.sample_interval_summary') into depth_res

    script:
        """
        java -Xmx15g -jar ${params.gatk_jar_path} \
            -T DepthOfCoverage \
            -R ${params.ref_sequence} \
            -L ${params.depth_target} \
            -I sample.clean.bam \
            -ct 1 -ct 10 -ct 20 -ct 30 -ct 50 -ct 100 -ct 200 -ct 1000 \
            -o sample.target.basedepth
        """
}

process UnifiedGenotyper_snp{
    label 'gatk'
    //gatk=3.6
    publishDir "${params.output}/snp/", mode: 'copy'
    input:
        set sample_name,file('sample.clean.bam') from bam_res0
    output:
        set sample_name,file('*.snp.vcf') into snp_vcf_res
    script:
        """
        java -Xmx15g -jar ${params.gatk_jar_path} \
            -T UnifiedGenotyper \
            -R ${params.ref_sequence} \
            -I sample.clean.bam \
            -glm SNP \
            -D ${params.dbsnp_137_hg19_vcf} \
            -o ${sample_name}.snp.vcf \
            -stand_call_conf 30 \
            -baqGOP 30 \
            -L ${params.gatk_snp_target} \
            -nct ${params.gatk_cpu} \
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
    label 'gatk'
    //gatk=3.6
    publishDir "${params.output}/indel/", mode: 'copy'
    input:
        set sample_name,file('sample.clean.bam') from bam_res2
    output:
        set sample_name,file('*.indel.vcf') into indel_vcf_res
    script:
        """
        java -Xmx15g -jar ${params.gatk_jar_path} \
            -T UnifiedGenotyper \
            -R ${params.ref_sequence} \
            -I sample.clean.bam \
            -glm INDEL \
            -D ${params.dbsnp_137_hg19_vcf} \
            -o ${sample_name}.indel.vcf \
            -stand_call_conf 30 \
            -baqGOP 30 \
            -L ${params.gatk_indel_target} \
            -nct ${params.gatk_cpu} \
            -U ALLOW_SEQ_DICT_INCOMPATIBILITY -A VariantType -A QualByDepth \
            -A HaplotypeScore -A BaseQualityRankSumTest \
            -A MappingQualityRankSumTest -A ReadPosRankSumTest \
            -A FisherStrand -A DepthPerAlleleBySample \
            -A ClippingRankSumTest
        """
}

process genotype{
    publishDir "${params.output}", mode: 'copy'
    stageInMode "copy"
    
    input:
        set sample_name,file(snp_vcf) from snp_vcf_res
        set sample_name_1,file(indel_vcf) from indel_vcf_res
        file('sample.target.basedepth.sample_interval_summary') from depth_res

    output:
        file("${sample_name}.geno") into genotype_res

    script:
        """
        Rscript ${workflow.projectDir}/scripts/dgsnpgenotype.R --args \
        $snp_vcf \
        $indel_vcf \
        sample.target.basedepth.sample_interval_summary \
        ./${sample_name} \
        ${params.snpbed} \
        ${params.indelbed}
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
}
