INput=$1
. $INput
. $(dirname $(realpath $0))/env.Variable

mkdir -p $WorkSH
cp $INput $SampIN

cat >$SampMK<<EOF
OUT=$OUT
SID=$SID
UID=$UID
WORKFLOW_DIR=$WORKFLOW_DIR

CATEGORY="fastqclink.${ID}"
CORES=7
MEMORY=20000
DISK=10
${Fawfq1} ${Fawfq2} : ${R1} ${R2}
    sh $WORKFLOW_DIR/FastqLink $SampIN

CATEGORY="cutadapt.${ID}"
CORES=7
MEMORY=20000
DISK=10
${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz : ${Fawfq1} ${Fawfq2}
    sh $WORKFLOW_DIR/CUTAdapt $SampIN

CATEGORY="fastqc.${ID}"
CORES=7
MEMORY=20000
DISK=10
${FASTQCpre}_raw_1_fastqc.html ${FASTQCpre}_trimmed_2_fastqc.html : $Fawfq1 $Fawfq2 ${CUTADAPTpre}_trimmed_1.fq.gz ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/FastQC $SampIN

CATEGORY="trim_galore.${ID}"
CORES=7
MEMORY=20000
DISK=10
${TRIMGALpre}_val_1.fq.gz  ${TRIMGALpre}_val_2.fq.gz : ${Fawfq1} ${Fawfq2}
    sh $WORKFLOW_DIR/TRIMGalore $SampIN

CATEGORY="bwa_mem.${ID}"
CORES=7
MEMORY=20000
DISK=10
${BAMpre}.bam : ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/BWAmem $SampIN ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz

CATEGORY="bwa.trans.${ID}"
CORES=7
MEMORY=20000
DISK=10
${BAMpre}.sort.bam $BAMmkdup ${BAMmkdup}.bai : ${BAMpre}.bam
    sh $WORKFLOW_DIR/BAMtrans $SampIN ${BAMpre}.bam

CATEGORY="bam.unmap.filter.${ID}"
CORES=7
MEMORY=20000
DISK=10
${BAMpre}.sort.markdup.unmapped.bam $BAMfilter ${BAMfilter}.bai : $BAMmkdup
    sh $WORKFLOW_DIR/BAMfilter $SampIN $BAMmkdup

CATEGORY="bam.stat.${ID}"
CORES=7
MEMORY=20000
DISK=10
$BAMstatout : $BAMmkdup
    sh $WORKFLOW_DIR/BAMstat $SampIN $BAMmkdup

CATEGORY="bam.deepstat.${ID}"
CORES=7
MEMORY=20000
DISK=10
$DEEPhmapgz : $BAMfilter
    sh $WORKFLOW_DIR/DEEPTstate $SampIN $BAMfilter

CATEGORY="callpeak.macs2.${ID}"
CORES=7
MEMORY=20000
DISK=10
${MACS2cppre}_peaks.narrowPeak ${MACS2cppre}_summits.bed : $BAMfilter 
    sh $WORKFLOW_DIR/MACS2_CallPeak  $SampIN  $BAMfilter 

CATEGORY="callpeak.genrich.${ID}"
CORES=7
MEMORY=20000
DISK=10
${GENRICHpre}.genrich.narrowPeak : $BAMfilter 
    sh $WORKFLOW_DIR/GENRICH_CallPeak $SampIN $BAMfilter 

CATEGORY="callpeak.annotate.homer.${ID}"
CORES=7
MEMORY=20000
DISK=10
${HOMERpre}_peaks.annotatePeaks.txt : ${MACS2cppre}_peaks.narrowPeak
    sh $WORKFLOW_DIR/HOMER_Anno $SampIN ${MACS2cppre}_peaks.narrowPeak

EOF

#$MAKEFLOW -T sge $SampMK
