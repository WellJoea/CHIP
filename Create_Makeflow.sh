INpath=$1
SID=$2
Module=$3

Lanemkflow()
{
    INp=$1
    OUp=$2
    cat $INp >> $OUp
    cat>>$OUp<<EOF


CATEGORY="fastqclink.${AID}"
CORES=7
MEMORY=30000
DISK=10
${RawFQpre}_r1.fq.gz ${RawFQpre}_r2.fq.gz : ${R1} ${R2}
    sh $WORKFLOW_DIR/FastqLink $INp

CATEGORY="cutadapt.trim.${AID}"
CORES=12
MEMORY=30000
DISK=10
${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz
    sh $WORKFLOW_DIR/CUTAdapt $INp

CATEGORY="fastqc.${AID}"
CORES=7
MEMORY=30000
DISK=10
${FASTQCpre}_r1_fastqc.html ${FASTQCpre}_trimmed_2_fastqc.html : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz ${CUTADAPTpre}_trimmed_1.fq.gz ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/FastQC $INp

CATEGORY="trim_galore.${AID}"
CORES=7
MEMORY=30000
DISK=10
${TRIMGALpre}_val_1.fq.gz  ${TRIMGALpre}_val_2.fq.gz : ${RawFQpre}_r1.fq.gz  ${RawFQpre}_r2.fq.gz
    sh $WORKFLOW_DIR/TRIMGalore $INp

CATEGORY="bwa_mem.${AID}"
CORES=7
MEMORY=20000
DISK=10
${BAMout}/${AID}.bam : ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz
    sh $WORKFLOW_DIR/BWAmem $INp ${CUTADAPTpre}_trimmed_1.fq.gz  ${CUTADAPTpre}_trimmed_2.fq.gz


EOF
}

Repmkflow()
{
    INp=$1
    OUp=$2
    cat $INp >> $OUp
    cat>>$OUp<<EOF


CATEGORY="bwa.trans.${TID}"
CORES=7
MEMORY=20000
DISK=10
$BAMmkdup : ${BAMreplane[@]}
    sh $WORKFLOW_DIR/BAMtrans $INp

CATEGORY="bam.filter.${TID}"
CORES=7
MEMORY=20000
DISK=10
$BAMfilter ${BAMfilter}.bai : $BAMmkdup
    sh $WORKFLOW_DIR/BAMfilter $INp $BAMmkdup

CATEGORY="bam.stat.${TID}"
CORES=7
MEMORY=20000
DISK=10
$BAMstatout : $BAMmkdup
    sh $WORKFLOW_DIR/BAMstat $INp $BAMmkdup

CATEGORY="bam.deepstat.${TID}"
CORES=7
MEMORY=20000
DISK=10
${DEEPTpre}.heatmap.gz ${DEEPTpre}.Heatmap1sortedRegions.bed : $BAMfilter $BMAControl
    sh $WORKFLOW_DIR/DEEPTstate $INp $BAMfilter $BMAControl
${BIGWIGpre}.scaled.bigWig : $BAMfilter
    sh $WORKFLOW_DIR/BIGWIG  $INp $BAMfilter

CATEGORY="callpeak.macs2.${TID}"
CORES=7
MEMORY=20000
DISK=10
${MACS2cppre}_peaks.narrowPeak ${MACS2cppre}_summits.bed : $BAMfilter $BMAControl
    sh $WORKFLOW_DIR/MACS2_CallPeak  $INp  $BAMfilter $BMAControl
${MACS2cppre}_peaks.narrowPeak.fCounts.txt : ${MACS2cppre}_peaks.narrowPeak $BAMfilter
    sh $WORKFLOW_DIR/FRIP $INp ${MACS2cppre}_peaks.narrowPeak $BAMfilter

CATEGORY="callpeak.genrich.${TID}"
CORES=7
MEMORY=20000
DISK=10
${GENRICHpre}_genrich.peaks.narrowPeak : $BAMfilter $BMAControl
    sh $WORKFLOW_DIR/GENRICH_CallPeak $INp $BAMfilter $BMAControl
${GENRICHpre}_genrich.peaks.narrowPeak.fCounts.txt : ${GENRICHpre}_genrich.peaks.narrowPeak $BAMfilter
    sh $WORKFLOW_DIR/FRIP $INp ${GENRICHpre}_genrich.peaks.narrowPeak $BAMfilter

CATEGORY="callpeak.annotate.homer.${TID}"
CORES=7
MEMORY=20000
DISK=10
${HOMERpre}_peaks.annotatePeaks.txt : ${MACS2cppre}_peaks.narrowPeak
    sh $WORKFLOW_DIR/HOMER_Anno $INp ${MACS2cppre}_peaks.narrowPeak
${HOMERpre}/homerResults.html : ${MACS2cppre}_peaks.narrowPeak
    sh $WORKFLOW_DIR/HOMER_Motifs $INp ${MACS2cppre}_peaks.narrowPeak

CATEGORY="ataqv.qc.${TID}"
CORES=7
MEMORY=20000
DISK=10
${ATAQCpre}.ataqv.json.gz : $BAMfilter ${MACS2cppre}_peaks.narrowPeak
    sh $WORKFLOW_DIR/ATAQV $INp $BAMfilter ${MACS2cppre}_peaks.narrowPeak
$ATAQCout/$SID/index.html : ${ATAQCrep[@]}
    sh $WORKFLOW_DIR/MKARV $INp
${IDRpre}.idrValues.txt : ${IDRPeak[@]}
    sh $WORKFLOW_DIR/IDR $INp


EOF
}


INSf=( $INpath/${SID}__D*${Module}.input )
. ${INSf[0]}
WorkSH=$AUT/WorkShell
SampMK=$WorkSH/${SID}_${Module}.mf

if [ -f $SampMK ];then
    rm $SampMK
fi
mkdir -p $WorkSH

#########each lane
for i in ${INSf[@]};do
    . $i
    . $WORKFLOW_DIR/env.Variable
    cp  $i $WorkSH/${AID}_${Module}.input
    Lanemkflow $WorkSH/${AID}_${Module}.input $SampMK
done

#########each rep
for j in ${TIDs[@]};do
    Tidf=( $WorkSH/${j}*${Module}.input )
    . ${Tidf[0]}
    . $WORKFLOW_DIR/env.Variable
    Repmkflow  ${Tidf[0]} $SampMK
done

