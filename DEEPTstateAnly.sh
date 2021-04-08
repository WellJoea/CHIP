#. $(dirname $(realpath $0))/env.Configure

WD=$1
IN=/datd/zhouwei/02production/20210223_0957
OU=$WD/MergePlot

Blacklist=/share/home/share/Repository/GenomeDB/Blacklist/mm10-blacklist.v2.bed
GENE_GTF=/share/home/share/Repository/GenomeDB/Reference/Mus_Musculus/ENSEMBL/Mus_musculus.GRCm38.100.gtf
mkdir -p  $OU

echo "`date +%Y%m%d-%H:%M:%S`: Deeptools multiple plot"
AID=(ND1-cor-1w-IgG ND1-cor-1w-K27me3 ND1-cor-1w-ND1)
ABAM=( ${IN}/ND1-cor-1w-IgG/CHIP/BAM_Align/ND1-cor-1w-IgG__D1.sort.markdup.remove2820.bam 
        ${IN}/ND1-cor-1w-K27me3/CHIP/BAM_Align/ND1-cor-1w-K27me3__D1.sort.markdup.remove2820.bam 
        ${IN}/ND1-cor-1w-ND1/CHIP/BAM_Align/ND1-cor-1w-ND1__D1.sort.markdup.remove2820.bam 
    )

CaseID=(ND1-cor-1w-IgG ND1-cor-1w-K27me3 ND1-cor-1w-ND1)
CaseBIG=( ${IN}/ND1-cor-1w-IgG/CHIP/DEEPT_Stat/ND1-cor-1w-IgG__D1.SeqDepthNorm.bw
          ${IN}/ND1-cor-1w-K27me3/CHIP/DEEPT_Stat/ND1-cor-1w-K27me3__D1.SeqDepthNorm.bw
          ${IN}/ND1-cor-1w-ND1/CHIP/DEEPT_Stat/ND1-cor-1w-ND1__D1.SeqDepthNorm.bw
)

getgenepred()
{
    gtfToGenePred -genePredExt -geneNameAsName2 $GENE_GTF $OU/GRCm38.100.genepred
    awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' $OU/GRCm38.100.genepred >  $OU/GRCm38.100.genepred.bed
    #plotFingerprint \
    # -b testFiles/*bam \
    #--labels H3K27me3 H3K4me1 H3K4me3 H3K9me3 input \
    #--minMappingQuality 30 --skipZeros \
    #--region 19 --numberOfSamples 50000 \
    #-T "Fingerprints of different samples"  \
    #--plotFile fingerprints.png \
    #--outRawCounts fingerprints.tab

    #plotFingerprint \
    #    -p 10 \
    #    --bamfiles $TBAM \
    #    --plotFile ${DEEPTpre}.plotFingerprint.pdf \
    #    --outRawCounts ${DEEPTpre}.plotFingerprint.raw.txt \
    #    --outQualityMetrics ${DEEPTpre}.plotFingerprint.qcmetrics.txt \
    #    --skipZeros
}


multisum()
{
    #multiBamSummary bins \
    #    --bamfiles ${ABAM[@]} \
    #    --labels ${AID[@]} \
    #    --minMappingQuality 10 \
    #    -p 20 \
    #    --blackListFileName $Blacklist \
    #    -out $OU/Bam.readCounts.npz \
    #    --outRawCounts $OU/Bam.readCounts.tab

    multiBigwigSummary bins \
        --bwfiles ${CaseBIG[@]} \
        --labels  ${CaseID[@]} \
        -p 20 \
        --blackListFileName $Blacklist \
        -out $OU/Bw.readCounts.npz \
        --outRawCounts $OU/Bw.readCounts.tab
}

plotsumm()
{
    RC=$1
    PH=${RC/.npz/}
    plotCorrelation \
        -in $RC  \
        -o  ${PH}.Correlation.pearson.scatter.pdf \
        --outFileCorMatrix ${PH}.Correlation.pearson.scatter.tab \
        --plotTitle  Correlation.pearson \
        --skipZeros \
        --colorMap RdBu \
        --corMethod pearson \
        --whatToPlot scatterplot 

    plotCorrelation \
        -in $RC \
        -o  ${PH}.Correlation.spearman.heatmap.pdf \
        --outFileCorMatrix ${PH}.Correlation.spearman.heatmap.tab \
        --plotTitle  Correlation.spearman \
        --skipZeros \
        --colorMap RdBu \
        --corMethod spearman \
        --plotNumbers \
        --whatToPlot heatmap 

    plotPCA \
        -in $RC \
        -o ${PH}.PCA.pdf
}

peakCom()
{
    computeMatrix scale-regions \
        -p 40 \
        -S ${CaseBIG[@]} \
        --samplesLabel ${CaseID[@]} \
        -R $GENE_GTF \
        --startLabel TSS \
        --endLabel TES \
        --beforeRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
        --skipZeros \
        --outFileName $OU/scale.regions.mat.gz \
        --outFileNameMatrix  $OU/scale.regions.IndividualValues.tab \
        --outFileSortedRegions $OU/scale.regions.HeatmapsortedRegions.bed \
        --blackListFileName $Blacklist

    computeMatrix reference-point \
        -p 40 \
        -R $GENE_GTF \
        -S ${CaseBIG[@]} \
        --samplesLabel ${CaseID[@]} \
        --referencePoint TSS \
        --beforeRegionStartLength 3000 \
        --afterRegionStartLength 3000 \
        --skipZeros \
        --outFileName $OU/ref.point.mat.gz \
        --outFileNameMatrix  $OU/ref.point.IndividualValues.tab \
        --outFileSortedRegions $OU/ref.point.HeatmapsortedRegions.bed \
        --blackListFileName $Blacklist 

}

plotpeak()
{
    MT=$1
    PH=${MT/.mat.gz/}

    plotProfile \
        --matrixFile $MT \
        --dpi 300 \
        --plotFileFormat pdf \
        --outFileName ${PH}.plotProfile.pdf \
        --outFileNameData ${PH}.plotProfile.tab
    ##--outFileSortedRegions ${DEEPTpre}.plotProfile.bed \

    plotHeatmap \
        --matrixFile $MT \
        --outFileName ${PH}.plotHeatmap.pdf \
        --outFileNameMatrix ${PH}.plotHeatmap.mat.tab 
    #    --colorMap YlGn \
    #    --kmeans 4  \
    #    --heatmapHeight 40 \
    #    --heatmapWidth 10
    #--whatToShow 'heatmap and colorbar' \
    #--outFileSortedRegions ${DEEPTpre}.plotHeatmap.bed \
    # --colorMap viridis \
    #--zMin -3 --zMax 3 --kmeans 4 \
}

multisum
plotsumm $OU/Bw.readCounts.npz
#plotsumm $OU/Bam.readCounts.npz
peakCom
plotpeak $OU/scale.regions.mat.gz 
plotpeak $OU/ref.point.mat.gz 

pygeno()
{
    make_tracks_file --trackFiles ${CaseBIG[@]} -o $OU/tracks.ini

    #gtfToGenePred

    pyGenomeTracks \
        --tracks $OU/tracks.ini \
        --outFileName $OU/all.ref.pdf \
        --BED /datd/zhouwei/02production/20210106_15140/MergePlot/hg38.ref.bed
}
#pygeno

<<COMMENT
/share/home/share/Pipeline/13SCATAC/ATACPipe1/ATACtools.py mergepeak \
    -B  /data/zhouwei/02production/20200721_17170/3T3-ATAC-sc4/ATAC/MACS2_CP/3T3-ATAC-sc4__D1_peaks.narrowPeak,/data/zhouwei/02production/20200721_17170/3T3-ATAC-sc5/ATAC/MACS2_CP/3T3-ATAC-sc5__D1_peaks.narrowPeak,/data/zhouwei/02production/20200721_17170/3T3-ATAC-sc6/ATAC/MACS2_CP/3T3-ATAC-sc6__D1_peaks.narrowPeak,/data/zhouwei/02production/20200721_17170/3T3-ATAC-sc7/ATAC/MACS2_CP/3T3-ATAC-sc7__D1_peaks.narrowPeak \
    -o ./

#grep "#sequence-region"  /share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gff3 |awk 'BEGIN{FS=" ";OFS="\t"}''{print $2,$3,$4}' > ref.bed
/data/zhouwei/02production/20200721_17170/3T3-ATAC-sc4/ATAC/WorkShell/get.merge.ref.bed.py
make_tracks_file \
    --trackFiles /data/zhouwei/02production/20200721_17170/3T3-ATAC-sc4/ATAC/DEEPTStat/3T3-ATAC-sc4__D1.SeqDepthNorm.rpkm.bw \
                 /data/zhouwei/02production/20200721_17170/3T3-ATAC-sc5/ATAC/DEEPTStat/3T3-ATAC-sc5__D1.SeqDepthNorm.rpkm.bw \
                 /data/zhouwei/02production/20200721_17170/3T3-ATAC-sc6/ATAC/DEEPTStat/3T3-ATAC-sc6__D1.SeqDepthNorm.rpkm.bw \
                 /data/zhouwei/02production/20200721_17170/3T3-ATAC-sc7/ATAC/DEEPTStat/3T3-ATAC-sc7__D1.SeqDepthNorm.rpkm.bw \
                 /share/home/share/Repository/GenomeDB/Reference/Mus_Musculus/ENSEMBL/Mus_musculus.GRCm38.100.gtf \
                 -o tracks.ini

pyGenomeTracks \
    --tracks tracks.ini \
    --outFileName all.nice_image.pdf \
    --BED ref.bed

pyGenomeTracks \
    --tracks tracks.ini \
    --outFileName all.ref.pdf \
    --BED merge.peaks_REF.bed

pyGenomeTracks \
    --tracks tracks.ini \
    --outFileName all.ref.pdf \
    --BED merge.peaks_REF.bed
COMMENT
