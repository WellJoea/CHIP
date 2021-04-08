. $1
. $(dirname $(realpath $0))/env.Configure
. $(dirname $(realpath $0))/env.Variable

mkdir -p $SPLADDERout
echo "`date +%Y%m%d-%H:%M:%S`: spladder build"



STAR_gtf=/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.100.gtf
STAR_ref=/share/home/share/Repository/GenomeDB/Reference/Homo_Sapiens/ENSEMBL/Homo_sapiens.GRCh38.dna.primary_assembly.fa
INBAM=/data/zhouwei/02production/20200914_1800/U2OS-b3-R11_L1_Q811605/SS2/HISAT2/U2OS-b3-R11_L1_Q811605__D1__L1.hisat2.sort.bam
INBAM=/data/zhouwei/02production/20200914_1800/U2OS-b3-R11_L1_Q811605/SS2/STAR/U2OS-b3-R11_L1_Q811605__D1.Aligned.sortedByCoord.sort.bam
SPLADDER=/share/home/share/software/Python-3.8.3/bin/spladder

python3 $SPLADDER build \
    -o /data/zhouwei/02production/20200914_1800/U2OS-b3-R11_L1_Q811605/SS2/Spladder3 \
    -l aa.log \
    --sparse-bam \
    -a $STAR_gtf \
    -b $INBAM \
    --labels U2OS-b3-R11_L1_Q811605 \
    --merge-strat single \
    --extract-ase \
    --quantify-graph \
    --output-txt \
    --output-gff3 \
    --output-struc \
    --output-struc-conf \
    --output-bed \
    --output-conf-bed
#    --pyproc
