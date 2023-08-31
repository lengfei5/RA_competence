#####################################
# this script is to do motifs analysis 
# using MEME suite
# the meme 4.12.0 is loaded from the module (current version)
#  
####################################
ml bedtools/2.25.0-foss-2018b
ml load meme/5.1.1-foss-2018b-python-3.6.6

#MEME_path="/groups/cochella/jiwang/local/meme/bin/"
cwd=`pwd`;

resDir="/groups/tanaka/People/current/jiwang/projects/RA_competence/results/scRNAseq_R13547_10x_mNT_20220813/motif_analysis"
input_peaks=${resDir}/peaks_bed
out_fasta=${resDir}/seqs_fasta
out_ame=${resDir}/ame
out_meme=${resDir}/meme_chip

mkdir -p $out_ame
mkdir -p $out_meme
mkdir -p $out_fasta

## prepare the background sequence .fa file
#cd /groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens
#bedtools getfasta -fi Mus_musculus.GRCm38.dna.toplevel.fa -bed mm10_background_2000.random.promoters.bed -name -fo mm10_background_2000.random.promoters.fa

#nb_cores=2;
pwms="/groups/tanaka/People/current/jiwang/Databases/motifs_TFs/JASPAR2022/JASPAR2022_CORE_vertebrates_nonRedundant.meme"
genome="/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.dna.toplevel.fa"
bg="/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/mm10_background_2000.random.promoters.fa"
#primMotif="/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/tbx_motif.meme"
pval=0.0001;

#for bb in `ls $input_peaks/*.bed `
#do 
    #echo $bb;
#    filename=`basename $bb`
#    filename="${filename%.*}"
    #echo $filename
#    if [ ! -e "${DIR}/${filename}.fa"  ]; then
#	bedtools getfasta -fi $genome -bed $bb -name -fo ${out_fasta}/${filename}.fa 
#    fi
#done

for ff in `ls ${out_fasta}/*.fa`
do 
    filename=`basename $ff`
    filename="${filename%.*}"
    echo $filename
    
    #meme-chip -db $pwms $ff -oc ${out_meme}/${filename} -centrimo-local;   
    #ame --oc ${out_ame}/${filename} --control $bg --scoring totalhits --method fisher $ff $pwms;
    ame --oc ${out_ame}/${filename} --control --shuffle--  --scoring avg --method fisher $ff $pwms;

    #break
done
