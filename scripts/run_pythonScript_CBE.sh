mkdir -p logs

inputDir="/groups/tanaka/People/current/jiwang/projects/RA_competence/images_data/CellPose_Embryos"
echo ${inputDir}

samples=`ls ${inputDir}/*_blur_cp_masks_cleaned.tif`

for f in ${samples};
do
    echo "---"
    echo "sample -- $f"
    fname=`basename ${f}`
    fname=${fname//"_blur_cp_masks_cleaned.tif"/}
    echo $fname
    script="logs/run_py_${fname}.sh"
    
    cat <<EOF > $script
#!/usr/bin/bash	
#SBATCH --export=ALL	
#SBATCH --qos=medium
#SBATCH --time=0-16:00:00
#SBATCH --partition=c
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 
#SBATCH -o logs/runPy_${fname}.out
#SBATCH -e logs/runpy_${fname}.err
#SBATCH --job-name runPy_${fname}


python collect_features_formMask.py $f

    
EOF
    cat $script
    sbatch $script
    #break
    
done


