runs=($(ls /scratch/$USER/Dengue/folders_for_fit))
i=0
args=("$@")
for i in "${!runs[@]}"; do
    name=${runs[i]}
    mkdir /scratch/$USER/Dengue/out/output/$name
    mkdir /scratch/$USER/Dengue/out/results/$name
    mkdir /scratch/$USER/Dengue/out/log/$name
    mkdir /scratch/$USER/Dengue/out/stats/$name
    sbatch \
        --output=/scratch/$USER/Dengue/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=${args[i]} \
        /scratch/$USER/Dengue/helpers/run.slurm
done
