runs=($(ls /scratch/$USER/folders_for_fit))
i=0
args=("$@")
for i in "${!runs[@]}"; do
    name=${runs[i]}
    mkdir /scratch/$USER/out/output/$name
    mkdir /scratch/$USER/out/results/$name
    mkdir /scratch/$USER/out/log/$name
    mkdir /scratch/$USER/out/stats/$name
    sbatch \
        --output=/scratch/$USER/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=${args[i]} \
        /scratch/$USER/helpers/run.slurm
done
