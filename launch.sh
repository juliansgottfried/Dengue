runs=($(ls /scratch/$USER/Dengue/folders_for_fit))

for i in "${!runs[@]}"; do
    name=${runs[i]}
    read -p "Enter time limit (in minutes) for $name: " run_time 
    mkdir /scratch/$USER/Dengue/out/output/$name
    mkdir /scratch/$USER/Dengue/out/results/$name
    mkdir /scratch/$USER/Dengue/out/log/$name
    mkdir /scratch/$USER/Dengue/out/stats/$name
    sbatch \
        --output=/scratch/$USER/Dengue/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
        /scratch/$USER/Dengue/helpers/run.slurm
done
