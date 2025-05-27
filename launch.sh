runs=($(ls /scratch/$USER/Dengue/folders_for_fit/$USER))

times=($(cat /scratch/$USER/Dengue/run_times.txt))

for i in "${!runs[@]}"; do
    name=${runs[i]}
    time=${times[i]}
    mkdir /scratch/$USER/Dengue/out/output/$name
    mkdir /scratch/$USER/Dengue/out/results/$name
    mkdir /scratch/$USER/Dengue/out/log/$name
    mkdir /scratch/$USER/Dengue/out/stats/$name
    sbatch \
        --output=/scratch/$USER/Dengue/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$time \
        /scratch/$USER/Dengue/helpers/run.slurm
done
