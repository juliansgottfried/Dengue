runs=($(ls folders_for_fit))

for i in "${!runs[@]}"; do
    name=${runs[i]}
    read -p "Enter time limit (in minutes) for $name: " run_time 
    mkdir out/output/$name
    mkdir out/results/$name
    mkdir out/log/$name
    mkdir out/stats/$name
    sbatch \
        --output=out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
        helpers/run.slurm
done
