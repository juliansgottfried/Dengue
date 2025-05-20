runs=($(ls folders_for_fit))
i=0
args=("$@")
for i in "${!runs[@]}"; do
    name=${runs[i]}
    mkdir out/output/$name
    mkdir out/results/$name
    mkdir out/log/$name
    mkdir out/stats/$name
    mkdir results/$name
    sbatch \
        --output=out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=${args[i]} \
        run.slurm
done
