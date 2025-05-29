path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

for i in "${!runs[@]}"; do
    name=${runs[i]}
    read -p "Enter time limit (in minutes) for $name: " run_time 
    mkdir $path/out/output/$name
    mkdir $path/out/results/$name
    mkdir $path/out/log/$name
    mkdir $path/out/stats/$name
    sbatch \
        --output=$path/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
        $path/helpers/run.slurm
done
