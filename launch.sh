path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

read -p $'Number of initial parameter combinations:\n(MUST be a multiple of 50, SHOULD be a multiple of 500) ' nseq
read -p $'Time limit (in minutes) for each run:\n(90 is a good first guess) ' run_time
read -p $'Number of mif2 refinement runs:\n(0, 1, or 2) ' n_refine
read -p "Panel fit? (y or n) " isPanel
read -p "Specific fit name (e.g., run_07_17_bb or n) " fit_name
read -p "Simulate MLE (e.g., y or n) " sim_mle

head -n 1 $path/helpers/email.sh > $path/helpers/line1.sh
echo "n_refine=$n_refine" > $path/helpers/line2.sh
echo "nseq=$nseq" > $path/helpers/line3.sh
echo "isPanel=$isPanel" > $path/helpers/line4.sh
cat $path/helpers/line*.sh > $path/helpers/email.sh
rm $path/helpers/line*.sh

source $path/helpers/email.sh

if [ $fit_name != "n" ]; then
    echo "running $fit_name"
    name=$fit_name
    mkdir $path/out/output/$name
    mkdir $path/out/results/$name
    mkdir $path/out/results_long/$name
    mkdir $path/out/log/$name
    mkdir $path/out/traces/$name
    mkdir $path/out/stats/$name

    sbatch \
        --output=$path/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
    --mail-user=$email \
        $path/helpers/run.slurm
else
   echo "running all fits in folders_for_fit"

    for i in "${!runs[@]}"; do
        name=${runs[i]}
        mkdir $path/out/output/$name
        mkdir $path/out/results/$name
        mkdir $path/out/results_long/$name
        mkdir $path/out/log/$name
        mkdir $path/out/traces/$name
        mkdir $path/out/stats/$name

        sbatch \
            --output=$path/out/output/$name/slurm-%A_%a.out \
            --job-name=$name \
            --time=$run_time \
            --mail-user=$email \
            $path/helpers/run.slurm
            
    done
fi