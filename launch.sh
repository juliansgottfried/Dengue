path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

read -p $'Number of initial parameter combinations:\n(MUST be a multiple of 50, SHOULD be a multiple of 500) ' nseq
read -p $'Time limit (in minutes) for each run:\n(90 is a good first guess, or 360 for panel fits) ' run_time
read -p $'Number of mif2 refinement runs:\n(0, 1, or 2) ' n_refine
read -p $'Panel fit?\n(y or n) ' isPanel
read -p $'Launch only one fit from folders_for_fit, or all:\n(one or all) ' which_fit
fit_name=""
if [ $which_fit == "one" ]; then
	read -p $'fit_name:\n(e.g. run_01_01_a) ' fit_name
fi

echo "n_refine=$n_refine" > $path/tmp_line1.sh
echo "nseq=$nseq" > $path/tmp_line2.sh
echo "isPanel=$isPanel" > $path/tmp_line3.sh
cat $path/tmp_line*.sh > $path/hyperparams.sh
rm $path/tmp_line*.sh

source $path/helpers/user_info.sh

if [ $which_fit == "one" ]; then
    echo "Running $fit_name"

    mv $path/hyperparams.sh $path/folders_for_fit/$name
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
   echo "Running all fits"

   cp $path/hyperparams.sh $path/folders_for_fit/$name
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
    rm $path/hyperparams.sh
    done
fi
