path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

read -p $'Number of initial parameter combinations:\n(MUST be a multiple of 50, SHOULD be a multiple of 500) ' nseq
read -p "Time limit (in minutes) for each run: " run_time
read -p "Number of mif2 refinement runs: " n_refine

head -n 1 email.sh > line1.sh
echo "n_refine=$n_refine" > line2.sh
echo "nseq=$nseq" > line3.sh
cat line*.sh > email.sh
rm line*.sh

source $path/helpers/email.sh

for i in "${!runs[@]}"; do
    name=${runs[i]}
    mkdir $path/out/output/$name
    mkdir $path/out/results/$name
    mkdir $path/out/log/$name
    mkdir $path/out/traces/$name
    mkdir $path/out/stats/$name

#    line_count=$(cat $path/folders_for_fit/$name/pars.csv | wc -l | tr -d ' ')
#    line_count_adj=$(expr $line_count - 1)
#    array_size=$(expr $nseq / 10)

    sbatch \
        --output=$path/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
	--mail-user=$email \
        $path/helpers/run.slurm
done
