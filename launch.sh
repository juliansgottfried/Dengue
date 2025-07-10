path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

read -p "Time limit (in minutes) for each run: " run_time
read -p "Number of mif2 refinements: " n_refine

head -n 1 email.sh > tmp1.sh
echo "n_refine=$n_refine" > tmp2.sh
cat tmp1.sh tmp2.sh > email.sh
rm tmp1.sh tmp2.sh

source $path/helpers/email.sh

for i in "${!runs[@]}"; do
    name=${runs[i]}
    mkdir $path/out/output/$name
    mkdir $path/out/results/$name
    mkdir $path/out/log/$name
    mkdir $path/out/traces/$name
    mkdir $path/out/stats/$name

    line_count=$(cat $path/folders_for_fit/$name/pars.csv | wc -l | tr -d ' ')
    line_count_adj=$(expr $line_count - 1)
    array_size=$(expr $line_count_adj / 10)

    sbatch \
        --output=$path/out/output/$name/slurm-%A_%a.out \
        --job-name=$name \
        --time=$run_time \
	--array=1-$array_size \
	--mail-user=$email \
        $path/helpers/run.slurm
done
