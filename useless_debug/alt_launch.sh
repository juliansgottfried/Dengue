path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

nseq=1
run_time=5
n_refine=2
isPanel="y"
fit_name="n"

head -n 1 $path/helpers/email.sh > $path/helpers/line1.sh
echo "n_refine=$n_refine" > $path/helpers/line2.sh
echo "nseq=$nseq" > $path/helpers/line3.sh
echo "isPanel=$isPanel" > $path/helpers/line4.sh
cat $path/helpers/line*.sh > $path/helpers/email.sh
rm $path/helpers/line*.sh

source $path/helpers/email.sh

for i in "${!runs[@]}"; do
    name=${runs[i]}
#    sbatch \
#        --output=$path/out/output/$name/slurm-%A_%a.out \
#        --job-name=$name \
#        --time=$run_time \
#        --mail-user=$email \
#        $path/helpers/run.slurm

	module unload r
	module unload gcc
	module load r/gcc/4.3.1

	Rscript $path/debug/panel_fit.R \
    	1 \
    	1 \
    	1 \
    	"test" \
    	2 \
    	1
done
