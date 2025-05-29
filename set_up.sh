path=/scratch/$LOGNAME/Dengue 

mkdir $path/out
mkdir $path/out/log
mkdir $path/out/results
mkdir $path/out/stats
mkdir $path/out/output

mkdir $path/folders_for_fit

read -p "Enter email: " email_address
echo "#SBATCH --mail-user=$email_address" > $path/helpers/parts/tmp.txt
cat $path/helpers/parts/top.txt $path/helpers/parts/tmp.txt $path/helpers/parts/bottom.txt > $path/helpers/run.slurm
rm $path/helpers/parts/tmp.txt

echo ""
echo "Activating renv..."

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript $path/helpers/init_renv.R
