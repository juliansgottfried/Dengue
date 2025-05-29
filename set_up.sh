mkdir out
mkdir out/log
mkdir out/results
mkdir out/stats
mkdir out/output

mkdir folders_for_fit

read -p "Enter email: " email_address
echo "#SBATCH --mail-user=$email_address" > helpers/parts/tmp.txt
cat helpers/parts/top.txt helpers/parts/tmp.txt helpers/parts/bottom.txt > helpers/run.slurm
rm helpers/parts/tmp.txt

echo ""
echo "Activating renv..."

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript helpers/init_renv.R
