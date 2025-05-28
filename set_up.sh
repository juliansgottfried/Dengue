mkdir /scratch/$USER/Dengue/out
mkdir /scratch/$USER/Dengue/out/log
mkdir /scratch/$USER/Dengue/out/results
mkdir /scratch/$USER/Dengue/out/stats
mkdir /scratch/$USER/Dengue/out/output

mkdir /scratch/$USER/Dengue/folders_for_fit

read -p "Enter email: " email_address
echo "#SBATCH --mail-user=$email_address" > /scratch/$USER/Dengue/helpers/parts/tmp.txt
cat /scratch/$USER/Dengue/helpers/parts/top.txt /scratch/$USER/Dengue/helpers/parts/tmp.txt /scratch/$USER/Dengue/helpers/parts/bottom.txt > /scratch/$USER/Dengu>
rm /scratch/$USER/Dengue/helpers/parts/tmp.txt

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript /scratch/$USER/Dengue/helpers/init_renv.R
