path=/scratch/$LOGNAME/Dengue 

mkdir $path/out
mkdir $path/out/log
mkdir $path/out/results
mkdir $path/out/results_long
mkdir $path/out/stats
mkdir $path/out/traces
mkdir $path/out/output

mkdir $path/folders_for_fit

read -p "Enter email: " email

echo "email=$email" > $path/helpers/email.sh

echo ""
echo "Activating renv..."

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript $path/helpers/init_renv.R
