path=/scratch/$LOGNAME/Dengue 

bash $path/helpers/get_times.sh

module unload r
module unload gcc
module load r/gcc/4.3.1
Rscript $path/helpers/collate.R

rm $path/times.txt
