path=/scratch/$LOGNAME/Dengue 

source $path/helpers/email.sh

bash $path/helpers/get_times.sh

module unload r
module unload gcc
module load r/gcc/4.3.1
Rscript $path/helpers/collate.R $isPanel

rm $path/times.txt
