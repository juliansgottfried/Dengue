path=/scratch/$LOGNAME/Dengue 
module unload r
module unload gcc
module load r/gcc/4.3.1
Rscript $path/helpers/collate.R
