mkdir /scratch/$USER/out
mkdir /scratch/$USER/out/log
mkdir /scratch/$USER/out/results
mkdir /scratch/$USER/out/stats
mkdir /scratch/$USER/out/output
mkdir /scratch/$USER/results
mkdir /scratch/$USER/folders_for_fit

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript /scratch/$USER/helpers/init_renv.R
