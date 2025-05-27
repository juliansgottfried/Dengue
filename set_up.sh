mkdir /scratch/$USER/Dengue/out
mkdir /scratch/$USER/Dengue/out/log
mkdir /scratch/$USER/Dengue/out/results
mkdir /scratch/$USER/Dengue/out/stats
mkdir /scratch/$USER/Dengue/out/output

touch .gitignore
true > .gitignore
echo .gitignore >> .gitignore
echo out >> .gitignore
echo folders_for_fit/* >> .gitignore
echo \!folders_for_fit/$USER >> .gitignore

module unload r
module unload gcc
module load r/gcc/4.3.1

Rscript /scratch/$USER/Dengue/helpers/init_renv.R
