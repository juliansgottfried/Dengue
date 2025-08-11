path=/scratch/$LOGNAME/Dengue 

source $path/helpers/user_info.sh

bash $path/helpers/get_times.sh

runs=($(ls $path/folders_for_fit))

touch $path/fit_info.csv
for i in "${!runs[@]}"; do
    name=${runs[i]}
    source $path/folders_for_fit/$name/hyperparams.sh
	echo $name,$fitType >> $path/fit_info.csv
done

module unload r
module unload gcc
module load r/gcc/4.3.1
Rscript $path/helpers/collate.R

rm $path/times.txt $path/fit_info.csv