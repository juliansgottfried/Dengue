path=/scratch/$LOGNAME/Dengue 

runs=($(ls $path/folders_for_fit))

touch $path/times.txt

for i in "${!runs[@]}"; do
    cat $path/times.txt > $path/tmp1.txt
    sacct -X --start now-1weeks --name=${runs[i]} --format=Elapsed,JobName | \
        grep -Eo "([0-9]{2}:{0,1}){3}" | \
        awk -v var="${runs[i]}" '{$2 = var; print}' > \
        $path/tmp2.txt
    cat $path/tmp*.txt > $path/times.txt
done

rm $path/tmp*.txt
