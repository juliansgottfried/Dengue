source path_vars.sh

if [ ! -d $local_path/folders_for_fit ]; then
        mkdir $local_path/folders_for_fit
fi

IFS=$'\n' a=($(cat $local_path/run_names.txt))
for i in $(seq ${#a[*]}); do
    [[ ${a[$i-1]} = $name ]] && echo "${a[$i]}"
done

for run in "${a[@]}"; do
        cp -r $store_path/$run $local_path/folders_for_fit
done
