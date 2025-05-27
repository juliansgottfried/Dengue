source path_vars.sh

if [ ! -d $local_path/folders_for_fit ]; then
        mkdir $local_path/folders_for_fit
fi
if [ ! -d $local_path/folders_for_fit/$user_name ]; then
        mkdir $local_path/folders_for_fit/$user_name
fi

names=($(cat $local_path/run_names.txt))

for name in "${names[@]}"; do
	cp -r $store_path/$name $local_path/folders_for_fit/$user_name
done
