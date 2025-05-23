source path_vars.sh

if [ ! -d $local_path/folders_for_fit ]; then
        mkdir $local_path/folders_for_fit
fi

args=("$@")
for arg in "${args[@]}"; do
	cp $store_path/$arg $local_path/folders_for_fit
done
