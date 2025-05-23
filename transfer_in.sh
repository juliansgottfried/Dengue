local_path=/Users/juliangottfried/Desktop/dengue_github/Dengue
store_path=/Users/juliangottfried/Desktop/NYU/thailand/cluster_runs

if [ ! -d $local_path/folders_for_fit ]; then
	mkdir $local_path/folders_for_fit
fi

args=("$@")
for arg in "${args[@]}"; do
	cp $store_path/$arg $local_path/folders_for_fit
done
