local_path=/Users/juliangottfried/Desktop/dengue_github/Dengue
save_path=/Users/juliangottfried/Desktop/NYU/thailand/cluster_runs

args=("$@")
for arg in "${args[@]}"; do
	cp $save_path/$arg $local_path/folders_for_fit
done
