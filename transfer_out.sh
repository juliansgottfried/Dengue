local_path=/Users/juliangottfried/Desktop/dengue_github/Dengue
store_path=/Users/juliangottfried/Desktop/NYU/thailand/cluster_runs

runs=($(ls $local_path/folders_for_fit))
for name in "${runs[@]}"; do
	cp $local_path/folders_for_fit/$name/results.csv \
		$local_path/folders_for_fit/$name/stats.csv \
		$store_path/$name	
done
