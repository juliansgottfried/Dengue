source path_vars.sh

runs=($(ls $local_path/folders_for_fit))
for name in "${runs[@]}"; do
	cp $local_path/folders_for_fit/$name/results.csv \
		$local_path/folders_for_fit/$name/stats.csv \
		$store_path/$name	
done
