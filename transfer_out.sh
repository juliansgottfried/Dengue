source path_vars.sh

runs=($(ls $local_path/folders_for_fit/$user_name))
for name in "${runs[@]}"; do
	cp $local_path/folders_for_fit/$user_name/$name/results.csv \
		$local_path/folders_for_fit/$user_name/$name/stats.csv \
		$store_path/$name
done
