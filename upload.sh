source path_vars.sh

#if [ ! -d $local_path/folders_for_fit ]; then
#        mkdir $local_path/folders_for_fit
#fi

names=($(cat $local_path/run_names.txt))

for name in "${names[@]}"; do
	ssh-keygen -R dtn.hpc.nyu.edu
	scp -o "StrictHostKeyChecking no" $store_path/$name jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit
#	cp -r $store_path/$name $local_path/folders_for_fit
done
