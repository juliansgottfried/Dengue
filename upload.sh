names=($(cat run_names.txt))

ssh-keygen -R dtn.hpc.nyu.edu
for name in "${names[@]}"; do
	scp -r -o "StrictHostKeyChecking no" $name jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit/
done
