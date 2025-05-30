mkdir tmp
ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit/ tmp/

paths=(tmp/*)
for path in "${paths[@]}"; do
	name=$(echo $path | awk -F/ '{print $NF}')
        mv path/results.csv path/stats.csv $name
done

rm -r tmp

ssh-keygen -R dtn.hpc.nyu.edu
