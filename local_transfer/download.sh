source user_name.sh

mkdir tmp
ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" $user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/folders_for_fit/* tmp/

paths=(tmp/*)
for path in "${paths[@]}"; do
	name=$(echo $path | awk -F/ '{print $NF}')
        mv $path/results.csv $path/stats.csv $name
done

rm -r tmp
