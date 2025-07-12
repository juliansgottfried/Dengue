source username.sh

mkdir tmp
ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" \
	$user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/folders_for_fit/* \
	$user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/summary.csv \
	tmp/

mv tmp/summary.csv .

paths=(tmp/*)
for path in "${paths[@]}"; do
	name=$(echo $path | awk -F/ '{print $NF}')
        mv $path/results.csv $path/stats.csv $path/traces.csv $name
done

rm -r tmp
