source helpers/username.sh

mkdir tmp
ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" \
	$user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/folders_for_fit/* \
	tmp/

mv tmp/summary.csv tmp_summary.csv

paths=(tmp/*)
for path in "${paths[@]}"; do
	name=$(echo $path | awk -F/ '{print $NF}')
        mv $path/* fitting_folders/$name
done

rm -r tmp

Rscript helpers/bind.R
rm tmp_summary.csv
