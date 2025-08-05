source helpers/user_info.sh

mkdir tmp
ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" \
	$user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/folders_for_fit/* \
	tmp/

mv tmp/summary.csv tmp_summary.csv

paths=(tmp/*)
touch fit_names.csv
for path in "${paths[@]}"; do
	name=$(echo $path | awk -F/ '{print $NF}')
	mv $path/* folders_for_fit/$name
	echo $name >> fit_names.csv
done

rm -r tmp

Rscript helpers/analyze.R
rm tmp_summary.csv
rm fit_names.csv
