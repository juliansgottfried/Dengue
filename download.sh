names=($(cat run_names.txt))

mkdir tmp

ssh-keygen -R dtn.hpc.nyu.edu
for name in "${names[@]}"; do
	scp -r -o "StrictHostKeyChecking no" jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit/$name/ tmp/
	mv tmp/results.csv tmp/stats.csv $name
	rm tmp/*
done

rm -r tmp
