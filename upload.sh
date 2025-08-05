source helpers/username.sh

touch tmp.txt

read -p $'Upload sequence of runs or pick specific runs:\n(seq or spec) ' upload_config
echo ""

if [ $upload_config == "seq" ]; then
	flagger="n"
	while [[ "$flagger" != "y" ]] ; do
		rm tmp.txt
		touch tmp.txt

		read -p "Date of run (format \"m_d\"): " run_date
		echo ""

		read -p "Enter first run ID (lowercase letter): " id_start
		echo ""
		read -p "Enter last run ID (lowercase letter): " id_end
		echo ""

		for run_id in $(eval echo "{$id_start..$id_end}"); do
			echo run_$run_date\_$run_id >> tmp.txt
		done
		echo ""

		echo "Printing folder names:"
		echo ""

		cat tmp.txt
		echo ""

		read -p "Are these folder names correct? (y/n): " flagger
		echo ""
	done 
else
	flagger="n"
	while [[ "$flagger" != "y" ]] ; do
		rm tmp.txt
		touch tmp.txt

		read -p "Date of run (format \"m_d\"): " run_date
		echo ""

		echo "Enter run IDs one at a time. Enter \"done\" when finished."
		echo ""

		run_id="NA"
		while [[ "$run_id" != "done" ]] ; do
			read -p "Enter ID (lowercase letter): " run_id
			echo run_$run_date\_$run_id >> tmp.txt
			echo ""
		done

		echo "Printing folder names:"
		echo ""

		cat tmp.txt
		echo ""

		read -p "Are these folder names correct? (y/n): " flagger
		echo ""
	done 
fi

names=($(cat tmp.txt))

rm tmp.txt

mkdir tmp

for name in "${names[@]}"; do
	cp -r folders_for_fit/$name tmp
done

ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" tmp/* $user_name@dtn.hpc.nyu.edu:/scratch/$user_name/Dengue/folders_for_fit/

rm -r tmp
