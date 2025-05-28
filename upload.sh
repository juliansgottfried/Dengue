touch tmp.txt

flagger="n"
while [[ "$flagger" != "y" ]] ; do
	rm tmp.txt
	touch tmp.txt

	read -p "Date of run (format \"m_d\"): " run_date

	echo ""
	echo "Below, enter run IDs (lowercase letters)"
	echo ""
	echo "Enter \"0\" once all run IDs are entered"
	echo ""
	i=1
	read -p "Enter ID #$i: " run_id
	while [[ "$run_id" != "0" ]] ; do
		echo run_$run_date\_$run_id >> tmp.txt
		((i++))
		read -p "Enter ID #$i: " run_id
	done

	cat tmp.txt
	read -p "Are these fit names correct? (y/n): " flagger
done

names=($(cat tmp.txt))

for name in "${names[@]}"; do
	ssh-keygen -R dtn.hpc.nyu.edu
        scp -r -o "StrictHostKeyChecking no" $name jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit/
done

rm tmp.txt
