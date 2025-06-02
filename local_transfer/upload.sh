touch tmp.txt

flagger="n"
while [[ "$flagger" != "y" ]] ; do
	rm tmp.txt
	touch tmp.txt

	read -p "Date of run (format \"m_d\"): " run_date

	echo ""
	echo "Below, enter run IDs (lowercase letters)"
	echo ""
	echo "Enter \"0\" once all run IDs are in"
	echo ""
	i=1
	read -p "ID #$i: " run_id
	while [[ "$run_id" != "0" ]] ; do
		echo run_$run_date\_$run_id >> tmp.txt
		((i++))
		read -p "ID #$i: " run_id
	done

	echo ""
        echo "Printing folder names:"
        echo ""
	cat tmp.txt
	echo ""
	read -p "Are these folder names correct? (y/n): " flagger
	echo ""
done

names=($(cat tmp.txt))

rm tmp.txt

mkdir tmp

for name in "${names[@]}"; do
	cp -r $name tmp
done

ssh-keygen -R dtn.hpc.nyu.edu
scp -r -o "StrictHostKeyChecking no" tmp/* jg8461@dtn.hpc.nyu.edu:/scratch/jg8461/Dengue/folders_for_fit/

rm -r tmp
