TO SET UP ENVIRONMENT

On the cluster:

Place this repo into your scratch home directory so that the path to the directory is "scratch/$USER/Dengue/". Renv should activate itself automatically. Run "bash set_up.sh" to create folder structure and complete renv activation; you may be asked to grant permission for some action.

On your local machine:

Place this repo somewhere on your computer. Make a secondary directory somewhere else to store fitting files and results. In the file "path_vars.sh", edit the "local_path" and "store_path" variables to be the path names to the repo and the secondary directory, respectively.

TO RUN

1. Within the secondary directory on your computer, create a folder for each fitting run. Name each folder "run_[month]\_[day]\_[x]", where "[x]" is a lowercase letter corresponding to the run. Inside each folder, place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv" ("pars.csv" should have exactly 500 rows).

2. In the local repo, run "bash transfer_in.sh run_1 run_2 run_3 ..." where "run_i" is the name of the ith run folder. Push the repo.

3. In the cluster, pull the repo.

4. Run "bash launch.sh t_1 t_2 t_3 ..." where "t_i" is the time limit in minutes for the ith run folder, in lexicographic order by folder name. 90 minutes is a good first guess!

5. Once every run is finished, run "bash collate.sh". Push the repo.

6. On your local machine, pull the repo. Run "bash transfer_out.sh". Each run folder in your secondary directory should now contain "results.csv" and "stats.csv". Verify that these files are sound.

7. On the cluster, run "bash clean.sh" and press ENTER to confirm. Push the repo and pull on your local machine.