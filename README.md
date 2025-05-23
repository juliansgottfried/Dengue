TO SET UP ENVIRONMENT

On the cluster:

Place the repo into your scratch home directory so that the path to the directory is "scratch/$USER/Dengue/". Renv should activate itself automatically. Run "bash set_up.sh" to create folder structure and complete renv activation; you may be asked to grant permission for some action.

Next, on your local machine:

Place the repo somewhere on your computer. Make a secondary directory somewhere else to store fitting files and results. In the file "path_vars.sh", edit the "local_path" and "store_path" variables to be the path names to the repo and the secondary directory, respectively.

TO RUN

1. Within the secondary directory on your computer, create a folder for each fitting run. Name each folder "run_[month]\_[day]\_[x]", where "[x]" is a letter corresponding to the run. Inside each folder, place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv" (this last file should have exactly 500 rows).

2. In your local repo, paste the names of the newly-created folders in the file "run_names.txt", each name on a new line.

3. Run "bash transfer_in.sh". Push the repo.

4. On the cluster, pull the repo. In the file "run_times.txt", write time limits in minutes for each run, each on a new line, in alphabetical order by folder name. 90 minutes is a good initial limit.

5. Run "bash launch.sh".

6. Once every run is finished, run "bash collate.sh". Push the repo.

7. On your local machine, pull the repo. Run "bash transfer_out.sh". Each run folder in your secondary directory should now contain "results.csv" and "stats.csv". Verify that these files are complete before proceeding with the last step.

8. On the cluster, run "bash clean.sh" and press "ENTER" to confirm. Push the repo and pull on your local machine to synchronize.
