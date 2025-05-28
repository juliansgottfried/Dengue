TO SET UP ENVIRONMENT

On the cluster:

Place the repo into your scratch home directory so that the path to the directory is "scratch/$USER/Dengue/". Renv should activate itself automatically. Run "bash set_up.sh". You will be asked for your email, and you may be asked to grant permission for some action.

Next, on your local machine:

Make a directory to store fitting files and results. Place the files "download.sh" and "upload.sh" inside this directory.

TO RUN

1. Within the directory on your computer, create a folder for each fitting run. Name each folder "run_[month]\_[day]\_[x]", where "[x]" is a letter corresponding to the run. Inside each folder, place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv" (this last file should have exactly 500 rows).

2. Run "bash upload.sh". You will be prompted for your cluster password once for each run folder (sorry about this).

4. On the cluster, pull the repo. In the file "run_times.txt", write time limits in minutes for each run, each on a new line, in alphabetical order by folder name. 90 minutes is a good initial limit.

5. Run "bash launch.sh".

6. Once every run is finished, run "bash collate.sh". Push the repo.

7. On your local machine, pull the repo. Run "bash transfer_out.sh". Each run folder in your secondary directory should now contain "results.csv" and "stats.csv". Verify that these files are complete before proceeding with the last step.

8. On the cluster, run "bash clean.sh" and press "ENTER" to confirm. Push the repo and pull on your local machine to synchronize.
