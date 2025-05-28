TO SET UP ENVIRONMENT

On the cluster:

In your scratch home directory, run "git clone https://github.com/juliansgottfried/Dengue.git". Run "bash Dengue/set_up.sh". You will be asked for your email, and you may be prompted by renv.

On your computer:

Make a directory to store fitting files. Download the files "download.sh" and "upload.sh" and place them inside the directory.

TO RUN

1. Within the directory on your computer, create a folder for each fitting run. Name each folder "run_month\_day\_x", where "x" is a lowercase letter corresponding to the run.

2. Inside each folder, place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv" (this last file should have exactly 500 rows).

3. From within the directory on your computer, run "bash upload.sh". Answer the prompts. You will need your cluster password.

4. On the cluster, run "bash launch.sh". You will be asked to provide time limits for each run.

5. Once all the runs are finished, run "bash collate.sh".

6. From within the directory on your computer, run "bash download.sh". You will need your cluster password. Each run folder in your directory should now contain "results.csv" and "stats.csv". Verify that these files are complete before proceeding with the final step.

7. On the cluster, run "bash clean.sh" and press "ENTER" to confirm.
