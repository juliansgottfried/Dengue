TO SET UP ENVIRONMENTS

On the cluster:

1. In your scratch home directory run "git clone https://github.com/juliansgottfried/Dengue.git". Change directory to the repo (all subsequent cluster commands will be made from within the repo).

2. Run "bash set_up.sh". You will be prompted for your email.

Next, on your computer:

1. Make a directory to store fitting folders.

2. Download the files inside the repo folder "local_transfer" and place them inside your new directory.

3. Edit the "user_name" variable in the file "username.sh" to be your HPC username.

TO RUN

1. Within the directory on your computer, create a folder for each fitting run. Name each folder "run_month\_day\_x", where "x" is a lowercase letter identifying each run. Sequences of runs should have contiguous identifying letters.

2. Inside each folder, place a data file "dataset.csv" and an object file "object.R". See the files within "example_files" for minimal examples.

3. Run "bash upload.sh". Answer the prompts. You will need your cluster password.

4. On the cluster, run "bash launch.sh". Answer the prompts. You will be emailed when each run begins and ends.

5. Once all the runs are finished, run "bash collate.sh" on the cluster.

6. On your computer, run "bash download.sh". You will need your cluster password.

7. Each run folder in your directory should now contain "results.csv", "stats.csv", "traces.csv", "pars.csv", and "plot.png". Additionally, your directory should now contain the file "summary.csv". Verify that these files are complete before proceeding with the final step.

8. On the cluster, run "bash clean.sh" and press "ENTER" to confirm.
