## To set up environments

### On the cluster

Run the following lines:
```
cd /scratch/$LOGNAME
git clone https://github.com/juliansgottfried/Dengue.git
cd Dengue
bash set_up.sh
```
You will be prompted for your email. All subsequent cluster commands should be made from within the Dengue repo.

### On your computer

Make a new directory. From within the directory, run the following line:
```
bash <(curl -s https://raw.githubusercontent.com/juliansgottfried/Dengue/refs/heads/main/local_transfer/helpers/set_up.sh)
```
You will be prompted for your HPC username. All subsequent computer commands should be made from within this directory.

## To run fitting

1. Within the folder "fitting_folders" in your computer directory, create a folder for each fitting run. Name each folder "run_month\_day\_x", where "x" is a lowercase letter identifying each run. Sequences of runs should have contiguous identifying letters.

2. Inside each folder, place a data file "dataset.csv" and an object file "object.R". See the files within "example_files" for minimal examples.

3. Run 'bash upload.sh'. Answer the prompts. You will need your cluster password.

4. On the cluster, run 'bash launch.sh'. Answer the prompts. You will be emailed when each run begins and ends.

5. Once all the runs are finished, run 'bash collate.sh' on the cluster.

6. On your computer, run 'bash download.sh'. You will need your cluster password.

7. Each fitting folder in your directory should now contain "results.csv", "stats.csv", "traces.csv", "pars.csv", and "plot.png". Additionally, your main directory should now contain the file "summary.csv". Verify that these files are complete before proceeding with the final step.

8. On the cluster, run 'bash clean.sh' and press 'ENTER' to confirm.
