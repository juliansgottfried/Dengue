## To set up directories

### On the cluster

Run the following lines:
```
cd /scratch/$LOGNAME
git clone https://github.com/juliansgottfried/Dengue.git
cd Dengue
bash cluster_set_up.sh
```
You will be prompted for your email. All subsequent cluster commands should be made from within the Dengue repo.

### On your computer

Change directory to where you want your local copy of the repo.
Run the following lines:
```
git clone https://github.com/juliansgottfried/Dengue.git
cd Dengue
bash local_set_up.sh
```
You will be prompted for your HPC username. All subsequent computer commands should be made from within the Dengue repo.

## To run fitting

**On your computer:**

1. Create a folder for each run, placing them within `folders_for_fit`. Name the folders `run_[month]_[day]_[x]`, where `[x]` is a lowercase letter identifying the run. Sequences of runs should have contiguous letters.

2. Inside each folder place a data file `dataset.csv` and an object file `object.R`. See `example_files` for specifications.

3. Run `bash upload.sh`. Answer the prompts. You will need your HPC password.

**On the cluster:**

4. Run `bash launch.sh`. Answer the prompts. You will be emailed when each run begins and ends.

5. Once all the runs are finished, run `bash collate.sh`.

**On your computer:**

6. Run `bash download.sh`. You will need your HPC password. Each folder within `folders_for_fit` should now contain `results.csv`, `stats.csv`, `traces.csv`, `pars.csv`, and, if a panel fit was performed, `results_long.csv`. Additionally, your repo should now contain the updated file `summary.csv`.

7. Open `analyze.R` and execute to generate simulation files and plots.

8. Verify that all files are complete before proceeding with the final step.

**On the cluster:**

9. Run `bash clean.sh` and press `Enter` to confirm.
