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

**On your computer:**

1. Create a folder for each run, placing them within `fitting_folders`. Name the folders `run_[month]_[day]_[x]`, where `[x]` is a lowercase letter identifying the run. Sequences of runs should have contiguous identifying letters.

2. Inside each folder, place a data file `dataset.csv` and an object file `object.R`. See the directory `example_files`.

3. Run `bash upload.sh`. Answer the prompts. You will need your HPC password.

**On the cluster:**

4. Run `bash launch.sh`. Answer the prompts. You will be emailed when each run begins and ends.

5. Once all the runs are finished, run `bash collate.sh`.

**On your computer:**

6. Run `bash download.sh`. You will need your HPC password. Each folder within `fitting_folders` should now contain `results.csv`, `stats.csv`, `traces.csv`, `pars.csv`, and `plot.png`. Additionally, your main directory should now contain the updated file `summary.csv`.

7. Verify that these files are complete before proceeding with the final step.

**On the cluster:**

8. Run `bash clean.sh` and press `ENTER` to confirm.
