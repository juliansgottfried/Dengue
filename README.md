TO SET UP ENVIRONMENT:

Place all files in repo into scratch home directory "/scratch/$USER/". Renv should activate automatically.

Run "bash set_up.sh". Renv may ask for permission for some action.

TO RUN:

Set-up:

Create folder(s) inside "folders_for_fit" called "run_[month]\_[day]\_[x]" where "[x]" is a lowercase letter corresponding to the run.
Inside each folder place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv".

The file "pars.csv" should have exactly 500 rows.

The folder "current_files" contains samples of the above files.

Run:

Run "bash launch.sh t1 t2 t3 ..." where "t\_i" is the time limit in minutes corresponding to the ith folder in "folders_for_fit", in order. 90 minutes is a good first guess.

Gather results:

Run "bash collate.sh". This places "results.csv" and "stats.csv" into "results/[run]" for each run folder in "folders_for_fit".

Transfer results:

Use Globus or some other method to transfer "results.csv" and "stats.csv" to your local machine.

Re-set directory:

Run "bash clean.sh" and press ENTER to confirm.
