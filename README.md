TO SET UP ENVIRONMENT:

Working in the cluster, place the repo into your scratch home directory "/scratch/$USER/". Renv should activate automatically.

Run "bash set_up.sh". Renv may ask for permission for some action.

TO RUN:

Set-up:

Working in your local machine, create a folder for each fit. Name the folders "run_[month]\_[day]\_[x]", where "[x]" is a lowercase letter corresponding to the run.
Inside each folder, place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv".
Place the folders inside "folders_for_fit" in your local copy of the repo, and transfer to the cluster repo.

Note -- the file "pars.csv" should have exactly 500 rows.

Run:

In the cluster, run "bash launch.sh t\_1 t\_2 t\_3 ..." where "t\_i" is the time limit in minutes corresponding to the ith folder in "folders_for_fit", in lexicographic order. 90 minutes is a good first guess.

Save results:

In the cluster, run "bash collate.sh". This places "results.csv" and "stats.csv" into the corresponding folder under "folders_for_fit".
Transfer each folder back to your local copy of the repo. Copy the results elsewhere in your machine for further work, because the repo will be re-set.

Re-set directory:

In the cluster, run "bash clean.sh" and press ENTER to confirm. This removes all input and output files.
