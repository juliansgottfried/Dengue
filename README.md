TO SET UP ENVIRONMENT:

Run "bash make_folders.sh".


TO RUN:

Set-up:

Create folder(s) inside "folders_for_fit" called "run_[month]\_[day]\_[x]" where "[x]" is a lowercase letter corresponding to the run.
Inside each folder place the data file "bk_df.csv", the object "object.R", and the parameter combinations "pars.csv".
The file "pars.csv" should have exactly 500 rows.
The folder "current_files" contains samples of the above files.


Run:

Run command "bash launch.sh a b c ..." where "a b c ..." are time limits in minutes corresponding to each folder in "folders_for_fit", in order. 90 minutes is a good guess.


Gather results:

Run "bash collate.sh". This places "results.csv" and "stats.csv" into "results/[folder]" where "[folder]" corresponds to each folder in "folders_for_fit".


Transfer results:

Use Globus or Git or ssh to transfer "results.csv" and "stats.csv" to your local computer.


Re-set directory:

Run "bash clean.sh" and press ENTER to confirm.
