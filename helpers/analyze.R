library(tidyverse)

source("helpers/local_helpers.R")

tmp <- read_csv("tmp_summary.csv")
if (file.exists("summary.csv")) {
	main <- read_csv("summary.csv")
	main <- rbind(main[,2:ncol(main)],tmp)
} else {
	main <- tmp
}
write_csv(main,"summary.csv")

fit_names <- unlist(read_csv("fit_names.csv",col_names=F))

for (name in fit_names) {
    #name <- "run_7_17_j"
    path <- paste0("folders_for_fit/",name,"/")
    print("Making trace plot")
    plot_trace(path,NULL)
    print("Making stats plot")
    plot_stats(path)
    
    suppressWarnings(rm(loc_key))
    source(paste0(path,"object.R"))
    isPanel=exists("loc_key")
    
    if (isPanel) {
        print("Simulating")
        sim_cases = simulate_panel_mle(path,T,T)[["sim_cases"]]
        print("Making simulation plot")
        panel_sim_plot(sim_cases,path,NULL,c("Northern","Northeastern","Central","Eastern","Southern"))
    } else {
        print("Simulating")
        sim_cases = simulate_mle(path,T,T)[["sim_cases"]]
        print("Making simulation plot")
        sim_plot(sim_cases,path,NULL)
    }
}
