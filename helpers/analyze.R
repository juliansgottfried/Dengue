library(tidyverse)

source("helpers/local_helpers.R")

fit_names <- unlist(read_csv("new_fits.csv",col_names=F))

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
