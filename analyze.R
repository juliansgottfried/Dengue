library(tidyverse)

# Parameter names to plot in the trace plots. If NULL, only loglik plotted.
plot_pars = NULL

# Ordering of spatial units when plotting from a panel fit.
ordering = c("Northern","Northeastern","Central","Eastern","Southern")

source("helpers/local_helpers.R")

fit_names <- unlist(read_csv("new_fits.csv",col_names=F))

for (name in fit_names) {
    path <- paste0("folders_for_fit/",name,"/")
    print("Making trace plot")
    plot_trace(path,plot_pars ,F)
    print("Making stats plot")
    plot_stats(path,F)
    
    suppressWarnings(rm(loc_key))
    source(paste0(path,"object.R"))
    isPanel=exists("loc_key")
    
    if (isPanel) {
        print("Simulating")
        sim_cases = simulate_panel_mle(path,T,T)[["sim_cases"]]
        print("Making simulation plot")
        panel_sim_plot(sim_cases,path,NULL,ordering)
    } else {
        print("Simulating")
        sim_cases = simulate_mle(path,T,T)[["sim_cases"]]
        print("Making simulation plot")
        sim_plot(sim_cases,path,NULL)
    }
}
