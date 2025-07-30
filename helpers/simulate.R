oldw <- getOption("warn")
options(warn = -1)

suppressMessages(library(tidyverse))
suppressMessages(library(ggdist))
suppressMessages(library(ggnewscale))
suppressMessages(library(panelPomp))

args      <- commandArgs(trailingOnly=TRUE)
isPanel   <- args[1]=="y"

log_name  <- Sys.getenv("LOGNAME")
path_name <- paste0("/scratch/", log_name, "/Dengue")

source(paste0(path_name, "/helpers/helper_functions.R"))
source(paste0(path_name, "/helpers/plot_functions.R"))

simulate_and_forecast <- function(fitting_path, save_sims=TRUE, save_filtered=TRUE)
{
    print("(saving simulation plot)")
    if (isPanel) {
        ##### this needs to be corrected as the new simulate.R does not pull mle before. #####
    	make_panel_plot(fitting_folder_path, mle, T)

    } else{
    	if (!file.exists(paste0(fitting_folder_path, "sim_cases.csv"))){
    		sim_df_list  <- simulate_mle(fitting_folder_path, save_sims=save_sims, save_filtered=save_filtered)
    		sim_cases_df <- sim_df_list$sim_cases
    	}else{
    		sim_cases_df <- read_csv(paste0(fitting_folder_path, "sim_cases.csv"), show_col_types=F)
    	}
    	if (file.exists(paste0(fitting_folder_path, "dataset_test.csv")) & !file.exists(paste0(fitting_folder_path, "forecast.csv"))){
    		fcast_df <- forecast_mle(fitting_folder_path)
    	}
    	if (!file.exists(paste0(fitting_folder_path, "plot.png"))){
    		make_plot(sim_cases_df, fitting_folder_path, F)
    	}
    }
}

names <- list.files(paste0(path_name, "/folders_for_fit"), full.names=F)
names <- names[names!= "summary.csv"]

for (name in names)
{
    fitting_folder_path <- paste0(path_name, "/folders_for_fit/", name, "/")
    simulate_and_forecast(fitting_folder_path)
}


