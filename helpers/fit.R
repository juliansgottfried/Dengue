library(tidyverse)
library(pomp)
library(doParallel)
library(doRNG)

seed <- 9087235
set.seed(seed)

args <- commandArgs(trailingOnly=TRUE)

n_array  <- as.numeric(args[1])
array_id <- as.numeric(args[2])
n_cores  <- as.numeric(args[3])
fit_name <- args[4]
n_refine <- as.numeric(args[5])
nseq     <- as.numeric(args[6])

log_name  <- Sys.getenv("LOGNAME")
repo_path <- paste0("/scratch/", log_name, "/Dengue")
path      <- paste0(repo_path, "/folders_for_fit/", fit_name, "/")

source(paste0(repo_path, "/helpers/helper_functions.R"))
source(paste0(path, "object.R"))

df        <- read_csv(paste0(path, "dataset.csv"))
po        <- construct_pomp(path, df)
pars_path <- paste0(path, "pars.csv")

if (!file.exists(pars_path)) {
	param_bounds <- param_bounds %>% as.data.frame
	init_vals <- sobol_design(
        	lower=unlist(param_bounds[1,]),
        	upper=unlist(param_bounds[2,]),
        	nseq=nseq)
	write_csv(init_vals,paste0(path,"pars.csv"))	
} else init_vals <- read_csv(pars_path)

est_pars  <- par_names[apply(init_vals,2,\(.) diff(range(.)))!=0]
len       <- nseq/n_array
init_vals <- init_vals[((array_id-1)*len+1):(array_id*len),]

non_ivp   <- c(0.02, 0.01,  0.001)
ivp       <- c(0.01, 0.005, 0.001)

lapply(1:3,function(order) {
    paste0(paste0("rdd",order),
           "<-rw_sd(",
           paste0(est_pars,ifelse(str_sub(est_pars,-2,-1)!="_0",
                                  paste0("=",non_ivp[order]),
                                  paste0("=ivp(",ivp[order],")")),
                  collapse=","),")") %>%
        str2expression %>%
        eval.parent })

run_fitting(po         = po,
            n_cores    = n_cores,
            parameters = init_vals,
            seed_num   = seed,
            rdd1=rdd1,rdd2=rdd2,rdd3=rdd3,
            n_refine    = n_refine,
            Np1=1000, Np2=2000, Nmif=50,
            result_path = paste0(repo_path, "/out/results/", fit_name, "/", as.character(array_id), ".csv"),
            log_path    = paste0(repo_path, "/out/log/",     fit_name, "/", as.character(array_id), ".txt"),
	        traces_path = paste0(repo_path, "/out/traces/",  fit_name, "/", as.character(array_id), ".csv"),
            stats_path  = paste0(repo_path, "/out/stats/",   fit_name, "/", as.character(array_id), ".csv"))
 