library(panelPomp)
library(tidyverse)
library(doParallel)
library(doRNG)

seed <- 9087235
set.seed(seed)

args <- commandArgs(trailingOnly=TRUE)

n_array <- as.numeric(args[1])
array_id <- as.numeric(args[2])
n_cores <- as.numeric(args[3])
fit_name <- args[4]
n_refine <- as.numeric(args[5])
nseq <- as.numeric(args[6])

log_name  <- Sys.getenv("LOGNAME")
repo_path <- paste0("/scratch/",log_name,"/Dengue")
path      <- paste0(repo_path,"/folders_for_fit/",fit_name,"/")

source(paste0(repo_path,"/debug/helper_functions.R"))
source(paste0(path,"object.R"))

ppo <- construct_panel_pomp(path,nseq)
init_vals <- read_csv(paste0(path,"pars.csv"),show_col_types=FALSE)

len <- nseq/n_array
init_vals <- init_vals[((array_id-1)*len+1):(array_id*len),]

non_ivp <- c(0.02,0.01,0.001)
ivp     <- c(0.01,0.005,0.001)

cut_names <- unlist(lapply(str_split(names(init_vals),"\\["),\(.).[1]))
est_pars <- unique(cut_names[apply(init_vals,2,\(.)diff(range(.)))!=0])
lapply(1:3,function(i) {
    paste0(paste0("rdd",i),
           "<-rw_sd(",
           paste0(est_pars,ifelse(str_sub(est_pars,-2,-1)!="_0",
                                  paste0("=",non_ivp[i]),
                                  paste0("=ivp(",ivp[i],")")),
                  collapse=","),")") %>%
        str2expression %>%
        eval.parent
})

run_panel_fitting(po=ppo,
            n_cores=n_cores,
            parameters=init_vals,
            seed_num=seed,
            rdd1=rdd1,rdd2=rdd2,rdd3=rdd3,
	    n_refine=n_refine,
            Np1=5,Np2=5,Nmif=5,
            resultw_path=paste0(repo_path,"/debug/resultw.csv"),
            resultl_path=paste0(repo_path,"/debug/resultl.csv"),
            log_path=paste0(repo_path,"/debug/keeplog.txt"),
	    traces_path=paste0(repo_path,"/debug/traces.csv"),
            stats_path=paste0(repo_path,"/debug/stats.csv"))
