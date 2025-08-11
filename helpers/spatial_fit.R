library(spatPomp)
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

#n_array <- 1
#array_id <- 1
#n_cores <- detectCores()
#fit_name <- "run_8_6_a"
#n_refine <- 0
#nseq <- n_cores

log_name  <- Sys.getenv("LOGNAME")
repo_path <- paste0("/scratch/",log_name,"/Dengue")
path      <- paste0(repo_path,"/folders_for_fit/",fit_name,"/")

#repo_path <- "~/Desktop/Dengue"
#path <- paste0(repo_path,"/folders_for_fit/",fit_name,"/")

source(paste0(repo_path,"/helpers/cluster_helpers.R"))
source(paste0(path,"object.R"))

spo <- construct_spatpomp(path)

pars_path <- paste0(path,"pars.csv")

if (!file.exists(pars_path)) {
    param_df <- param_bounds %>% as.data.frame
    init_vals <- sobol_design(
        lower=unlist(param_df[1,]),
        upper=unlist(param_df[2,]),
        nseq=nseq)
    shared_init <- init_vals[,paste0(shared_pars,"1")]
    shared_init <- map2(1:ncol(shared_init),colnames(shared_init),\(x,y) {
        duplicated <- suppressMessages(bind_cols(replicate(U,shared_init[,x])))
        names(duplicated) <- par_list[[str_sub(y,1,-2)]]
        duplicated
        }) %>% bind_cols()
    specific_init <- init_vals[,unname(unlist(par_list[c(specific_pars)]))]
    init_vals <- bind_cols(specific_init,shared_init)
    write_csv(init_vals,paste0(path,"pars.csv"))	
} else init_vals <- read_csv(pars_path)

len <- nseq/n_array
init_vals <- init_vals[((array_id-1)*len+1):(array_id*len),]

non_ivp <- c(0.02,0.01,0.001)
ivp     <- c(0.01,0.005,0.001)

est_pars <- colnames(init_vals)[apply(init_vals,2,\(.)diff(range(.)))!=0]
est_ivps <- est_pars[unlist(lapply(str_split(est_pars,"_"),\(.) .[1]))%in%str_remove(ivp_pars,"_")]
lapply(1:3,function(i) {
    paste0(paste0("rdd",i),
           "<-rw_sd(",
           paste0(est_pars,ifelse(est_pars%in%est_ivps,
                                  paste0("=ivp(",ivp[i],")"),
                                  paste0("=",non_ivp[i])),
                  collapse=","),")") %>%
        str2expression %>%
        eval.parent
})

est_pars_basic <- paste0(unique(unlist(lapply(str_split(est_pars,"_"),\(.) .[1]))),"_")

est_shared <- shared_pars[shared_pars%in%est_pars_basic]
est_specific <- specific_pars[specific_pars%in%est_pars_basic]

run_spatial_fitting(po=spo,
            n_cores=n_cores,
            parameters=init_vals,
            unitParNames=est_specific,
            sharedParNames=est_shared,
            seed_num=seed,
            rdd1=rdd1,rdd2=rdd2,rdd3=rdd3,
	        n_refine=n_refine,
            Np1=1000,Np2=2000,Nbpf=50,
            block_size=2,
            resultw_path=paste0(repo_path,"/out/results/",fit_name,"/",as.character(array_id),".csv"),
            resultl_path=paste0(repo_path,"/out/results_long/",fit_name,"/",as.character(array_id),".csv"),
            log_path=paste0(repo_path,"/out/log/",fit_name,"/",as.character(array_id),".txt"),
            traces_path=paste0(repo_path,"/out/traces/",fit_name,"/",as.character(array_id),".csv"),
            stats_path=paste0(repo_path,"/out/stats/",fit_name,"/",as.character(array_id),".csv"))
