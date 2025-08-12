library(spatPomp)
library(tidyverse)
library(doParallel)
library(doRNG)

seed <- 9087235
set.seed(seed)

#args <- commandArgs(trailingOnly=TRUE)

#n_array <- as.numeric(args[1])
#array_id <- as.numeric(args[2])
#n_cores <- as.numeric(args[3])
#fit_name <- args[4]
#n_refine <- as.numeric(args[5])
#nseq <- as.numeric(args[6])

n_array <- 1
array_id <- 1
n_cores <- detectCores()
fit_name <- "run_8_6_a"
n_refine <- 1
nseq <- n_cores

#log_name  <- Sys.getenv("LOGNAME")
#repo_path <- paste0("/scratch/",log_name,"/Dengue")
#path      <- paste0(repo_path,"/folders_for_fit/",fit_name,"/")

repo_path <- "~/Desktop/Dengue"
path <- paste0(repo_path,"/folders_for_fit/",fit_name,"/")

source(paste0(repo_path,"/helpers/cluster_helpers.R"))

df <- read_csv(paste0(path,"dataset.csv"))
source(paste0(path,"object.R"))

U <- nrow(unique(df[,unit_name]))
par_list <- setNames(lapply(paste0(c(specific_pars,shared_pars),"_"),\(.) paste0(rep(.,U),1:U)),
                     c(specific_pars,shared_pars))
basic_par_list <- c(par_list[c(specific_pars)],
                    lapply(par_list[c(shared_pars)],\(.) .[1]))

get_transf <- function(pars) {unname(unlist(par_list[pars]))}
log_transf <- get_transf(log_transf)
logit_transf <- get_transf(logit_transf)
barycentric_transf <- get_transf(barycentric_transf)

par_trans <- parameter_trans(
    log = log_transf,
    logit = logit_transf,
    barycentric = barycentric_transf)

spo <- construct_spatpomp(path,df,
                          unname(unlist(par_list)),par_trans)

bounds_df <- unlist(lapply(names(param_bounds),\(.) {
    params <- basic_par_list[[.]]
    setNames(rep(list(param_bounds[[.]]),length(params)),
             basic_par_list[[.]])
}), recursive=F) |> as.data.frame()

pars_path <- paste0(path,"pars.csv")
if (!file.exists(pars_path)) {
    init_vals <- sobol_design(
        lower=unlist(bounds_df[1,]),
        upper=unlist(bounds_df[2,]),
        nseq=nseq)
    shared_init <- setNames(init_vals[,rep(paste0(shared_pars,"_1"),each=U)],
                            unname(unlist(par_list[shared_pars])))
    specific_init <- init_vals[,unname(unlist(par_list[c(specific_pars)]))]
    init_vals <- bind_cols(specific_init,shared_init)
    write_csv(init_vals,paste0(path,"pars.csv"))	
} else init_vals <- read_csv(pars_path)

len <- nseq/n_array
init_vals <- init_vals[((array_id-1)*len+1):(array_id*len),]

non_ivp <- c(0.005,0.00125,0.00125)
ivp     <- c(0.005,0.00125,0.00125)

est_pars <- names(bounds_df)[bounds_df[1,]-bounds_df[2,]!=0]
est_pars_cut <- unlist(lapply(str_split(est_pars,"_"),\(.) .[1]))
est_ivps <- est_pars[est_pars_cut%in%ivp_pars]
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

est_shared <- shared_pars[shared_pars%in%est_pars_cut]
est_specific <- specific_pars[specific_pars%in%est_pars_cut]

run_spatial_fitting(po=spo,
            n_cores=n_cores,
            parameters=init_vals,
            unitParNames=paste0(est_specific,"_"),
            sharedParNames=paste0(est_shared,"_"),
            par_names=c(specific_pars,shared_pars),
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
