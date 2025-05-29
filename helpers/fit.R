library(tidyverse)
library(pomp)
library(ggridges)
library(doParallel)
library(doRNG)
library(ggtext)

log_name <- Sys.getenv("LOGNAME")
path_name<-paste0("/scratch/",log_name,"/Dengue")

args <- commandArgs(trailingOnly=TRUE)
n_array <- as.numeric(args[1])
array_id <- as.numeric(args[2])
n_cores <- as.numeric(args[3])
fit_name <- args[4]
path <- paste0(path_name,"/folders_for_fit/",fit_name,"/")

source(paste0(path_name,"/helpers/helper_function.R"))

seed <- 9087235
set.seed(seed)

source(paste0(path,"object.R"))
init_vals <- read_csv(paste0(path,"pars.csv"))

init_vals <- init_vals[1:min(500,nrow(init_vals)),par_names]

est_pars <- par_names[apply(init_vals,2,function(x) {max(x)-min(x)})!=0]
fixed_pars <- par_names[!par_names%in%est_pars]

bk <- read_csv(paste0(path,"bk_df.csv"),show_col_types=FALSE)

covariates <- covariate_table(
    bk %>% select(all_of(covars)),
    times = "time")
  
t_extrap <- with(bk, c(2 * time[1] - time[2], time))
covariates <- repair_lookup_table(covariates, t_extrap)

pomped <- pomp(
    data = bk %>% select(all_of(data_vars)) %>% na.omit,
    time = "time",
    t0 = with(bk, 2*time[1]-time[2]),
    covar = covariates,
    rprocess = euler(step.fun = rproc, delta.t = 1/365),
    rmeasure = rmeas,
    dmeasure = dmeas,
    rinit = rinit,
    paramnames = par_names,
    partrans = parameter_trans(
        log = log_transf,
        logit = logit_transf,
        barycentric = barycentric_transf
        ),
    accumvars = accum_names,
    statenames = c(accum_names,state_names),
    verbose = TRUE)

len <- nrow(init_vals)/n_array
init_vals <- init_vals[((array_id-1)*len+1):(array_id*len),]

non_ivp <- c(0.02,0.01,0.001)
ivp <- c(0.01,0.005,0.001)
lapply(1:3,function(order) {
    paste0(paste0("rdd",order),
           "<-rw_sd(",
           paste0(est_pars,ifelse(str_sub(est_pars,-2,-1)!="_0",
                                  paste0("=",non_ivp[order]),
                                  paste0("=ivp(",ivp[order],")")),
                  collapse=","),")") %>%
        str2expression %>%
        eval.parent
})

run_fitting(po=pomped,
            n_cores=n_cores,
            parameters=init_vals,
            seed_num=seed,
            rdd1=rdd1,rdd2=rdd2,rdd3=rdd3,
            result_path=paste0(path_name,"/out/results/",fit_name,"/",as.character(array_id),".csv"),
            log_path=paste0(path_name,"/out/log/",fit_name,"/",as.character(array_id),".txt"),
            stats_path=paste0(path_name,"/out/stats/",fit_name,"/",as.character(array_id),".csv"))
