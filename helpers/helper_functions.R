construct_pomp <- function(path, df) {

	source(paste0(path, "object.R"))
	covariates <- df %>% select(-all_of(obs_vars))
	covariates <- covariate_table(covariates, times = "time")
	t_extrap   <- with(df, c(2 * time[1] - time[2], time))
	covariates <- repair_lookup_table(covariates, t_extrap)

	po <- pomp(
		data       = df %>% select(time, all_of(obs_vars)) %>% na.omit,
		times      = "time",
		t0         = with(df, 2*time[1]-time[2]),
		covar      = covariates,
		rprocess   = euler(step.fun = rproc, delta.t = 1/365),
		rmeasure   = rmeas,
		dmeasure   = dmeas,
		rinit      = rinit,
		paramnames = par_names,
		partrans   = parameter_trans(
			log    = log_transf,
			logit  = logit_transf,
			barycentric = barycentric_transf),
		accumvars  = accum_names,
		statenames = c(accum_names,state_names),
		verbose = F)

	return(po)
}

construct_panel_pomp <- function(path, nseq) {
    source(paste0(path,"object.R"))

    df <- read_csv(paste0(path,"dataset.csv"),show_col_types=FALSE)
    df[,loc_key] <- str_remove_all(unlist(df[,loc_key])," ")

    locs    <- unname(unlist(df[,loc_key]))
    keys    <- unique(locs)
    tmp_ppo <- lapply(keys,function(x)
    {
        tmp_df = df[locs==x,] %>% select(-all_of(loc_key))
        po <- construct_pomp(path, tmp_df)
    })

    names(tmp_ppo) <- keys

    pars_path <- paste0(path,"pars.csv")
    if (!file.exists(pars_path)) {
        param_bounds <- param_bounds %>% as.data.frame
        init_vals <- runif_panel_design(
            lower=unlist(param_bounds[1,]),
            upper=unlist(param_bounds[2,]),
            nseq=nseq,
            specific_names=specific_names,
            unit_names=keys)
        write_csv(init_vals,pars_path)
    } else init_vals <- read_csv(pars_path,show_col_types=FALSE)
    
    ppo <- panelPomp(tmp_ppo,params=unlist(init_vals[1,]))
    return(ppo)
}

run_fitting <- function(
        po, n_cores, parameters,
        seed_num, rdd1, rdd2, rdd3,
	n_refine,
	Np1,Np2,Nmif,
        result_path,
        log_path,
	traces_path,
        stats_path) {
    
    cl <- parallel::makeCluster(n_cores)
    registerDoParallel(cl)
    registerDoRNG(seed = as.integer(round(abs(seed_num))) + 1234)
    
    ## for each parameter row, run mif
    r1 <- foreach::foreach(
        i = seq_len(nrow(parameters)),
        .packages = c("pomp", "dplyr", "readr")
    ) %dopar% {
        param <- as.numeric(parameters[i, ])
        names(param) <- colnames(parameters)

        cat(paste("Starting iteration", i, "\n"),
		file = log_path,
                append = TRUE)

	rdds = list(rdd1,rdd2,rdd3)

        mifout <- tryCatch(po |>
                               mif2(Np = Np1,
                                    Nmif = Nmif,
                                    cooling.type = "geometric",
                                    cooling.fraction.50 = 0.5,
                                    params = param,
                                    rw.sd = rdds[[1]]),
                           error = function(e) e)

	traces <- data.frame(mifout@traces,iter=0:Nmif,run=1)

	if (n_refine > 0) {
		for (j in 1:n_refine) {
			mifout <- mifout |> mif2(rw.sd=rdds[[j+1]])
			traces_j <- data.frame(mifout@traces,iter=0:Nmif,run=j+1)
			traces <- bind_rows(traces,traces_j)
		}
	}
	if (file.exists(traces_path)) {
            read_csv(traces_path) %>% bind_rows(traces) %>% write_csv(traces_path)
        } else write_csv(traces,traces_path)

        stats <- data.frame(cond=mifout@cond.logLik,
                            eff=mifout@eff.sample.size,
                            time=po@times)
        if (file.exists(stats_path)) {
            read_csv(stats_path) %>% bind_rows(stats) %>% write_csv(stats_path)
        } else write_csv(stats,stats_path)

        result <- c(rep(NA, length(param) + 4))
        names(result) <- c("sample",
                           colnames(parameters),
                           "loglik",
                           "loglik.se",
                           "flag")
        result[1] <- i
        if (length(coef(mifout)) > 0) {
            loglik_mif <- tryCatch(replicate(n = 10,
                                             logLik(pfilter(po,
                                                            params = coef(mifout),
                                                            Np = Np2))),
                                   error = function(e) e)

            if (is.numeric(loglik_mif)) {
                bl <- logmeanexp(loglik_mif, se = TRUE)
                loglik_mif_est <- bl[1]
                loglik_mif_se <- bl[2]
                cat(paste(i, loglik_mif_est, "\n"), file = log_path, append = TRUE)
                result[length(result)] <- 2
            }
            if (is.numeric(loglik_mif)) {
                par_out <- coef(mifout)
                result[(length(param) + 2):(length(param) + 3)] <- bl
                result[(2):(length(param) + 1)] <- par_out
                result[length(result)] <- 1
            }
        } else {
            result[length(result)] <- 3
            cat(paste(i, "failed", "\n"), file = log_path, append = TRUE)
        }
        result
    }
    r1 <- r1 |>
        bind_rows() |>
        remove_missing()
    write.table(r1,
                result_path,
                append = TRUE,
                col.names = !file.exists(result_path),
                row.names = FALSE, sep = ",")
    
    stopCluster(cl)
    return(r1)
}

run_panel_fitting <- function(
        po, n_cores, parameters,
        seed_num, rdd1, rdd2, rdd3,
        n_refine,
	Np1,Np2,Nmif,
        resultw_path,
        resultl_path,
        log_path,
        traces_path,
        stats_path) {
    
    cl <- parallel::makeCluster(n_cores)
    registerDoParallel(cl)
    registerDoRNG(seed = as.integer(round(abs(seed_num))) + 1234)
    
    ## for each parameter row, run mif
    rs <- foreach::foreach(
        i = seq_len(nrow(parameters)),
        .packages = c("panelPomp", "dplyr", "readr")
    ) %dopar% {
        param <- as.numeric(parameters[i, ])
        names(param) <- colnames(parameters)
        
        cat(paste("Starting iteration", i, "\n"),
            file = log_path,
            append = TRUE)
        
        rdds = list(rdd1,rdd2,rdd3)
        keys <- names(po@unit_objects)
        
        mifout <- tryCatch(po |>
                               mif2(Np = Np1,
                                    Nmif = Nmif,
                                    cooling.type = "geometric",
                                    cooling.fraction.50 = 0.5,
                                    start = param,
                                    rw.sd = rdds[[1]]),
                           error = function(e) e)
        
        get_trace <- function(i) {
            lapply(keys,\(.) {
                data.frame(mifout@unit_objects[[.]]@traces,iter=0:Nmif,run=i,unit=.)
                }) |> bind_rows()
        }

        traces <- get_trace(1)

	if (n_refine > 0) {
                for (j in 1:n_refine) {
                        mifout <- mifout |> mif2(rw.sd=rdds[[j+1]])
                        traces_j <- get_trace(j+1)
                        traces <- bind_rows(traces,traces_j)
                }
        }

        if (file.exists(traces_path)) {
            read_csv(traces_path) %>% bind_rows(traces) %>% write_csv(traces_path)
        } else write_csv(traces,traces_path)
        
        stats <- lapply(keys,\(.) {
            data.frame(cond=mifout@unit_objects[[.]]@cond.logLik,
                       eff=mifout@unit_objects[[.]]@eff.sample.size,
                       time=mifout@unit_objects[[.]]@times,
                       unit=.)
            }) |> bind_rows()

        if (file.exists(stats_path)) {
            read_csv(stats_path) %>% bind_rows(stats) %>% write_csv(stats_path)
        } else write_csv(stats,stats_path)
        
        resultw <- c(rep(NA, length(param) + 4))
        names(resultw) <- c("sample",colnames(parameters),
                           "loglik","loglik.se","flag")
        
        unique_pars <- c(names(shared(po)),rownames(specific(po)))
        resultl <- matrix(rep(NA, length(keys)*(length(unique_pars) + 5)),nrow=length(keys))
        colnames(resultl) <- c("sample","unit",unique_pars,
                               "loglik","loglik.se","flag")
        resultl <- data.frame(resultl)
        resultl$unit <- keys
        
        resultw["sample"] <- i
        resultl$sample <- i
        
        if (length(coef(mifout)) > 0) {
            loglik_mif <- tryCatch(replicate(n = 10,
                                             logLik(pfilter(po,
                                                            params = coef(mifout),
                                                            Np = Np2))),
                                   error = function(e) e)
            
            if (is.numeric(loglik_mif)) {
                bl <- logmeanexp(loglik_mif, se = TRUE)
                names(bl) <- c("loglik","loglik.se")
                loglik_mif_est <- bl[1]
                cat(paste(i, loglik_mif_est, "\n"), file = log_path, append = TRUE)
                
                resultw["flag"] <- 2
                resultl$flag <- 2
            }
            if (is.numeric(loglik_mif)) {
                resultw[names(bl)] <- bl
                resultl[,names(bl)] <- t(matrix(bl))[rep(1,length(keys)),]
                
                resultw["flag"] <- 1
                resultl$flag <- 1
                
                par_out <- coef(mifout)
                resultw[names(par_out)] <- par_out
                
                par_out <- specific(mifout) |>
                    t() |>
                    data.frame() |> 
                    bind_cols(data.frame(t(shared(mifout))))
                
                resultl[,colnames(par_out)] <- par_out
            }
        } else {
            resultw["flag"] <- 3
            resultl$flag <- 3
            cat(paste(i, "failed", "\n"), file = log_path, append = TRUE)
        }
        list(resultw,resultl)
    }
    rw <- lapply(rs,\(. ).[[1]]) |> 
        bind_rows() |>
        remove_missing()
    rl <- lapply(rs,\(. ).[[2]]) |> 
        bind_rows() |>
        remove_missing()
    
    write.table(rw, 
                resultw_path,
                append = TRUE,
                col.names = !file.exists(resultw_path),
                row.names = FALSE, sep = ",")
    write.table(rl, 
                resultl_path,
                append = TRUE,
                col.names = !file.exists(resultl_path),
                row.names = FALSE, sep = ",")

    stopCluster(cl)
}

simulate_mle <- function(path, save_sims=TRUE, save_filtered = TRUE) {
    df   <- read_csv(paste0(path, "dataset.csv"), show_col_types=FALSE)
    po   <- construct_pomp(path, df)
    mle  <- read_csv(paste0(path, "results.csv"), show_col_types=FALSE) %>% head(1)%>% select(all_of(par_names))

	coef(po, names(mle)) <- mle
	sims      <- po %>% simulate(nsim  = 1000,
				                include.data = TRUE,
				                format       = "data.frame")

    sim_states_df        <- sims %>% select(time, .id, all_of(c(state_names, accum_names)))
    sim_states_df$type   <- "sim_states"

	set_0        <- function(x) (ifelse(is.na(x),0,x))
	sim_cases_df <- sims %>% select(time,.id, all_of(obs_vars)) %>% mutate_at(obs_vars, set_0)

    if (save_filtered==T){

        filt_sims_states_df      <- pfilter(po, save.states="filter", Np=1000) %>% saved_states(format="data.frame")
        filt_sims_states_df      <- filt_sims_states_df %>% pivot_wider(names_from="name")
        filt_sims_states_df$type <- "filter_states"
        states_df                <- rbind(sim_states_df, filt_sims_states_df)

    } else {
        states_df <- sim_states_df
    }

    if (save_sims==T) {
        write_csv(sim_cases_df, paste0(path, "sim_cases.csv"))
        write_csv(states_df, paste0(path, "sim_states.csv"))
    }
    return(list(states=states_df, sim_cases=sim_cases_df))
}

## - ## - ##
forecast_mle <- function(path)
{
    source(paste0(path, "object.R"))
    fcast_df <- read_csv(paste0(path, "dataset_test.csv"), show_col_types=FALSE)
    po_fcast <- construct_pomp(path, fcast_df)

    if (file.exists(paste0(path, "sim_states.csv"))){
        states_df <- read_csv(paste0(path, "sim_states.csv"), show_col_types=FALSE)
        states_df <- states_df[states_df$type == 'filter_states', ]

        if (nrow(states_df) == 0) {
            states_df <- simulate_mle(path, save_sims=F, save_filtered=T)
            states_df <- states_df[states_df$type == 'filter_states', ]
        }

    } else {
        states_df <- simulate_mle(path, save_sims=F, save_filtered=T)
        states_df <- states_df[states_df$type == 'filter_states', ]
    }

    x0             <- states_df[states_df$time == tail(sort(unique(states_df$time)), 1), ]
    x0             <- x0 %>% select(.id, all_of(c(state_names)))
    x0             <- x0 %>% pivot_longer(-.id) %>% spread(.id, value) %>% column_to_rownames("name") %>% as.matrix()

    mle                <- read_csv(paste0(path, "results.csv"), show_col_types=FALSE) %>% head(1)%>% select(all_of(par_names))
    x0                 <- x0 / colSums(x0)
    rownames(x0)       <- paste0(rownames(x0), "_0")
    pp                 <- parmat(unlist(mle %>% select(all_of(par_names))), ncol(x0))
    pp[rownames(x0), ] <- x0

    timezero(po_fcast) <- tail(sort(unique(states_df$time)), 1)
    forecast           <- po_fcast %>% simulate(params=pp, format="data.frame")
    fcast_df           <- forecast %>% select(time, .id, all_of(c(state_names, accum_names, obs_vars)))
    path_to_save_fcast <- paste0(path, "forecast.csv")
    fcast_df$type      <- "forecast"

    write.csv(fcast_df, path_to_save_fcast, row.names = FALSE)
    return(fcast_df)
}
