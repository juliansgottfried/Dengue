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

construct_panelpomp <- function(path, nseq) {
    source(paste0(path,"object.R"))

    df <- read_csv(paste0(path, "dataset.csv")) |> filter(train)
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

construct_spatpomp <- function(path,df,
                               basic_par_names,
                               par_trans) {
    source(paste0(path,"object.R"))
    spo <- spatPomp(
        data = df %>% select(all_of(c(time_name,unit_name,obs_name))) %>% na.omit,
        covar = df %>% select(all_of(c(time_name,unit_name,unit_covar,shared_covar))),
        times = time_name,
        units = unit_name,
        t0 = df$time[1]-1/365,
        unit_statenames = c(state_names,accum_names),
        unit_accumvars = accum_names,
        paramnames = basic_par_names,
        rinit = rinit,
        rprocess = euler(rproc, delta.t=1/365),
        dunit_measure = dunit_meas,
        runit_measure = runit_meas,
        dmeasure = dmeas,
        rmeasure = rmeas,
        globals = global_vals,
        partrans = par_trans
    )
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
                               mif2(Np                  = Np1,
                                    Nmif                = Nmif,
                                    cooling.type        = "geometric",
                                    cooling.fraction.50 = 0.5,
                                    params              = param,
                                    rw.sd               = rdds[[1]]),
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
	    Np1, Np2, Nmif,
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
            cbind(as.data.frame(traces(mifout)),
                  data.frame(iter=0:Nmif,
                             run=i))
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

    suppressWarnings(write.table(rw, 
                resultw_path,
                append = TRUE,
                col.names = !file.exists(resultw_path),
                row.names = FALSE, sep = ","))
    suppressWarnings(write.table(rl, 
                resultl_path,
                append = TRUE,
                col.names = !file.exists(resultl_path),
                row.names = FALSE, sep = ","))

    stopCluster(cl)
}

run_spatial_fitting <- function(po,n_cores,parameters,
                                unitParNames,sharedParNames,
                                par_names,
                                seed_num,
                                rdd1,rdd2,rdd3,
                                n_refine,Np1,Np2,Nbpf,
                                block_size=2,
                                resultw_path,resultl_path,
                                log_path,traces_path,stats_path
                                ) {
    
    #po=spo
    #parameters=init_vals
    #unitParNames=paste0(est_specific,"_")
    #sharedParNames=paste0(est_shared,"_")
    #par_names=c(specific_pars,shared_pars)
    #seed_num=seed
    #Np1=5
    #Np2=5
    #Nbpf=1
    #block_size=2
    #resultw_path=paste0(path,"results.csv")
    #resultl_path=paste0(path,"results_long.csv")
    #log_path=paste0(path,"log.txt")
    #traces_path=paste0(path,"traces.csv")
    #stats_path=paste0(path,"stats.csv")
    
    cl <- parallel::makeCluster(n_cores)
    registerDoParallel(cl)
    registerDoRNG(seed = as.integer(round(abs(seed_num))) + 1234)
    
    ## for each parameter row, run ibpf
    rs <- foreach::foreach(
        i = seq_len(nrow(parameters)),
        .packages = c("spatPomp", "dplyr", "readr", "stringr")
    ) %dopar% {
        #i = 1
        param <- as.numeric(parameters[i, ])
        names(param) <- colnames(parameters)
        
        #coef(po) <- param
        #sim <- po %>% simulate(nsim=10,format="data.frame")
        #bpfs <- po %>% bpfilter(Np=10,block_size=2,save_states=T)
        #bpfs@loglik
        
        cat(paste("Starting iteration", i, "\n"),
            file = log_path,
            append = TRUE)
        
        rdds = list(rdd1,rdd2,rdd3)
        
        bpfout <- tryCatch(po |>
                               ibpf(Np = Np1,
                                    Nbpf = Nbpf,
                                    cooling.type = "geometric",
                                    cooling.fraction.50 = 0.5,
                                    params = param,
                                    unitParNames = unitParNames,
                                    sharedParNames = sharedParNames,
                                    spat_regression = 0.1,
                                    block_size = block_size,
                                    rw.sd = rdds[[1]]),
                           error = function(e) e)
        
        get_trace <- function(i) {
            cbind(as.data.frame(bpfout@traces),data.frame(iter=0:Nbpf,run=i))
        }
        
        traces <- get_trace(1)
        
        if (n_refine > 0) {
            for (j in 1:n_refine) {
                bpfout <- bpfout |> ibpf(rw.sd = rdds[[j+1]],
                                         unitParNames = unitParNames,
                                         sharedParNames = sharedParNames)
                traces_j <- get_trace(j+1)
                traces <- bind_rows(traces,traces_j)
            }
        }
        
        if (file.exists(traces_path)) {
            read_csv(traces_path) %>% bind_rows(traces) %>% write_csv(traces_path)
        } else write_csv(traces,traces_path)
        
        stats <- data.frame(cond=bpfout@cond.loglik,time=bpfout@times)
        
        if (file.exists(stats_path)) {
            read_csv(stats_path) %>% bind_rows(stats) %>% write_csv(stats_path)
        } else write_csv(stats,stats_path)
        
        resultw <- c(rep(NA, length(param) + 4))
        names(resultw) <- c("sample",colnames(parameters),
                            "loglik","loglik.se","flag")
        
        keys <- po@unit_names
        resultl <- matrix(rep(NA, length(keys)*(length(par_names) + 5)),nrow=length(keys))
        colnames(resultl) <- c("sample","unit",par_names,
                               "loglik","loglik.se","flag")
        resultl <- data.frame(resultl)
        resultl$unit <- keys
        
        resultw["sample"] <- i
        resultl$sample <- i
        
        if (length(coef(bpfout)) > 0) {
            loglik_bpf <- tryCatch(replicate(n = 10,
                                             logLik(bpfilter(po,
                                                             params = coef(bpfout),
                                                             Np = Np2,
                                                             block_size = 2))),
                                   error = function(e) e)
            
            if (is.numeric(loglik_bpf)) {
                bl <- logmeanexp(loglik_bpf, se = TRUE)
                names(bl) <- c("loglik","loglik.se")
                loglik_bpf_est <- bl[1]
                cat(paste(i, loglik_bpf_est, "\n"), file = log_path, append = TRUE)
                
                resultw["flag"] <- 2
                resultl$flag <- 2
            }
            if (is.numeric(loglik_bpf)) {
                resultw[names(bl)] <- bl
                resultl[,names(bl)] <- t(matrix(bl))[rep(1,length(keys)),]
                
                resultw["flag"] <- 1
                resultl$flag <- 1
                
                par_out <- coef(bpfout)
                resultw[names(par_out)] <- par_out
                
                ordered_par_names <- unique(unlist(lapply(str_split(names(par_out),"_"),\(.) .[1])))
                
                listed <- setNames(rep(list(1:length(keys)),length(par_out)/length(keys)),ordered_par_names)
                for (i in 1:length(listed)) {
                    listed[[i]] <- listed[[i]]+5*(i-1)
                    listed[[i]] <- unname(par_out[listed[[i]]])
                }
                par_out <- data.frame(listed)
                
                resultl[,colnames(par_out)] <- par_out
            }
        } else {
            resultw["flag"] <- 3
            resultl$flag <- 3
            cat(paste(i, "failed", "\n"), file = log_path, append = TRUE)
        }
        list(resultw,resultl)
    }
    
    rw <- lapply(rs,\(.) .[[1]]) |> 
        bind_rows() |>
        remove_missing()
    rl <- lapply(rs,\(.) .[[2]]) |> 
        bind_rows() |>
        remove_missing()
    
    suppressWarnings(write.table(rw, 
                                 resultw_path,
                                 append = TRUE,
                                 col.names = !file.exists(resultw_path),
                                 row.names = FALSE, sep = ","))
    suppressWarnings(write.table(rl, 
                                 resultl_path,
                                 append = TRUE,
                                 col.names = !file.exists(resultl_path),
                                 row.names = FALSE, sep = ","))
    
    stopCluster(cl)
}
