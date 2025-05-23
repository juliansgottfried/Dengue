run_fitting <- function(
        po, n_cores, parameters,
        seed_num, rdd1, rdd2, rdd3,
        start_index,
        result_path,
        log_path,
        stats_path) {
    
    cl <- parallel::makeCluster(n_cores)
    registerDoParallel(cl)
    registerDoRNG(seed = as.integer(round(abs(seed_num))) + 1234)
    
    ## for each parameter row, run mif
    r1 <- foreach::foreach(
        i = seq_len(nrow(parameters)),
        .packages = c("pomp", "dplyr")
    ) %dopar% {
        ## collect parameter vector
        param <- as.numeric(parameters[i, ])
        names(param) <- colnames(parameters)
        
        cat(paste("Starting iteration",
                  start_index + i - 1, "\n"),
            file = log_path,
            append = TRUE)
            
        mifout <- tryCatch(po |>
                               mif2(Np = 1000,
                                    Nmif = 50,
                                    cooling.type = "geometric",
                                    cooling.fraction.50 = 0.5,
                                    params = param,
                                    rw.sd = rdd1) |>
                               mif2(rw.sd = rdd2) |>
                               mif2(rw.sd = rdd3),
                           error = function(e) e)

        stats <- data.frame(cond=mifout@cond.logLik,
                            eff=mifout@eff.sample.size,
                            time=po@times)
        if (file.exists(stats_path)) {
            read.csv(stats_path) %>% bind_rows(stats) %>% write.csv(stats_path)
        } else write.csv(stats,stats_path)
            
        result <- c(rep(NA, length(param) + 4))
        names(result) <- c("sample",
                           colnames(parameters),
                           "loglik",
                           "loglik.se",
                           "flag")
        result[1] <- start_index + i - 1
        if (length(coef(mifout)) > 0) {
            loglik_mif <- tryCatch(replicate(n = 10,
                                             logLik(pfilter(po,
                                                            params = coef(mifout),
                                                            Np = 2000))),
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
