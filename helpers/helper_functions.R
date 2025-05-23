#' Run fitting process
#'
#' This function runs the fitting process using the provided parameters.
#'
#' @param po The pomp object representing the model.
#' @param n_cores The number of CPU cores to use for parallel processing.
#' @param parameters A data frame containing the parameter values to be used in the fitting process.
#' @param seed_num The seed number for random number generation.
#' @param rdd The random walk standard deviation.
#' @param allow_parallel A logical value indicating whether to allow parallel processing.
#' @param start_index The starting index for the iterations.
#' @param output_to_file A logical value indicating whether to output the results to a file.
#' @param output_result The file path for the output results.
#' @param output_log The file path for the output log.
#' @param output_format The output format for the results (for now only "data.frame").
#'
#' @return A data frame containing the fitting results.
run_fitting <- function(
        po, n_cores, parameters,
        seed_num, rdd1, rdd2, rdd3, allow_parallel,
        start_index,
        output_to_file, output_result,
        output_log, output_format,
        stats_path) {
    if (allow_parallel) {
        ## if parallel register cluster
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
            
            if (output_to_file) {
                cat(paste(
                    "Starting iteration",
                    start_index + i - 1, "\n"
                ), file = output_log, append = TRUE)
            }
            
            mifout <- tryCatch(po |>
                                   mif2(
                                       Np = 1000, Nmif = 50,
                                       cooling.type = "geometric",
                                       cooling.fraction.50 = 0.5,
                                       params = param,
                                       rw.sd = rdd1
                                   ) |>
				   mif2(
				       rw.sd = rdd2
				   ) |>
				   mif2(
				       rw.sd = rdd3
				   ), error = function(e) e)

            stats <- data.frame(cond=mifout@cond.logLik,
				eff=mifout@eff.sample.size,
				time=po@times)
            if (file.exists(stats_path)) {
                read.csv(stats_path) %>% bind_rows(stats) %>% write.csv(stats_path)
            } else write.csv(stats,stats_path)
            
            result <- c(rep(NA, length(param) + 4))
            names(result) <- c(
                "sample", colnames(parameters),
                "loglik", "loglik.se", "flag"
            )
            result[1] <- start_index + i - 1
            if (length(coef(mifout)) > 0) {
                loglik_mif <- tryCatch(
                    replicate(
                        n = 10,
                        logLik(pfilter(po,
                                       params = coef(mifout),
                                       Np = 2000
                        ))
                    ),
                    error = function(e) e
                )
                
                
                if (is.numeric(loglik_mif)) {
                    bl <- logmeanexp(loglik_mif, se = TRUE)
                    loglik_mif_est <- bl[1]
                    loglik_mif_se <- bl[2]
                    if (output_to_file) {
                        cat(paste(i, loglik_mif_est, "\n"), file = output_log, append = TRUE)
                    }
                    result[length(result)] <- 2
                }
                ### SAVE OUTPUT
                if (is.numeric(loglik_mif)) {
                    par_out <- coef(mifout)
                    result[(length(param) + 2):(length(param) + 3)] <- bl
                    result[(2):(length(param) + 1)] <- par_out
                    result[length(result)] <- 1
                }
            } else {
                result[length(result)] <- 3
                if (output_to_file) {
                    cat(paste(i, "failed", "\n"), file = output_log, append = TRUE)
                }
            }
            result
        }
        if (output_format == "data.frame") {
            r1 <- r1 |>
                bind_rows() |>
                remove_missing()
            if (output_to_file) {
                write.table(r1, output_result,
                            append = TRUE,
                            col.names = !file.exists(output_result),
                            row.names = FALSE, sep = ","
                )
            }
        }
        stopCluster(cl)
    } else {
        ## for each parameter row, run mif
        r1 <- foreach::foreach(
            i = seq_len(nrow(parameters)),
            .packages = c("pomp", "dplyr")
        ) %do% {
            ## collect parameter vector
            param <- as.numeric(parameters[i, ])
            names(param) <- colnames(parameters)
            
            if (output_to_file) {
                cat(paste(
                    "Starting iteration",
                    start_index + i - 1, "\n"
                ), file = output_log, append = TRUE)
            }
            mifout <- tryCatch(po |>
                                   mif2(
                                       Np = 1000, Nmif = 50,
                                       cooling.type = "geometric",
                                       cooling.fraction.50 = 0.5,
                                       params = param,
                                       rw.sd = rdd1
                                   ), error = function(e) e)
            
            result <- c(rep(NA, length(param) + 4))
            names(result) <- c(
                "sample", colnames(parameters),
                "loglik", "loglik.se", "flag"
            )
            result[1] <- start_index + i - 1
            if (length(coef(mifout)) > 0) {
                loglik_mif <- tryCatch(
                    replicate(
                        n = 10,
                        logLik(pfilter(po,
                                       params = coef(mifout),
                                       Np = 1000
                        ))
                    ),
                    error = function(e) e
                )
                
                
                if (is.numeric(loglik_mif)) {
                    bl <- logmeanexp(loglik_mif, se = TRUE)
                    loglik_mif_est <- bl[1]
                    loglik_mif_se <- bl[2]
                    if (output_to_file) {
                        cat(paste(i, loglik_mif_est, "\n"), file = output_log, append = TRUE)
                    }
                    result[length(result)] <- 2
                }
                ### SAVE OUTPUT
                if (is.numeric(loglik_mif)) {
                    par_out <- coef(mifout)
                    result[(length(param) + 2):(length(param) + 3)] <- bl
                    result[(2):(length(param) + 1)] <- par_out
                    result[length(result)] <- 1
                }
            } else {
                result[length(result)] <- 3
                if (output_to_file) {
                    cat(paste(i, "failed", "\n"), file = output_log, append = TRUE)
                }
            }
            result
        }
        if (output_format == "data.frame") {
            r1 <- r1 |>
                bind_rows() |>
                remove_missing()
            if (output_to_file) {
                write.table(r1, output_result,
                            append = TRUE,
                            col.names = !file.exists(output_result),
                            row.names = FALSE, sep = ","
                )
            }
        }
    }
    return(r1)
}
