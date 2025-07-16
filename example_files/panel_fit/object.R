rproc <- Csnippet("
    S += (dpopdt + pop*delta - beta*S*I/N - delta*S)*dt;
    I += (beta*S*I/N - gamma*I - delta*I)*dt;
    R += (gamma*I - delta*R)*dt;
    
    C += rho*(beta*S*I/N)*dt;
    
    if (S < 0.0) S = 0.0;
	if (I < 0.0) I = 0.0;
	if (R < 0.0) R = 0.0;
")

rinit <- Csnippet("
    double m = pop/(S_0+I_0+R_0);
                    
    S = nearbyint(m*S_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    
    C = 0;
")

dmeas <- Csnippet("
    static double tol = 1.0e-18;

    double size = 1.0/sigma/sigma;
    double prob = size/(C+size);
    
    lik = dnbinom(cases,size,prob+tol,0)+tol;
    if (give_log) lik = log(lik);
")

rmeas <- Csnippet("
    double size = 1.0/sigma/sigma;
    double prob = size/(C+size);
    
    cases = rnbinom(size,prob);
")

par_names <- c("delta","beta","gamma",
               "rho",
               "S_0","I_0","R_0",
               "sigma")

specific_names <- c("rho",
                    "S_0","I_0","R_0")

state_names <- c("S","I","R")

accum_names <- c("C")

log_transf <- c("delta","beta","gamma",
                "sigma")

logit_transf <- c("rho",
                  "S_0","I_0","R_0")

barycentric_transf <- c()

obs_vars <- "cases" 

loc_key <- "loc"

aggregate_key <- "aggregate"

param_bounds <- list(
    delta=c(1/60,1/60),
    beta=c(5,20),
    gamma=c(365/30,365/7),
    rho=c(0.001,0.05),
    S_0=c(0.5,1),
    I_0=c(0,0.01),
    R_0=c(0,0.4),
    sigma=c(0.01,0.1))
