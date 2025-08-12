time_name <- "time"
unit_name <- "loc"
obs_name <- "cases"

unit_covar <- c("pop","dpopdt")
shared_covar <- c()

state_names <- c("S","I","R")
accum_names <- c("C")

ivp_pars <- paste0(state_names,"0_")
specific_pars <- c("rho_","sigma_",ivp_pars)
shared_pars <- c("delta_","beta_","gamma_")

log_transf <- c("delta_","beta_","gamma_","sigma_")

logit_transf <- c("rho_",ivp_pars)

barycentric_transf <- NULL

param_bounds <- list(
    delta_=c(1/60,1/60),
    beta_=c(5,20),
    gamma_=c(365/30,365/7),
    rho_=c(0.001,0.05),
    S0_=c(0.5,1),
    I0_=c(0,0.01),
    R0_=c(0,0.4),
    sigma_=c(0.01,0.1))

global_vals <-  spatPomp_Csnippet("
    const double D[4][4] = {
        {0, 1, 1, 1},
        {1, 0, 1, 1},
        {1, 1, 0, 1},
        {1, 1, 1, 0}
    };
")

rinit <- spatPomp_Csnippet(
    unit_statenames = c(state_names,accum_names),
    unit_covarnames = unit_covar,
    unit_paramnames = paste0(c(specific_pars,shared_pars),"_"),
    code = "
        for (int u=0; u<U; u++) {
            double m = pop[u]/(S0_[u]+I0_[u]+R0_[u]);
            
            S[u] = nearbyint(m*S0_[u]);
            I[u] = nearbyint(m*I0_[u]);
            R[u] = nearbyint(m*R0_[u]);
            
            C[u] = 0;
        }
")

rproc <- spatPomp_Csnippet(
    unit_statenames = c(state_names,accum_names),
    unit_covarnames = unit_covar,
    unit_paramnames = paste0(c(specific_pars,shared_pars),"_"),
    code = "
        for (int u=0; u<U; u++) {
            S[u] = S[u]>0 ? floor(S[u]) : 0;
            I[u] = I[u]>0 ? floor(I[u]) : 0;
            R[u] = R[u]>0 ? floor(R[u]) : 0;
            
            S[u] += (dpopdt[u] + pop[u]*delta_[0] - beta_[0]*S[u]*I[u]/pop[u] - delta_[0]*S[u])*dt;
            I += (beta_[0]*S[u]*I[u]/pop[u] - gamma_[0]*I[u] - delta_[0]*I[u])*dt;
            R += (gamma_[0]*I[u] - delta_[0]*R[u])*dt;
            
            C += rho_[u]*(beta_[0]*S[u]*I[u]/pop[u])*dt;
        }
")

dunit_meas <- spatPomp_Csnippet(
    code = "
        double *sigma_ = &sigma_1;
        static double tol = 1.0e-18;
    
        double size = 1.0/sigma_[u]/sigma_[u];
        double prob = size/(C+size);
        
        lik = dnbinom(cases,size,prob+tol,0)+tol;
        if (give_log) lik = log(lik);
")

runit_meas <- spatPomp_Csnippet(
    code = "
        double *sigma_ = &sigma_1;
        double size = 1.0/sigma_[u]/sigma_[u];
        double prob = size/(C+size);
        
        cases = rnbinom(size,prob);
")

dmeas <- spatPomp_Csnippet(
    code = "
        double *cases = &cases1;
        double *C = &C1;
        double *sigma_ = &sigma_1;
        
        static double tol = 1.0e-18;
        double size;
        double prob;
        int u;
        
        lik = 0;
        for (u = 0; u < U; u++) {
            size = 1.0/sigma_[u]/sigma_[u];
            prob = size/(C[u]+size);
            lik += dnbinom(cases[u],size,prob+tol,1)+tol;
        }
        if (!give_log) lik = exp(lik);
")

rmeas <- spatPomp_Csnippet( 
    code = "
        double *cases = &cases1;
        double *C = &C1;
        double *sigma_ = &sigma_1;

        double size;
        double prob;
        int u;
        
        for (u = 0; u < U; u++) {
            size = 1.0/sigma_[u]/sigma_[u];
            prob = size/(C[u]+size);
            cases[u] = rnbinom(size,prob);
        }
")
