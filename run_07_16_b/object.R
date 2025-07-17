rproc <- Csnippet("
    double noise;
    double k, beta, foi;
    int sign;
    double dpop, births;
    double rate[8], trans[8];

    noise = rgammawn(sigma,dt);

    k    = theta*exp(season6*beta_eff_cov*(bONI*oni + bRONI*roni));
    beta = k*beta_eff;

    foi  = (beta_out+beta*I/pop)*(noise/dt);

    sign = (0<dpopdt)-(dpopdt<0);
    dpop = sign*rpois(sign*dpopdt*dt);

    births = rbinom(nearbyint(pop), 1-exp(-delta*dt));

    // Exit from S
    rate[0] = foi;
    rate[1] = delta;

    // Exit from E
    rate[2] = muEI;
    rate[3] = delta;

    // Exit from I
    rate[4] = muIR;
    rate[5] = delta;

    // Exit from R
    rate[6] = muRS;
    rate[7] = delta;

    reulermultinom(2, nearbyint(S), &rate[0],dt, &trans[0]);
    reulermultinom(2, nearbyint(E), &rate[2],dt, &trans[2]);
    reulermultinom(2, nearbyint(I), &rate[4],dt, &trans[4]);
    reulermultinom(2, nearbyint(R), &rate[6],dt, &trans[6]);

    S += dpop + births + trans[6] - trans[0] - trans[1];
    E += trans[0] - trans[2] - trans[3];
    I += trans[2] - trans[4] - trans[5];
    R += trans[4] - trans[6] - trans[7];

    C += rho*trans[2];

    if (S < 0.0) S=0.0;
	if (E < 0.0) E=0.0;
	if (I < 0.0) I=0.0;
	if (R < 0.0) R=0.0;
")

rinit <- Csnippet("
    double m = pop/(S_0+E_0+I_0+R_0);

    S = nearbyint(m*S_0);
    E = nearbyint(m*E_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);

    C = 0;
")

dmeas <- Csnippet("
    double size = 1.0/sig_obs/sig_obs;
    double prob = size/(C+size);
    static double tol = 1.0e-18;
    lik = dnbinom(cases,size,prob+tol,0)+tol;
    if (give_log) lik = log(lik);
")

rmeas <- Csnippet("
    double size = 1.0/sig_obs/sig_obs;
    double prob = size/(C+size);
    cases = rnbinom(size,prob);
")

par_names <- c("theta",
               "beta_out",
               "sigma",
               "bONI", "bRONI",
               "delta",
               "muEI", "muIR", "muRS",
               "tau",
               "rho",
               "S_0", "E_0","I_0","R_0",
               "sig_obs")

state_names <- c("S", "E", "I", "R")
accum_names <- c("C")

log_transf <- c("theta",
                "delta",
                "beta_out",
                "sigma",
                "muEI", "muIR", "muRS",
                "tau",
                "K_0","F_0",
                "sig_obs")

logit_transf       <- c("rho")
barycentric_transf <- c("S_0","E_0","I_0","R_0")

param_bounds <- list(
            theta    = c(10, 40),
            beta_out = c(0.001, 0.01),
            sigma    = c(0.001, 0.1),
            bONI     = c(0, 3),
            bRONI    = c(0, 3),
            delta    = c(1/60, 1/60),
            muEI     = c(365/30, 365/7),
            muIR     = c(365/30, 365/7),
            muRS     = c(0.5, 2),
            rho      = c(0.01, 0.05),
            S_0      = c(0.6, 0.99),
            E_0      = c(0, 0.01),
            I_0      = c(0, 0.01),
            R_0      = c(0.05, 0.3),
            sig_obs  = c(0.05, 0.5))

data_vars <- c("time", "cases")
obs_vars  <- "cases"
covars    <- c("time", "pop", "dpopdt", "beta_eff", "oni", "roni", "beta_eff_cov")