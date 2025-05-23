rproc <- Csnippet("
    double noise;
    double k, foi;
    double dpop, births;
    double rate[8], trans[8];

    noise = rgammawn(sigma,dt);

    k = (theta+exp(bNino*nino+bNina*nina))*phi;
    foi = (beta_out+k*beta_eff*pow(I,alpha)/pop)*(noise/dt);

    dpop = rpois(dpopdt*dt);

    births = rbinom(nearbyint(pop),1-exp(-delta*dt));

    // Exit from S
    rate[0] = F;
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

    reulermultinom(2,nearbyint(S),&rate[0],dt,&trans[0]);
    reulermultinom(2,nearbyint(E),&rate[2],dt,&trans[2]);
    reulermultinom(2,nearbyint(I),&rate[4],dt,&trans[4]);
    reulermultinom(2,nearbyint(R),&rate[6],dt,&trans[6]);

    S += dpop + births + trans[6] - trans[0] - trans[1];
    E += trans[0] - trans[2] - trans[3];
    I += trans[2] - trans[4] - trans[5];
    R += trans[4] - trans[6] - trans[7];
    
    K += ((foi-K)/(tau/2.0))*dt;
    F += ((K-F)/(tau/2.0))*dt;
    
    C += rho*trans[2];
    
    if (S < 0.0) S=0.0;
	if (E < 0.0) E=0.0;
	if (I < 0.0) I=0.0;
	if (R < 0.0) R=0.0;
	if (K < 0.0) K=0.0;
	if (F < 0.0) F=0.0;
    
    if (t<2000.01) {
        // Rprintf(\"%lg %lg %lg %lg %lg \\n\",S,E,I,R,C);
        // Rprintf(\"%lg %lg \\n\",births, dpop);
    }
")

rinit <- Csnippet("
    double m = pop/(S_0+E_0+I_0+R_0);
                    
    S = nearbyint(m*S_0);
    E = nearbyint(m*E_0);
    I = nearbyint(m*I_0);
    R = nearbyint(m*R_0);
    
    K = K_0;
    F = F_0;
  
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

par_names <- c("sigma",
               "theta","phi",
               "bNina","bNino",
               "beta_out","alpha",
               "delta",
               "muEI","muIR","muRS",
               "tau",
               "rho",
               "S_0","E_0","I_0","R_0",
               "K_0","F_0",
               "sig_obs")

state_names <- c("S","E","I","R",
                 "K","F")

accum_names <- c("C")

log_transf <- c("sigma",
                "phi",
                "beta_out","alpha",
                "delta",
                "muEI","muIR","muRS",
                "tau",
                "K_0","F_0",
                "sig_obs")

logit_transf <- c("rho")

barycentric_transf <- c("S_0","E_0","I_0","R_0")

params <- c(sigma=0.005,
            theta=3,phi=200,
            bNina=-0.5,bNino=2,
            beta_out=0.01,alpha=0.9,
            delta=1/60,
            muEI=365/10,muIR=365/10,muRS=0.5,
            tau=0.01,
            rho=0.001,
            S_0=0.5,E_0=0.01,I_0=0.1,R_0=0.5,
            K_0=0.1,F_0=0.1,
            sig_obs=0.05)

data_vars <- c("time","cases")
covars <- c("time","pop","dpopdt","beta_eff","nino","nina")
