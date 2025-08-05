plot_trace <- function(path,
                       plot_pars,
                       plot_lik=T,
                       plot_mean=F) {
    
    #plot_pars <- c("S_0","muEI","muRS","I_0","R_0","muRS")
    if (plot_lik) plot_pars <- c(plot_pars,"loglik")
    
    plot_rows <- 2
    plot_cols = ceiling(length(plot_pars)/plot_rows)
    
    traces <- read_csv(paste0(path,"traces.csv")) %>% 
        arrange(run,iter) %>% 
        na.omit() %>% suppressMessages()
    
    Nmif <- max(traces$iter)+1
    Nseq <- nrow(traces[traces$iter==0 & traces$run==1,])
    Nrun <- max(traces$run)
    
    traces$i=traces$iter+Nmif*(traces$run-1)
    traces$sample <- rep(1:Nseq,Nmif*Nrun)
    traceplot <- traces %>%
        select(all_of(plot_pars),sample,i) %>% 
        pivot_longer(all_of(plot_pars)) %>% 
        group_by(i,name) %>% 
        summarize(mean=mean(value),
                  median=median(value),
                  upper=quantile(value,0.975),
                  lower=quantile(value,0.025)) %>%
        suppressMessages() %>% 
        ggplot(aes(x=i))+
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)+
        geom_line(aes(y=median),linewidth=0.5,color="#4f4f4f")+
        theme_classic()+
        labs(x="Iteration",y="Value")+
        facet_wrap(~name,scales="free_y",nrow=plot_rows)
    
    if (plot_mean) traceplot <- traceplot+geom_line(aes(y=mean),linewidth=0.5,color="#ff1239")

    ggsave(filename="traceplot.png",
           plot = traceplot,
           path = path,
           width = plot_cols*3,
           height = 6,
           units = "in") %>% suppressWarnings()
    
    #traces %>% 
    #    select(value=muIR,loglik) %>% 
    #    ggplot(aes(x=log(value),y=loglik))+
    #    geom_point()+
    #    geom_vline(aes(xintercept=log(median(value))))+
    #    geom_vline(aes(xintercept=log(mean(value))),color="#ff1239")+
    #    theme_classic()
}

plot_stats <- function(path,plot_mean=F) {
    stats <- read_csv(paste0(path,"stats.csv"))
    statplot <- stats %>% 
        pivot_longer(cols=c(cond,eff)) %>% 
        group_by(name,time) %>% 
        summarize(mean=mean(value),
                  median=median(value),
                  upper=quantile(value,0.975),
                  lower=quantile(value,0.025)) %>% 
        ggplot(aes(x=time))+
        geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)+
        geom_line(aes(y=median),linewidth=0.5,color="#4f4f4f")+
        facet_wrap(~name,
                   labeller=as_labeller(c(cond="conditional loglik",
                                        eff="effective sample size")),
                   scales="free_y")+
        theme_classic()+
        theme(text=element_text(size=14))+
        labs(x="Time",y="Value")
    
    if (plot_mean) statplot <- statplot+geom_line(aes(y=mean),linewidth=0.5,color="#ff1239")

    ggsave(filename="statplot.png",
           plot = statplot,
           path = path,
           width = 12,
           height = 6,
           units = "in") %>% suppressWarnings()
}
