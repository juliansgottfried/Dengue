simulate_mle <- function(path, save_sims=TRUE, save_filtered = TRUE) {

    df   <- read_csv(paste0(path, "dataset.csv"), show_col_types=FALSE) |> 
        suppressMessages() |>
        filter(train)
    
    po   <- construct_pomp(path, df)
    mle  <- read_csv(paste0(path, "results.csv"), show_col_types=FALSE) %>% 
        suppressMessages() |>
        head(1) %>% select(all_of(par_names))

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

simulate_panel_mle <- function(path, save_sims=TRUE, save_filtered = TRUE) {

    source(paste0(path, "object.R"))
    df <- read_csv(paste0(path, "dataset.csv"), show_col_types=FALSE) |> 
        suppressMessages() |>
        filter(train)
    
    results <- read_csv(paste0(path,"results.csv")) |> 
        suppressMessages()
    
    df[,loc_key]       <- str_remove_all(unlist(df[,loc_key])," ")
    df[,aggregate_key] <- str_remove_all(unlist(df[,aggregate_key])," ")
    
    po                 <- construct_panel_pomp(path,500)
    
    mle <- results %>% select(all_of(names(coef(po)))) %>% slice_head(n=1)
    
    coef(po,names(mle)) <- mle
    
    set_0 <- function(x) (ifelse(is.na(x),0,x))
    keys  <- names(po@unit_objects)
    
    agg_matches <- df |>
        select(all_of(c(loc_key,aggregate_key))) |>
        distinct() |>
        pull(all_of(aggregate_key))
    names(agg_matches) <- keys
    
    sims <- lapply(keys, \(.) {
        obj <- po@unit_objects[[.]]
        coef(obj) <-  c(po@shared,po@specific[,.])
        obj |> simulate(nsim=1000,
                        include.data=TRUE,
                        format="data.frame") |>
            mutate(unit=.,
                   agg=agg_matches[.]) |>
            mutate_at(obs_vars,set_0)
    }) |> bind_rows()

    sim_states_df        <- sims %>% select(time, .id, unit, agg, all_of(c(state_names, accum_names)))
    sim_states_df$type   <- "sim_states"
    
    sim_cases_df <- sims %>% select(time, .id, unit, agg, all_of(obs_vars))
    
    if (save_filtered==T){
        filt_sims_states_df <- lapply(keys, \(.) {
            obj <- po@unit_objects[[.]]
            coef(obj) <-  c(po@shared,po@specific[,.])
            obj |> pfilter(Np=1000,save.states="filter") |>
                saved_states(format="data.frame") |>
                mutate(unit=.,
                       agg=agg_matches[.])
            }) |> bind_rows()
        
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

forecast_mle <- function(path)
{
    source(paste0(path, "object.R"))
    fcast_df <- read_csv(paste0(path, "dataset.csv"), show_col_types=FALSE) |> filter(!train)
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

sim_plot <- function(sim_cases_dataframe, path_to_save_fig, enso_df) {
	sims_cases <- sim_cases_dataframe %>%
                    mutate(.id=ifelse(.id=="data", "data", "sim")) 

	if (length(enso_df)!=0) {
		enso_bounds <- enso_df %>%
    			filter(!is.na(cases)) %>%
    			mutate(enso=nino+nina) %>%
    			select(time,enso) %>%
    			mutate(upper=sims_cases %>%
    				filter(.id!="data") %>%
    				group_by(time) %>%
    				summarize(upper=quantile(cases, 0.95)) %>%
    				pull(upper))
	}

	fill_colors  <- c("black", "#40d6ed")
	color_colors <- c("black", "#16a4ba")
	labels       <- c("Data","Simulation")
	bg.color     <- "#e6e6e6"
	bg.color     <- "white"

	plotter <- ggplot() +
		ggdist::stat_lineribbon(data=sims_cases,
					mapping   = aes(x=time,
						y     = cases,
						color = .id,
						fill. = .id),
					.width    = c(0.1,0.9),
					alpha     = 0.65,
					linewidth = 0.8) +
		labs(x = "Year",
			y  = "Cases",
			title = "Dengue in Thailand") +
		scale_fill_manual(values=fill_colors,
				labels=labels) +
		scale_color_manual(values=color_colors,
				labels=labels) +
		guides(color="none",
			fill=guide_legend(title=""))

	if (length(enso_df)!=0) {
		plotter <- plotter   +
		    new_scale_fill() +
    			geom_rect(data=enso_bounds,
    				mapping=aes(xmin=time,
    					xmax = time+1/12,
    					ymin = upper,
    					ymax = Inf,
    					fill = enso),
    				alpha    = 0.45,
    				inherit.aes = F)+
    			scale_fill_gradient2(
                    low  = "#f2c51f",
					mid  = bg.color,
                    high = "#f54997",
                    name = "ENSO")
	}

	plotter <- plotter +
	    theme_classic()+
		theme(legend.position="bottom",
			legend.title=element_text(size=10,vjust=0.8),
			text=element_text(size = 14),
			legend.background = element_rect(fill=bg.color,color=bg.color),
			panel.background  = element_rect(fill=bg.color,color=bg.color),
			plot.background   = element_rect(fill=bg.color,color=bg.color))+
		ylim(0, 5000)

	ggsave(filename="simplot.png",
		plot = plotter,
		path = path_to_save_fig,
		width = 12,
		height = 6,
		units = "in") %>% suppressWarnings()
}

panel_sim_plot <- function(sim_cases_dataframe, path_to_save_fig, enso_df, ordering) {
    sims_cases <- sim_cases_dataframe %>%
        mutate(.id=ifelse(.id=="data", "data", "sim")) 

    if (length(enso_df)!=0) {
        enso_bounds <- enso_df |>
            filter(!is.na(cases)) |>
            mutate(enso=nino+nina) |>
            select(time,enso) |>
            distinct() |>
            mutate(upper=sims_cases |>
                       filter(.id!="data") |>
                       group_by(time)|>
                       summarize(upper=quantile(cases,0.95))|>
                       pull(upper))
    }

    fill_colors  <- c("black","#40d6ed")
    color_colors <- c("black","#16a4ba")
    
    labels   <- c("Data","Simulation")
    bg.color <- "#e6e6e6"
    bg.color <- "white"

    sims_cases$agg <- factor(sims_cases$agg,levels=ordering)

    plotter <- ggplot() +
        ggdist::stat_lineribbon(data=sims_cases,
                                mapping=aes(x=time,
                                            y=cases,
                                            color=.id,
                                            fill=.id),
                                .width = c(0.1,0.9),
                                alpha=0.65,
                                linewidth=0.8)+
        labs(x="Year",
             y="Cases",
             title="Dengue in Thailand")+
        scale_fill_manual(values=fill_colors,
                          labels=labels)+
        scale_color_manual(values=color_colors,
                           labels=labels)+
        guides(color="none",
               fill=guide_legend(title=""))
    
    if (length(enso_df)!=0) {
        plotter <- plotter + 
            new_scale_fill()+
            geom_rect(data=enso_bounds,
                      mapping=aes(xmin=time,
                                  xmax=time+1/12,
                                  ymin=upper,
                                  ymax=Inf,
                                  fill=enso),
                      alpha=0.45,
                      inherit.aes=F)+
            scale_fill_gradient2(low="#f2c51f",
                                 mid=bg.color,
                                 high="#f54997",
                                 name="ENSO")
    }

    plotter <- plotter +
        theme_classic()+
        theme(legend.position="bottom",
              legend.title=element_text(size=10,vjust=0.8),
              text=element_text(size = 14),
              legend.background=element_rect(fill=bg.color,color=bg.color),
              panel.background=element_rect(fill=bg.color,color=bg.color),
              plot.background=element_rect(fill=bg.color,color=bg.color))+
        facet_wrap("agg",ncol=1,scales="free_y")

    ggsave(filename="simplot.png",
           plot = plotter,
           path = path,
           width = 12,
           height = 3+3*length(unique(agg_matches)),
           units = "in") %>% suppressWarnings()
}

plot_trace <- function(path,
                       plot_pars,
                       plot_mean=F) {
    
    #plot_pars <- c("S_0","muEI","muRS","I_0","R_0","muRS")
    if (length(plot_pars)==0) {
        plot_pars <- "loglik"
    } else plot_pars <- c(plot_pars,"loglik")
    
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
           width = plot_cols*6,
           height = 6,
           units = "in") %>% suppressMessages()
    
    #traces %>% 
    #    select(value=muIR,loglik) %>% 
    #    ggplot(aes(x=log(value),y=loglik))+
    #    geom_point()+
    #    geom_vline(aes(xintercept=log(median(value))))+
    #    geom_vline(aes(xintercept=log(mean(value))),color="#ff1239")+
    #    theme_classic()
}

plot_stats <- function(path,plot_mean=F) {
    stats <- read_csv(paste0(path,"stats.csv")) %>% suppressMessages()
    statplot <- stats %>% 
        pivot_longer(cols=c(cond,eff)) %>% 
        group_by(name,time) %>% 
        summarize(mean=mean(value),
                  median=median(value),
                  upper=quantile(value,0.975),
                  lower=quantile(value,0.025)) %>%
        suppressMessages() %>% 
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
           units = "in") %>% suppressMessages()
}
