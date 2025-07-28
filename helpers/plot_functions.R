

make_plot <- function(sim_cases_data.frame, path_to_save_fig, enso) {

    set_0      <- function(x) (ifelse(is.na(x),0,x))
	sims_cases <- sim_cases_data.frame %>%
                    mutate(.id=ifelse(.id=="data", "data", "sim")) %>%
		            mutate_at(obs_vars, set_0)

	if (enso) {
		enso_bounds <- df %>%
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

	if (enso) {
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

	ggsave(filename="plot.png",
		plot = plotter,
		path = path_to_save_fig,
		width = 12,
		height = 6,
		units = "in") %>% suppressWarnings()
}

make_panel_plot <- function(path, mle, enso) {

    source(paste0(path, "object.R"))
    df <- read_csv(paste0(path, "dataset.csv"), show_col_types=FALSE)

    df[,loc_key]       <- str_remove_all(unlist(df[,loc_key])," ")
    df[,aggregate_key] <- str_remove_all(unlist(df[,aggregate_key])," ")

    po                 <- construct_panel_pomp(path,df,NA)
    po                 <- construct_panel_pomp(path,df,NA)

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
            select(time,.id,all_of(obs_vars)) |>
            mutate(.id=ifelse(.id=="data","data","sim"),
                   unit=.,
                   agg=agg_matches[.]) |>
            mutate_at(obs_vars,set_0)
        }) |> bind_rows()

    if (enso) {
        enso_bounds <- df |>
            filter(!is.na(cases)) |>
            mutate(enso=nino+nina) |>
            select(time,enso) |>
            distinct() |>
            mutate(upper=sims |>
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

    plotter <- ggplot() +
        ggdist::stat_lineribbon(data=sims,
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
    
    if (enso) {
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
        ylim(0,5000)+
        facet_wrap(aggregate_key,ncol=1)

    ggsave(filename="plot.png",
           plot = plotter,
           path = path,
           width = 12,
           height = 3+3*length(unique(region_matches)),
           units = "in") %>% suppressWarnings()
}