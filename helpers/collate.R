oldw <- getOption("warn")
options(warn = -1)

suppressMessages(library(tidyverse))
suppressMessages(library(ggdist))
suppressMessages(library(ggnewscale))
suppressMessages(library(panelPomp))

isPanel <- args[1]

log_name  <- Sys.getenv("LOGNAME")
path_name <- paste0("/scratch/",log_name,"/Dengue")

source(paste0(path_name, "/helpers/helper_functions.R"))

time_df = read.table(paste0(path_name,"/times.txt"), sep=" ") |> 
    t() |> 
    na.omit() |> 
    t() |> 
    data.frame()
colnames(time_df) = c("time","run")

time_df=time_df |>
    mutate(time  = hms(time)) |>
    mutate(hours = hour(time),minutes=minute(time)) |>
    mutate(time  = hours*60+minutes) |>
    select(-c(hours,minutes)) |>
    group_by(run) |>
    summarize(time=max(time)) |>
    ungroup()

process <- function(path, par_names) {
	read_csv(path) %>%
		select(contains(c("sample", "unit", par_names,"loglik","loglik.se","flag",
			"time","cond","eff",
			"iter","run"))) %>%
		mutate_all(as.numeric) %>% suppressMessages()
	}

result_type <- c("results", "traces", "stats")
if (isPanel) result_type <- c("results_long", result_type)

summary <- lapply(result_type,function(type) {
    print(paste0(toupper(type)))

    paths<-list.files(paste0(path_name,"/out/",type),full.names=T)
    names<-list.files(paste0(path_name,"/out/",type),full.names=F)

    summary <- map2(paths,names,function(path,name) {
        print(paste0(name))

	fitting_folder_path <- paste0(path_name,"/folders_for_fit/",name,"/")
	source(paste0(fitting_folder_path, "object.R"))

        files <- list.files(path,full.names=T)
        accum <- process(files[1], par_names)

        for (i in 2:length(files)) {
		add   <- process(files[i], par_names)
		accum <- bind_rows(accum,add)
	}
	
	summary=NULL
	if (type=="results") {
		
		accum['ll_lowIQ'] <- accum$loglik - 0.6745 * accum$loglik.se
		accum    		  <- accum %>% arrange(-ll_lowIQ)

		init_vals <- read_csv(paste0(fitting_folder_path, "pars.csv"), show_col_types=F)
		k      	  <- sum(apply(init_vals,2,\(.) diff(range(.)))!=0)
	
		maxlik       <- accum %>% pull(loglik)  %>% head(1)
		maxlik_lowIQ <- accum %>% pull(ll_lowIQ) %>% head(1)

		aic       <- 2*(k-maxlik)
		aic_lowIQ <- 2*(k-maxlik_lowIQ)

		time    <- time_df %>% filter(run==name) %>% pull(time)
		summary <- c(run=name,time=time,loglik=maxlik,loglik_lowIQ=maxlik_lowIQ,k=k,aic=aic,aic_lowIQ=aic_lowIQ)

		mle <- accum %>% slice(1)
		
		print("(saving simulation plot)")

		if (isPanel) {
			make_panel_plot(fitting_folder_path,mle,T)
		} else make_plot(fitting_folder_path, mle, F)
	}
	
        write_csv(accum,paste0(path_name,"/folders_for_fit/",name,"/",type,".csv"))
	return(summary)
    })
    return(summary)
})

summary_data <- bind_rows(summary[[1]])
write_csv(summary_data,paste0(path_name,"/folders_for_fit/summary.csv"))

options(warn = oldw)
