oldw <- getOption("warn")
options(warn = -1)

suppressMessages(library(tidyverse))
suppressMessages(library(pomp))

log_name <- Sys.getenv("LOGNAME")
path_name <- paste0("/scratch/",log_name,"/Dengue")

time_df = read.table(paste0(path_name,"/times.txt"),sep=" ") |> 
    t() |> 
    na.omit() |> 
    t() |> 
    data.frame()
colnames(time_df)=c("time","run")

time_df=time_df |>
    mutate(time=hms(time)) |>
    mutate(hours=hour(time),minutes=minute(time)) |>
    mutate(time=hours*60+minutes) |>
    select(-c(hours,minutes)) |>
    group_by(run) |>
    summarize(time=max(time)) |>
    ungroup()

process<-function(path,par_names) {
	read_csv(path) %>%
		select(matches(c("sample",par_names,"loglik","loglik.se","flag",
			"time","cond","eff",
			"iter","run"))) %>%
		mutate_all(as.numeric) %>% suppressMessages()
	}

result_type <- c("results","traces","stats")

summary <- lapply(result_type,function(type) {
    print(paste0(toupper(type)))

    paths<-list.files(paste0(path_name,"/out/",type),full.names=T)
    names<-list.files(paste0(path_name,"/out/",type),full.names=F)

    summary <- map2(paths,names,function(path,name) {
        print(paste0(name))

	source(paste0(path_name,"/folders_for_fit/",name,"/object.R"))

        files<-list.files(path,full.names=T)
        accum<-process(files[1],par_names)

        for (i in 2:length(files)) {
		add<-process(files[i],par_names)
		accum<-bind_rows(accum,add)
	}
	
	summary=NULL
	if (type=="results") {
		init_vals <- read_csv(paste0(path_name,"/folders_for_fit/",name,"/pars.csv"),show_col_types=F)
		k <- length(par_names[apply(init_vals,2,function(x) {max(x)-min(x)})!=0])
		maxlik <- accum %>% arrange(-loglik) %>% pull(loglik) %>% head(1)
		aic <- 2*(k-maxlik)
		time <- time_df %>% filter(run==name) %>% pull(time)
		
		summary <- c(run=name,time=time,loglik=maxlik,k=k,aic=aic)
	}
	
        write_csv(accum,paste0(path_name,"/folders_for_fit/",name,"/",type,".csv"))
	return(summary)
    })
    return(summary)
})

summary_data <- bind_rows(summary)
write_csv(summary_data,paste0(path_name,"/folders_for_fit/summary.csv"))

options(warn = oldw)
