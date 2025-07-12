library(tidyverse)

log_name <- Sys.getenv("LOGNAME")
path_name <- paste0("/scratch/",log_name,"/Dengue")

time_df = read.table(paste0(path_name,"/times.txt"),sep=" ") |> 
    t() |> 
    na.omit() |> 
    t() |> 
    data.frame()
colnames(time_df)=c("time","name")

time_df=time_df |>
    mutate(time=hms(time)) |>
    mutate(hours=hour(time),minutes=minute(time)) |>
    mutate(time=hours*60+minutes) |>
    select(-c(hours,minutes)) |>
    group_by(name) |>
    summarize(time=max(time)) |>
    ungroup()

process<-function(path) {
	read_csv(path) %>%
		select(matches(c("sample",par_names,"loglik","loglik.se","flag",
			"time","cond","eff",
			"iter","run"))) %>%
		mutate_all(as.numeric)
	}

summary_data <- c(rep(NA,5))
names(summary_data) <- c("run","time","loglik","k","aic")

source(paste0(path,"object.R"))
pars_path <- paste0(path,"pars.csv")
init_vals <- read_csv(pars_path)
est_pars <- par_names[apply(init_vals,2,function(x) {max(x)-min(x)})!=0]

result_type <- c("results","traces","stats")

lapply(result_type,function(type) {
    print(paste0(toupper(type)))

    paths<-list.files(paste0(path_name,"/out/",type),full.names=T)
    names<-list.files(paste0(path_name,"/out/",type),full.names=F)

    map2(paths,names,function(path,name) {
        print(paste0(name))

        files<-list.files(path,full.names=T)
        accum<-process(files[1])

        for (i in 2:length(files)) {
		add<-process(files[i])
		accum<-bind_rows(accum,add)
	}

	if (type=="results") {
		source(paste0(path,"object.R"))
		init_vals <- read_csv(paste0(path,"pars.csv"))

		k <- length(par_names[apply(init_vals,2,function(x) {max(x)-min(x)})!=0])
		maxlik <- accum %>% arrange(-loglik) %>% pull(loglik) %>% head(1)
		aic <- 2*(k-maxlik)
		time <- time_df %>% filter(name==name) %>% pull(time)
		
		current <- c(name,time,maxlik,k,aic)
		names(current) <- names(summary_data)

		summary_data <- bind_rows(summary_data,current)	
		summary_data <- summary_data %>% na.omit(cols=run)
	}

        write_csv(accum,paste0(path_name,"/folders_for_fit/",name,"/",type,".csv"))
    })
})

write_csv(summary_data,paste0(path_name,"/summary.csv"))
