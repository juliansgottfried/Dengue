library(tidyverse)

log_name<-Sys.getenv("LOGNAME")
path_name<-paste0("/scratch/",log_name,"/Dengue")

data<-c("results","traces","stats")

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

lapply(data,function(type) {
    print(paste0(toupper(type)))

    paths<-list.files(paste0(path_name,"/out/",type),full.names=T)
    names<-list.files(paste0(path_name,"/out/",type),full.names=F)

    map2(paths,names,function(path,name) {
        print(paste0(name))

	source(paste0(path,"object.R"))
	
	init_vals <- read_csv(paste0(path,"pars.csv"))
	k <- length(par_names[apply(init_vals,2,function(x) {max(x)-min(x)})!=0])

        files<-list.files(path,full.names=T)
        accum<-process(files[1])

        for (i in 2:length(files)) {
		add<-process(files[i])
		accum<-bind_rows(accum,add)
	}

	if (type=="results") {
		maxlik <- accum %>% arrange(-loglik) %>% pull(loglik) %>% head(1)
		aic <- 2*(k-maxlik)
		current <- c(name,NA,maxlik,k,aic)
		names(current) <- names(summary_data)
		summary_data <- bind_rows(summary_data,current)	
	}

	summary_data <- summary_data %>% na.omit(cols=run)
	write_csv(summary_data,paste0(path_name,"/folders_for_fit/summary.csv"))

        write_csv(accum,paste0(path_name,"/folders_for_fit/",name,"/",type,".csv"))
    })
})
