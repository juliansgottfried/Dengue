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

lapply(data,function(type) {
    print(paste0(toupper(type)))

    paths<-list.files(paste0(path_name,"/out/",type),full.names=T)
    names<-list.files(paste0(path_name,"/out/",type),full.names=F)

    map2(paths,names,function(path,name) {
        print(paste0(name))

	source(paste0(path,"object.R"))

        files<-list.files(path,full.names=T)
        accum<-process(files[1])

        for (i in 2:length(files)) {
		add<-process(files[i])
		accum<-bind_rows(accum,add)
	}

        write.csv(accum,paste0(path_name,"/folders_for_fit/",name,"/",type,".csv"),row.names=F)
    })
})
