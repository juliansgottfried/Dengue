library(dplyr)
library(purrr)

print("Ignoring params containing 'X'")

data<-c("results","stats")

lapply(data,function(type) {
    print(paste0("data = ",type))
    paths<-list.files(paste0("../out/",type),full.names=T)
    names<-list.files(paste0("../out/",type),full.names=F)

    map2(paths,names,function(path,name) {
        print(paste0("name = ",name))
        files<-list.files(path,full.names=T)
        accum<-read.csv(files[1])%>%mutate_all(as.numeric)
        for (i in 2:length(files)) {
		add<-read.csv(files[i])%>%mutate_all(as.numeric)
		accum<-bind_rows(accum,add)
	}
	accum<-accum%>%select(-contains("X",ignore.case=F))
        write.csv(accum,paste0("../results/",name,"/",type,".csv"),row.names=F)
        print("Success")
    })
})
