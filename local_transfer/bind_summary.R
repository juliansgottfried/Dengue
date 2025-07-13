tmp <- read.csv("tmp_summary.csv")

if (file.exists("summary.csv")) {
	main <- read.csv("summary.csv")
	main <- rbind(main,tmp)
} else {
	main <- tmp
}

write.csv(main,"summary.csv")
