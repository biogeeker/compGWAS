args<-commandArgs(trailingOnly=TRUE)
library(foreach)
library(doParallel)

CHI2 <- function(i) {
	if (args[6]=="CDS") {
		freq <- matrix(as.numeric(chi2d[i,c(2,3,4,5)]),2,2)
	} else if (args[6]=="SNP") {
		freq <- matrix(as.numeric(chi2d[i,c(3,4,5,6)]),2,2)
	} else if (args[6]=="nonCDS") {
		freq <- matrix(as.numeric(chi2d[i,c(4,5,6,7)]),2,2)
	}
	chiT <- chisq.test(freq,correct=TRUE)
	res_freq[i,1:5] <- c(as.numeric(freq[1,1]),as.numeric(freq[2,1]),as.numeric(freq[1,2]),as.numeric(freq[2,2]),as.numeric(chiT$statistic))
	freq[freq==0] <- 1
	Odds <- ((freq[1,1] * freq[2,2]) / (freq[1,2] * freq[2,1]))
	res_freq[i,6:7]=c(as.numeric(Odds),(chiT$p.value))
	res_freq[i,]
}

if (args[6]=="CDS") {
	chi2d <- read.table(args[1],header=FALSE,sep="\t",as.is=TRUE);
	chi2d <- chi2d[order(chi2d[,1]),];
	INDEX <- chi2d[,1];
} else if (args[6]=="SNP") {
	chi2d <- read.table(args[1],header=FALSE,sep="\t",as.is=TRUE);
	chi2d <- chi2d[order(chi2d[,1],chi2d[,2]),];
	INDEX <- paste(chi2d[,1],chi2d[,2],sep="_!T!E!S!T!_");
	Ctig <- matrix(unlist(strsplit(as.vector(INDEX),"_!T!E!S!T!_")),2,length(INDEX))[1,];
	LOC <- matrix(unlist(strsplit(as.vector(INDEX),"_!T!E!S!T!_")),2,length(INDEX))[2,];
} else if (args[6]=="nonCDS") {
	chi2d <- read.table(args[1],header=FALSE,sep="\t",as.is=TRUE);
	chi2d <- chi2d[order(chi2d[,1],chi2d[,2],chi2d[,3]),];
	INDEX <- paste(chi2d[,1],chi2d[,2],sep="_!T!E!S!T!_");
	Ctig <- matrix(unlist(strsplit(as.vector(INDEX),"_!T!E!S!T!_")),2,length(INDEX))[1,];
	LOC <- matrix(unlist(strsplit(as.vector(INDEX),"_!T!E!S!T!_")),2,length(INDEX))[2,];
}

rowNum <- nrow(chi2d);
res_freq <- matrix(0,rowNum,7);
rownames(res_freq) <- INDEX;
colnames(res_freq) <- c(args[2],args[3],args[4],args[5],"CHI_squared","OR","CHI2_Pvalue");

cl <- makeCluster(as.numeric(args[7]), type="FORK");
registerDoParallel(cl);
res <- foreach(i=1:rowNum,
              .combine=rbind,
              .export = c("chi2d","res_freq")) %dopar% {
       CHI2(i)
}
stopImplicitCluster()

res_all <- data.frame(res,stringsAsFactors=FALSE);

if (args[6]=="CDS") {
	Loc_data <- data.frame(Loc=INDEX,stringsAsFactors=FALSE);
	res_all <- cbind(Loc_data,res_all);
	for (i in 2:ncol(res_all)) {res_all[,i] <- signif(as.numeric(res_all[,i]),3)}
} else if (args[6]=="SNP") {
	Loc_data <- data.frame(RefName=Ctig,Loc=LOC,stringsAsFactors=FALSE);
	res_all <- cbind(Loc_data,res_all);
	for (i in 3:ncol(res_all)) {res_all[,i] <- signif(as.numeric(res_all[,i]),3)}
} else if (args[6]=="nonCDS") {
	Loc_data <- data.frame(RefName=Ctig,Loc=LOC,Sub=chi2d[,3],stringsAsFactors=FALSE);
	res_all <- cbind(Loc_data,res_all);
	for (i in 4:ncol(res_all)) {res_all[,i] <- signif(as.numeric(res_all[,i]),3)}

}

options(max.print=.Machine$integer.max,width=9999,digits = 2)

write.table(res_all,args[8],row.names = FALSE,col.names=TRUE,quote = FALSE, sep = "\t");


