
nPop <- 3;
sZero <- "ab";
sOne <- "Am";
bDrawOnly <- FALSE;
sOutPrefix <- "27_73_unimodal";
sOutFile <- "states.out.txt";
sGraphPDF <- "states.out.pdf";
nGraphW <- 14; #in inches
nGraphH <- 7; 
sFilePrefix <- paste(sOutPrefix, "/psi/", sep="");
sFileSuffix <- "_markers.txt";
arrGenerations <- c(1,2, seq(10,1000,10));

#These will be passed to phenotype php
sPheOut <- paste(sFilePrefix , "/phenosum.txt",sep="");
nSex <- 1;
nColToSum <- 7;
nPheLowB <- 0;
nPheUpB <- 16;
nPheStep <- 1;

sPHPScript <- paste("php sumphenotypedistr2R.php -d \"", 
				sFilePrefix,
				"\" -o\"",
				sPheOut,
				"\" -p",
				nPop,
				" -s",
				nSex,
				" -c",
				nColToSum,
				" -l",
				nPheLowB,
				" -u",
				nPheUpB,
				" -t",
				nPheStep
				,sep="");

if (!bDrawOnly) {
cat("Analysing marker output ...\n");
for(nGen in arrGenerations) {
	cat("Gen ", nGen, "\n");
	markers <- read.table(paste(sFilePrefix , "/Gen", nGen,  sFileSuffix , sep = "" ), header=FALSE, as.is=TRUE);
	pop3markers <- markers[markers[,2]==nPop,]
	pop3markers[pop3markers==sZero] <- 0;
	pop3markers[pop3markers==sOne] <- 1;
	pop3nobound <- pop3markers[,pop3markers[1,]!=-1];

	hap1 <- pop3nobound[seq(1, nrow(pop3nobound)-1, 2),];
	hap2 <- pop3nobound[seq(2, nrow(pop3nobound)  ,2),];

	hap1 <- data.matrix(hap1[,seq(6, ncol(hap1))]) ;
	hap2 <- data.matrix(hap2[,seq(6, ncol(hap2))]);

	rownames(hap1) <- seq(1, nrow(hap1));
	rownames(hap2) <- seq(1, nrow(hap2));
	genotypes <- hap1 + hap2;

	#nHybridIndex <- (sum(hap1) + sum(hap2) ) / (nrow(hap1) * ncol(hap1) * 2);
	nHybridIndex <- (rowSums(genotypes)) / (ncol(genotypes) * 2);
	
	nP <- (colSums(hap1) + colSums(hap2)) / ( nrow(hap1)*2 );
	nQ <- 1 - nP;

	hetgenotypes <- genotypes
	hetgenotypes[hetgenotypes!=1] <-0;

	nObservedHetFreq <- colSums(hetgenotypes) / nrow(hetgenotypes);

	nFis <- 1 - nObservedHetFreq / (2 * nP * nQ);
	
	nFis[is.nan(nFis)] <- 1; #fixed, treated as "structured"
	
	oFrmToWrite <- data.frame(nGen , mean(nHybridIndex), mean(nFis), data.frame(t(as.vector(nFis))));
	write.table( oFrmToWrite, paste(sFilePrefix, "/", sOutFile,sep=""), append=TRUE, eol="\n", sep="\t", row.names = FALSE, col.names = FALSE);

	oFrmToWrite <- data.frame(nGen , mean(nHybridIndex), data.frame(t(as.vector(nHybridIndex))));
	write.table( oFrmToWrite, paste(sFilePrefix, "/hybridindex_", sOutFile,sep=""), append=TRUE, eol="\n", sep="\t", row.names = FALSE, col.names = FALSE);
	
	#write.table( paste(nGen, nHybridIndex, mean(nFis), sep = "\t"  ), sOutFile, append=TRUE, eol="", sep="" , row.names = FALSE, col.names = FALSE );
	#write.table( as.vector(nFis), sOutFile, append=TRUE, eol="\t", sep="\t", row.names = FALSE, col.names = FALSE);
}
}

cat("Plotting Fis ...\n");
plotdata <- read.table(paste(sFilePrefix, "/", sOutFile,sep=""),header=FALSE);
pdf(paste(sFilePrefix, "/","Fis_",sGraphPDF,sep="") , width=nGraphW, height=nGraphH);
oplot <- boxplot(t(plotdata[,seq(4,ncol(plotdata))]), names=plotdata[,1], range=0, xlab="Generation", ylab="FIS", lwd=1,ylim=c(-0.1,1))
lines(oplot$stats[3,], lwd=2, col="blue",lty=1)
lines(plotdata[,3], lwd=2, col="red",lty=1)
dev.off();

cat("Plotting Ancestry ...\n");
plotdata <- read.table( paste(sFilePrefix, "/hybridindex_", sOutFile,sep=""),header=FALSE, fill=TRUE);
pdf(paste(sFilePrefix, "/","ancestry_",sGraphPDF,sep="") , width=nGraphW, height=nGraphH);
oplot <- boxplot(t(plotdata[,seq(3,ncol(plotdata))]), names=plotdata[,1], range=0, xlab="Generation", ylab="Ancestry", lwd=1,ylim=c(0,1))
lines(plotdata[,2], lwd=2, col="red", lty=1)
dev.off();

cat("Plotting Phenotype ...\n");
system(sPHPScript); 
dat <- read.csv(sPheOut, header=TRUE, row.names=1, sep="\t", check.names=FALSE);
pdf(paste(sFilePrefix, "/","phenotype_col_",nColToSum,"_",sGraphPDF,sep="") , width=nGraphW, height=nGraphH);
image(as.numeric(rownames(dat)), as.numeric(colnames(dat)), sqrt(as.matrix(dat)[,seq(1,length(dat))]), col = rainbow(256, start=.6, end=.1, alpha = 1) , xlab="Generation", ylab="Phenotype Value");

dev.off();
