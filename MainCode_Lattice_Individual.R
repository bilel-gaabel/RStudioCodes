
rm(list=ls())
library(lme4)


# Working directory     (change "\" by "/")
dir_root <- "D:/2016/BIGM/Students/Sara"

# file input name
file_name <- "all analyses 1 annee.csv"

setwd(dir_root)
datos <- read.csv(file_name)
head(datos)
tail(datos)
str(datos)

datos$Rep <- factor(datos$Rep)
datos$block <- factor(datos$block)
datos$Geno <- factor(datos$Geno)

boxplot(area_mm2~Location,data=datos,las=2)
boxplot(area_mm2~Geno,data=datos,las=2)
boxplot(area_mm2~Rep+Location,data=datos,las=2)

envirs <- levels(datos$Location)
envirs
  

MEANSALL <- NULL

	
for(i in 1:length(envirs)){
	envir <- envirs[i]

	# Select just those line belonging to the selected site
	datos2 <- droplevels(subset(datos, datos$Location==envirs[i]))

	statis <- data.frame(matrix(NA,nrow=7,ncol=1))
	rownames(statis) <- c("n Reps","Error Var","Genotypic Var","Heritability","Grand Mean","LSD","CV")
	
	
	  
	 #  Compute BLUPs   (Entry as random effect)		
	fm <- lmer(area_mm2~(1|Geno)+(1|Rep)+(1|block:Rep),data=datos2)

	fixEffect <- as.matrix(fixef(fm))
	randEffect <- as.matrix(ranef(fm))
	interc <- fixEffect[1,1]
	EntryEffect <- randEffect[["Geno",1]]
    
	BLUP <- interc+EntryEffect
	colnames(BLUP) <- "BLUP"	
       
	#  Compute statistics
	varcorr <- VarCorr(fm)
	varErr <- attr(varcorr,'sc')^2
	varG <- as.vector(varcorr$'Geno')
	#varGE <- as.vector(varcorr$'Entry:Loc')
	#varLoc <- as.vector(varcorr$'Loc')
	nRep <- max(as.numeric(as.character(datos2$Rep)))
	#nLoc <- length(unique(as.character(datos.tmp$Loc)))
	LSD <- 1.96*sqrt(varErr)

	#  Compute BLUE (Entry as fixed effect)
	fm <- lmer(area_mm2~0+Geno+(1|Rep)+(1|block:Rep),data=datos2)
				  
	BLUE <- fixef(fm)
	names(BLUE) <- substr(names(BLUE),5,nchar(names(BLUE)))
	if(sum(rownames(BLUP)!=names(BLUE))>0)  stop("names BLUE and names BLUP dont match")

	results <- cbind(Entry=rownames(BLUP),BLUP,BLUE)
	rownames(results) <- NULL
	results$Location <- envirs[i]

	# Calculate h2 and fill the table of statistics
	media <- mean(BLUE,na.rm=TRUE)
	CV <- 100*sqrt(varErr)/abs(media)

	h2 <- varG/(varG + varErr/nRep)
	statis[,i] <- round(c(nRep,varErr,varG,h2,media,LSD,CV),4)
  
MEANSALL <- rbind(MEANSALL, results)
}
