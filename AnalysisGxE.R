
rm(list=ls())
library(lme4)

# Working directory     (change "\" by "/")
dir_root <- "C:/Users/bilel/Desktop/R workspace sana/R Code From Zakaria"

# file input name
file_name <- "all analyses 1 annee.csv"
setwd(dir_root)
datos <- read.csv(file_name)
head(datos)
tail(datos)
str(datos)
sapply(datos, is.numeric)
#######################################
datos$Rep <- factor(datos$Rep)
datos$block <- factor(datos$block)
datos$Geno <- factor(datos$Geno)

boxplot(area_mm2~Location,data=datos,las=2)

envirs <- levels(datos$Location)
envirs

statis <- matrix(NA,nrow=9,ncol=1)
rownames(statis) <- c("n Locs","n Reps","Error Var","Genotypic Var","GenxEnv Var","Heritability","Grand Mean","LSD","CV")


#  Compute BLUP   (Entry as random effect)		
fm <- lmer(area_mm2~-1+Location+(1|Geno)+(1|Rep:Location)+(1|block:Rep:Location)+(1|Location:Geno),data=datos)

fixEffect <- as.matrix(fixef(fm))
rownames(fixEffect) <- envirs
 
randEffect <- as.matrix(ranef(fm))
#interc <- fixEffect[1,1]
EntryEffect <- randEffect[["Location:Geno",1]]
EntryEffect$intercept <- rep(fixEffect,each=24)

BLUP <- EntryEffect$intercept+EntryEffect$`(Intercept)`
BLUP <- as.matrix(BLUP, ncol = 1)
colnames(BLUP) <- "BLUP"

#  Compute statistics
varcorr <- VarCorr(fm)
varErr <- attr(varcorr,'sc')^2
varG <- as.vector(varcorr$'Geno')
varGE <- as.vector(varcorr$'Location:Geno')
#varLoc <- as.vector(varcorr$'Location')
nRep <- max(as.numeric(as.character(datos$Rep)))
nLoc <- length(unique(as.character(datos$Loc)))
LSD <- 1.96*sqrt(varErr)

#  Compute BLUE (Entry as fixed effect)
fm <- lmer(area_mm2~-1+Geno:Location+(1|Rep:Location)+(1|block:Rep:Location),data=datos)

BLUE <- fixef(fm)
names(BLUE) <- substr(names(BLUE),6,nchar(names(BLUE)))
if(sum(rownames(BLUP)!=names(BLUE))>0)  stop("names BLUE and names BLUP dont match")

results <- cbind(Entry=rownames(BLUP),BLUP,BLUE)
rownames(results) <- rownames(EntryEffect)
results

write.csv(results, 'blue_blup_GxE.csv')

media <- mean(BLUE,na.rm=TRUE)
CV <- 100*sqrt(varErr)/abs(media)
h2 <- varG/(varG + varGE/nLoc + varErr/(nRep*nLoc))
statis[,1] <- round(c(nLoc,nRep,varErr,varG,varGE,h2,media,LSD,CV),4)
statis2 <- as.data.frame(statis)
statis2
write.csv(statis, 'statisticsGxE.csv')
