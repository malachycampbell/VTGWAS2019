###########
##RR BLUP##
###########

genos <- read.table("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day1Ex/NumericGenos.txt", sep = "\t", header = T)
genos <- genos[order(genos$chr, genos$bp) ,]
phenos <- read.csv("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day1Ex/Phenotypes.csv")

MAP <- genos[1:3]
genos <- genos[4:ncol(genos)]
genos <- as.matrix(genos)

#Check if is a matrix 
is.matrix(genos)

#Set the row names of the matrix as the marker names
row.names(genos) <- MAP$rs

#Run RR-BLUP
rrblup.mod <- mixed.solve(y = phenos$Y, Z = t(genos), K = NULL, X = NULL, method="REML", SE=F)

SNPeff <- rrblup.mod$u

#Create dataframe for manhattan plot
rrBLUP.mrkeffs <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(SNPeff))

############################
###Manhattan Plot Function##
############################

#Create manhattan plot for marker effects
pdf("~/Desktop/RRBLUP.pdf", h = 5, w = 7)
manhattan.Beta(rrBLUP.mrkeffs, Title = "rrBLUP")
dev.off()

########
##BGLR##
########

##Loading and formatting data
genos <- read.table("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day1Ex/NumericGenos.txt", sep = "\t", header = T)
genos <- genos[order(genos$chr, genos$bp) ,]
phenos <- read.csv("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day1Ex/Phenotypes.csv")

MAP <- genos[1:3]
genos <- genos[4:ncol(genos)]
genos <- as.matrix(genos)

#Check if is a matrix 
is.matrix(genos)

#Set the row names of the matrix as the marker names
row.names(genos) <- MAP$rs

genos <- t(genos + 1)

#######
##BRR##
#######

ETA <- list(list(X = genos, model = "BRR"))

Y <- scale(phenos$Y, center = T, scale = F)

#Running BGLR
BRRmodel <- BGLR(y = Y, ETA = ETA, 
                 nIter = 1500, burnIn = 150, 
                 thin = 5, 
                 saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BRRmodel.1500.150_")

#Plot trace plots for BGLR
resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BRRmodel.1500.150_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BRRmodel$burnIn/BRRmodel$thin, col = "black", lty = 3)

str(BRRmodel)

BRRmodel.mrkeff <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BRRmodel$ETA[[1]]$b))


pdf("~/Desktop/BRR.pdf", h = 5, w = 7)
manhattan.Beta(BRRmodel.mrkeff, Title = "BRR")
dev.off()


##Second run with more iterations
Y <- scale(phenos$Y, center = T, scale = F)

BRRmodel2 <- BGLR(y = Y, ETA = ETA, 
                 nIter = 12000, burnIn = 2000, 
                 thin = 5, saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BRRmodel.12000.2000_")

#Trace plots for second run
resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BRRmodel.12000.2000_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BRRmodel$burnIn/BRRmodel$thin, col = "black", lty = 3)

#Manhattan plot
BRRmodel.mrkeff2 <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BRRmodel2$ETA[[1]]$b))

pdf("~/Desktop/BRR2.pdf", h = 5, w = 7)
manhattan.Beta(BRRmodel.mrkeff2, Title = "BRR")
dev.off()

#########
##BayesA#
#########

ETA <- list(list(X = genos, model = "BayesA"))

Y <- scale(phenos$Y, center = T, scale = F)

BayesAmodel <- BGLR(y = Y, ETA = ETA, nIter = 12000, burnIn = 2000, thin = 5, saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesAmodel.12000.2000_")

resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesAmodel.12000.2000_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BayesAmodel$burnIn/BayesAmodel$thin, col = "black", lty = 3)

BayesAmodel.mrkeff <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BayesAmodel$ETA[[1]]$b))

pdf("~/Desktop/BayesA.pdf", h = 5, w = 7)
manhattan.Beta(BayesAmodel.mrkeff, Title = "BayesA")
dev.off()

###########
##Bayes B##
###########

ETA <- list(list(X = genos, model = "BayesB"))
Y <- scale(phenos$Y, center = T, scale = F)

BayesBmodel <- BGLR(y = Y, ETA = ETA, 
                    nIter = 12000, burnIn = 2000, 
                    thin = 5, saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesBmodel.12000.2000_")

resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesBmodel.12000.2000_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BayesBmodel$burnIn/BayesBmodel$thin, col = "black", lty = 3)

BayesBmodel.mrkeff <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BayesBmodel$ETA[[1]]$b))

pdf("~/Desktop/BayesB.pdf", h = 5, w = 7)
manhattan.Beta(BayesBmodel.mrkeff, Title = "BayesB")
dev.off()

###########
##Bayes B##
###########

ETA3 <- list(list(X = genos, model = "BayesB", counts = 8, probIn = 0.368))
Y <- scale(phenos$Y, center = T, scale = F)

BayesBmodel3 <- BGLR(y = Y, ETA = ETA3, 
                    nIter = 12000, burnIn = 2000, 
                    thin = 5, saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesBmodel.12000.2000_")

resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesBmodel.12000.2000_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BayesBmodel$burnIn/BayesBmodel$thin, col = "black", lty = 3)

BayesBmodel.mrkeff3 <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BayesBmodel3$ETA[[1]]$b))

pdf("~/Desktop/BayesB3.pdf", h = 5, w = 7)
manhattan.Beta(BayesBmodel.mrkeff3, Title = "BayesB")
dev.off()

##########
##BayesC##
##########

ETA <- list(list(X = genos, model = "BayesC"))

Y <- scale(phenos$Y, center = T, scale = F)

BayesCmodel <- BGLR(y = Y, ETA = ETA, 
                    nIter = 12000, burnIn = 2000, 
                    saveAt = "~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesCmodel.12000.2000_")

resVar <- scan("~/Documents/Dropbox/Work/Presentations/VT_GWAS_workshop/Day2Ex/BGLR/BayesCmodel.12000.2000_varE.dat")
plot(resVar, type = 'o', col = "red", cex = 0.5, ylab = expression(var[e]))
abline(v = BayesCmodel$burnIn/BayesCmodel$thin, col = "black", lty = 3)

BayesCmodel.mrkeff <- data.frame(CHR = MAP$chr, BP = MAP$bp, Beta = abs(BayesCmodel$ETA[[1]]$b))

pdf("~/Desktop/BayesC.pdf", h = 5, w = 7)
manhattan.Beta(BayesCmodel.mrkeff, Title = "BayesC")
dev.off()
