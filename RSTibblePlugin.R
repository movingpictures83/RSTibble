library(Matrix)
library(glmnet)
library(tidyverse)
library(bestNormalize)

library(tidyverse)
library(reshape2)
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }


x <<- read.table(paste(pfix, parameters["x", 2], sep="/"), row.names=1, header=TRUE, sep="\t")
y <<- read.table(paste(pfix, parameters["y", 2], sep="/"), row.names=1, header=TRUE, sep="\t")
p.fac <<- read.table(paste(pfix, parameters["pfac", 2], sep="/"), row.names=1, header=TRUE, sep="\t")
}

run <- function() {}

output <- function(outputfile) {
(BNobject <- bestNormalize(as.matrix(y), quiet=T))
fit <- glmnet(as.matrix(x), BNobject$x.t, penalty.factor=as.numeric(unlist(p.fac)), family="gaussian", alpha=1)
y.t <- BNobject$x.t
coef <- as.matrix(coef(fit))
colnames(coef) <- apply(abs(fit$beta), 2, sum) #set L1 norm as the header
coef <- coef[,round(as.numeric(colnames(coef)), digits=2)<=0.75]
coef <- coef[rowSums(coef)!=0,]
coef <- coef[order(coef[,dim(coef)[2]], decreasing=T),]
coef <- coef[-1,]
x <- as_tibble(x, rownames="RefSeq")
y <- as_tibble(y, rownames="RefSeq")
fit$dev
write.table(coef, paste(outputfile, "coef.txt", sep="/"), sep="\t", row.names=TRUE, col.names=NA, quote=F)

df <- as_tibble(melt(coef))
colnames(df) <- c('RBP', 'norm', 'coefficient')
df %>% ggplot(aes(x=norm, y=coefficient, group=RBP, color=RBP)) + geom_line(size=1.5) + theme_bw(12) +
        xlab("L1 Norm") + ylab("Coefficient") +
        theme(text = element_text(size=20), panel.grid.minor = element_blank(), panel.grid.major = element_blank())
ggsave(paste(outputfile, "coef.png", sep="/"), width=4.5, height=3)
}
