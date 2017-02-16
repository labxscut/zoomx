#!/usr/bin/env Rscript
library("RColorBrewer")
library("ggplot2")

setwd("~/ws/10x6")
files = list.files(pattern="mcv$"); data = lapply(files, "read.table")
gridSize = c(100000,100000000,10000,10000000,1000000)
averageFragLength = c(24444)
fragNumber = c(22650985)
coveredGenome = rep(3200000000,5) #length hg19
#coveredGenome = rep(3137161264,5) #length hg19
#of frags per 10K, 100K, 1M, 10M, 100M
cov = lapply(seq_along(data), FUN=function(x){ head(data[[x]]); data[[x]][,4][data[[x]][,4]!=0] } )
#cov = lapply(seq_along(data), FUN=function(x){ head(data[[x]]); data[[x]][,4][data[[x]][,4]!=0]*averageFragLength/gridSize[x] } )

labels = files
coveragePerBlock = round(fragNumber / (coveredGenome/gridSize))  # divide by 2 because of half overlap sliding windows
titleText = element_text(size=22, family="Arial", face="bold")
labelText = element_text(size=22, family="Arial")
plotMargin = unit(c(1.5,1.5,1.5,1.5), "cm")
panelMargin = unit(100, "lines")

plotDenCov = function(i, cov, files, coveragePerBlock){
  Coverage = cov[[i]]; prefix=files[[i]]; cpb=coveragePerBlock[[i]]
  cat("plotting", i, "-th figure", paste(prefix, "png", sep="."), "\n")
  xlim=c(0,2*cpb)
  ylim=c(0, 1.25*max(dpois(seq(1,2*cpb),lambda=cpb)))
  cat("theory: Poisson lamba=", cpb, "\n")
  cat("xlim:", xlim, "\n")
  cat("ylim:", ylim, "\n")
  png(paste(prefix, "png", sep="."), h=750, w=750, pointsize=20)
  p0 = qplot(Coverage, geom = 'blank') +   
    geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +
    #stat_function(fun = function(x, lambda) { dpois(round(x), lambda=lambda) }, aes(colour = 'Poisson'), args=list(lambda=cpb)) +                       
    geom_histogram(aes(y = ..density..), alpha = 0.4) +
    geom_line(aes(x = seq(min(xlim), round(max(xlim))), y = dpois(seq(min(xlim), round(max(xlim))), lambda=cpb), colour = 'Poisson')) +
    scale_colour_manual(name = '', values = c('red', 'blue')) + 
    xlim(xlim) + ylim(ylim) + ggtitle("Genome Molecule Coverage (Density)") + ylab("Density") +
    theme(legend.position = c(0.85, 0.85), plot.title=titleText, text=labelText, plot.margin=plotMargin, panel.margin=panelMargin)
  print(p0)
  dev.off()
}

plotDenCov(3,cov,files,coveragePerBlock)

lapply(seq_along(cov), "plotDenCov", cov, files, coveragePerBlock)

#quit()

library(MASS)

testDenCov = function(i, cov, files, coveragePerBlock) {
  cv = cov[[i]]; prefix=files[[i]]; cpb=coveragePerBlock[[i]]
  my.mle<-fitdistr(cv, densfun="poisson")
  BIC(my.mle)
  cv_select = seq(round(0.5*median(cv)),round(1.5*median(cv))) #exclude outliers
  tcv = cv[ cv %in% cv_select ]
  tb = as.data.frame(table(tcv))
  tb$tcv = as.integer(as.character(tb$tcv))
  tb$Freq = as.integer(as.character(tb$Freq))
  tcv_max = max(tb$tcv)
  tcv_min = min(tb$tcv)
  for(j in seq(tcv_min, tcv_max)) if(! j %in% tb$tcv) tb=rbind(tb, data.frame(tcv=j,Freq=0))
  tb$Expect = dpois(tb$tcv, median(tcv)) # or mean of (tcv)
  tb$P = tb$Expect/(sum(tb$Expect))
  chisq.test(tb$Freq, p=tb$P, simulate.p.value=TRUE)
  cat("testing", i, "-th figure", paste(prefix, "png", sep="."), "\n")
  #ks.test(cv+rnorm(length(cv))/1000, "ppois", mean(cv))
  cat("theory: Poisson lamba=", cpb, "\n")
  cat("xlim:", xlim, "\n")
  cat("ylim:", ylim, "\n")
}


# plotDenCov = function(i, cov, files, coveragePerBlock){
#   cols <- brewer.pal(3, "Set1")
#   cv = cov[[i]]; prefix=files[[i]]; cpb=coveragePerBlock[[i]]
#   cat("plotting",paste(prefix, "png", sep="."),"\n")
#   png(paste(prefix, "png", sep="."), h=1000, w=1000, pointsize=20)
#   hp = hist(cv, breaks=100, ylim=c(0,0.01), xlim=c(0,2*cpb), prob=T, col=cols[1])            # prob=TRUE for probabilities not counts
#   print(hp)
#   poix = seq(min(hp$breaks),max(hp$breaks)) #hp$breaks
#   poiy = dpois(poix,lambda=cpb) #cpbp coveragePerBlock
#   #lines(density(cv), col=cols[2])             # add a density estimate with defaults
#   lines(density(cv, adjust=2), lty="dotted", ylim=hp$ylim*2, col=cols[2])   # add another "smoother" density
#   #points(poix, poiy, col=cols[3])
#   lines(poix, poiy, lty="dotted", col=cols[3])   # add another "smoother" density
#   # Add a legend using the nice sample labeles rather than the full filenames.
#   #legend("topright", legend=labels, col=cols, lty=1, lwd=4)
#   dev.off()
# }
# 
# lapply(seq_along(cov), "plotDenCov", cov, files, coveragePerBlock)

#ggplot2
# create some data to work with

x = rnorm(1000);

# overlay histogram, empirical density and normal density
p0 = qplot(x, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
  theme(legend.position = c(0.85, 0.85))

print(p0)