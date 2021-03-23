## Set Working Directory
## Make sure data file is located in working directory

#install.packages(c("devtools","Hmisc", "corrplot", "dplyr", "vegan", "reshape2"))

## Turn on packages:
library(ggplot2)
library(devtools)
library(dplyr)


## import data via csv (make sure data is in working directory):
data.BM=read.csv("2018_Beach_Monitoring.csv", header = T)
## check the dimensions of dataset
## Data should have 517 rows, 22 columns
dim(data.BM)
## Look at data.BM dataset, notice there are null or blanks
## remove null data/black cells, and dataset should have 40 observations:
data.BM<-na.omit(data.BM)
dim(data.BM) #Check dimensions

## Check the columns------
## Microcystin is identified as a factor
## Convert this to numeric
data.BM$Microcystin=as.numeric(paste(data.BM$Microcystin))


## Is the data normally distributed? ------
## null hypothesis: data is normally distributed, 
## if p-value <0.05, then data is not normally distributed
## In such cases, certain statistical methods may not be appropriate
## and can lead to inaccurate interpretations
shapiro.test(data.BM$TDFe)
shapiro.test(data.BM$Microcystin)

## We can also test normality using K-S test
ks.test(data.BM$TDFe, "pnorm", 1, 2) #mean =1, 2= std dev
## Can also visually check for normal distribution
qqPlot(data.BM$Microcystin) #you can also try other variables


## Scale data to max and min-----
## Data values vary in range (from 0 to 100+)
## The function to normalize data is (x - min(x))/(max(x) - min(x))
## First, normalize env't dataset using apply() and converting to a dataframe
data_norm<-as.data.frame(apply(data.BM[, 9:19], 2, 
                                 function(x)(x - min(x))/(max(x)-min(x))))
### Research Question-----
## What environmental conditions associate with microcystin?

## Compute p-value and Spearman correlation-------
## Often, Spearman correlation is used for environmental samples
## Because environmental samples are not normally distributed
## Pearson correlation is used when dataset is normally distributed
library(Hmisc) #use this package for correlations
res.corr<-rcorr(as.matrix(data_norm[,c(1:11)], type = "spearman"))
res.corr #See results

## Graphing correlation plot------
## Extract p-value from results:
pval<-res.corr$P
## Extract Spearman correlation coefficient from results:
rval<-res.corr$r

## Plot correlation with p-value:
library(corrplot) #use this package to plot correlation
corr.plot<-corrplot(rval, type = "upper", order ="FPC",
                     p.mat = pval, sig.level = 0.05,
                     tl.col = "black", tl.srt = 45, tl.cex = 1.2,
                     cl.cex = 1.2, insig = "blank")
## tl = text label, cl = color label, col = columns, .cex magnifies the text size
## p-mat = p-value, type ??corrplot for more info

## plot correlation using ggplot
## Note: this is a bit more complicated with coding 
## but feel free to explore
#install.packages("reshape2")
library(reshape2)
res.corr2<-round(cor(data_norm[,c(1:11)], method = "spearman"),2) #round 2 decimal places
melted_corr<-melt(res.corr2) #rthe melt() will earrange correlation data for plotting
P1<-ggplot(data = melted_corr, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + scale_fill_gradient(low = "blue", high = "red")
P1

## Note that, a correlation matrix has redundant information. 
## We'll use the function below to grab the upper portion of
## matric and set half of it to NA
 get_upper_tri<-function(res.corr2){
    res.corr2[lower.tri(res.corr2)]<- NA
    return(res.corr2)
    }
upper_tri<-get_upper_tri(res.corr2) #grabs the upper portion of matrix

## heatmap of upper portion of matrix
P2<-ggplot(data = melt(upper_tri, na.rm = TRUE), aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile(colour = "white") + 
  scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                      midpoint = 0, limit = c(-1,1), space = "Lab",
                      name = "Spearman Rank Correlation") + #correlation coefficients range from -1 to 1.
  theme(panel.grid.major = element_blank (),panel.grid.minor = element_blank(),
              panel.background = element_rect (fill = "white"),text = element_text(size = 14),
              axis.line = element_line(colour = "black"),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 14), 
              axis.text.y = element_text(size = 14),
              axis.text.x = element_text(angle = 45, size = 14, hjust = 1)) +
  coord_fixed() #ensures that one unit on x-axis is same as y-axis

P2

## Add r value to heatmap
P3<-P2 +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4)
P3


##Research Question ------
##What environmental conditions associate with microcystin?
## NMDS uses rank orders to observe how species or composition
## change from one community to the next
## k can be increased if stress is too high
library(vegan)
mds.out=metaMDS(data_norm[,1:11],distance="bray", k=3)

## Put NMDS output into a dataframe------
data.sm<-as.data.frame(scores(mds.out))
data.sm

## Add columns from data.BM to mds.out dataframe
data.sm$Date=as.character(data.BM$Date)
data.sm$Site<-as.character(data.BM$Site) #alternative way to add 
data.sm$Month=as.character(data.BM$Month)

## Plot NMDS using ggplot-----
P4<-ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Month)) + 
  geom_point(size=3) + xlim(-0.6, 0.6) + ylim(-0.6, 0.6)+
  geom_hline(yintercept=0.0, colour="grey", lty=2)+
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  stat_ellipse() +
  theme(legend.background = element_rect(fill="white", size=0.3, 
                                         linetype="solid", colour="black"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.key = element_rect(fill = "white"))
P4

## adding arrows
env.monitoring=envfit(mds.out, data.BM[,9:19], perm=1001)
env.monitoring

arrow.monitoring=as.data.frame(scores(env.monitoring,display="vectors"))
arrow.monitoring
arrow.monitoring=cbind(arrow.monitoring,Species=rownames(arrow.monitoring))

##Plotting with arrows
P5<-ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2)) + 
  geom_point(aes(colour=Month), size=3) + 
  xlim(-0.8, 0.8) + ylim(-0.8, 0.8) +
  geom_hline(yintercept=0.0, colour="grey", lty=2) +
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  scale_color_brewer(palette="Set1") + theme(legend.position = c(0.9, 0.7)) +
  theme(legend.background = element_rect(fill="white", size=0.3, 
                                         linetype="solid", colour="black"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=14), axis.title = element_text(size = 14),
        legend.key = element_rect(fill = "white")) +
  geom_segment(data=arrow.monitoring, aes(x=0, xend=2*NMDS1, y=0, yend=2*NMDS2), 
               arrow=arrow(length=unit(0.25,"cm")))+
  geom_text(data=arrow.monitoring, aes(x=2.2*NMDS1, y=2.2*NMDS2, label=Species), size=5)
P5

## MRPP & ANOSIM

envt=as.matrix(data_norm[,2:11])
Microcystin=as.numeric(data_norm[,1])


mrpp(envt, group=Microcystin, distance="bray")
anosim(envt, grouping=Microcystin, permutations=999, distance="bray")
 
