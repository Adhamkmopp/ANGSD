# ANGSD :: Bioinformatics Project I,
## <u> Extending the ANGSD summary statistics, 2017, University of Copenhagen </u> 
### Course Description

This was a course in learning C++ through extending the [ANGDS](https://github.com/ANGSD/angsd); an open source multi-threaded program for generating summary statistics of NGS data. The only modifications are in abcSum where summary statistics were extended to include distributions of raw base count and ratio by depth, base quality scores by depth and allowing the user to set a maximum depth beyond the default.

### Learning Outcomes:
1. Constructor and destructor logic and functionality
2. Object/Class construction, methods and super class inheritence
3. Pointers/Pointers-to-pointers construction, free store dealocation, decimal types and overflow

Official ANGDS repo: [ANGDS](https://github.com/ANGSD/angsd)

ANGSD Publication in BMC: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4

Aside from abcSum, R and ggplot were used to automate (as possible) the plot building process

Files with the appropriate extension are loaded to in their respective dataframes 

Base quality score by depth :: qsDep
Base quality score by forward position:: qsPosi
Base quality score by backward position :: qsIsop
Base quality counts :: qsBase
Depth counts :: sDep

```

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)


temp = list.files(pattern="*.qsDep")
ancient = temp[3]

qsDep<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t", colClasses="character"))

temp = list.files(pattern="*.qsPosi")
ancient = temp[3]

qsPosi<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t",colClasses="character"))

temp = list.files(pattern="*.qsIsop")
ancient = temp[3]

qsIsop<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t",colClasses="character"))
temp = list.files(pattern="*.qsBase")
ancient = temp[3]

qsBase<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t"))
qsMapq<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t"))
temp = list.files(pattern="*.sDepth")
ancient = temp[3]

sDepth<-  lapply(setNames(temp, make.names(temp)), function(x) 
  read.table(x, sep="\t",colClasses="character"))
subject.names <-gsub(temp, pattern=".sDepth$", replacement="")
glist <- sapply(subject.names,function(x) NULL)
```


```
drawfunction <- function(individual.names, posi,isop,baseqs,depthqs, depthbase){
  
  
  
  for(subject in 1:length(individual.names)){ 
    ################################ POSI ################################
    x<-as.data.frame(posi[subject])
    x<-cbind(x,0:99)
    for(i in 1:3){
      x[,i]<-as.numeric(x[,i])
    }
    colnames(x) <-c("Position", "Count","Quality")
    g.posi <- ggplot(x[x$Count>0,], aes( group=Position, x=Position, y = Quality, weight=Count)) +
      geom_boxplot() 
    
    counts <- table(x$Position[x$Count>0])
    counts <-as.vector(counts)
    medians<-ggplot_build(g.posi)
    medians <-medians$data[[1]]$middle
    medians<-rep(medians, counts)
    x <-cbind(x[x$Count>0,],medians)
    g.posi <- ggplot(x[x$Count>0,], aes( group=Position, x=Position, y = Quality, weight=Count, fill=medians)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_gradient(low="#05D9F6",
                          high="#5011D1",space = "Lab" )+
      labs(title= "Base Quality Score by Position (Forward)", subtitle=names(glist[subject]), x="Position of Read", y="Quality Score of Base") +
      theme_bw()
    glist[[subject]][1] <- list(g.posi)
    
    
    ################################ ISOP ################################
    x<-as.data.frame(isop[subject])
    x<-cbind(x,0:99)
    for(i in 1:3){
      x[,i]<-as.numeric(x[,i])
    }
    colnames(x) <-c("Position", "Count","Quality")
    g.isop <- ggplot(x[x$Count>0,], aes( group=Position, x=Position, y = Quality, weight=Count)) +
      geom_boxplot()
    counts <- table(x$Position[x$Count>0])
    counts <-as.vector(counts)
    medians<-ggplot_build(g.isop)
    medians <-medians$data[[1]]$middle
    medians<-rep(medians, counts)
    x <-cbind(x[x$Count>0,],medians)
    
    g.isop <- ggplot(x[x$Count>0,], aes( group=Position, x=Position, y = Quality, weight=Count, fill=medians)) +
      geom_boxplot(outlier.shape = NA) +
      scale_fill_gradient(low="#05D9F6",
                          high="#5011D1",space = "Lab" )+
      labs(title= "Base Quality Score by Position(Backward)", x="Position of Read", y="Quality Score of Base", subtitle=names(glist[subject]))+theme_bw()
    glist[[subject]][2]  <- list(g.isop)
    
    ################################ qsBASE ################################
    x<-as.data.frame(baseqs[subject])
    x<-cbind(x,0:99)
    for(i in 2:3){
      x[,i]<-as.numeric(x[,i])
    }
    colnames(x) <-c("Base", "Count","Quality")
    x<-x[x$Count>0,]
    
    
    x$Base <- factor(x$Base , levels = c("A", "T","C","G","N"))
    ratios= x %>% 
      group_by(base=x$Base) %>% 
      summarise(freq = sum(Count))%>%
      mutate(ratio = freq / sum(freq))
    
    ratios=format(round(ratios$ratio, 2), nsmall = 2)
    g.baseqs <- ggplot(x, aes( group=Base, x=Base, y = Quality, weight=Count,fill=Base)) +
      geom_violin()+
      geom_boxplot(outlier.shape = NA, width=0.065)+
      labs(title= "Base Quality Score", x="Base", y= "Quality Score of Base", subtitle=names(glist[subject]))+
      scale_fill_discrete(name="Base",
                          breaks=c("A", "T", "C", "G","N"),
                          label = c(paste("A:",ratios[1]),paste("T:",ratios[2]),paste("C:",ratios[3]),paste("G:",ratios[4]),paste("N:",ratios[5])))+
      theme(legend.background = element_rect(fill=alpha("blue",0.2),
                                             size=10, linetype="solid", 
                                             colour ="darkblue"),
            panel.background = element_rect(fill = "white", colour = "grey50"))+
      theme_bw()
    
    glist[[subject]][3]  <- list(g.baseqs)
    
    ################################ sDEPTH ################################
    x<-as.data.frame(depthbase[subject])
    for(i in 1:9){
      x[,i]<-as.numeric(x[,i])
    }
    x.c<- x[,c(1,2,3,4,5)]
    x.r <-x[,c(1,6,7,8,9)]
    colnames(x.c) <- c("Depth", "A","C","G", "T")
    colnames(x.r) <- c("Depth", "A","C","G", "T")
    
    temp <- with(x.c, colSums(x.c[Depth>50, 1:5]))
    x.c<- rbind(subset(x.c, Depth<50),temp)
    x.c$Depth[length(x.c$Depth)] = 52
    onedeep.c <- x.c %>%
      gather(Base, Count, c("A","C","G","T"), na.rm = TRUE)
    onedeep.r <- x.r %>%
      gather(Base, Count, c("A","C","G","T"), na.rm = TRUE)
    onedeep.c$Base <- factor(onedeep.c$Base , levels = c("A", "T","C","G"))
    onedeep.r$Base <- factor(onedeep.r$Base , levels = c("A", "T","C","G"))
    max.depth <- max(onedeep.c$Depth[onedeep.c$Count!=0 ])
    g.depthbase.c <- ggplot(onedeep.c[onedeep.c$Depth<max.depth+1,], aes(Depth, Count)) +   
      geom_bar(aes(fill = Base), position = "dodge", stat="identity")+
      annotate("text",x=52.5, y=max(onedeep.c$Count)/2, size=1,label="italic(Depth>50)", parse=TRUE)+
      geom_vline(xintercept = 50,linetype="dotted")+
      labs(title= "Base Count by Depth", x="Depth of Reads", y= "Base Count", subtitle=names(glist[subject]))+theme_bw()
    
    max.depth <- max(onedeep.r$Depth[onedeep.r$Count!=0 ])
    if(subject==1){
      max.depth=2500
    }
    g.depthbase.r <- ggplot(onedeep.r[onedeep.r$Depth<max.depth,]) + 
      geom_bar(aes(y =Count, x=Depth, fill= Base),stat="identity")+
      labs(title= "Base Ratio by Depth", x="Depth of Reads", y="Base Ratio", size=20, subtitle=names(glist[subject]))+theme_bw()
    
    glist[[subject]][4]  <- list(g.depthbase.c)
    glist[[subject]][5]  <- list(g.depthbase.r)
    
    ################################ qsDEP ################################
    
    x<-as.data.frame(depthqs[subject])
    x<-cbind(x,0:99)
    for(i in 1:3){
      x[,i]<-as.numeric(x[,i])
    }
    if(subject==2){
      x<-x[1:20000,]
      x <-cbind(x, rep(1:4, each=5000))
      colnames(x) <-c("Depth", "Count","Quality", "lock")
      
      x <-x[x$Count>0,]
      
      
      maxDepth = max(x$Depth)
      remainder =maxDepth%%1
      x$cut<- cut(x[,1], breaks=seq(from = 0, to = maxDepth+remainder
                                    , by = 1))
      g.dep <- ggplot(x, aes( x=cut, y = Quality, weight=Count, group=cut)) +
        geom_boxplot(outlier.shape = NA)
      z<- aggregate(data.frame(count = x$cut), list(value = x$cut), length)
      
      medians<-ggplot_build(g.dep)
      medians <-medians$data[[1]]$middle
      medians<-rep(medians,z$count)
      x <-cbind(x,medians)
      
      g.dep <- ggplot(x, aes( x=cut, y = Quality, weight=Count, group=cut, fill=medians)) +
        geom_boxplot(outlier.shape = NA)+
        scale_fill_gradient(low="#05D9F6",
                            high="#5011D1",space = "Lab" )+
        coord_cartesian(ylim = c(0, 60))+
        labs(title= "Base Quality Score by Depth", x="Depth of Reads", y="Quality Score of Base", subtitle=names(glist[subject]))+
        theme_bw()+
        facet_wrap(~lock,scales="free", ncol=1)+
        theme(axis.text.x=element_text(angle=90,hjust=1, size=6))
      
      glist[[subject]][7] <- list(g.dep)
    }else{
      x<-x[1:20000,]
      x <-cbind(x, rep(1:4, each=5000))
      colnames(x) <-c("Depth", "Count","Quality", "lock")
      x <-x[x$Count>0,]
      maxDepth = max(x$Depth)
      remainder=maxDepth%%1
      x$cut<- cut(x[,1], breaks=seq(from = 0, to = maxDepth+remainder
                                    , by = 1))
      
      g.dep <- ggplot(x, aes( x=cut, y = Quality, weight=Count, group=cut)) +
        geom_boxplot(outlier.shape = NA)
      
      z<- aggregate(data.frame(count = x$cut), list(value = x$cut), length)
      
      medians<-ggplot_build(g.dep)
      medians <-medians$data[[1]]$middle
      medians<-rep(medians,z$count)
      x <-cbind(x,medians)
      
      g.dep <- ggplot(x, aes( x=cut, y = Quality, weight=Count, group=cut, fill=medians)) +
        geom_boxplot(outlier.shape = NA)+
        scale_fill_gradient(low="#05D9F6",
                            high="#5011D1",space = "Lab" )+
        coord_cartesian(ylim = c(0, 60))+
        labs(title= "Base Quality Score by Depth", x="Depth of Reads", y="Quality Score of Base", subtitle=names(glist[subject]))+
        theme_bw()+
        facet_wrap(~lock,scales="free", ncol=1)+
        theme(axis.text.x=element_text(angle=90,hjust=1, size=6))
      
      glist[[subject]][7] <- list(g.dep)
    }
    z <- x %>% 
      group_by(Depth) %>% 
      summarise(Count = sum(Count))
    zz<-cumsum(z$Count)/sum(z$Count)
    zzz<-(z$Depth[zz>0.9][1])
    
    maxmid <-(max(x$Quality[x$Count>0 & x$Depth<(zzz)]))/2
    g.cake <- ggplot(x[x$Depth<(zzz+50),], aes(x=Depth,y=Count, fill=Quality)) + 
      geom_bar(stat="identity", colour="black")+
      labs(title= "Depth Count Distribution", x="Depth of Reads", y="Counts of Depths", subtitle=names(glist[subject]))+
      geom_vline(xintercept = zzz,linetype="dotted")+
      scale_fill_gradient2(low = "red", mid = "white",
                           high = "deepskyblue4", midpoint =25, space = "Lab")+
      xlim(0, zzz+50)+
      theme_bw()
    glist[[subject]][6] <- list(g.cake)
    
    
    
    
  }
  return(glist)
  
}

p.list <- drawfunction(glist,qsPosi,qsIsop,qsBase,qsDep,sDepth)

```
