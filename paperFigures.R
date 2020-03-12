########### Figures ###############
library(dplyr)   ## for %>% and general data frame manipulation
library(ggplot2) ## for plotting!
library(igraph)
library(MASS)
library(graphics)
library(ggpubr)

##############################################
# approx. to cumulative distribution binomial
# D function relative entropy
##################################
D_re <- function(a,p)
{
  a * log(a/p) + (1-a)*log((1-a)/(1-p))
}
upper_D <- function(k,n,p)
{
  exp(-n * D_re(k/n,p))
}
lower_D <- function(k,n,p)
{
  (1/sqrt(2*n)) * exp(-n * D_re(k/n,p))
}

############################################
hoeffding <- function(k,n,p)
{
  exp((-2/n * (n*p - k)^2))
}

pan.loss <- function(n,k,s)
{
  ((n-k)/(n+k*s))^n
}
kimura <- function(p0,Nalpha)
{
  num <- 1-exp(-2 * Nalpha * p0 )
  den <- 1 - exp(-2 * Nalpha)
  
  num/den
}

kimura.time <- function(p0,N)
{
  -4*N*((p0*log(p0)+((1-p0)*log(1-p0))))
}

k2.time <- function(p,Ne)
{
  -1/p*(4*Ne*(1-p)*log(1-p))
}
############################################
#   NETWORK PROPERTIES
############################################

# Mean path length - averaged over all 100 spaces
#
# For each network in folder, calculate path length
# calculate mean of all (100).
# Return mean + sd
#
mean.graph.path.length <- function(folder="Spacefixed/networks/pan64/")
{
  f <- list.files(path=folder)

  res <- matrix(nrow=length(f),ncol=2)
  
  for (i in 1:length(f))
  {
    g = as.matrix(read.table(file=paste(folder,f[i],sep="")))
    diag(g) <- 0
    g <- graph_from_adjacency_matrix(g,mode="undirected",weighted=TRUE)
    d.g <- distances(g)
    res[i,1] <- mean(colSums(d.g))/100
    res[i,2] <- mean_distance(g) #sd(rowSums(d.g))
  }
  
  c(mean(res[,1]),mean(res[,2]))
}
all.path.lengths <- function()
{
  res <- matrix(nrow=6,ncol=2)
  res[1,] <- mean.graph.path.length(folder="Spacefixed/networks/pan64/")
  res[2,] <- mean.graph.path.length(folder="Spacefixed/networks/lattice64/")
  res[3,] <- mean.graph.path.length(folder="Spacefixed/networks/ring64/")
  res[4,] <- mean.graph.path.length(folder="Spacefixed/networks/sf64p1/")
  res[5,] <- mean.graph.path.length(folder="Spacefixed/networks/sf64p2/")
  res[6,] <- mean.graph.path.length(folder="Spacefixed/networks/star64/")
  res
}


#############################################
#   NETWORKS FIGURE
##############################################
plot.network <- function(graph="Spacefixed/networks/pan64/1.txt",
                         layout=layout_with_fr,...)
{
  g = as.matrix(read.table(file=graph))
  diag(g) <- 0
  op <- par(mar=c(0,0,0,0))
  g <- graph_from_adjacency_matrix(g,mode="undirected",weighted=TRUE)
  plot(g,layout=layout,vertex.size=3.5,vertex.label=NA,
       vertex.color="black",...)
  par(op)
}
plot.all.spaces <- function()
{
  t.x = -0.05
  t.y = -1.15
  t.cex = 2.2
  op <- par(mar=c(0,0,0,0))

  graphics::layout(mat=matrix(data=c(1,2,3,4,5,6),ncol=3,byrow=T))
  plot.network(graph="Spacefixed/networks/pan64/1.txt")
#               layout=layout_in_circle)
  labs=expression(paste("Panmictic"))
  text(t.x,t.y,labels=labs,cex=t.cex)
  plot.network(graph="Spacefixed/networks/lattice64/1.txt")
  labs=expression(paste("Lattice"))
  text(t.x,t.y,labels=labs,cex=t.cex)
  plot.network(graph="Spacefixed/networks/ring64/1.txt",
               layout=layout_in_circle,edge.curved=TRUE)
  labs=expression(paste("Ring"))
  text(t.x,t.y,labels=labs,cex=t.cex)
  plot.network(graph="Spacefixed/networks/sf64p1/3.txt")
  labs=expression(paste("Scale-Free (p=1)"))
  text(t.x,t.y,labels=labs,cex=t.cex)
  plot.network(graph="Spacefixed/networks/sf64p2/9.txt")
  labs=expression(paste("Scale-Free (p=2)"))
  text(t.x,t.y,labels=labs,cex=t.cex)
  plot.network(graph="Spacefixed/networks/star64/2.txt")
  labs=expression(paste("Star"))
  text(t.x,t.y,labels=labs,cex=t.cex)
}
jpg.all.spaces <- function()
{
  jpeg(filename="figures/all.networks.jpg",
       units="cm",
       width=30,height=25,res=400,
       quality=100)
  plot.all.spaces()
  dev.off()
}

###
## JUST LOAD, COMBINE and save
###
load.SUPP.figvaryN.data <- function()
{
  q1 <- read.csv("Spacefixed/output/panfig4VaryNSUPP.10000.P1.csv",
                 header=T) 
  
  q2 <- read.csv("Spacefixed/output/panfig4VaryNSUPP.10000.P2.csv",
                 header=T) 
  q3 <- read.csv("Spacefixed/output/panfig4VaryNSUPP.10000.P3.csv",
                 header=T) 
 # q4 <- read.csv("Spacefixed/output/panfig4VaryNSUPP.1000.P3.csv",
#               header=T) 
  
  res <- rbind(q1,q2,q3)
  write.csv(res,"Spacefixed/output/panfig4VaryNSUPP.ALL.csv",
            row.names=FALSE)
  
  
}
plot.figvaryN.SUPP <- function(
  in.fname="Spacefixed/output/panfig4VaryNSUPP.ALL.csv",
                          titlestr = expression(paste(p[0]," = 0.1, N", alpha," = 8, N",beta," = 4, Nc = 1/4")),
                          type=1, 
                          ylab.expn=expression(paste(mu,"(",p[0],")")),
                          kimura.approx=kimura(0.3,4))
  
{
  res <- load.fig4.fix.gg(in.fname,type=type)
  res$N <- factor(res$N)
  res %>% group_by(N,q) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=q,y=stat)) +
    ggtitle(titlestr) +
    geom_line(aes(group=N,colour=N,x=q,y=stat)) +
    #   geom_point(aes(group=N,x=q,y=stat,colour=N)) +
    geom_ribbon(aes(group=N,ymin=stat-ci, 
                    ymax=stat+ci,fill=N),
                colour=NA,
                alpha=0.25) +
    #ylim(0.74,0.86) +   # Include this for S1 figure
    ylim(0.68,0.93) +   # Include this for S1 figure
    xlab(expression(q[0])) + 
    ylab(ylab.expn) +
    xlim(0.0,1.0) + 
    geom_hline(yintercept=kimura(0.1,8), linetype="dashed", 
               color = "red", size=1)
}

plot.figONEvaryN.SUPP <- function(
  in.fname="Spacefixed/output/panfig4VaryNSUPP.ALL.csv",
  titlestr = 
    expression(paste(p[0]," = 0.1, ",q[0]," = 0.2, N", alpha," = 8, N",
                     beta," = 4, Nc = 1/4")),
  type=1, 
  ylab.expn=expression(paste(mu,"(",p[0],")"))
)
  
{
  res <- load.fig4.fix.gg(in.fname,type=type)
  
  # and now just pull out the data for q0 = 0.1
  
  res <- subset(res,q==0.2)
  
  #  res$N <- factor(res$N)
  res %>% 
    group_by(N) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=N,y=stat)) + 
    geom_errorbar(mapping=aes(x=N, ymin=stat-ci, 
                              ymax=stat+ci), 
                  width=0.2, size=1, color="black") + 
    geom_point(aes(x=N,y=stat),size=3,shape=21,fill="white") +
    ggtitle(titlestr)+
    scale_x_continuous("Population Size N (log axis)", 
                       breaks=c(16,32,64,128,256,512), 
                       labels=c(16,32,64,128,256,512))+
    #xlab("Population Size")+
    ylab(ylab.expn) +
    coord_trans(x = "log2") +
    ylim(0.75,0.85)
}

################## VARY POPULATION SIZE ##################
# Try putting two plots on same figure
##########################################################

plot.figvaryN <- function(in.fname="Spacefixed/output/panfig4VaryN.10000.csv",
                           titlestr = expression(paste(p[0]," = 0.3, N", alpha," = 4, N",beta," = 8, Nc = 1/4")),
                           type=1, 
                           ylab.expn=expression(paste(mu,"(",p[0],")")),
                           kimura.approx=kimura(0.3,4))
  
{
  res <- load.fig4.fix.gg(in.fname,type=type)
  res$N <- factor(res$N)
  res %>% group_by(N,q) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=q,y=stat)) +
    ggtitle(titlestr) +
    geom_line(aes(group=N,colour=N,x=q,y=stat)) +
    #   geom_point(aes(group=N,x=q,y=stat,colour=N)) +
    geom_ribbon(aes(group=N,ymin=stat-ci, 
                    ymax=stat+ci,fill=N),
                colour=NA,
                alpha=0.25) +
    #ylim(0.74,0.86) +   # Include this for S1 figure
    ylim(0.68,0.93) +   # Include this for S1 figure
    xlab(expression(q[0])) + 
    ylab(ylab.expn) +
    xlim(0.0,1.0) + 
    geom_hline(yintercept=kimura.approx, linetype="dashed", 
               color = "red", size=1)
}

plot.figONEvaryN <- function(
  in.fname="Spacefixed/output/panfig4VaryN.10000.csv",
  titlestr = 
    expression(paste(p[0]," = 0.3, ",q[0]," = 0.1, N", alpha," = 4, N",
                     beta," = 8, Nc = 1/4")),
  type=1, 
  ylab.expn=expression(paste(mu,"(",p[0],")"))
)
  
{
  res <- load.fig4.fix.gg(in.fname,type=type)
  
  # and now just pull out the data for q0 = 0.1
  
  res <- subset(res,q==0.1)
  
  #  res$N <- factor(res$N)
  res %>% 
    group_by(N) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=N,y=stat)) + 
    geom_errorbar(mapping=aes(x=N, ymin=stat-ci, 
                              ymax=stat+ci), 
                  width=0.2, size=1, color="black") + 
    geom_point(aes(x=N,y=stat),size=3,shape=21,fill="white") +
    ggtitle(titlestr)+
    scale_x_continuous("Population Size N (log axis)", 
                       breaks=c(16,32,64,128,256,512), 
                       labels=c(16,32,64,128,256,512))+
    #xlab("Population Size")+
    ylab(ylab.expn) +
    coord_trans(x = "log2") +
    ylim(0.675,0.775)
}
plot.combined.figvaryN.SUPP <- function()
{
  # See: 
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/#create-some-plots
  
  f1 <- plot.figvaryN.SUPP()
  f2 <- plot.figONEvaryN.SUPP()
  ggarrange(f1, f2,  
            labels = c("A", "B"),
            ncol = 2, nrow = 1)
}
plot.combined.figvaryN <- function()
{
# See: 
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/#create-some-plots
  
  f1 <- plot.figvaryN()
  f2 <- plot.figONEvaryN()
  ggarrange(f1, f2,  
            labels = c("A", "B"),
            ncol = 2, nrow = 1)
}
jpg.combined.figN <- function()
{
  plot.combined.figvaryN()
ggsave("figures/combined.figvaryN.jpg",dpi=400,
       units="cm",
       width=22,height=10)

}
jpg.combined.figN.SUPP <- function()
{
  plot.combined.figvaryN.SUPP()
  ggsave("figures/combined.figvaryN.SUPP.jpg",dpi=400,
         units="cm",
         width=22,height=10)
  
}
#########################################
# Load the data for prob. fixation A
#
##########################################
load.fig1.for.gg <- function(in.fname)
{
  res <- read.csv(in.fname,header=TRUE)
  res$fix <- ifelse(res$fixAa==1,1,0)
  res$NiB <- res$Nibeta
  res$Nibeta <- "0"
  nib <- sort(unique(res$NiB))
  nibtext <- c("0","2","4","8","16","32") 
  for (index in 1:length(nib))
  {
    res$Nibeta[which(res$NiB==nib[index])] <- nibtext[index]
  }
  res
}

load.fig1.ttf <- function(in.fname)
{
  res <- read.csv(in.fname,header=TRUE)
  res <- subset(res,fixAa==1)
  res$Nibeta <- factor(res$Nibeta)
  res
}

### #######################################################
###  Figure 1 HR paper - but with different spaces
###
### USING N = 64 
###########################################################
fig.probFixA.spaces <- function() # Hard wire
{
  res <- load.fig1.for.gg("Spacefixed/output/pan64probFixA.csv")
  res$Model <- "Panmictic"
  res2 <- load.fig1.for.gg("Spacefixed/output/ring64probFixA.csv")
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/lattice64probFixA.csv")
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
#  res2 <- load.fig1.for.gg("Spacefixed/output/sf64p1probFixA.csv")
  res2 <- load.fig1.for.gg("Spacefixed/output/rnd64probFixA.4edge.csv")
  
    res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/sf64p2probFixA.csv")
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/star64probFixA.csv")
  res2$Model <- "Star"
  res <- rbind(res,res2)
  

  titlestr <- expression(paste("N = 64, ",p[0]," = 0.3,", " N", alpha," = 4, Nc = 0"))
  
  # Reorder for drawing legend
  res$Nibeta <- factor(res$Nibeta, 
              levels=c("0","2","4","8","16","32"), 
              labels=c("0","2","4","8","16","32"))
  
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                         "Scale Free (p=2)","Star"),
          #      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
          labels=c("Panmictic","Lattice","Ring","Random 4 edge",
                                  
                  "Scale-Free (p=2)","Star"))
  
  Nibeta.txt <- expression(paste("N",beta))
  
  res %>% group_by(Model,NiB,Nibeta,q) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=q,y=stat)) +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    geom_line(aes(group=NiB,colour=Nibeta,x=q,y=stat)) +
    geom_point(aes(group=NiB,x=q,y=stat,colour=Nibeta)) +
    geom_ribbon(aes(group=NiB,ymin=stat-ci, 
                     ymax=stat+ci,fill=Nibeta),
                 colour=NA,
                 alpha=0.25) +
    labs(colour=Nibeta.txt,fill=Nibeta.txt) +
    scale_x_continuous(expression(q[0]),
                       breaks=c(0,0.25,0.5,0.75,1),
                       labels=c(0,0.25,0.5,0.75,1),
                       limits=c(0,1)) +
    
    ylab(expression(paste(mu,"(",p[0],")"))) + 
    geom_hline(yintercept=kimura(0.3,4), linetype="dashed", 
               color = "red", size=1)

}

fig.probFixA.spaces.jpeg <- function()
{

  fig.probFixA.spaces()
  ggsave("figures/fig.probFixA.space.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}

fig.probFixA.spaces.ttf <- function() # Hard wire
{
  res <- load.fig1.ttf("Spacefixed/output/pan64probFixA.csv")
  res$Model <- "Panmictic"
  res2 <- load.fig1.ttf("Spacefixed/output/ring64probFixA.csv")
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig1.ttf("Spacefixed/output/star64probFixA.csv")
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  #res2 <- load.fig1.ttf("Spacefixed/output/sf64p1probFixA.csv")
  
  res2 <- load.fig1.ttf("Spacefixed/output/rnd64probFixA.csv")
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig1.ttf("Spacefixed/output/sf64p2probFixA.csv")
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  
  res2 <- load.fig1.ttf("Spacefixed/output/lattice64probFixA.csv")
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  
  titlestr <- expression(paste("N = 64, ",p[0]," = 0.3,", " N", alpha," = 4, Nc = 0"))
  
  # Reorder for drawing legend
  res$Nibeta <- factor(res$Nibeta, 
                       levels=c("0","2","4","8","16","32"), 
                       labels=c("0","2","4","8","16","32"))
  
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)","Star"),
                #      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
 
                               labels=c("Panmictic","Lattice","Ring","Poission Mean 10",
                                        "Scale-Free (p=2)","Star"))
  
  Nibeta.txt <- expression(paste("N",beta))
  
#  ylab.txt <- expression(paste("Mean time to fixation. Axis "," ",log[10]," scale",sep=""))
  ylab.txt <- expression(paste("Mean time to fixation"))
  
  res %>% group_by(Model,Nibeta,q) %>%
    summarise(stat=mean(gensAa),ci=qt(0.975, n() - 1) * sd(gensAa) / sqrt(n())) %>%
    ggplot(aes(x=q,y=stat)) +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    
   #coord_trans(y = "log10") +
    
#    scale_y_log10(breaks=c(25,50,100,150,300),
#                  labels=c(25,50,100,150,300),
#                  limits=c(20,300)) +
    
    scale_x_continuous(expression(q[0]),
                       breaks=c(0,0.25,0.5,0.75,1),
                       labels=c(0,0.25,0.5,0.75,1),
                       limits=c(0,1)) +
    geom_line(aes(group=Nibeta,colour=Nibeta,x=q,y=stat)) +
    geom_point(aes(group=Nibeta,x=q,y=stat,colour=Nibeta)) +
    geom_ribbon(aes(group=Nibeta,ymin=stat-ci, 
                    ymax=stat+ci,fill=Nibeta),
                colour=NA,
                alpha=0.25) +
    labs(colour=Nibeta.txt,fill=Nibeta.txt) +
    xlab(expression(q[0])) + ylab(ylab.txt) 
}


fig.probFixA.spaces.ttf.jpeg <- function()
{
  
  fig.probFixA.spaces.ttf()
  ggsave("figures/fig.probFixA.spaces.ttf.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}
####################################
################# FIGURE 2 (HR)
####################################

load.fig2.for.gg <- function(in.fname,type)
{
  res <- read.csv(in.fname,header=TRUE)
  if (type==1) res$fix <- ifelse(res$fixAa==1,1,0)
  if (type==2) res$fix <- ifelse(res$fixAa==2,1,0)
  if (type==3) res$fix <- ifelse(res$fixBb==1,1,0)
  if (type==4) res$fix <- ifelse(res$fixBb==2,1,0)
  res$Nc <- factor(res$Nc)
  res
}

fig2.probA.varyB.Nc.spaces <- function(type=1,
                       ylabel=expression(paste(mu,"(",p[0],")"))) # Hard wire
{
  res <- load.fig2.for.gg("Spacefixed/output/pan64FixAvaryBNc.csv",type)
  res$Model <- "Panmictic"
  res2 <- load.fig2.for.gg("Spacefixed/output/lattice64FixAvaryBNc.csv",type)
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/ring64FixAvaryBNc.csv",type)
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/sf64p1FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/sf64p2FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/star64FixAvaryBNc.csv",type)
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, ", " N", alpha," = 8"))
  xlabel <- expression(paste("N",beta))
  
  
  res$Nc <- factor(res$Nc,
                   levels=c(0,0.25,1,4,8,16),
                   labels=c("0","1/4","1","4","8","16")) 
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)","Star"),
                      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
                               "Scale-Free (p=2)","Star"))
  
  res %>% group_by(Model,Nc,Nibeta) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=Nibeta,y=stat)) +
    coord_trans(x = "log2") +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    geom_line(aes(group=Nc,colour=Nc,x=Nibeta,y=stat)) +
    geom_point(aes(group=Nc,x=Nibeta,y=stat,colour=Nc)) +
    geom_ribbon(aes(group=Nc,ymin=stat-ci, 
                    ymax=stat+ci,fill=Nc),
                colour=NA,
                alpha=0.25) +
    scale_x_continuous(xlabel,
                       breaks=c(1,2,4,8,16,32),
                       labels=c(1,2,4,8,16,32),
                       limits=c(1,34)) + 
    geom_hline(yintercept=kimura(0.1,8), linetype="dashed", 
               color = "red", size=1) + 
    
    ylab(ylabel) 
}
fig2.probAvaryB.Nc.spaces.jpeg <- function()
{
  fig2.probA.varyB.Nc.spaces()
  ggsave("figures/fig.probFixA.VaryB.Nc.spaces.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}
fig2.SUPP.probA.varyB.Nc.spaces <- function(type=1,
                                       ylabel=expression(paste(mu,"(",p[0],")"))) # Hard wire
{
  res <- load.fig2.for.gg("Spacefixed/output/pan64SUPPFixAvaryBNc.csv",type)
  res$Model <- "Panmictic"
  res2 <- load.fig2.for.gg("Spacefixed/output/lattice64SUPPFixAvaryBNc.csv",type)
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/ring64SUPPFixAvaryBNc.csv",type)
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/sf64SUPPp1FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/sf64SUPPp2FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  res2 <- load.fig2.for.gg("Spacefixed/output/star64SUPPFixAvaryBNc.csv",type)
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.3, ", " N", alpha," = 4"))
  xlabel <- expression(paste("N",beta))
  
  
  res$Nc <- factor(res$Nc,
                   levels=c(0,0.25,1,4,8,16),
                   labels=c("0","1/4","1","4","8","16")) 
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)","Star"),
                      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
                               "Scale-Free (p=2)","Star"))
  
  res %>% group_by(Model,Nc,Nibeta) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=Nibeta,y=stat)) +
    coord_trans(x = "log2") +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    geom_line(aes(group=Nc,colour=Nc,x=Nibeta,y=stat)) +
    geom_point(aes(group=Nc,x=Nibeta,y=stat,colour=Nc)) +
    geom_ribbon(aes(group=Nc,ymin=stat-ci, 
                    ymax=stat+ci,fill=Nc),
                colour=NA,
                alpha=0.25) +
    scale_x_continuous(xlabel,
                       breaks=c(1,2,4,8,16,32),
                       labels=c(1,2,4,8,16,32),
                       limits=c(1,34)) + 
    geom_hline(yintercept=kimura(0.3,4), linetype="dashed", 
               color = "red", size=1) + 
    
    ylab(ylabel) 
}
fig2.SUPP.probAvaryB.Nc.spaces.jpeg <- function()
{
  fig2.SUPP.probA.varyB.Nc.spaces()
  ggsave("figures/fig.SUPP.probFixA.VaryB.Nc.spaces.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}

load.fig2.ttf.for.gg <- function(in.fname,type=1)
{
  res <- read.csv(in.fname,header=TRUE)
  if (type==1)
  {
    res <- subset(res,fixAa==1)
    res$ttf <- res$gensAa
  }
  if (type==2)
  {
    res <- subset(res,fixAa==2)
    res$ttf <- res$gensAa
  }
  if (type==3)
  {
    res <- subset(res,fixBb==1)
    res$ttf <- res$gensBb
  }
  if (type==4)
  {
    res <- subset(res,fixBb==2)
    res$ttf <- res$gensBb
  }
  res
}

fig2.ttfA.varB.Nc.spaces <- function(type=1,ylabel="Mean time to fixation allele A") # Hard wire
{
  res <- load.fig2.ttf.for.gg("Spacefixed/output/pan64FixAvaryBNc.csv",type)
  res$Model <- "Panmictic"
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/lattice64FixAvaryBNc.csv",type)
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/ring64FixAvaryBNc.csv",type)
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/sf64p1FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/sf64p2FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/star64FixAvaryBNc.csv",type)
  res2$Model <- "Star"
  res <- rbind(res,res2)

  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, "," N",alpha," = 8"))
  xlabel <- expression(paste("N",beta))
  
  # Reorder for drawing legend
  
  res$Nc <- factor(res$Nc)
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)","Star"),
                      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
                               "Scale-Free (p=2)","Star"))
  
  res %>% group_by(Model,Nc,Nibeta) %>%
    summarise(stat=mean(ttf),ci=qt(0.975, n() - 1) * sd(ttf) / sqrt(n())) %>%
    ggplot(aes(x=Nibeta,y=stat)) +
    coord_trans(x = "log2") +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    geom_line(aes(group=Nc,colour=Nc,x=Nibeta,y=stat)) +
    geom_point(aes(group=Nc,x=Nibeta,y=stat,colour=Nc)) +
    geom_ribbon(aes(group=Nc,ymin=stat-ci, 
                    ymax=stat+ci,fill=Nc),
                alpha=0.25) +
    scale_x_continuous(xlabel,
                       breaks=c(1,2,4,8,16,32),
                       labels=c(1,2,4,8,16,32),
                       limits=c(1,34)) + 
    
    ylab(ylabel)
    #geom_hline(yintercept=kimura.time(0.1,64), linetype="dashed", 
    #           color = "red", size=1) 
    
}

fig2.ttfFixAvaryB.Nc.spaces.jpeg <- function()
{
  fig2.ttfA.varB.Nc.spaces()
  ggsave("figures/fig.ttf.FixA.VaryB.Nc.spaces.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}

fig2.SUPP.ttfA.varB.Nc.spaces <- function(type=1,ylabel="Mean time to fixation allele A") # Hard wire
{
  res <- load.fig2.ttf.for.gg("Spacefixed/output/pan64SUPPFixAvaryBNc.csv",type)
  res$Model <- "Panmictic"
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/lattice64SUPPFixAvaryBNc.csv",type)
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/ring64SUPPFixAvaryBNc.csv",type)
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/sf64SUPPp1FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  res2 <- load.fig2.ttf.for.gg("Spacefixed/output/sf64SUPPp2FixAvaryBNc.csv",type)
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.3, "," N",alpha," = 4"))
  xlabel <- expression(paste("N",beta))
  
  # Reorder for drawing legend
  
  res$Nc <- factor(res$Nc)
  
  # Try reordering the Model term for drawing facet order
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)","Star"),
                      labels=c("Panmictic","Lattice","Ring","Scale-Free (p=1)",
                               "Scale-Free (p=2)","Star"))
  
  res %>% group_by(Model,Nc,Nibeta) %>%
    summarise(stat=mean(ttf),ci=qt(0.975, n() - 1) * sd(ttf) / sqrt(n())) %>%
    ggplot(aes(x=Nibeta,y=stat)) +
    coord_trans(x = "log2") +
    facet_wrap(~Model) +
    ggtitle(titlestr) +
    geom_line(aes(group=Nc,colour=Nc,x=Nibeta,y=stat)) +
    geom_point(aes(group=Nc,x=Nibeta,y=stat,colour=Nc)) +
    geom_ribbon(aes(group=Nc,ymin=stat-ci, 
                    ymax=stat+ci,fill=Nc),
                alpha=0.25) +
    scale_x_continuous(xlabel,
                       breaks=c(1,2,4,8,16,32),
                       labels=c(1,2,4,8,16,32),
                       limits=c(1,34)) + 
    
    ylab(ylabel) 
#    geom_hline(yintercept=kimura.time(0.1,64), linetype="dashed", 
#               color = "red", size=1) 
  
}

fig2.SUPP.ttfFixAvaryB.Nc.spaces.jpeg <- function()
{
  fig2.SUPP.ttfA.varB.Nc.spaces()
  ggsave("figures/fig.SUPP.ttf.FixA.VaryB.Nc.spaces.jpg",dpi=400,
         units="cm",
         width=15,height=13)
}


############### HAPLOID UTILITIES ###########################
########################
# buildtable.Nibeta.ttf
##########################################################
# TTF - for gametes must take the maximum of the time
# Builds the complete table for a specific NBeta
# for AA, Ab, etc. time to fixation
##########################################################
buildtable.Nibeta.ttf <- function(res.all)
{
  res <- res.all
  res <- subset(res,(res$fixAa==1 & res$fixBb==1)) # AB
  res$ttf <- apply(res[,c("gensAa","gensBb")],1,max)
  res$type <- "AB"
  
  res2 <- res.all
  res2 <- subset(res2,(res2$fixAa==1 & res2$fixBb==2)) # Ab
  res2$ttf <- apply(res2[,c("gensAa","gensBb")],1,max)
  res2$type <- "Ab"
  res <- rbind(res,res2)
  
  res2 <- res.all
  res2 <- subset(res2,(res2$fixAa==2 & res2$fixBb==1)) # aB
  res2$ttf <- apply(res2[,c("gensAa","gensBb")],1,max)
  res2$type <- "aB"
  res <- rbind(res,res2)
  
  res2 <- res.all
  res2 <- subset(res2,(res2$fixAa==2 & res2$fixBb==2)) # ab
  res2$ttf <- apply(res2[,c("gensAa","gensBb")],1,max)
  res2$type <- "ab"
  res <- rbind(res,res2)
  
  res$Nc <- factor(res$Nc)
  res
}

#####################################################
# Note res.all here is just the res.all table
# but selected for a specific Nibeta value
####################################################
buildtable.Nibeta <- function(res.all)
{
  res <- res.all
  res$fix <- 0
  res[which((res$fixAa==1) & (res$fixBb==1)),"fix"] <- 1
  res$type <- "AB"
  
  res2 <- res.all
  res2$fix <- 0
  res2[which((res2$fixAa==1) & (res2$fixBb==2)),"fix"] <- 1
  res2$type <- "Ab"
  res <- rbind(res,res2)
  
  res2 <- res.all
  res2$fix <- 0
  res2[which((res2$fixAa==2) & (res2$fixBb==1)),"fix"] <- 1
  res2$type <- "aB"
  res <- rbind(res,res2)
  
  res2 <- res.all
  res2$fix <- 0
  res2[which((res2$fixAa==2) & (res2$fixBb==2)),"fix"] <- 1
  res2$type <- "ab"
  res <- rbind(res,res2)
  
  res
}
#################################################################


fig.HaplotypeFix.VaryB.Nc <- function()
{
  title=expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, N",
                         alpha," = 8, N",beta,"(panels) 2, 4, 7, 8"))
  
  res.all <- read.csv("Spacefixed/output/pan64FixAvaryBNc.csv",header=T)
  
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  
  # and keep just some Nibeta
  
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  
  res <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res <- rbind(res,
                 buildtable.Nibeta(
                   res.all[which(res.all$Nibeta==beta),]))
  }
  res$Model <- "Panmictic"
  
  # res.all <- read.csv("Spacefixed/output/lattice64FixAvaryBNc.csv",header=T)
  # # First use fixAa and fixBb to create the 4 possible datasets
  # res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # # and keep just some Nibeta
  # nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  # res.all <- res.all[nikeep,]
  # res2 <- NULL
  # for (beta in unique(res.all$Nibeta))
  # {
  #   res2 <- rbind(res2,
  #                buildtable.Nibeta(
  #                  res.all[which(res.all$Nibeta==beta),]))
  # }
  # res2$Model <- "Lattice"
  # res <- rbind(res,res2)
  # 
  # res.all <- read.csv("Spacefixed/output/ring64FixAvaryBNc.csv",header=T)
  # # First use fixAa and fixBb to create the 4 possible datasets
  # res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # # and keep just some Nibeta
  # nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  # res.all <- res.all[nikeep,]
  # res2 <- NULL
  # for (beta in unique(res.all$Nibeta))
  # {
  #   res2 <- rbind(res2,
  #                 buildtable.Nibeta(
  #                   res.all[which(res.all$Nibeta==beta),]))
  # }
  # res2$Model <- "Ring"
  # res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/sf64p1FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/sf64p2FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/star64FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  
  res$Nc <- factor(res$Nc, 
                   levels=rev(c("16","8","4","1","0.5","0.25","0.0625","0.001")),
                   labels=rev(c("16","8","4","1","1/2","1/4","1/16","0")))
  
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Scale Free (p=1)",
                               "Scale Free (p=2)", "Star"),
                      labels=c("Panmictic","SF (p=1)",
                               "SF (p=2)", "Star"))
  
#  levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
#           "Scale Free (p=2)", "Star"),
#  labels=c("Panmictic","Lattice","Ring","SF (p=1)",
#           "SF (p=2)", "Star"))

  
  res$type <- factor(res$type)
  res$type <- factor(res$type,
                     levels=c("AB","Ab","aB","ab"),
                     labels=c("AB","Ab","aB","ab"))
  
  res %>% group_by(Model,Nibeta,Nc,type) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=Nc,y=stat)) +
    ggtitle(title) +
    facet_grid(Model ~ Nibeta) +
    coord_trans(x = "log2") +
    geom_line(aes(group=type,colour=type,x=Nc,y=stat)) + 
    geom_ribbon(aes(group=type,ymin=stat-ci, 
                    ymax=stat+ci,fill=type),
                colour=NA,
                alpha=0.25)+
    # ylim(c(0,1)) + 
    ylab("Haplotype Fixation Probability")
}

fig.HaplotypeFix.VaryB.Nc.jpeg <- function()
{
  fig.HaplotypeFix.VaryB.Nc()
  ggsave("figures/fig.HaplotypeFix.VaryB.Nc.jpg",dpi=400,
         units="cm",
         width=30,height=25)
}

################ haplotype ttf ####################


fig.HaplotypeFix.VaryB.Nc.ttf <- function()
{
  title=expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, N",
                         alpha," = 8, N",beta,"(panels) 2, 4, 7, 8"))
  
  res.all <- read.csv("Spacefixed/output/pan64FixAvaryBNc.csv",header=T)
  
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  
  # and keep just some Nibeta
  
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  
  res <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res <- rbind(res,
                 buildtable.Nibeta.ttf(
                   res.all[which(res.all$Nibeta==beta),]))
  }
  res$Model <- "Panmictic"
  
  res.all <- read.csv("Spacefixed/output/lattice64FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Lattice"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/ring64FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/sf64p1FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Scale Free (p=1)"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/sf64p2FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/star64FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  
  res$Nc <- factor(res$Nc, 
                   levels=rev(c("16","8","4","1","0.5","0.25","0.0625","0.001")),
                   labels=rev(c("16","8","4","1","1/2","1/4","1/16","0")))
  
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Lattice","Ring","Scale Free (p=1)",
                               "Scale Free (p=2)", "Star"),
                      labels=c("Panmictic","Lattice","Ring","SF (p=1)",
                               "SF (p=2)", "Star"))
  
  res$type <- factor(res$type)
  res$type <- factor(res$type,
                     levels=c("AB","Ab","aB","ab"),
                     labels=c("AB","Ab","aB","ab"))
  
  res %>% group_by(Model,Nibeta,Nc,type) %>%
    summarise(stat=mean(ttf),ci=qt(0.975, n() - 1) * sd(ttf) / sqrt(n())) %>%
    ggplot(aes(x=Nc,y=stat)) +
    ggtitle(title) +
    facet_grid(Model ~ Nibeta) +
    coord_trans(x = "log2") +
    scale_y_log10("Time to Fixation",
                  breaks=c(10,100,200,400,800),
                  labels=c(10,100,200,400,800),
                  limits=c(5,450)) +
    
    geom_line(aes(group=type,colour=type,x=Nc,y=stat)) + 
    geom_ribbon(aes(group=type,ymin=stat-ci, 
                    ymax=stat+ci,fill=type),
                colour=NA,
                alpha=0.25)+
    # ylim(c(0,1)) + 
    ylab("Time to Fixation")
}

fig.HaplotypeFix.VaryB.Nc.jpeg.ttf <- function()
{
  fig.HaplotypeFix.VaryB.Nc.ttf()
  ggsave("figures/fig.HaplotypeFix.VaryB.Nc.ttf.jpg",dpi=400,
         units="cm",
         width=30,height=30)
}

fig.REDUCED.HaplotypeFix.VaryB.Nc.ttf <- function()
{
  title=expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, N",
                         alpha," = 8, N",beta,"(panels) 2, 4, 7, 8"))
  
  res.all <- read.csv("Spacefixed/output/pan64FixAvaryBNc.csv",header=T)
  
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  
  # and keep just some Nibeta
  
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  
  res <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res <- rbind(res,
                 buildtable.Nibeta.ttf(
                   res.all[which(res.all$Nibeta==beta),]))
  }
  res$Model <- "Panmictic"
  
  res.all <- read.csv("Spacefixed/output/sf64p1FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SF1"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/sf64p2FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Scale Free (p=2)"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/star64FixAvaryBNc.csv",header=T)
  # First use fixAa and fixBb to create the 4 possible datasets
  res.all$Nc[which(res.all$Nc==0.0)] <- 0.001  # do this first....
  # and keep just some Nibeta
  nikeep <- which(res.all$Nibeta==c(2,4,7,8))
  res.all <- res.all[nikeep,]
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star"
  res <- rbind(res,res2)
  
  
  res$Nc <- factor(res$Nc, 
                   levels=rev(c("16","8","4","1","0.5","0.25","0.0625","0.001")),
                   labels=rev(c("16","8","4","1","1/2","1/4","1/16","0")))
  
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","SF1",
                               "Scale Free (p=2)", "Star"),
                      labels=c("Panmictic","SF (p=1)",
                               "SF (p=2)", "Star"))
  
  res$type <- factor(res$type)
  res$type <- factor(res$type,
                     levels=c("AB","Ab","aB","ab"),
                     labels=c("AB","Ab","aB","ab"))
  
  res %>% group_by(Model,Nibeta,Nc,type) %>%
    summarise(stat=mean(ttf),ci=qt(0.975, n() - 1) * sd(ttf) / sqrt(n())) %>%
    ggplot(aes(x=Nc,y=stat)) +
    ggtitle(title) +
    facet_grid(Model ~ Nibeta) +
    coord_trans(x = "log2") +
    scale_y_log10("Time to Fixation",
                  breaks=c(10,100,200,400,800),
                  labels=c(10,100,200,400,800),
                  limits=c(5,450)) +
    
    geom_line(aes(group=type,colour=type,x=Nc,y=stat)) + 
    geom_ribbon(aes(group=type,ymin=stat-ci, 
                    ymax=stat+ci,fill=type),
                colour=NA,
                alpha=0.25)+
    # ylim(c(0,1)) + 
    ylab("Time to Fixation")
}
fig.REDUCED.HaplotypeFix.VaryB.Nc.jpeg.ttf <- function()
{
  fig.REDUCED.HaplotypeFix.VaryB.Nc.ttf()
  ggsave("figures/fig.REDUCED.HaplotypeFix.VaryB.Nc.ttf.jpg",dpi=400,
         units="cm",
         width=30,height=25)
}


##################### STEADY STATE TESTS #######################
################################################################

fig.SSW.test1.FixA <- function()
{
  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, "," N",alpha," = 4, ",
                               " N",beta," = 4"))
  
  res <- load.fig1.for.gg("Spacefixed/output/SSW/pan64.FixA.csv")
  res$Model <- "Panmictic"
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/ring64.FixA.csv")
  res2$Model <- "Ring"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/sf64.FixA.csv")
  res2$Model <- "SF"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/star64.FixA.csv")
  res2$Model <- "Star"
  res <- rbind(res,res2)
  res$Group <- "HR"
  
  # Added group for SS without setting diag = 0
  #
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSpan64.T1.FixA.csv")
  res2$Model <- "Panmictic (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSring64.T1.FixA.csv")
  res2$Model <- "Ring (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSsf64p2.T1.FixA.csv")
  res2$Model <- "SF (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSstar64.T1.FixA.csv")
  res2$Model <- "Star (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  ######
  # Added group for SS without setting diag = 0 but GLOBAL p1
  #
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSGpan64.T1.FixA.csv")
  res2$Model <- "Panmictic (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSGring64.T1.FixA.csv")
  res2$Model <- "Ring (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSGsf64p2.T1.FixA.csv")
  res2$Model <- "SF (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSGstar64.T1.FixA.csv")
  res2$Model <- "Star (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  ######
  
    
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWpan64.T1.FixA.csv")
  res2$Model <- "Panmictic (WL)"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWring64.T1.FixA.csv")
  res2$Model <- "Ring (WL)"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWsf64p2.T1.FixA.csv")
  res2$Model <- "SF (WL)"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWstar64.T1.FixA.csv")
  res2$Model <- "Star (WL)"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWGpan64.T1.FixA.csv")
  res2$Model <- "Panmictic (WG)"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWGring64.T1.FixA.csv")
  res2$Model <- "Ring (WG)"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWGstar64.T1.FixA.csv")
  res2$Model <- "Star (WG)"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  
  res2 <- load.fig1.for.gg("Spacefixed/output/SSW/SSWGsf64p2.T1.FixA.csv")
  res2$Model <- "sf"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Ring","SF","Star",
                               "Panmictic (SS)", "Ring (SS)",
                               "SF (SS)","Star (SS)",
                               "Panmictic (SG)", "Ring (SG)",
                               "SF (SG)","Star (SG)",
                               "Panmictic (WL)", "Ring (WL)",
                               "SF (WL)","Star (WL)",
                               "Panmictic (WG)", "Ring (WG)",
                               "sf","Star (WG)"),
                      labels=c("FC", "R", 
                               "SF","S",
                               "FC-SS","R-SS",
                               "SF-SS","S-SS",
                               "FC-SG","R-SG",
                               "SF-SG","S-SG",
                               "FC-WL","R-WL",
                               "SF-WL","S-WL",
                               "FC-WG","R-WG",
                               "SF-WG","S-WG"))

  res$Group <- factor(res$Group,
                      levels=c("HR","SS","SG","WL","WG"),
                      labels=c("HR","SS","SG","WL","WG"))
    
  res %>% 
    group_by(Model,Group) %>%
    summarise(stat=mean(fix),ci=qt(0.975, n() - 1) * sd(fix) / sqrt(n())) %>%
    ggplot(aes(x=Model,y=stat)) + 
    geom_errorbar(mapping=aes(x=Model, ymin=stat-ci, 
                              ymax=stat+ci), 
                  width=0.2, size=1, color="black") + 
    geom_point(aes(x=Model,y=stat,color=Group),size=3) +
    ggtitle(titlestr) +
    ylab(expression(paste(mu,"(",p[0],")"))) +
    ylim(0.1,0.5)
  
}
fig.SSW.test1.FixA.ttf <- function()
{
  titlestr <- expression(paste("N = 64, ", p[0]," = ",q[0]," = 0.1, "," N",alpha," = 4, ",
                               " N",beta," = 4"))
  
  res.all <- read.csv("Spacefixed/output/SSW/pan64.FixA.csv",header=T)
  res <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res <- rbind(res,
                 buildtable.Nibeta.ttf(
                   res.all[which(res.all$Nibeta==beta),]))
  }
  res$Model <- "Panmictic"
  res$Group <- "HR"
  
  res.all <- read.csv("Spacefixed/output/SSW/ring64.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Ring"
  res2$Group <- "HR"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/sf64.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SF"
  res2$Group <- "HR"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/star64.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                 buildtable.Nibeta.ttf(
                   res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star"
  res2$Group <- "HR"
  res <- rbind(res,res2)
  res$ttf <- res$ttf * unique(res$N) # so it gets changed back after
  
  ##### SS ###########
  res.all <- read.csv("Spacefixed/output/SSW/SSpan64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Panmictic (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSring64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Ring (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSsf64p2.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SF (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSstar64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star (SS)"
  res2$Group <- "SS"
  res <- rbind(res,res2)
  #####################
  
  ##### SSG ###########
  res.all <- read.csv("Spacefixed/output/SSW/SSGpan64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Panmictic (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSGring64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Ring (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSGsf64p2.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SF (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSGstar64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "Star (SG)"
  res2$Group <- "SG"
  res <- rbind(res,res2)
  ########## END SG ###########
  
  res.all <- read.csv("Spacefixed/output/SSW/SSWpan64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "PanWL"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSWring64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "RingWL"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  
  res <- rbind(res,res2)

  res.all <- read.csv("Spacefixed/output/SSW/SSWsf64p2.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SFWL"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  
  res.all <- read.csv("Spacefixed/output/SSW/SSWstar64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "StarWL"
  res2$Group <- "WL"
  res <- rbind(res,res2)
  
  
  res.all <- read.csv("Spacefixed/output/SSW/SSWGpan64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "PanWG"
  res2$Group <- "WG"
  res <- rbind(res,res2)

  res.all <- read.csv("Spacefixed/output/SSW/SSWGring64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "RingWG"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  
    
  res.all <- read.csv("Spacefixed/output/SSW/SSWGstar64.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "StarWG"
  res2$Group <- "WG"
  res <- rbind(res,res2)

  res.all <- read.csv("Spacefixed/output/SSW/SSWGsf64p2.T1.FixA.csv",header=T)
  res2 <- NULL
  for (beta in unique(res.all$Nibeta))
  {
    res2 <- rbind(res2,
                  buildtable.Nibeta.ttf(
                    res.all[which(res.all$Nibeta==beta),]))
  }
  res2$Model <- "SFWG"
  res2$Group <- "WG"
  res <- rbind(res,res2)
  
    
  res$ttf <- res$ttf/unique(res$N) # Normalise by population
  
  #####################################################
  res$Model <- factor(res$Model,
                      levels=c("Panmictic","Ring","SF","Star",
                               "Panmictic (SS)", "Ring (SS)",
                               "SF (SS)","Star (SS)",
                               "Panmictic (SG)", "Ring (SG)",
                               "SF (SG)","Star (SG)",
                               "PanWL", "RingWL",
                               "SFWL","StarWL",
                               "PanWG", "RingWG",
                               "SFWG","StarWG"),
                      labels=c("FC", "R", 
                               "SF","S",
                               "FC-SS","R-SS",
                               "SF-SS","S-SS",
                               "FC-SG","R-SG",
                               "SF-SG","S-SG",
                               "FC-WL","R-WL",
                               "SF-WL","S-WL",
                               "FC-WG","R-WG",
                               "SF-WG","S-WG"))
  
  ####################################################
  res$Group <- factor(res$Group,
                      levels=c("HR", "SS","SG", "WL","WG"),
                      labels=c("HR","SS","SG","WL","WG"))
  

  ylab.txt <- expression(paste("Mean time to fixation (/N). Axis "," ",log[10]," scale",sep=""))
  
  res %>% 
    group_by(Model,Group) %>%
    summarise(stat=mean(ttf),ci=qt(0.975, n() - 1) * sd(ttf) / sqrt(n())) %>%
    ggplot(aes(x=Model,y=stat)) + 
    geom_errorbar(mapping=aes(x=Model, ymin=stat-ci, 
                              ymax=stat+ci), 
                  width=0.2, size=1, color="black") + 
    geom_point(aes(x=Model,y=stat,color=Group),size=3) +
    ggtitle(titlestr) +
    scale_y_log10(ylab.txt,
                  breaks=c(0.1,1,10,100,500,1000,5000,10000),
                  labels=c(0.1, 1,10,100,500,1000,5000, 10000),
                  limits=c(10,5000))
  
}
jpg.SSW.Models <- function()
{
  f1 <- fig.SSW.test1.FixA()
  ggsave("figures/fig.SSW.FixA.jpg",dpi=400,
         units="cm",
         width=30,height=15)
  f2 <- fig.SSW.test1.FixA.ttf()
  ggsave("figures/fig.SSW.FixA.TTF.jpg",dpi=400,
         units="cm",
         width=30,height=15)
}
plot.combined.SSW.Models <- function()
{
  f1 <- fig.SSW.test1.FixA()
  f2 <- fig.SSW.test1.FixA.ttf()
  ggarrange(f1, f2,  
            labels = c("A", "B"),label.x=c(1,1),
            label.y=c(1,1),
            ncol = 2, nrow = 1)
}
