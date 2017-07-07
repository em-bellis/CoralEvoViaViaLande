library(MASS)
library(car)
library(sfsmisc)
library(ggplot2)
library(grid)

##choose parameter values
#env1 = hot/stress env; env2 = control/refugia env
gen=30  #total number of generations to run
gen2100=20  #number of generations completed by the year 2100
theta = c(0.1,0.1) #optimum trait values
zbar1.start=0.371 #mean trait value in stressful env.
zbar2.start=0.106 #mean trait value in refugia env.
G11=0.033 #additive genetic variance in stressful env.
G22=0.001 
P11=0.052 #phenotypic variance ...
P22=0.009
omega1=3*P11 #width of the fitness function
omega2=3*P22
r12=c(0.999, 0.75, 0.375, 0, -0.375, -0.75, -0.999)   #genetic correlation

##determine qs (proportion population inhabiting stressful environment) for first gen2100 generations
#rcp 8.5/6 (linear increase from 0.001 to 0.98 over gen2100 generations)
delta.q <- (0.98-0.001)/gen2100
curr.q <- 0.001
for (k in 1:gen2100) {
	new.q <- curr.q[k] + delta.q
	curr.q <- c(curr.q, new.q)
}
qs <- c(0.001,curr.q,rep(0.98,(gen-gen2100-2)))
j.2 <- seq(2, gen+1, 1)
qs.var <- cbind.data.frame(j.2,qs)

##calculate evolutionary trajectories
df.T <- NULL
for(i in 1:length(r12)){ ##this loop controls different values of genetic correlation, creating a trajectory for each
  
  ##create matrices from parameter values, except for G
  omega.P = matrix(c(omega1+P11, 0, 0, omega2+P22), 2)
  inverse = solve(omega.P)
  zbar.mat=matrix(c(zbar1.start, zbar2.start),1)
  start=zbar.mat[1,]
  zbar=zbar.mat[1,]
  G12 <- r12[i]*sqrt(G11*G22)
  G=matrix(c(G11,G12,G12,G22),2)
  prod = G %*% inverse
  P=NULL #clear P for each iteration; P stores each set of points for one trajectory sans the first row

  delta=c(0,0)
  for (j in 2:gen) {
  	  q1=qs.var$qs[j]
  	  q=matrix(c(q1,(1-q1)),2)
	  delta = prod %*% (q*(theta - zbar))
	  zbar=zbar+delta
	  P.1=cbind(j, r12[i],t(zbar))
    P=rbind(P,P.1)
  }
  T.1 <- cbind(1, r12[i],t(start)) #create rows to id gen and r
  T <- rbind(T.1,P)  #add start position to P
  T <- data.frame(T)
  df.T <- rbind(df.T, T) #df.T combines P for each 
}

## create plot
interval = seq(from=10, to=length(df.T$j), by=10)  ##arrowheads every how many generations? Total # of generations run needs to be divisible by this number!
p <- ggplot(df.T, aes(x=V3, y=V4, group=V2, color=V2)) + geom_line(size=1.25)+ theme_bw()+ scale_colour_gradientn(name="Genetic\ncorrelation", colours = c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")) + xlab("Larval mortality,\nStressful Environment")+ ylab("Larval mortality,\nRefugia Environment")  
for (i in 1:length(interval)) {
    p <- p + geom_segment(x= df.T[interval[i]-1,3], y=df.T[interval[i]-1,4], xend=df.T[interval[i],3], yend=df.T[interval[i],4], colour="black", arrow=arrow(type="closed", length=unit(0.2,"cm")))
}
p +  annotate("text", x = theta[1], y = theta[2], label = "+", size=15, color="black")+ xlim(c(0,0.5)) + ylim(c(0,0.5))+theme(legend.position=c(0.85,0.75))+geom_vline(xintercept=0.24, color="grey50", linetype=2)+geom_vline(xintercept=0.305, color="grey20", linetype=2)+annotate("text", x=0.265, y=0.4, label="+50%", color="grey50") +annotate("text", x=0.33, y=0.4, label="+25%", color="grey20")

## return values of interest (% distance evolved to high temp. optimum)
#evolved value in stressful is column V3; in refugia is V4

#subset(df.T, df.T$j==gen2100+1)  ##output for all genetic correlation values
evolved.zbar <- (subset(df.T, df.T$j==gen2100+1)$V3)[4]  ##evolved phenotypic val. in stressful environment when genetic correlation = 0
distance.to.evolve <- theta[1]-zbar1.start
percent.evolved <- (evolved.zbar - zbar1.start)/distance.to.evolve

percent.evolved
