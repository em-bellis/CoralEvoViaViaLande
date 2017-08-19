library(MASS)
library(car)
library(sfsmisc)
library(ggplot2)
library(grid)

##choose fixed parameter values
##env1 = hot/stress env; env2 = control/refugia env
gen=100  #total number of generations to run
theta = c(0.0015,0.0015)
zbar1.start=-0.028
zbar2.start=-0.031
#G11=0.001197
G22=0.006
P11=0.013
P22=0.03
#omega1=P11
#omega2=P22
r12=c(0.999, 0.75, 0.375, 0, -0.375, -0.75, -0.999)   #genetic correlation

###function
traject <- function(omega1, omega2, G11, gen2050, gen2100) {
	##rcp 2.6 (linear increase from 0.001 to 0.98 until 2050, then remain at 0.98 for remaining generations)
	delta.q <- (0.98-0.001)/gen2050
	curr.q <- 0.001
	for (k in 1:gen2050) {
		new.q <- curr.q[k] + delta.q
		curr.q <- c(curr.q, new.q)
	}
	qs <- c(0.001,curr.q,rep(0.98,(gen-gen2050-2)))
	j.2 <- seq(2, gen+1, 1)
	qs.var <- cbind.data.frame(j.2,qs)
	
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
	#evolved phenotypic value in stressful is column V3; in refugia is V4
	#subset(df.T, df.T$j==gen2100+1)  ##output for all genetic correlation values
	evolved.zbar <- (subset(df.T, df.T$j==gen2100+1)$V3)[4]  ##evolved phenotypic val. in stressful environment when genetic correlation = 0
	distance.to.evolve <- theta[1]-zbar1.start
	percent.evolved <- (evolved.zbar - zbar1.start)/distance.to.evolve
	
	return(percent.evolved)
}

gen2050s.to.run <- c(2,4,6,8,10,12)
gen2100s.to.run <- c(5,10,15,20,25,30)
#omega1.to.run <- c(0*P11, 3*P11, 8*P11, 15*P11, 24*P11, 35*P11, 48*P11, 63*P11, 80*P11, 99*P11) #corresponds to adaptive landscape 1..10x wider than character distribution 
#omega2.to.run <- c(0*P22, 3*P22, 8*P22, 15*P22, 24*P22, 35*P22, 48*P22, 63*P22, 80*P22, 99*P22) #corresponds to adaptive landscape 1..10x wider than character distribution
#h2.to.run <- c(0.9*P11, 0.8*P11, 0.7*P11, 0.6*P11, 0.5*P11, 0.4*P11, 0.3*P11, 0.2*P11, 0.1*P11, 0*P11) #corresponds to narrow-sense heritability of 0, 0.1, 0.2 ... 0.9
omega1.to.run <- c(0*P11, 1*P11, 2*P11, 3*P11, 4*P11, 5*P11, 6*P11, 7*P11, 8*P11, 9*P11, 10*P11, 11*P11, 12*P11, 13*P11, 14*P11, 15*P11, 16*P11, 17*P11, 18*P11, 19*P11, 20*P11)
omega2.to.run <- c(0*P22, 1*P22, 2*P22, 3*P22, 4*P22, 5*P22, 6*P22, 7*P22, 8*P22, 9*P22, 10*P22, 11*P22, 12*P22, 13*P22, 14*P22, 15*P22, 16*P22, 17*P22, 18*P22, 19*P22, 20*P22) #corresponds to adaptive landscape 1..10x wider than character distribution 
h2.to.run <- c(0.95*P11, 0.9*P11, 0.85*P11, 0.8*P11, 0.75*P11, 0.7*P11, 0.65*P11, 0.6*P11, 0.55*P11, 0.5*P11, 0.45*P11, 0.4*P11, 0.35*P11, 0.3*P11, 0.25*P11, 0.2*P11, 0.15*P11, 0.1*P11, 0.05*P11, 0.0*P11) #corresponds to narrow-sense heritability of 0, 0.1, 0.2 ... 0.9


df.U <- NULL
for (l in 1:length(gen2050s.to.run)) {
	for (m in 1:length(omega1.to.run)) {
		for (n in 1:length(h2.to.run)) {
			corr0 <- traject(omega1.to.run[m],omega2.to.run[m],h2.to.run[n],gen2050s.to.run[l],gen2100s.to.run[l])
			U <- cbind(omega1.to.run[m]/P11,omega2.to.run[m]/P22,h2.to.run[n]/P11,gen2100s.to.run[l], corr0)
			df.U <- rbind(df.U, U)
		}
	}
}

df.U <- as.data.frame(df.U)
colnames(df.U) <- c("WidthSS_env1","WidthSS_env2","h2_env1","GensBy2100","PercentEvolved")

##data frame is ready, now to make heatmap
library(RColorBrewer)
library(pheatmap)
library(tidyr)
library(viridis

#display.brewer.all()
col.pal <- inferno(20)
#fixInNamespace("draw_colnames","pheatmap")  ##change hjust to 0.5 and rot to 0 to rotate x axis labels

##gen 5
gen.df <- subset(df.U, GensBy2100==5)
gen.mat <- as.matrix(xtabs(PercentEvolved~h2_env1+WidthSS_env1, data=gen.df))
gen.mat2 <- apply(gen.mat, 2, rev)

pheatmap(gen.mat2, cluster_row=F, cluster_cols = F, color=col.pal, legend = F, fontsize_row=5, fontsize_col = 5, labels_col=c(1,2,3,4,5,6,7,8,9,10), filename="/Users/Weissem/Desktop/CRA_draft/Past/Past5.pdf", width=1, height=1)

##gen 10
gen.df <- subset(df.U, GensBy2100==10)
gen.mat <- as.matrix(xtabs(PercentEvolved~h2_env1+WidthSS_env1, data=gen.df))
gen.mat2 <- apply(gen.mat, 2, rev)

pheatmap(gen.mat2, cluster_row=F, cluster_cols = F, color=col.pal, legend = F, fontsize_row=5, fontsize_col = 5, labels_col=c(1,2,3,4,5,6,7,8,9,10), filename="/Users/Weissem/Desktop/CRA_draft/Past/Past10.pdf", width=1, height=1)

##gen 15
gen.df <- subset(df.U, GensBy2100==15)
gen.mat <- as.matrix(xtabs(PercentEvolved~h2_env1+WidthSS_env1, data=gen.df))
gen.mat2 <- apply(gen.mat, 2, rev)

pheatmap(gen.mat2, cluster_row=F, cluster_cols = F, color=col.pal, legend = F, fontsize_row=5, fontsize_col = 5, labels_col=c(1,2,3,4,5,6,7,8,9,10), width=1, height=1)

##gen 20
gen.df <- subset(df.U, GensBy2100==20)
gen.mat <- as.matrix(xtabs(PercentEvolved~h2_env1+WidthSS_env1, data=gen.df))
gen.mat2 <- apply(gen.mat, 2, rev)

pheatmap(gen.mat2, cluster_row=F, cluster_cols = F, color=col.pal, legend = F, fontsize_row=5, fontsize_col = 5, labels_col=c(1,2,3,4,5,6,7,8,9,10), width=1, height=1)

#etc...
