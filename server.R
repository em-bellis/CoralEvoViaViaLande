library(shiny)

library(MASS)

library(car)

library(sfsmisc)

library(ggplot2)

library(gridExtra)




r12=c(0.999, 0.75, 0.375, 0, -0.375, -0.75, -0.999)

df.T <- NULL

shinyServer(function(input, output) {

  evolvedResp <- reactiveValues()

  output$evoPlot <- renderPlot({

    df.T <- NULL
    totalgen = as.numeric(input$gen)

    for(i in 1:length(r12)){ ##this loop controls different values of genetic correlation, creating a trajectory for each
	
      ##create matrices from parameter values

      omega1inP = input$omega1*input$P11
      
      omega2inP = input$omega2*input$P22
      
      omega.P = matrix(c(omega1inP+input$P11, 0, 0, omega2inP+input$P22), 2)

      inverse = solve(omega.P)

      zbar.mat=matrix(c(input$zbar1.start, input$zbar2.start),1)

      start=zbar.mat[1,]

      zbar=zbar.mat[1,]

      theta = c(input$theta1,input$theta2)

      G12 <- r12[i]*sqrt(input$G11*input$G22)

      G=matrix(c(input$G11,G12,G12,input$G22),2)

      prod = G %*% inverse

      P=NULL #clear P for each iteration; P stores each set of points for one trajectory sans the first row

      ##q matirx
      
      delta.q <- (0.98-0.001)/totalgen

	  curr.q <- 0.001

	  for (k in 1:totalgen) {

		new.q <- curr.q[k] + delta.q

		curr.q <- c(curr.q, new.q)

	  }

	  qs <- c(0.001,curr.q)

      j.2 <- seq(1, totalgen+2, 1)

     qs.var <- cbind.data.frame(j.2,qs)

      

      delta=c(0,0)

      for (j in 2:totalgen+2) {

        q1 = qs.var$qs[j]
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

    evolvedResp$a <- df.T$V3[totalgen]

    evolvedResp$b <- df.T$V4[totalgen]

    evolvedResp$c <- df.T$V3[totalgen*7]

    evolvedResp$d <- df.T$V4[totalgen*7]

    

    ### create plot

    interval = NULL

    for (h in 0:(length(r12)-1)) {

    	interval.x = seq(from=5+(totalgen*h), to=totalgen+(totalgen*h), by=5)

  		interval = c(interval, interval.x)

 	}

    


    p <- ggplot(df.T, aes(x=V3, y=V4, group=V2, color=V2)) + geom_line(size=1.25)+ theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=18))+scale_colour_gradientn(name="Genetic correlation", colours = c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")) + xlab("Trait value\nStressful environment")+ ylab("Trait value\nRefugia") 

    for (i in 1:length(interval)) {

      p <- p + geom_segment(x= df.T[interval[i]-1,3], y=df.T[interval[i]-1,4], xend=df.T[interval[i],3], yend=df.T[interval[i],4], colour="black", arrow=arrow(type="closed", length=unit(0.2,"cm")))

    }

    p +  annotate("text", x = theta[1], y = theta[2], label = "+", size=15, color="black") 

  })

 

  output$text <- renderText({

    paste("(",round(evolvedResp$a,2),",", round(evolvedResp$b,2),").")

 })

 

 output$text2 <- renderText({

   paste("(",round(evolvedResp$c,2),",", round(evolvedResp$d,2),").")

 })

 

 output$text3 <- renderText({

   paste(83/as.numeric(input$gen), "years")

 })

})
