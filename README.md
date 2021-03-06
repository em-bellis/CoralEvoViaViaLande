# CoralEvoViaViaLande

This repository contains code associated with computer simulations and figures presented in X (submitted).  We simulated trajectories of coral thermal tolerance trait evolution based on a quantitative genetic model detailed by [Via & Lande (1985)](http://www.jstor.org/stable/2408649) which describes evolution of a single character expressed in two environments.  Our extension of this model projects trait evolution assuming a spatially heterogenous environment, where part of the population experiences 'refugia' conditions characterized by temperature regimes similar to the present day and part of the population experiences 'stressful' conditions characterized by mean temperatures between 1 and 4 degrees higher than current summer maxima.  We model two scenarios: 1) a linear decrease in the proportion of the population experiencing refugium temperatures until 2100 (RCP scenarios 6.0, 8.5) and 2) a linear decrease in the proportion of the population experiencing refugium temperatures until 2050, at which point the majority of the population is present in the stressful environment (RCP scenario 2.6).  Thus, the model incorporates the direct response to selection in the stressful environment and an indirect response to selection experienced by individuals in the refugium environment.  Included are scripts for:

#####Evolutionary trajectory simulation:
Original code to calculate evolutionary trajectories was shared by Steve J. Arnold.  `evoTraj_noX11_2.6.R` and `evoTraj_noX11_6.0or8.5.R`are modified versions that do not use X11 and incorporate variation in the proportion of the population in each environment (q) over time.  They depict simulated evolutionary trajectories for trait evolution in two environments, as shown in Figure S1.

If you want to run the `evoTraj_noX11_6.0or8.5.R` script via a user-friendly web app and try changing different parameter values yourself, there is a Shiny web app currently hosted at https://embellis.shinyapps.io/CoralEvoViaViaLande based on the `ui.R` and `server.R` files.

#####Matrix of simulation results:
`SimMat_6.0or8.5.R` and `SimMat_2.6.R` essentially run `evoTraj_noX11_6.0or8.5.R` or `evoTraj_noX11_2.6.R` for several different generation times, selection regimes, and values of heritability and then plot the simulated change in population mean trait value by the year 2100 as a heatmap.  For these particular simulations, we assume a genetic correlation of 0, though other values can be explored by changing the line `evolved.zbar <- (subset(df.T, df.T$j==gen2100+1)$V3)[4]` to pick other components of the vector.

Note: code for an earlier version of this model (which does not implement variable q over time) is available in [another github repository](https://github.com/em-bellis/bleachingevolution).

#####Running an interactive web app for evolutionary trajectory simulations:
`ui.R` and `server.R` contain code needed to create a user-friendly web app to run simulations via [Shiny for R](http://shiny.rstudio.com).  There is a 25 active hour monthly limit for my shinyapps.io account, so if the site is not working, it may be necessary to launch the app locally.  To do this, you will need [RStudio](https://www.rstudio.com), the 'shiny' package for R, and the `ui.R` and `server.R` files. Once everything is downloaded and installed, run the following commands from the RStudio terminal, replacing `~/path/to/shinyappfolder` with the path to the folder that contains `ui.R` and `server.R`:
```
> library(shiny)
> runApp("~/path/to/shinyappfolder")
```
You may also need to install the following packages: MASS, car, sfsmisc, ggplot2, and gridExtra.
