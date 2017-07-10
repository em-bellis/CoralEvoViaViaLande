library(shiny)




# Define UI for miles per gallon application

shinyUI(fluidPage(

  

  # Application title

  titlePanel("Coral Trait Evolution, via Via & Lande (1985)"),

  

  fluidRow(

    column(3, 

      h3("Model Inputs"),

      wellPanel(
    
      selectInput("gen", 

                 label = "Number of generations by 2100", choices = list("5" = 5, "10" = 10, "15"=15, "20" = 20, "25"= 25, "30"=30), selected = 20),
  	  
  	  numericInput("G11", 

             
              label = "Additive genetic variance, stressful environment (G11)", value = 0.033, min=0),    


      numericInput("G22", 

                 label = "Additive genetic variance, refugia (G22)", value = 0.001, min=0),

      numericInput("P11", 

                 label = "Phenotypic variance, stressful environment (P11)", value = 0.052, min=0),

      numericInput("P22", 

                 label = "Phenotypic variance, refugia (P22)", value = 0.009, min=0),

      numericInput("zbar1.start", 

                 label = "Mean trait value of initial population, stressful environment", value = 0.37),

      numericInput("zbar2.start", 

                 label = "Mean trait value of initial population, refugia", value = 0.11),

      numericInput("theta1", 

                 label = "Optimal trait value, stressful environment", value = 0.1),

      numericInput("theta2", 

                 label = "Optimal trait value, refugia", value = 0.1),

      numericInput("omega1", 

                 label = "Strength of stabilizing selection relative to P11, stressful environment (smaller=stronger)", value = 3, min=0, max=99),

      numericInput("omega2", 

                 label = "Strength of stabilizing selection relative to P22, refugia (smaller=stronger)", value = 3, min=0, max=99)

    )),

    column(6,

      plotOutput("evoPlot"), 

      

      p(""),

      p("This plot simulates evolution of coral thermal tolerance traits until the year 2100 based on a quantitative genetic model that describes phenotypic evolution of a single character expressed in two environments (Via & Lande 1985).  It takes into account direct selection acting on the trait in the stressful environment, as well as a correlated response to selection in the refugium environment.  We assume that the proportion of the population inhabiting refugia decreases linearly over time from close to 100% in 2017 to almost 0% in 2100.  The joint phenotypic optimum is indicated by the '+' and arrowheads are plotted every 5 generations.  Note that choosing a value of '3' for the selection strength corresponds to very strong selection with an adaptive landscape 2x as wide as the trait distribution, a value of '99' corresponds to very weak selection with an adaptive landscape 10x the width of the trait distribution, etc."),
      
      p("Additional information regarding the model parameters can be found in X (submitted).  Please contact weissem[at]science[dot]oregonstate[dot]edu with suggestions, questions, or comments."),

      p("X (submitted). Title."),

      p("Via, S. & R. Lande (1985). Genotype-environment interaction and the evolution of phenotypic plasticity. Evolution 39: 505-512.")

    ),

    column (3,

      h3("Model Outputs"),

      wellPanel(

        p("If genetic correlation = +0.99, evolved bleaching response is estimated as"),

        verbatimTextOutput("text"),      

        p("If genetic correlation = -0.99, evolved bleaching response is estimated as"),

        verbatimTextOutput("text2"),

        p("This corresponds to a generation time of "),

        verbatimTextOutput("text3")


      )

    )  

  )

))

