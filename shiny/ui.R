library(shiny)

shinyUI(
  fluidPage(
  #----------------------------------------------
  
  titlePanel("2D IR spectroscopy"),
  
    sidebarLayout(
      
      sidebarPanel(
        p("Please choose the value of the following parameters"),
        sliderInput("axRange", label = h3("X and Y axes range"), min = 950, 
                    max = 1100, value = c(1000, 1065)
                    )
      ),
      
      mainPanel(
        plotOutput("plot2dir",width="100%",height="800px")
      )
      
      
      
    )
  
  #----------------------------------------------
  )
)
