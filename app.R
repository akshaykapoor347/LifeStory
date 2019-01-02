library(shiny)
library(datasets)
library(ggplot2)
library(png)
library(corrplot)
library(protViz)
library(LinearRegressionMDE)
library(dplyr)
library(magick)

#setwd('C:/Users/Akshay/Downloads/Data mining LifeStory final')

all_details <- read.csv('alldetailsdata9-10-18(1).csv')

health <- read.csv("health and consumer data(1).csv", na.strings=c(""), header=T, stringsAsFactors = FALSE)

all_details$Peptide <- as.character(all_details$Peptide)
all_details$Modification <- as.character(all_details$Modification)

#peptide analyses
all_details$parent_ion_mass <- parentIonMass(all_details$Peptide)

#mass of modification
all_details$mod_mass <- all_details$Observed.mass..Da. - all_details$parent_ion_mass
min_mod_mass <- min(all_details$mod_mass)
max_mod_mass <- max(all_details$mod_mass)
median_mod_mass <- median(all_details$mod_mass)

a1 <- filter(all_details, mod_mass == min_mod_mass) #min modification
b1 <- filter(all_details, mod_mass == max_mod_mass) #max modification
c1 <- filter(all_details, mod_mass == median_mod_mass) #median modification
peptides <- c(a1$Peptide, b1$Peptide, c1$Peptide)
pim<-parentIonMass(peptides)
fi<-fragmentIon(peptides)

### Health
Tampon <- table(health$Tampons.)
Flow <- table(health$Flow)
Birth <- table(health$If.yes.to.birth.control..please.list.method.and.dose)
Soap <- table(health$Bath.Soap.)

###img


colnames(all_details) <- c("Patient ID","Peptide","Modification","Observed_mass_Da","Observed_RT_min","Observed.M.Z.","Charge","Matched.1st.Gen.Primary.Ions","X_Matched_1st_Gen_Primary.Ions","Assigned_intensity","parent_ion_mass","mod_mass")

a <- colnames(all_details)
a <- a[-c(1,2,3)]

x <- cor(all_details[4:12], method = "pearson")
corrplot(x, type = "lower", order = "hclust", method = "color", insig = "p-value", title = "X")

p.mat <- cor.mtest(x)$p
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))


plotType <- function(x, type) {
  switch(type,
         A = corrplot(x, type = "lower", order = "hclust", method = "color", insig = "p-value", title = "X")
         ,
         B = corrplot(x, method = "color", col = col(200),
                      type = "lower", order = "hclust", number.cex = .7,
                      addCoef.col = "black", # Add coefficient of correlation
                      tl.col = "black", tl.srt = 90, # Text label color and rotation
                      # Combine with significance
                      p.mat = p.mat, sig.level = 0.01, insig = "blank", 
                      # hide correlation coefficient on the principal diagonal
                      diag = FALSE)
  )
}

plotType_p <- function(x,type) {
  switch(type,
         A = barplot(Tampon, main="Hygine", 
                     xlab="Tampons", ylab = "Count", col = "lightcoral" )
         ,
         B = barplot(Flow, main="Flow", 
                     xlab="Flow", ylab = "Count", col = "lightcoral")
         ,
         C =  
           barplot(Birth, main="Birth Control Types", 
                   xlab="method/brand", ylab = "Count", col = "lightcoral")
         ,
         D = 
           barplot(Soap, main="Soap brands", 
                   xlab="Soap brands", ylab = "Count", col = "lightcoral")
         
  )
}

imag <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)



# Define UI for application that draws a histogram
ui <- navbarPage("LifeStory Health",
           tabPanel("Image Comparision",sidebarLayout(      
             
             
             numericInput('c1', 'First image', 1,
                          min = 1, max = 23),
             numericInput('c2', 'Second image', 2,
                          min = 1, max = 23)
             
           ), mainPanel(
             plotOutput("picPlot",  width = "100%")
             
             
             
           )),
           tabPanel("Histogram",
                    sidebarLayout(      
                      
                      # Define the sidebar with one input
                      sidebarPanel(
                        selectInput("region", "Attributes:", 
                                    choices=a),
                        hr(),
                        helpText("Univariate Analysis")
                      ),
                      
                      # Create a spot for the barplot
                      mainPanel(
                        plotOutput("phonePlot")  
                        
                      )
                      
                    )
           ),
           tabPanel("Boxplot",
                    sidebarLayout(      
                      
                      # Define the sidebar with one input
                      sidebarPanel(
                        selectInput("region1", "Attributes:", 
                                    choices=a),
                        hr(),
                        helpText("Univariate Analysis")
                      ),
                      
                      # Create a spot for the barplot
                      mainPanel(
                        plotOutput("boxPlot")  
                      )
                      
                    ) ), tabPanel("Corrplot", radioButtons("pType", "Choose plot type:",
                                                           list("A", "B")),
                                  
                                  mainPanel(
                                    plotOutput("boxPlot1",  width = "100%") ) 
                    ), tabPanel("Peakplot", 
                                mainPanel(
                                  plotOutput("boxPlot2") ) 
                    ),
           tabPanel("Sampleplot", 
                    mainPanel(
                      plotOutput("boxPlot3",  width = "100%") ) 
           ), tabPanel("Menstrual Flow", radioButtons("pType1", "Choose plot type:",
                                                      list("A", "B", "C", "D")),
                       
                       mainPanel(
                         plotOutput("boxPlot4") ) 
           )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  # Fill in the spot we created for a plot
  output$phonePlot <- renderPlot({
    
    hist(all_details[,input$region]*1000,
         main=input$region,
         ylab="Frequency",
         xlab=input$region,
         col="pink",
         border="red")
  })
  
  output$boxPlot <- renderPlot({
    
    boxplot(all_details[,input$region1]*1000,
            main=input$region,
            ylab="Frequency",
            xlab=input$region,
            col="pink",
            border="red")
  })
  
  output$boxPlot1 <- renderPlot({
    
    plotType(x, input$pType)
  }, height = 700, width = 1400 )
  
  output$boxPlot2 <- renderPlot({
    
    peakplot(a1$Peptide, msms[[1]])
  }, height = 700, width = 1400 )
  
  output$boxPlot3 <- renderPlot({
    par(mfrow=c(3,1));
    
    
    for (i in 1:length(peptides)){
      plot(0,0,
           xlab='m/Z',
           ylab='',
           xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
           ylim=c(0,1),
           type='n',
           axes=FALSE,
           sub=paste( pim[i], "Da"));
      box()
      axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))
      pepSeq<-strsplit(peptides[i],"")
      axis(3,fi[i][[1]]$b,pepSeq[[1]])
      
      abline(v=fi[i][[1]]$b, col='red',lwd=2)
      abline(v=fi[i][[1]]$c, col='orange')
      abline(v=fi[i][[1]]$y, col='blue',lwd=2)
      abline(v=fi[i][[1]]$z, col='cyan')
      
    }
  }, height = 700, width = 1400 )
  
  output$boxPlot4 <- renderPlot({
    
    plotType_p(x,input$pType1)
  })
  
  output$picPlot <- renderPlot({
    A <- toString(input$c1)
    B <- toString(input$c2)
    A <- image_read(A)
    B <- image_read(B)
    A <- c(A,B)
    A <- image_scale(A, "5000x2500")
    plot(image_flatten(A, 'Add'))
    
  }, height = 700, width = 1400 )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

