#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(igraph)
library(scales)
library(lubridate)
library(reshape2)


# Define User Interface for the application 
ui <- fluidPage( #fluidPage 

    # Application title
    titlePanel("Sanjay's shiny app"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
            
        sidebarPanel(
                
                fileInput("file1", "Choose CSV File",
                          multiple = TRUE,
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")),
                # Horizontal line ----
                tags$hr(),
                
                # Input: Select number of rows to display ----
                radioButtons("disp", "Display data",
                             choices = c(Head = "head",
                                         Structure = "str"),
                             selected = "head"),
                
                # Horizontal line ----
                tags$hr(),
                
                 sliderInput("timewind",
                         "Time window (seconds)",
                         min = 0,
                         max = 50,
                         value = 0),
                
                #actionButton("netbutton", "Submit", class = "btn-success")
                
                textInput("timecomp",
                          "Time windows to compare",
                          placeholder = "Enter numbers separated by commas"),
                
                downloadButton("downloadData", "Download")
        ),

        # Show a plot of the generated distribution
        mainPanel(
                
                # Output: Data file ----
                tableOutput("headdata"),
                
                verbatimTextOutput("strdata"),
                
                # Output: Network plot ----
                plotOutput("netplot"),
                
                # Output: Individual network strength
                tableOutput("strength")
                
        )
    )
)


#########
#########
#########
#########

# Define server logic (this is where you define the relationship between inputs and outputs)
server <- function(input, output) {

        
        output$headdata <- renderTable({
                
                # input$file1 will be NULL initially. After the user selects
                # and uploads a file, head of that data file by default,
                # or all rows if selected, will be shown.
                
                req(input$file1)
                
                df <- read.csv(input$file1$datapath)#,
                               # header = input$header,
                               # sep = input$sep,
                               # quote = input$quote)
                
                if(input$disp == "head") {
                        return(head(df))
                }
                
        })
        
        output$strdata <- renderPrint({
                req(input$file1)
                
                df <- read.csv(input$file1$datapath)#,
                # header = input$header,
                # sep = input$sep,
                # quote = input$quote)
                
                if(input$disp == "str") {
                        return(str(df))
                }
        })
        
        
        
        output$netplot <- renderPlot({
                
                if (input$timewind != 0) { #need this if statement to make sure a plot only shows up if slider is not on zero
                #input$netbutton #this adds dependency on the action button (button needs to be clicked before this code is run again)
                
                req(input$file1)
                
                df <- readr::read_csv(input$file1$datapath)
                
                df$hour <- lubridate::hour(df$Time)
                df$minute <- lubridate::minute(df$Time)
                df$second <- lubridate::second(df$Time)
                
                # custom rounding function
                mround <- function(x,base){
                        base*round(x/base)
                }
                
                # round seconds to nearest "timewind" seconds
                df$time.round <- mround(df$second, input$timewind) #add isolate() around the stuff on the right side of the arrow to make it so that the network only updates when the submit button is clicked
                
                # force two digits to display
                df$time.round <- sprintf("%02d", df$time.round)
                
                # add a key with hour, minute, second, and antennaID
                df$groupKEY <- paste(df$Date,
                                                   paste(df$hour, 
                                                         df$minute, 
                                                         df$time.round, 
                                                         sep=":"),
                                                   df$Antenna, 
                                                   sep="_")
                
                #count the number of times each individual is at each antenna at each time using the groupKey, then add whether each individual is "present"
                ct.penguinXkey <- df %>% 
                        group_by(Penguin, groupKEY) %>% 
                        summarize(count.pings=n())
                
                # add a column that is just "present" (1) at each time interval at each antenna
                ct.penguinXkey$present <- 1
                
                presence.penguinXkey <- unique(ct.penguinXkey[c("Penguin", "groupKEY", "present")])
                
                
                # Reshape the data into individual X group format
                indivXgrp <- reshape2::dcast(presence.penguinXkey, Penguin~groupKEY)
                
                # convert to matrix format using matrix.please function
                matrix.please <- function(x) {
                        m<-as.matrix(x[,-1])
                        rownames(m)<-x[,1]
                        m
                }
                
                indivXgrp.mx <- matrix.please(indivXgrp)
                indivXgrp.mx[is.na(indivXgrp.mx)] <- 0
                
                
                #project the data to show the number of times each individual is seen at the same time and place (antennaID) as others
                penguinXpenguin <- indivXgrp.mx %*% t(indivXgrp.mx)
                
                #convert to network object and plot
                game.net <- igraph::graph_from_adjacency_matrix(penguinXpenguin, "undirected", weighted=TRUE, diag=FALSE)
                
                # plot
                set.seed(42)
                plot(game.net,
                     
                     #format nodes
                     vertex.color="lightblue", #node fill color
                     vertex.frame.color="white", #node outline color
                     
                     #format edges
                     edge.color=alpha("black", 0.5), #edge color
                     edge.width=E(game.net)$weight*2, #width of edges
                     
                     #format labels
                     vertex.label.cex=1, #font size of labels
                     vertex.label.color="black", #font color
                     vertex.label.dist=2.5, #moves label away from node center
                )
                # add title
                title("Penguin Game: association = present at the same antenna in the same time window", cex.main=0.75)
                }        
        })
        
        
        
        net.summary.data <- reactive({
                
                text = gsub(" ", "", input$timecomp)
                split = strsplit(text, ",", fixed = FALSE)[[1]]
                guide <- as.numeric(split)
               
                
                net.summary <- data.frame(rounding.scheme = character(),
                                          penguin = character(),
                                          degree = numeric(),
                                          strength = numeric())
                
                req(input$file1)
                
                df <- readr::read_csv(input$file1$datapath)
                
                df$hour <- lubridate::hour(df$Time)
                df$minute <- lubridate::minute(df$Time)
                df$second <- lubridate::second(df$Time)
                
                mround <- function(x,base){
                        base*round(x/base)
                }
                
                matrix.please <- function(x) {
                        m<-as.matrix(x[,-1])
                        rownames(m)<-x[,1]
                        m
                }
                
                for (net in 1:length(guide)) { 
                        #net=1
                        loop.round <- guide[net] #what are we rounding by for this loop?
                        
                        #rename data for loop
                        loop.data <- df
                        
                        # split of each part of the time, rounding by whatever is in loop.round
                        loop.data$sec.rounded <- mround(loop.data$second, loop.round)
                        
                        # force two digits to display
                        loop.data$sec.rounded <- sprintf("%02d", loop.data$sec.rounded)
                        
                        # add a key with hour, minute, second, and antennaID
                        loop.data$groupKEY <- paste(loop.data$Date,
                                                    paste(loop.data$hour, 
                                                          loop.data$minute, 
                                                          loop.data$sec.rounded, 
                                                          sep=":"),
                                                    loop.data$Antenna, 
                                                    sep="_")
                        
                        # summarize
                        ct.penguinXkey <- loop.data %>% 
                                group_by(Penguin, groupKEY) %>% 
                                summarize(count.pings=n())
                        
                        # add a column that is just "present" (1) at each time interval at each antenna
                        ct.penguinXkey$present <- 1
                        
                        # extract all unique rows
                        presence.penguinXkey <- unique(ct.penguinXkey[c("Penguin", "groupKEY", "present")])
                        
                        #check
                        #head(presence.penguinXkey)
                        
                        #make penguin by group dataframe then convert to a matrix
                        indivXgrp <- reshape2::dcast(presence.penguinXkey, Penguin~groupKEY)
                        indivXgrp.mx <- matrix.please(indivXgrp) # make into matrix
                        indivXgrp.mx[is.na(indivXgrp.mx)] <- 0  # set any NA cells to 0
                        
                        #multiple matrix by transpose to get penguin by penguin matrix
                        penguinXpenguin <- indivXgrp.mx %*% t(indivXgrp.mx)
                        
                        #convert to network object
                        game.net <- igraph::graph_from_adjacency_matrix(penguinXpenguin, "undirected", weighted=TRUE, diag=FALSE)
                        
                        
                        # summarize network structure
                        
                        #find the density
                        edge_density(game.net)
                        
                        #find degree
                        degree <- degree(game.net)
                        
                        #find strength
                        strength <- strength(game.net)
                        
                        #save penguin degree and strength then write to main dataframe
                        penguinID <- names(strength) #make labels for penguins
                        
                        rounding.scheme <- rep(loop.round, length(penguinID)) #make labels for rounding scheme
                        
                        combo <- cbind.data.frame(penguinID, rounding.scheme, degree, strength) #assemble this loop's data
                        rownames(combo) <- NULL #get rid of row names
                        
                        net.summary <- rbind.data.frame(net.summary, combo) # save this loop's data to the main dataframe
                        
                }
                net.summary # this becomes the object that is assigned to net.summary.data
                #return(net.summary)
                # cast data into wider format
                #net.summary.wide <- reshape2::dcast(net.summary, penguinID~rounding.scheme)
                
                #return(net.summary.wide)
        })
        
        
         output$strength <- renderTable({
                 if (input$timecomp != "") {net.summary.data()}
         })
    
        
        # Downloadable csv of data ----
        output$downloadData <- downloadHandler(
                filename = function() {
                        paste("net_summary", ".csv", sep = "")
                },
                content = function(file) {
                        write.csv(net.summary.data(), file, row.names = FALSE)
                }
        )
}

# Run the application 
shinyApp(ui = ui, server = server)
