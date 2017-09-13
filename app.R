library(shiny)
library(ggplot2)

server <- function(input, output, session) {
  
  # Listing possible options for variable controls
  choices <- reactiveValues( 
    xVar = c("Random spread","Alphabet"),
    yVar = c("AllelicDiv","GCContent","RatioCount"),
    cVar = c("Labelled loci","Category","Missing","Multiple copies","AvgLength")
  )
  
  
  observe({ # Dynamically updating possible options for variable controls
    if ( "Category" %in% colnames(df.all()) )
    { choices$xVar = c("Random spread","Alphabetical","Category") }
    
    if ( "VSitesNuc" %in% colnames(df.all()) )
    { choices$yVar = c("AllelicDiv","GCContent","RatioCount", "VSitesNuc", "VSitesAA", "RatioVS") }
    
    updateSelectInput( session = session, inputId = "xVar", choices = choices$xVar )
    updateSelectInput( session = session, inputId = "yVar", choices = choices$yVar )
    updateSelectInput( session = session, inputId = "cVar", choices = choices$cVar )
  })
  
  
  
  # reading the table in to one whole data frame, including popping out isolate count from header
  df.all <- reactive({
    inFile <- input$resultstable
    
    validate( need(input$resultstable != "", "Please upload a table on the left") )
    
    # getting the total number of isolates scanned from the header of the results file
    headtext <- scan(inFile$datapath, what = "complex", sep = "\n", n = 5)
    
    isocount$total = as.numeric( strsplit(headtext[4], " ")[[1]][length(strsplit(headtext[4], " ")[[1]])] )
    
    # putting in the whole table to "df.all"
    read.table(inFile$datapath, # open the file once it's in
               header=TRUE, sep="\t", # tab-separated with a header
               quote='"', skip = 5)
  })
  
  # subselecting the data frame depending on level of exclusion
  df <- reactive({
    df <- df.all()
    
    locount$total <- nrow(df) # the count of loci is the count of rows in the main table 
    
    if ( input$percexc != 0 ) # if not excluding any loci
    {
      isocount$cutoff <- floor((1-(input$percexc/100))*isocount$total) # cutoff = % * total of isolates, rounded down 
      
      # exclude from data frame all those which are worse off that the cutoff value
      df <- df[df$Missing < isocount$cutoff,]
      validate( need(nrow(df) != "0", "You removed all loci. Please adjust the slider on the left.") )
    }
    
    locount$filtered <- locount$total - nrow(df) 
    
    df$randomX <- sample(seq(from=0, to=1, by=.001), size = nrow(df), replace = TRUE)
    
    # Removing excluded points
    if ( !is.null(labeltag$removed) )
    {  df <- df[ ! df$Locus %in% labeltag$removed, ]  }
    
    return ( df )
    
  })
  
  # text output that describes total number of loci in data table and number excluded / remaining
  
  # establishing data frames with integer values 
  locount  <- reactiveValues( total = 0, filtered = 0, removed = 0, remain = 0 )
  isocount <- reactiveValues( total = 0, cutoff  = 0, remaining = 0 )
  
  output$exclusions <- renderText({
    locount$removed    <- locount$total - locount$remain - locount$filtered
    locount$remain     <- nrow(df()) # remaining loci are total minus the removed from above
    isocount$remaining <- isocount$total - isocount$cutoff
    
    line1 <- paste('Isolates reviewed:', isocount$total)
    line2 <- paste('Loci currently in analysis:', locount$remain)
    line3 <- paste(' ')
    line4 <- paste('Total loci:', locount$total)
    line5 <- paste('Cutoff of',input$percexc,'% filters out loci in fewer than',isocount$remaining,'isolates.')
    line6 <- paste('Loci filtered out:', locount$filtered)
    line7 <- paste('Loci removed manually:', locount$removed)
    
    
    paste (line1, line2, line3, line4, line5, line6, line7, sep="\n")
  })
  
  
  
  # selects the data that will be used for the graph, taking into account categ, yVar and zscore
  df.sel <- reactive ({
    
    df <- df() # making visible copies of the relevant large data table
    
    # make the x-variable into either the category or a random number between 0 and 1 for visibility
    if ( input$xVar == "Random spread" ) {
      x1 <- df$randomX }
    
    if ( input$xVar == "Alphabetical" ) {
      x1 <- as.numeric(as.factor( df$Locus )) }
    
    if ( input$xVar == "Category" ) {
      x1 <- df$Category }
    
    
    # make the y-variable whatever the selected variable is, with z-score scaling if necessary
    if ( input$yVar == "AllelicDiv" ) { y1 <- df$AllelicDiv }
    if ( input$yVar == "GCContent" )  { y1 <- df$GCContent }
    if ( input$yVar == "VSitesNuc" )  { y1 <- df$VSitesNuc }
    if ( input$yVar == "VSitesAA" )   { y1 <- df$VSitesAA }
    if ( input$yVar == "RatioCount" ) { y1 <- df$RatioCount }
    if ( input$yVar == "RatioVS" )    { y1 <- df$RatioVS }
    
    if ( input$logyVar == TRUE )  { y1 <- log ( y1 ) }
    if ( input$zscore  == TRUE )  { y1 <- scale( y1, center = TRUE, scale = TRUE) }
    
    # make the selected data frame
    df.sel <- data.frame( df$Locus, df$Missing, df$MultipleCopies, df$AvgLength, x1, y1)
    colnames(df.sel) <- c( "Locus", "Missing", "MultipleCopies", "AvgLength", "xsel", "ysel")
    
    # including Category if it was in the results file
    if ( "Category" %in% colnames(df()) == TRUE ) { df.sel$Category <- df$Category }
    if ( "Name" %in% colnames(df()) == TRUE )     { df.sel$Name     <- df$Name }
    
    return ( df.sel )
  })
  
  # output of the data table, with exclusions already removed, and labelled points on column
  output$data.table <- renderDataTable({
    df <- df()
    df$Label<- "No"
    df$Label[df$Locus %in% labeltag$list] <- "Yes"
    return ( df )
  })
  
  # plots distribution of Missing values across entire dataset
  output$distplot <- renderPlot({
    
    # make a copy of the table that will be used
    df <- df.all()
    
    p <- ggplot( df, aes(x=Missing) ) +
      geom_histogram( binwidth=(isocount$total/40) ) + # so that there are always 30 bins
      geom_vline( xintercept=ifelse(input$percexc,isocount$cutoff,isocount$total),
                  size=1, colour="red", linetype="dashed") +
      coord_cartesian( xlim = c( 0, isocount$total ) ) +
      theme_minimal() +
      ggtitle("Distribution of loci by number of isolates for which there was has no known sequence")
    
    return ( p )
  })
  
  # plots Missing vs Allelic Div
  output$corrplot <- renderPlot({
    df <- df()
    
    ggplot( df, aes (x = Missing, y = AllelicDiv) )+
      geom_point() +
      geom_smooth( span = .01*locount$total ) + # so that span scales with n
      coord_cartesian( xlim = c( 0, isocount$total ) ) +
      theme_minimal() +
      ggtitle("Regression of 'missing' designation and allelic diversity value per locus")
    
  })
  
  brange   <- reactiveValues ( xmin = NULL, ymin = NULL, xmax = NULL, ymax = NULL )
  labeltag <- reactiveValues ( list = NULL, removed = NULL, use = FALSE )
  
  observeEvent(input$resetlabel, {
    labeltag$list <- NULL
    labeltag$use <- FALSE
  })
  
  observeEvent(input$resetremovals, {
    labeltag$removed <- NULL
  })
  
  observeEvent(input$removelabel, {
    labeltag$removed <- append(labeltag$removed, labeltag$list)
  })
  
  
  
  # # if a double click happens, add the brushed points' Locus ID to labeltag$list
  observeEvent( { input$mp_dblclick }, {
    
    brush <- input$mp_brush
    
    if (!is.null(brush)) {
      brange$xmin <- brush$xmin
      brange$ymin <- brush$ymin
      brange$xmax <- brush$xmax
      brange$ymax <- brush$ymax
    } else {
      brange$xmin <- NULL
      brange$ymin <- NULL
      brange$xmax <- NULL
      brange$ymax <- NULL
    }
    
    df.sel <- df.sel()
    
    # revise the labeltag list into non-duplicated of appended
    if (input$xVar == "Category") {
      templist <- append ( labeltag$list,
                           as.character( subset ( df.sel,
                                                  Category %in% levels(df.sel$Category)[round(brange$xmin):round(brange$xmax)] &
                                                    ysel  >  brange$ymin &
                                                    ysel  <  brange$ymax )[,"Locus"])) }
    else {
      templist <- append ( labeltag$list,
                           as.character( subset (df.sel,
                                                 ysel  >  brange$ymin &
                                                   ysel  <  brange$ymax &
                                                   xsel  >  brange$xmin &
                                                   xsel  <  brange$xmax)[,"Locus"])) }
    
    
    labeltag$list <- templist[!duplicated(templist)]
    templist <- NULL
    
    # turn on the added label layer only if there are points to label, avoids error messages
    if (!is.null(labeltag$list))
    { labeltag$use <- TRUE }
  })
  
  # the main plot of the whole thing
  scatterPlot <- reactive({
    
    # make a copy of the table that will be used
    df.sel <- df.sel()
    
    # add the colour variable
    if ( input$cVar == "Labelled loci" ) {
      df.sel$csel <- ifelse(df.sel$Locus %in% labeltag$list, "Labelled", "Non-labelled")
      df.sel$csel <- ordered(df.sel$csel, levels = c("Non-labelled", "Labelled"))
      colourtitle <- "Labels"
    }
    
    if ( input$cVar == "Category" ) {
      df.sel$csel <- df.sel$Category
      colourtitle <- "Category of the locus" }
    
    if ( input$cVar == "Missing" ) {
      df.sel$csel <- df.sel$Missing / isocount$total
      colourtitle <- "Percentage of isolates where locus is not tagged" }
    
    if ( input$cVar == "Multiple copies" )    {
      df.sel$csel <- df.sel$MultipleCopies / isocount$total
      colourtitle <- "Percentage of isolates where locus has multiple alleles" }
    
    if ( input$cVar == "AvgLength" ) {
      df.sel$csel <- df.sel$AvgLength
      colourtitle <- "Average length of locus" }
    
    # the actual plot!
    p <- ggplot(df.sel) +
      geom_point( data = df.sel,
                  aes( x=xsel, y=ysel, colour = csel ),
                  position = position_jitter(width = ifelse(input$xVar=="Category",0.3,0)),
                  size=5, alpha=.5 ) +
      ggtitle("Genetic diversity of loci") +
      xlab("") +
      theme_minimal() +
      theme( axis.text.x  = element_text(size=14),
             axis.ticks   = element_line(0.5),
             plot.title   = element_text(face="bold"),
             axis.title.x = element_text(vjust=-.5, size=14),
             legend.position = "bottom", legend.key.width = unit(50, "pt") )
    
    # Average, running average or category averages
    if ( input$xVar == "Random spread" )
    { p <- p + geom_hline( yintercept = mean(df.sel$ysel),
                           size=1, colour="royalblue4", linetype="dashed" ) +
      theme( axis.text.x = element_blank() ) }
    
    if ( input$xVar == "Alphabetical" )
    { p <- p + geom_smooth( data = df.sel, aes( x=xsel, y=ysel ), method="loess", span=0.1, se=FALSE,
                            colour="royalblue4", size=1, linetype="dotted" ) }
    
    if ( input$xVar == "Category" )
    { p <- p + geom_boxplot( data = df.sel, aes( x=xsel, y=ysel, alpha=1 ),
                             outlier.shape = NA ) +
      theme(axis.text.x = element_text(face="bold", angle=45 )) }
    
    
    
    # If the colour scale is a numeric one (not categorical)
    if ( is.numeric(df.sel$csel) ) {
      p <- p + guides(size = FALSE, alpha = FALSE, fill = FALSE,
                      colour = guide_colorbar(title.position = "top") ) +
        scale_color_gradient2(       low="yellow",
                                     mid="red",
                                     high="blue",
                                     midpoint = median(df.sel$csel),
                                     guide_legend(title = colourtitle))
    }
    else {  p <- p + guides(size = FALSE, alpha = FALSE, fill = FALSE) + 
      scale_colour_discrete(guide_legend(title = colourtitle)) }
    
    # Display locus ID/Name for labelled points
    if ( labeltag$use )
    { p <- p +
      geom_text(  data = subset (df.sel, Locus %in% labeltag$list),
                  size=4, alpha=.8, vjust=-.5, angle = 30,
                  aes( x = xsel, y = ysel,
                       label = as.character(Locus) ) )  }
    
    
    # Significance lines
    if ( input$logyVar )
    { ysel <- df.sel$ysel
    p <- p + geom_hline( yintercept =  qnorm (0.05 / length(ysel)) * sd(ysel) + mean(ysel),
                         size = 1, color = "seagreen4", linetype = "dotted") + 
      geom_hline( yintercept = -qnorm (0.05 / length(ysel)) * sd(ysel) + mean(ysel),
                  size = 1, color = "seagreen4", linetype = "dotted")
    }
    
    
    # # putting in correct yVar label depending on variable chosen
    if ( input$yVar == "AllelicDiv" )
    { p <- p + ylab("Alleles per nucleotide") }
    
    if ( input$yVar == "GCContent" )
    { p <- p + ylab("Proportion of G or C bases in locus") }
    
    if ( input$yVar == "VSitesNuc" )
    { p <- p + ylab("Proportion of sites in nucleotide alignment which show any variation") }
    
    if ( input$yVar == "VSitesAA" )
    { p <- p + ylab("Proportion of sites in nucleotide alignment which show any variation") }
    
    if ( input$yVar == "RatioCount" )
    { p <- p + ylab("Ratio of unique nucleotide to unique amino acid sequences per locus") }
    
    if ( input$yVar == "RatioVS" )
    { p <- p + ylab("Ratio of variable sites in nucleotide to amino acid format per locus") }
    
    return ( p )
  })
  
  # what actually does the plot
  output$mainplot <- renderPlot({
    return ( scatterPlot() )
  })
  
  # prints the download button and knows what file to make available
  output$downloadPlot <- downloadHandler(
    filename = function()
    { "GUAVA-Plot.pdf" },
    content = function(file)
    { ggsave(file, width = 8, height = 8, plot = scatterPlot(), device = "pdf") }
  )
  
  output$downloadTable <- downloadHandler(
    filename = function()
    { "GUAVA-DataTable.txt" },
    content = function(file)
    { write.table( df(), file, sep = "\t", row.names = FALSE ) }
  )
  
  
}

ui <- fluidPage(
  titlePanel("GUAVA Visualiser"),
  
  # the whole thing is on the "Side Bar" style layout
  sidebarLayout(
    
    # the control panel on the "Side" of the layout
    sidebarPanel(
      
      h3("Control Panel"), 
      helpText( a("For examples and help, click here.",
                  href="https://github.com/sohauck/GUAVA",
                  target = "_blank") ),
      tags$hr(), 
      
      # First heading
      h4("Upload your GUAVA table"),
      
      fileInput('resultstable', 'Choose ResultsTable file',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      # slider from 0-100 to choose how strict to place tagging cut off
      sliderInput("percexc", "Exclude loci where percentage of isolates with allele designation is below this percentage:", 
                  min = 0, max = 100, value = 5, step = .1),
      
      # blurb of text explaining how many isolates, where cut off is, how many loci total and remaining, etc.
      tagAppendAttributes(textOutput("exclusions"), style="white-space:pre-wrap;"),
      
      # some white space
      tags$hr(),
      
      # Second heading
      h4("Plotting options"),
      
      # Horinzontal axis
      selectInput(inputId = "xVar",
                  label = "Variable plotted on the horizontal x-axis:",
                  choices = c("Random spread","Alphabetical"),
                  selected = "Random spread"),
      
      # Vertical axis
      selectInput(inputId = "yVar",
                  label = "Variable plotted on the vertical y-axis:",
                  choices = c("AllelicDiv","GCContent","RatioCount"),
                  selected = "AllelicDiv"),  
      
      # Colour
      selectInput(inputId = "cVar",
                  label = "Variable plotted on the point colours:",
                  choices = c("Labelled loci","Category","Missing","Multiple copies","AvgLength")),
      
      checkboxInput("zscore",  label = "Use z-scores for y-axis values", value = FALSE),
      checkboxInput("logyVar", label = "Use logarithm of y-axis values, and include Bonferroni-corrected (p < 0.05) significance threshholds", value = FALSE),
      
      actionButton('resetlabel',    'Reset labels'), 
      actionButton('resetremovals', 'Reset  removals'),
      actionButton('removelabel',   'Remove labelled points'),
      
      tags$hr(), # white space
      
      # Download options, either Scatter Plot or Table depending on current tab
      
      conditionalPanel(condition = "input.conditionaltab == 1",
                       downloadButton('downloadPlot', 'Download Plot')
      ),
      
      conditionalPanel(condition = "input.conditionaltab == 3",
                       downloadButton('downloadTable', 'Download Table')
      )
      
      
      
    ), # closes sidebar
    
    # Main display panel in the layout
    mainPanel(
      
      # Multiple tabs as display options
      tabsetPanel(id = "conditionaltab", type = "tabs",
                  
                  # tab #1: main plot, including double-click and brush info          
                  tabPanel("Scatter Plot", value = "1",
                           plotOutput(outputId = "mainplot", height = 700,
                                      dblclick = "mp_dblclick",
                                      brush = brushOpts(
                                        id = "mp_brush"))
                  ),
                  
                  # tab #2: two graphs on how to set the exclusion sliders
                  tabPanel("Excluding loci", value = "2",
                           helpText("Choose what percentage of isolates a locus can be marked 'missing' (due to no allele designation)",
                                    "in and still included in the analysis. In other words, a max 'missing' percentage.",
                                    "Use the slider on the left to choose a cutoff that excludes loci that have been",
                                    "measured as having low diversity (low y-values in the bottom plot) when they have",
                                    "high 'missing' counts (high x-values in the same). The upper plot shows the distribution",
                                    "of loci over the range of 0 to the total count of isolates, where the red-dotted line indicated that all",
                                    "loci to its right have been excluded from all other graphs and tables in this application."),
                           plotOutput(outputId = "distplot", height = 400),
                           plotOutput(outputId = "corrplot")),
                  
                  # tab #3: big interactive data table with all the data
                  tabPanel("Data Table", value = "3",
                           dataTableOutput("data.table"))
                  
      ) # closes tabsetPanel
      
    ) # closes main panel
    
  ) # closes side bar layout
  
) #closes fluidPage & shinyUI


shinyApp(ui = ui, server = server)