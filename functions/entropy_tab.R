###### ENTROPY #############
#########################################
entropy_tab <- tabPanel('Entropy',
         sidebarPanel(
           checkboxInput('tsflag',strong("Time series data"),FALSE),
           conditionalPanel(
             condition = "input.tsflag==true",
             textInput('entropy_timepoints',"Number of time points"),
             helpText("* Please order the sample columns in input file properly. Time series data of the same genotype should be put in adjacent columns.")
           ),
           radioButtons('entropy_graph_type',"Graph type:",
                        c('Bar chart','Line chart')),
           downloadButton("downloadentropy","Download as PNG"),
           conditionalPanel(
             condition = "input.tsflag==true",
             h5("Specify names of the genotypes"),
             uiOutput("expand_genonames_entropy")
           )
         ),
         mainPanel(
           h3('Shannon entropy'),
           plotlyOutput('entropy.plot')
         )
)
