noise_tab <- tabPanel('Noise',
         sidebarPanel(
           radioButtons('noise_situation',"Select desired noise plot between",choices = c('replicates'='a','genotypes (average of replicates)'='b','genotypes (no replicate)'='c')),
           conditionalPanel(
             condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
             textInput('noise_numOfRep',"Number of replicates",value=1),
             helpText("* Please order the sample columns in input file properly. Replicates of the same genotype should be put in adjacent columns.")
           ),
           conditionalPanel(
             condition = "input.noise_situation=='b'",
             uiOutput("noise_anchor_choices")
           ),
           conditionalPanel(
             condition = "input.noise_situation=='c'",
             selectInput('noise_anchor_c',"Anchor genotype",choices = "")
           ),
           radioButtons('noise_graph_type',"Graph type:",
                        c('Bar chart','Line chart')),
           downloadButton("downloadnoise","Download as PNG"),
           actionButton('noise_plot',"Plot",width='65px',style="color: #fff; background-color: #337ab7; border-color:#337ab7;float:right"),
           conditionalPanel(
             condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
             h5("Specify names of the genotypes"),
             uiOutput("expand_genonames_noise")
           )
         ),
         mainPanel(
           conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                            div(img(src="load.gif",width=240,height=180),h4("Processing ... Please wait"), style="text-align: center;")
           ), 
           conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                            plotlyOutput('noise.plot')
           )
         ))