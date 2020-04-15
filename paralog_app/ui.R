library(shiny)
library(shinythemes)
library(shinyjs)
library(shinycssloaders)

# css <- HTML(".pull-left{float: left !important;}
#               .pull-right{float: right !important;}")
# 
# js <- HTML("$(function(){
#         setTimeout(function(){
#            $('.dataTables_filter').addClass('pull-left');
#            $('.dataTables_length').addClass('pull-right');
#            }, 200);
#            });")

fluidPage(
  useShinyjs(),
  # tags$head(
  #   tags$style(css),
  #   tags$script(js)
  #   ),
  theme=shinytheme("cosmo"), # eg. lumen # https://rstudio.github.io/shinythemes/
    navbarPage(
      "PARALOG Annotator",
      tabPanel("Search",
               #h2("Missense Variant Annotation for Inherited Cardiac Conditions",align="center"),
               br(),
               sidebarLayout(
                  sidebarPanel(
                   # img(src = "paralogo2.png", width = "100%"),
                    id = "myapp",
                    
                      #dont need below as now have shinycssloaders
#                     tags$head(tags$style(type="text/css", "#loadmessage {
# position: fixed;
# top: 100px;
# left: 100px;
# width: 100%;
# padding: 5px 0px 5px 0px;
# text-align: center;
# font-weight: bold;
# font-size: 100%;
# color: #000000;
# background-color: #CCFF66;
# z-index: 105;
#              }")),
                    h3("Input your variant"),
                    br(),
                    radioButtons("format",label=NULL,
                                 # Here a new input method can be inserted eg. upload a file with variants
                                 # eg. choices = list("upload file"="upload",
                                 # fileInput("file", NULL,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))),
                                 choices = list("Paste variants"="paste", "Upload Variants"="upload"),
                                 selected = NULL),
                      #NOTE EXAMPLES BELOW NO LONGER WORK AS REAL DATA USES DIF BUILD
                      p("e.g. 1:115256528:T:G, 3:38592567:T:A, or X:70443591:G:A"),

                      conditionalPanel(
                        condition="input.format=='paste'",
                        textAreaInput("var",label=NULL,placeholder = "1:114713907:T:G")
                        ),
                      conditionalPanel(
                        condition="input.format=='upload'",
                        fileInput("file", NULL,accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv"))),
                    #dont need below as now have shinycssloaders
                    # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                    #                  tags$div("Loading...",id="loadmessage")),
                    actionButton("sumbit_button","Submit"),
                    actionButton("reset", "Reset form")
                      ),
                   mainPanel(
                     tabsetPanel(
                       id = "All_results",
                       type = "tabs",
                       tabPanel("Query Variant",
                                h4("If your variant(s) are already documented in ClinVar they will appear here"),
                                conditionalPanel(condition = "input.sumbit_button", withSpinner(dataTableOutput("known_clinvar")))
                                
                       ),
                       tabPanel("Known pathogenic variants in paralogous positions",
                                h4("Equivalent missense variant(s) identified by Paralogue Annotation"),
                                conditionalPanel(condition = "input.sumbit_button", withSpinner(dataTableOutput("paralog"))),
                                conditionalPanel("output.paralog",downloadButton("download","Download"))
                       )
                       
                       )
                     )
                  )),
      tabPanel("About",
      style = "width:80%; margin-right:auto; margin-left:auto", 
      includeHTML("about.html"), # This is an HTML page that is read in from dir 
      br()
      )
      
    )
  )
