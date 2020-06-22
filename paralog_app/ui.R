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
  theme=shinytheme("yeti"), # eg. cosmo # https://rstudio.github.io/shinythemes/
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
    navbarPage(
      #title = "PARALOG Annotator DEMO version 0.2.1",
      title = "PARALOG Annotator DEMO version 0.2.2",
      
      id = "navbar",
      tabPanel("Search",
               #h2("Missense Variant Annotation for Inherited Cardiac Conditions",align="center"),
               br(),
               #"test",
               sidebarLayout(
                  sidebarPanel(
                   # img(src = "paralogo2.png", width = "100%"),
                    id = "myapp",
                    width = 2,
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
                    h5("Genome build GRCh37"),
                    br(),
                    radioButtons("format",label=NULL, 
                                 # Here a new input method can be inserted eg. upload a file with variants
                                 # eg. choices = list("upload file"="upload",
                                 # fileInput("file", NULL,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))),
                                 #choices = list("Paste variants"="paste", "Upload Variants"="upload"),
                                 choiceNames = list(
                                   #h5("Paste Variants"),
                                   HTML("<div style='font-size:14px'>Paste Variants</div>"),
                                   #h5("Upload Variants")
                                   HTML("<div style='font-size:14px'>Upload Variants</div>")
                                 ),
                                 choiceValues = list(
                                   "paste", "upload"
                                 ),
                                 selected = NULL,
                                 width = "100%"),
                      #NOTE EXAMPLES BELOW NO LONGER WORK AS REAL DATA USES DIF BUILD
                      #HTML("e.g. <br>1:115256528:T:G<br>3:38592567:T:A<br>X:70443591:G:A<br>"),
                      HTML("e.g. <br>1:115256528:T:C<br>1:115256528:T:G<br>3:38592567:T:A<br>21:44592214:C:T<br>21:47421902:G:A<br>X:70443591:G:A<br>"),

                      #textOutput('text1'),
                      #tags$head(tags$style("#text1{color: red; font-size: 12px; }")), 
                      conditionalPanel(
                        condition="input.format=='paste'",
                        textAreaInput("var",label=NULL,placeholder = "Paste variants here...")
                        ),
                      conditionalPanel(
                        condition="input.format=='upload'",
                        fileInput("file", NULL,accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv",
                                      ".txt",
                                      ".vcf",
                                      ".gz"))),
                    #dont need below as now have shinycssloaders
                    # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                    #                  tags$div("Loading...",id="loadmessage")),
                    actionButton("submit_button","Submit"),
                    actionButton("reset_button", "Reset form")
                      ),
                   mainPanel(
                     width = 10,
                     tabsetPanel(
                       id = "All_results",
                       type = "tabs",
                       tabPanel(value = "tab1",
                                title = h4("Paralogue Annotation"),
                                h4("Equivalent missense variant(s) identified by Paralogue Annotation"),
                                #tags$head(tags$style("#paralog  {white-space: nowrap;  }")), #set nowrap for table column names
                                conditionalPanel(condition = "input.submit_button", withSpinner(dataTableOutput("paralog"))),
                                br(),
                                conditionalPanel("output.paralog",downloadButton("download_paralog","Download (.tsv)"),downloadButton("download_paralog_excel","Download (.xslx)")),
                                #conditionalPanel("output.paralog",downloadButton("download_paralog_excel","Download (.xslx)")),
                                br(),
                                br()
                       ),
                       tabPanel(value = "tab2",
                                title = h4("All Paralogous Positions"),
                                h4("Equivalent positions identified by Paralogue Annotation"),
                                #tags$head(tags$style("#paraloc  {white-space: nowrap;  }")), #set nowrap for table column names
                                conditionalPanel(condition = "input.submit_button", withSpinner(dataTableOutput("paraloc"))),
                                br(),
                                conditionalPanel("output.paraloc",downloadButton("download_paraloc","Download (.tsv"),downloadButton("download_paraloc_excel","Download (.xlsx)")),
                                #conditionalPanel("output.paraloc",downloadButton("download_paraloc_excel","Download (.xlsx)")),
                                br(),
                                br()
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
