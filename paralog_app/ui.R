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
    navbarPage(title = "PARALOG Annotator DEMO version 0.2.5", id = "navbar",selected = "tab1",
               
               # Main box search and description -----------------------------------------
               tabPanel(title = "Home", value = "tab1",
                        br(),br(),
                        h2("PARALOG Annotator", align = "center"),
                        br(),br(),
                        fluidRow(
                          column(width = 6, offset = 3, align = "center",
                                 wellPanel(
                                   #"Type a variant ID",
                                           style = "background-color: #333333;
                                         color: white;
                                         border-top-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-bottom: 0px;
                                         padding:5px"), 
                                
                                 wellPanel(br(),
                                           HTML( "Input a query variant bellow or click"),
                                           actionLink("link_to_tabpanel_b", "here"),
                                           HTML("to input a list or variants or upload a vcf file. All variants in Genome build GRCh37 coordinates"),
                                           br(),br(),
                                           textInput(inputId = "line", label = NULL),#, value = "clinvar"),
                                           HTML("e.g. <br>1-115256528-T-G<br>"),

                                           br(),
                                           actionButton(inputId ="search_button", label = "Search"
                                                        #, class = "btn-primary"
                                                        ), 
                                           style = "background-color: #ffffff;
                                         border-bottom-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-top: 0px")
                          ) # WellPanel
                        ), #Fluid row
                        fluidRow(column(width = 6, offset = 3, br(), br(), p("Paralogue Annotation utilizes information from evolutionarily related proteins, specifically paralogues, to help inform the clinical significance of missense variants associated with human diseases.",
                                                                             align = "center"), p(""), style = "background-color: #ffffff"),
                                 ),
                        fluidRow(column(12, align="center",br(), br(),
                                 #img(src='./data/Logo_for_Imperial_College_London.svg.png', align = "center"),
                                 #titlePanel(title=div(img(src="nhl.jpg")))
                                 
                                 tags$a(href='https://www.imperial.ac.uk/', target="_blank",tags$img(src='Logo_for_Imperial_College_London.svg.png',height='50',width='200' ))
                                 ))
                                 
               ),
               
               
               # Results - search "left" side -------------------------------------------------------------
               tabPanel(title = "Results", value = "tab2",
               br(),
               sidebarLayout(
                  sidebarPanel(
                    id = "tab2_search",
                    width = 2,

                    h3("Input your variants"),
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
                      HTML("e.g. <br>1-115256528-T-C<br>1-115256528-T-G<br>3-38592567-T-A<br>21-44592214-C-T<br>X-70443591-G-A<br>"),

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
                  
                  # Results - tables "right" side -------------------------------------------------------------
                  mainPanel(
                     width = 10,
                     tabsetPanel(
                       id = "All_results",
                       type = "tabs",
                       tabPanel(value = "tab1",
                                title = h4("Paralogue Annotation"),
                                h4("Equivalent missense variant(s) identified by Paralogue Annotation"),
                                #tags$head(tags$style("#paralog  {white-space: nowrap;  }")), #set nowrap for table column names
                                conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paralog"))),
                                                 #condition = "input.search_button", withSpinner(dataTableOutput("paralog"))),
                                br(),
                                #conditionalPanel("output.paralog",downloadButton("download_paralog","Download (.txt)"),downloadButton("download_paralog_excel","Download (.xslx)")),
                                #conditionalPanel("output.paralog",downloadButton("download_paralog_excel","Download (.xslx)")),
                                br(),
                                br()
                       ),
                       tabPanel(value = "tab2",
                                title = h4("All Paralogous Positions"),
                                h4("Equivalent positions identified by Paralogue Annotation"),
                                #tags$head(tags$style("#paraloc  {white-space: nowrap;  }")), #set nowrap for table column names
                                conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paraloc"))),
                                br(),
                               # conditionalPanel("output.paraloc",downloadButton("download_paraloc","Download (.txt"),downloadButton("download_paraloc_excel","Download (.xlsx)")),
                                #conditionalPanel("output.paraloc",downloadButton("download_paraloc_excel","Download (.xlsx)")),
                                br(),
                                br()
                                )
                       )
                     )
                  )),
      tabPanel(title = "About", value = "tab3",
      style = "width:80%; margin-right:auto; margin-left:auto", 
      #includeHTML("about.html"), # This is an HTML file that is read in from dir 
      includeMarkdown("about.md"), # This is an md file  page that is read in from dir 
      br()
      )
      
    )
  )
