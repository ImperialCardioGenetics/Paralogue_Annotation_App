suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(shinycssloaders)
  library(kableExtra)
}))



clinvar_version=20250623

fluidPage(
  useShinyjs(),

  #theme=shinytheme("yeti"), # eg. cosmo # https://rstudio.github.io/shinythemes/
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
    navbarPage(title = "PARALOG Annotator v0.4.4", id = "navbar",selected = "tab1",
               
               # Load awesome fonts stylesheet to be used in the app
               header=tags$head(tags$link(rel = "stylesheet", href = "https://maxcdn.bootstrapcdn.com/font-awesome/4.4.0/css/font-awesome.min.css")),
               
    # Main box search and description -----------------------------------------
    tabPanel(title = "Home", value = "tab1",
            br(),br(),
            h2("PARALOG Annotator", align = "center"),
            br(),br(),
            fluidRow(
              column(width = 6, offset = 3, align = "center",
                      wellPanel(style = "background-color: #333333;
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
                                HTML("to input a list or variants or upload a vcf file. All variants in Genome build GRCh37 coordinates."),
                                HTML(paste0("The ClinVar version used is ",clinvar_version,".")),
                                br(),br(),
                                textInput(inputId = "line", label = NULL),
                                HTML("e.g. <br>1-115256528-T-G<br>"),

                                br(),
                                actionButton(inputId ="search_button", label = "Search"
                                            ), 
                                style = "background-color: #ffffff;
                              border-bottom-color: #333333;
                              border-left-color: #333333;
                              border-right-color: #333333;
                              box-shadow: 3px 3px 3px #d8d8d8;
                              margin-top: 0px")
                      ) # WellPanel
            ), #Fluid row
            fluidRow(column(width = 6, offset = 3,
                            br(),
                            br(),
                            p("Paralogue Annotation utilizes information from evolutionarily related proteins, specifically paralogues, to help inform the clinical significance of missense variants associated with human diseases.",align = "center"), 
                            p(""), style = "background-color: #ffffff")
            ),
            
            # fluidRow(column(12, align="center",
            #                 br(),
            #                 # tags$p("You can download the data ",
            #                 #        tags$a("here", id = "download_file", href = "session/download/download_file", download = "paralog_data.txt.gz"))),
            #                 tags$p("You can download the data ",downloadLink("download_file", "here"))),
            #                 br()
            # ),
            fluidRow(column(12, align="center",
                            br(), 
                            br(),
                            tags$a(href='https://www.imperial.ac.uk/', target="_blank",tags$img(src='Logo_for_Imperial_College_London.svg.png',height='50',width='200' )),
                            br(),
                            br(),
                            br(),
                            tags$a(href='https://lms.mrc.ac.uk/', target="_blank",tags$img(src='UKRI_MRC_LMS_VERTICAL_RGB.svg',height='150',width='500' )))
            ),
            br(),
            br(),
            br()
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
                              choiceNames = list(
                                HTML("<div style='font-size:14px'>Paste Variants</div>"),
                                HTML("<div style='font-size:14px'>Upload Variants</div>")
                              ),
                              choiceValues = list(
                                "paste", "upload"
                              ),
                              selected = "paste",
                              width = "100%"),
                  HTML("e.g. <br>1-115256528-T-C<br>1-115256528-T-G<br>3-38592567-T-A<br>21-44592214-C-T<br>X-70443591-G-A<br>"),
                br(),
                h5("Maximum number of variants for upload: 50"),
                  conditionalPanel(
                    condition="input.format=='paste'",
                    textAreaInput("var",label=NULL,placeholder = "Paste variants here...")
                    ),
                  conditionalPanel(
                    condition="input.format=='upload'",
                    fileInput("file",label=NULL,accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv",
                                  ".txt",
                                  ".vcf",
                                  ".gz"))),
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
                            title = h4("Paralogous Annotations"),
                            br(),
                            HTML("Shown below are missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b> found at the equivalent amino acid residue of other members of the protein family by Paralogue Annotation."),
                            br(),
                            HTML(paste0("The <b>ClinVar</b> version used is ",clinvar_version,".")),
                            br(),
                            br(),
                            fluidRow(conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paralog")))),
                            br(),
                            br(),
                            br(),
                            br()
                            ),
                    tabPanel(value = "tab2",
                            title = h4("All Paralogous Positions"),
                            br(),
                            HTML("Shown below are all the equivalent amino acid positions across all members of the paralogous family"),
                            br(),
                            br(),
                            conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paraloc"))),
                            br(),
                            br(),
                            br()
                            ),
                    
                    tabPanel(value = "tab4",
                            title = h4("Homologous Pfam Annotations"),
                            br(),
                            HTML("Shown below are missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b> found at the equivalent amino acid residue of other homologous proteins that share a pfam protein domain."),
                            br(),
                            br(),
                            conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("homolog"))),
                            br(),
                            br(),
                            br()
                            )
                  )
              )
            )
    ),
      
      
    tabPanel(title = "Data", value = "tab3",
             style = "width:80%; margin-right:auto; margin-left:auto", 
             h2("Downloads"),
             br(),
             p(HTML("<b>1.</b> Download the list of all possible missense variants (855K) and their annotated gene that are predicted to be pathogenic by paralogous annotation using the ClinVar P/LP variants (version 20190114).")),
             p(HTML("<i>All data mapped to GRCh37 reference genome.</i>")),
             # br(),
             downloadButton("all_pos_mis_PA", label = HTML(sprintf("Download</button>")),class = "dt-button-clone"),
             # uiOutput("table_mis_PA"),
             br(),
             br(),
             br(),
             br(),
             p(HTML("<b>2.</b> Download all possible amino acid substitutions on all paralogous genes annotated by paralogous annotation for all possible missiense variants.")),
             p(HTML("<i>All data mapped to GRCh37 reference genome.</i>")),
             br(),
             uiOutput("table"),
             br()
    ),
      
      
    tabPanel(title = "About", value = "tab4",
               style = "width:80%; margin-right:auto; margin-left:auto", 
               withSpinner(htmlOutput("about")), # This is an html file page that is read in from dir 
               # withSpinner(uiOutput("about")), # This is an html file page that is read in from dir 
               br()
    )
  )
)
