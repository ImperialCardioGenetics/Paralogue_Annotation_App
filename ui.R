suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinyjs)
  library(shinycssloaders)
}))



fluidPage(
  useShinyjs(),
  
  #theme=shinytheme("yeti"), # eg. cosmo # https://rstudio.github.io/shinythemes/
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  navbarPage(
    title = "PARALOG Annotator v0.4.2",
    id = "navbar",
    selected = "tab1",
    
    
    # Main box search and description -----------------------------------------
    tabPanel(
      title = "Home",
      value = "tab1",
      br(),
      br(),
      h2("PARALOG Annotator", align = "center"),
      br(),
      br(),
      fluidRow(
        column(
          width = 6,
          offset = 3,
          align = "center",
          wellPanel(
            style = "background-color: #333333;
                                color: white;
                                border-top-color: #333333;
                                border-left-color: #333333;
                                border-right-color: #333333;
                                box-shadow: 3px 3px 3px #d8d8d8;
                                margin-bottom: 0px;
                                padding:5px"
          ),
          wellPanel(
            br(),
            HTML("Input a query variant bellow or click"),
            actionLink("link_to_tabpanel_b", "here"),
            HTML(
              "to input a list or variants or upload a vcf file. All variants in Genome build GRCh37 coordinates."
            ),
            br(),
            br(),
            textInput(inputId = "line", label = NULL),
            HTML("e.g. <br>1-115256528-T-G<br>"),
            
            br(),
            actionButton(inputId = "search_button", label = "Search"),
            style = "background-color: #ffffff;
                              border-bottom-color: #333333;
                              border-left-color: #333333;
                              border-right-color: #333333;
                              box-shadow: 3px 3px 3px #d8d8d8;
                              margin-top: 0px"
          )
        ) # WellPanel
      ),
      #Fluid row
      fluidRow(
        column(
          width = 6,
          offset = 3,
          br(),
          br(),
          p(
            "Paralogue Annotation utilizes information from evolutionarily related proteins, specifically paralogues, to help inform the clinical significance of missense variants associated with human diseases.",
            align = "center"
          ),
          p(""),
          style = "background-color: #ffffff"
        )
      ),
      fluidRow(
        column(
          12,
          align = "center",
          br(),
          br(),
          tags$a(
            href = 'https://www.imperial.ac.uk/',
            target = "_blank",
            tags$img(
              src = 'Logo_for_Imperial_College_London.svg.png',
              height = '50',
              width = '200'
            )
          ),
          br(),
          br(),
          br(),
          tags$a(
            href = 'https://lms.mrc.ac.uk/',
            target = "_blank",
            tags$img(
              src = 'UKRI_MRC_LMS_VERTICAL_RGB.svg',
              height = '150',
              width = '500'
            )
          )
        )
      ),
      br(),
      br(),
      br()
    ),
    
    
    # Results - search "left" side -------------------------------------------------------------
    tabPanel(
      title = "Results",
      value = "tab2",
      br(),
      sidebarLayout(
        sidebarPanel(
          id = "tab2_search",
          width = 2,
          h3("Input your variants"),
          h5("Genome build GRCh37"),
          br(),
          radioButtons(
            "format",
            label = NULL,
            # Here a new input method can be inserted eg. upload a file with variants
            # eg. choices = list("upload file"="upload",
            choiceNames = list(
              HTML("<div style='font-size:14px'>Paste Variants</div>"),
              HTML("<div style='font-size:14px'>Upload Variants</div>")
            ),
            choiceValues = list("paste", "upload"),
            selected = "paste",
            width = "100%"
          ),
          HTML(
            "e.g. <br>1-115256528-T-C<br>1-115256528-T-G<br>3-38592567-T-A<br>21-44592214-C-T<br>X-70443591-G-A<br>"
          ),
          br(),
          h5("Maximum number of variants for upload: 50"),
          conditionalPanel(
            condition = "input.format=='paste'",
            textAreaInput("var", label = NULL, placeholder = "Paste variants here...")
          ),
          conditionalPanel(condition = "input.format=='upload'",
                           fileInput(
                             "file",
                             label = NULL,
                             accept = c(
                               "text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv",
                               ".txt",
                               ".vcf",
                               ".gz"
                             )
                           )),
          actionButton("submit_button", "Submit"),
          actionButton("reset_button", "Reset form")
        ),
        
        # Results - tables "right" side -------------------------------------------------------------
        mainPanel(
          width = 10,
          tabsetPanel(
            id = "All_results",
            type = "tabs",
            tabPanel(
              value = "tab1",
              title = h4("Paralogous Annotations"),
              br(),
              HTML(
                "Shown below are missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b> found at the equivalent amino acid residue of other members of the protein family by Paralogue Annotation"
              ),
              br(),
              br(),
              fluidRow(
                conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paralog")))
              ),
              br(),
              br(),
              br(),
              br()
            ),
            tabPanel(
              value = "tab2",
              title = h4("All Paralogous Positions"),
              br(),
              HTML(
                "Shown below are all the equivalent amino acid positions across all members of the paralogous family"
              ),
              br(),
              br(),
              conditionalPanel(condition = "input.submit_button || input.search_button", withSpinner(dataTableOutput("paraloc"))),
              br(),
              br(),
              br()
            ),
            
            tabPanel(
              value = "tab4",
              title = h4("Homologous Pfam Annotations"),
              br(),
              HTML(
                "Shown below are missense variants annotated as <b>Pathogenic</b> or <b>Likely Pathogenic</b> in <b>ClinVar</b> found at the equivalent amino acid residue of other homologous proteins that share a pfam protein domain"
              ),
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
    
    
    tabPanel(
      title = "Data",
      value = "tab3",
      style = "width:80%; margin-right:auto; margin-left:auto",
      h2("Downloads"),
      p(
        "Download all possible amino acid substitutions on all paralogous genes annotated by PARALOG Annotator for all possible missiense variants"
      ),
      p("All data mapped to GRCh37 reference genome"),
      br(),
      dataTableOutput("table"),
      br()
    ),
    
    
    tabPanel(
      title = "About",
      value = "tab4",
      style = "width:80%; margin-right:auto; margin-left:auto",
      # withSpinner(htmlOutput("about")), # This is an html file page that is read in from dir
      HTML("
      <p>Paralogue Annotation utilizes information from evolutionarily related proteins, specifically paralogues, to help inform the clinical significance of missense variants associated with human diseases. 
      The original methodology and implementation of <a href= 'https://www.cardiodb.org/paralogue_annotation/'>Paralogue Annotation on arrhythmia syndrome genes'</a> 
      was published <a href= 'https://onlinelibrary.wiley.com/doi/full/10.1002/humu.22114'>here</a> and <a href= 'https://jmg.bmj.com/content/51/1/35'>here</a>. 
      This web app extends Paralogue Annotation exome-wide, using paralogues defined by <a href= 'https://www.ensembl.org/Help/View?id=137'>Ensembl gene trees</a> and pathogenic/likely pathogenic missense variants defined by <a href= 'https://www.ncbi.nlm.nih.gov/clinvar/'>ClinVar</a>.</p>
      <br>
      <p>This web app is currently being built using <a href= 'http://shiny.rstudio.com'>Shiny</a>, the source code is available at <a href= 'https://github.com/ImperialCardioGenetics/Paralogue_Annotation_App'>https://github.com/ImperialCardioGenetics/Paralogue_Annotation_App</a>.
      <br>
      <br>
      <br>
      <h2>Frequently Asked Questions (FAQ)</h2>
      <br>
      <b>Q. What genome build coordinates do my variants need to be in?</b>
      <b>A.</b> Currently only GRCh37 coordinates are supported. 
      We recommend using <a href= 'https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter'>Ensembl liftover service</a>, for coordinate conversions.
      <br>
      <br>
      <b>Q. What are the Para_z scores?</b>
      <br>
      <b>A.</b> The Para_z scores are a measure of paralogue conservation independently derived by <a href= 'https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00725-6'>Lal et al. (2020)</a>. 
      You may therefore find in your results that some Para_z scores do not agree with your expectations. 
      This is because the paralogue alignments used to generate the scores are different to the alignments used here. 
      The Para_z scores can thus be thought as a third-party confidence score of paralogue conservation across aligned positions.
      <br>
      <br>
      <b>Q. How can we use the conservation of Ref/Alt alleles to filter out results?</b>
      <br>
      <b>A.</b>That is not currently available in this version of the web app.
      <br>
      <br>
      <b>Q. What paralogue alignments do you use here?</b>
      <br>
      <b>A.</b> We utilize paralogue alignments at the protein level generated by Ensembl, which were obtained through <a href= 'http://www.ensembl.org/info/docs/api/compara/compara_tutorial.html'>Compara</a>.
      <br>
      <br>
      <b>Q. Why do the results for arrhythmia genes from the original <a href= 'https://www.cardiodb.org/paralogue_annotation/'>Paralogue Annotation</a> and here differ?</b>
      <br>
      <b>A.</b> This is mainly because the original Paralogue Annotation utilized T-COFFEE for the alignments, whereas Ensembl alignments are generated by CLUSTAL W instead. 
      Furthermore, variants from <a href= 'http://www.hgmd.cf.ac.uk/ac/index.php'>HGMD</a> were used instead of <a href= 'https://www.ncbi.nlm.nih.gov/clinvar/'>ClinVar</a>.
      <br>
      <br>
      <b>Q. What formats do my input variants have to be in?</b>
      <br>
      <b>A.</b> Currently variants have to be submitted using their chromosome, position, reference allele, 
      and alternate allele using any delimiter in the format of 'CHROM:POS:REF:ALT' with separate variants on newlines. 
      Alternatively we also accept <a href= 'https://www.internationalgenome.org/wiki/Analysis/vcf4.0/'>VCF</a> as well.
      <br>
      <br>
      <b>Q: I have previously identified a known paralogous variant to my variant of interest, but why does Paralogue Annotation not find it?</b>
      <br>
      <b>A:</b> We currently utilize Pathogenic and Likely pathogenic missense variants from ClinVar. 
      Your variant may not exist in ClinVar or may be considered Variant of Uncertain Significance (VUS). We do not at this time look at other databases, e.g. HGMD.
      <br>
      <br>
      <b>Q: Do you have all homologous pfam position available?</b>
      <br>
      <b>A:</b> We do not currently, but we are working towards providing this information in a near future update.
      <br>
      <br>
      <br>
      For more details on specific methods, code of how Paralogue Annotation functions, or any other questions please email <a href='mailto:nyl112@ic.ac.uk' class='email'>nyl112@ic.ac.uk</a>
      <br>
      <br>
      <br>
      This web app is a <b>work in progress,<b> final version may differ.
    ")
    )
  )
)
