library(shiny)
library(shinythemes)

fluidPage(theme=shinytheme("cosmo"), # eg. lumen # https://rstudio.github.io/shinythemes/
    navbarPage(
      "PARALOG Annotator",
      tabPanel("Search",
               #h2("Missense Variant Annotation for Inherited Cardiac Conditions",align="center"),
               br(),
               sidebarLayout(
                  sidebarPanel(img(src = "paralogo2.png", width = "100%"),
                      h3("Input your variant"),
                      br(),
                      radioButtons("format",label=NULL,
                                   # Here a new input method can be inserted eg. upload a file with variants
                                   # eg. choices = list("upload file"="upload",
                                   # fileInput("file", NULL,accept = c("text/csv","text/comma-separated-values,text/plain",".csv")))),
                                   choices = list("Choose position"="pick","Paste variants"="paste", "Upload Variants"="upload"),selected = NULL),
                      p("e.g. 1:114713907:T:G, 3:38551076:T:A, or X:71223741:G:A"),
                      conditionalPanel(condition="input.format=='pick'",
                      selectInput(inputId = "chr",
                                  label = "Chromosome:", 
                                  selected = "1",width = "80",multiple = F,selectize = F,
                                  choices = c(1:23)),
                      textInput(inputId = "pos",
                                label = "Position:",
                                width = "100",
                                placeholder = "",
                                c("114713907")),
                      selectInput(inputId = "ref",
                                  label = "Reference:", 
                                  selected = "T",width = "80",multiple = F,selectize = F,
                                  choices = c("A","G","T","C")),
                      selectInput(inputId = "alt",
                                  label = "Alternate:", 
                                  selected = "A",width = "80",multiple = F,selectize = F,
                                  choices = c("A","G","T","C"))),
                      conditionalPanel(condition="input.format=='paste'",textAreaInput("var",label=NULL,placeholder = "1:114713907:T:G")),
                      conditionalPanel(condition="input.format=='upload'",fileInput("file", NULL,accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")))
                      ),
                   mainPanel(h2("Missense Variant Paralogue Annotation",align="center"),
                            dataTableOutput("paralog"),
                            conditionalPanel("output.paralog",downloadButton("download","Download"))
               ))),
      tabPanel("About",
      style = "width:80%; margin-right:auto; margin-left:auto", 
      includeHTML("about.html"), # This is an HTML page that is read in from dir 
      br()
      )
      
    )
  )
