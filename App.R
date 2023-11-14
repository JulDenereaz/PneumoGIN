library(shiny) 
library(shinydashboard)
library(dplyr)
library(visNetwork)
library(igraph)
library(DT)
edges <- read.table("edges.txt", header=T, sep="\t")
nodes <- read.table("D39V_targets.txt", header=T, sep=";") %>%
  filter(id %in% edges$from | id %in% edges$to) %>%
  mutate(label=Gene.Name)
edges <- edges %>%
  mutate(color.color = case_when(
    Group == "Negative" ~ '#edae49',
    Group == "Positive" ~ '#00798c',
    TRUE ~ '#d1495b')
  ) %>%
  mutate(color.highlight = color.color) %>%
  mutate(label=paste0("", round(epsilonSum, 2))) %>%
  mutate(length = round((max(abs(epsilonSum))-abs(epsilonSum) + 1) * 25))

cogs <- read.table("cogs.txt", header=T)

url1 <- "https://veeninglab.com/pbd39v/?loc=D39V.0%3A"
url2 <- "&tracks=DNA%2Cgff_operons%2Cgff_genes%2Cgff_TSS_term%2Cgff_regulatory%2CsgRNAlib0mm_XUE&highlight="

##### UI #####
ui <- dashboardPage(
  dashboardHeader(title = "PneumoGIN"),
  dashboardSidebar(width="300px",
    selectizeInput('ngenes', 'Select Nodes', choices = NULL, multiple = T, options=list(placeholder="pbp2b...")),
    actionButton("reset", "Reset nodes selection", class='btn-warning'),
    sliderInput("epsilon", "Epsilon Threshold:", min=1, max = ceiling(max(edges$epsilonSum)), value=1, step=0.5),
    actionButton("update", "Update Network", icon = icon("arrows-rotate"),class = "btn-xl")
  ),
  dashboardBody(
    fluidRow(
      tabBox(
        width="100%",
        tabPanel(
          title="Network",
          span(textOutput("textUIOut"), style="color:red;font-size:25px"),
          shinycssloaders::withSpinner(
            visNetworkOutput("networkplot", height = "85vh", width="100%"),
            type=1,
            color.background ="white"
          )
        ),
        tabPanel(
          title="Data Table",
          DT::dataTableOutput("mytable")
        )
        
      )
    )
  )
)




##### Server function #####
server <- function(input, output, session) {
  allch <- c(unlist(strsplit(as.character(nodes$Gene.Name), ",")), nodes$id)
  updateSelectizeInput(session, 'ngenes',  choices = allch, server=T)
  observeEvent(input$reset, {
    updateSelectizeInput(session, 'ngenes',selected ="")
  })
  output$networkplot <- NULL
  output$textUIOut <- renderText("Please select at least one node")
  
  output$mytable = DT::renderDataTable({
    DT::datatable(edges[-which(names(edges) %in% c("length", "label", "color.highlight", "color.color"))], options = list(pageLength = 50))
  })
  
  observeEvent(input$update, {
    if(is.null(input$ngenes)) {
      output$textUIOut <- renderText("Please select at least one Node")
      output$networkplot <- NULL
      return()
    }
    output$networkplot <- renderVisNetwork({
      output$textUIOut <- NULL
      isolate(
        ed <- edges %>% 
          filter(
            from %in% nodes$id[grepl(paste0(input$ngenes, collapse="|"), nodes$Gene.Name)] |
            to %in% nodes$id[grepl(paste0(input$ngenes, collapse="|"), nodes$Gene.Name)] |
            from %in% input$ngenes |
            to %in% input$ngenes
          ) %>%
          filter(
            case_when(
              epsilonSum < 0 ~ epsilonSum < -input$epsilon,
              epsilonSum > 0 ~ epsilonSum > input$epsilon  
            )
            
          )
      )
      nd <- nodes %>% 
        filter(id %in% ed$from | id %in% ed$to)
      nd$title <- apply(nd, 1, function(c) {
        return(paste0(
          "<p><h4><u><b>", 
          c[[1]],
          "</b></u></h4>",
          "<b>D39V target(s):</b><br>",
          "<em>-",
          c[[3]],
          "</em>",
          "</b></u><br>",
          "<b>COGs:</b><br>",
          paste0("-", cogs$Desc[match(unique(unlist(strsplit(as.character(c[[4]]), ","))), cogs$COG)], collapse="<br>"),
          "<br>",
          "<b>PneumoBrowse link(s):</b><br>",
          paste0(
            "-",
            '<a target="_blank" href = "',
            url1,
            unlist(strsplit(as.character(c[[5]]), ",")),
            url2,
            '">',
            unlist(strsplit(as.character(c[[5]]), ",")),
            '</a>',
            collapse = "<br>"
          ),
          "</p>"))
      })
      visNetwork(nd, ed) %>%
        visOptions(highlightNearest = T, autoResize = T) %>%
        visEdges(width=5, length=length, shadow = list(enabled = TRUE)) %>%
        visPhysics(stabilization = T) %>%
        visIgraphLayout(randomSeed=123) %>%
        visGroups(groupname="S", color="grey") %>%
        visNodes(shadow=T, size = 22, font=list(size=30)) %>%
        addFontAwesome()

    })

  })
  
}



shinyApp(ui=ui, server=server)


