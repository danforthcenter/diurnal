library(shiny)
library(shinyBS)
library(ggplot2)
library(plotly)
library(jsonlite)
library(data.table)
library(DT)
library(stringr)
library(dqshiny)

# Define UI for application that draws a histogram
ui <- shiny::fluidPage(
    # Application title
  shiny::titlePanel("Diurnal Gene Expression"),
    # Sidebar with a slider input for number of bins 
    shiny::sidebarLayout(
      shiny::sidebarPanel(
          shiny::selectInput(
            "species",
            "Species:",
            choices = c(
              "arabidopsis_thaliana", "brachypodium_distachyon", "diurnal_data", 
              "glycine_max", "oryza_sativa", "populus_trichocarpa", "zea_mays"
            ),
            selected = "arabidopsis_thaliana"
          ),
          shiny::uiOutput("conditions_buttons"),
          # autocompletion relies on `dqshiny`, which is probably not maintained.
          # should anything go wrong with the autocompletion change to use the textInput
          # below. The rest should be set up to allow for multiple genes already
          # except for line 104 (making groups in the plot).
          #textInput("gene_names", "Enter Gene Name(s)", placeholder = "gene1, gene2, ..."),
          shiny::uiOutput("gene_name_autocomplete"),
          # shinyBS does not work well with fileInput,
          # so I'm attaching a tooltip to a separate header text
          tags$div(id = "staticText", strong("Optional File Input:")),
          shinyBS::bsTooltip(
            id = "staticText",
            title = "File should contain gene names on new lines or a comma separated list of gene names",
            placement = "top", trigger = "hover"
          ),
          shiny::fileInput("gene_file", "",
            placeholder = "Gene1, Gene2, ..., GeneN",
            ),
          shiny::numericInput("correlation", "Correlation Cutoff", value = 0.8, min = 0, max = 1)
        ),
        shiny::mainPanel(
          # Show a plot of the generated distribution
          plotly::plotlyOutput("mainplot"),
          # Allow Downloads
          shiny::downloadButton('download',"Download the data"),
          # show dataframe
          DT::DTOutput("datatable"),
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # set reactive variables for inputs
  species_choice <- shiny::reactive({input$species})

  # gene choice can come from the text entry or from a file upload. If both are present they should append.
  gene_choice <- shiny::reactive({
    gene_list <- list()
    if (shiny::isTruthy(input$gene_names)) {
      # note this is no longer necessary, but the rest is set up for a list and I might
      # switch back to textInput so I'll leave it for now.
      gene_list <- append(
        gene_list,
        lapply(strsplit(input$gene_names, ",")[[1]], function(g) {
          stringr::str_to_lower(trimws(g))
        })
      )
    }
    if (shiny::isTruthy(input$gene_file)) {
      genes_from_file <- # hacky way to do this to allow for several file types.
        lapply( # reads in data, forces it to be a list, flattens it to vector, rebuilds 1 level list.
          unlist( # this is just to allow for txt or csv with consistency in output
            as.list(
              data.table::fread(input$gene_file$datapath, header = FALSE)
            )
          ), identity)
      gene_list <- append(
        genes_from_file, gene_list
      )
    }
    return(gene_list)
  })
  condition_choice <- shiny::reactive({input$conditions_radio})
  cor_cutoff <- shiny::reactive({input$correlation})
  # make condition buttons based on species selection
  output$conditions_buttons <- shiny::renderUI({
    keys <- jsonlite::fromJSON("/mnt/species_to_condition.json")[[species_choice()]]
    shiny::tagList(
      shiny::radioButtons(
        "conditions_radio", "Condition: ",
        choices = keys, selected = keys[1], inline = TRUE)
    )
  })
  # should anything go wrong with autocompletion use textInput in UI side and remove this output
  output$gene_name_autocomplete <- shiny::renderUI({
    dqshiny::autocomplete_input("gene_names", "Gene Name:", dir(paste0("/mnt/", species_choice())),
      max_options = 50)
  })
  # make reactive dataframe
  df <- shiny::reactive({
    shiny::req(
      # require at least one way of identifying genes
      isTruthy(length(gene_choice()) >= 1),
      # require species/conditions
      isTruthy(input$species), isTruthy(input$conditions_radio)
    )
    genes <- gene_choice()
    full_dt <- data.table::rbindlist(lapply(genes, function(gene) {
      files = dir(path = paste0("/mnt/", species_choice(), "/", gene),
        pattern = paste0("^", condition_choice(), "__.*.csv"), full.names = TRUE)
      dt <- data.table::rbindlist(lapply(files, function(file) {
        d <- data.table::fread(file)
        d$probe <- gsub(".*__(.*).csv", "\\1", file)
        d$gene <- gene
        return(d)
      }))
      if ("correlation" %in% colnames(dt)){
        dt <- dt[is.na(correlation) | correlation > cor_cutoff(), ]
      }
      return(dt)
    }))
    return(full_dt)
  })
  # send dataframe to UI
  output$datatable <- DT::renderDataTable({
    df()
  })
  # main plot output
  output$mainplot <- plotly::renderPlotly({
    shiny::req(
      # require at least one way of identifying genes
      isTruthy(length(gene_choice()) >= 1),
      # require species/conditions
      isTruthy(input$species), isTruthy(input$conditions_radio)
    )
    df <- as.data.frame(df())
    modlist <- unique(na.omit(df$model))
    # if multiple genes were possible at once then this would need to be changed
    df[is.na(df$model), "model"] <- "Expression"
    df$model <- factor(df$model, ordered = TRUE,
      levels = c("Expression", modlist))
    line_layer <- ggplot2::geom_line(aes(linetype = probe))
    point_layer <- ggplot2::geom_point(aes(shape = probe))
    if (length(unique(df$probe)) == 1) {
      df$probe = ""
      line_layer <- ggplot2::geom_line()
      point_layer <- ggplot2::geom_point()
    }
    ggp <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = val,
      color = model, group = gene
    )) +
      line_layer +
    point_layer +
    ggplot2::labs(
      x = "Time of Day (h)", y = "Expression", linetype = "Probe", color = "Model",
      title = paste0(species_choice(), ": ", paste(gene_choice(), collapse = ", "))
    ) +
    ggplot2::scale_color_hue()+
    ggplot2::theme_bw()
    return(plotly::ggplotly(ggp))
  })

  output$download <- shiny::downloadHandler(
    filename = function(){
      paste0(species_choice(), "-",
        paste(gene_choice(), collapse = "_"), "-",
        condition_choice(), ".csv"
      )
    },
    content = function(fname){
      write.csv(df(), fname)
    }
  )
}

# Run the application 
shiny::shinyApp(ui = ui, server = server)
