#' Main function to start EpiShiny
#'
#' @description This application was built based on the python package
#' epiScanpy, the tool for single-cell epigenomic analysis (Danese et al. 2019).
#' The application also provide the basic analysis functions such as clustering
#' and plotting from the python package scanpy (Wolf, Angerer, and Theis 2018).
#' Therefore, this function supports the input for both scRNA-seq and scATAC-seq.
#' Note that The function requires no parameters and always returns NULL.
#'
#'
#' @import shiny
#' @import reticulate
#' @import plotly
#' @import stringr
#'
#' @export
#' @md
#'
#' @references
#'
#' Danese, Anna, Maria L. Richter, David S. Fischer, Fabian J. Theis, and
#' Maria Colomé-Tatché. 2019. “EpiScanpy: Integrated Single-Cell Epigenomic
#' Analysis.” Preprint. Bioinformatics. https://doi.org/10.1101/648097.
#'
#' Wolf, F. Alexander, Philipp Angerer, and Fabian J. Theis. 2018. “SCANPY:
#' Large-Scale Single-Cell Gene Expression Data Analysis.” Genome Biology 19
#' (1): 15. https://doi.org/10.1186/s13059-017-1382-0.


#'
#'
#' @examples
#' \dontrun{
#' library(epiview)
#' epiview()
#' }
#'

# read this https://mastering-shiny.org/scaling-packaging.html
# https://docs.rstudio.com/tutorials/user/using-python-with-rstudio-and-reticulate/
# https://rstudio.github.io/shinydashboard/
#
# create new virtual environment
# virtualenv .venv
# activate virtual environment
# source .venv/bin/activate


epishiny <- function(...){
    # Define UI for application that draws a histogram
    python_file_path <- file.path(system.file(package="epishiny"), "interface_anndata.py")
    source_python(python_file_path)

    adata <- NULL
    filename1 <- NULL
    selected_data <- NULL
    options(shiny.maxRequestSize=5000*1024^2)

    ui <- bootstrapPage(
        tabsetPanel(
            tabPanel("Workspace",
                sidebarLayout(
                     sidebarPanel(
                         fileInput("file1", "Choose h5ad file", accept = c(".h5ad",".csv",".rds")),
                         fileInput("workspacefile", "Load workspace", accept = c(".rda")),
                         downloadButton("saveWorkspace", "Save workspace"),
                     ),
                     mainPanel(
                         plotlyOutput("mainWorkspace")
                     )
                )
            ),
            tabPanel("Interactive plot",
                 sidebarLayout(
                     sidebarPanel(
                         uiOutput('obsm_menu'),
                         uiOutput('obs_menu'),
                         uiOutput('recal_menu'),
                         downloadButton("saveAllData", "Save All Data"),
                         downloadButton("saveSelectedData", "Save Selected Data")
                     ),
                     mainPanel(
                         plotlyOutput("mainInterative")
                     )
                 )
            ),
            tabPanel("Static plot",
                 sidebarLayout(
                     sidebarPanel(
                         downloadButton("saveFigure", "Save figure")
                     ),
                     mainPanel(
                         plotlyOutput("mainStatic")
                     )
                 )
            ),
            tabPanel("Pathway",
                 sidebarLayout(
                     sidebarPanel(
                         downloadButton("exportPathway", "Export")
                     ),
                     mainPanel(
                         plotlyOutput("mainPathway")
                     )
                 )
            )
        )
    )

    # Define server logic required to draw a histogram
    server <- function(input, output) {

        all_obsm = c('-','PCA 2D','PCA 3D','UMAP 2D','UMAP 3D','TSNE 2D','DRAW_GRAPH_FA 2D','DIFFMAP 2D','DIFFMAP 3D')

        output$obsm_menu = renderUI({
            if (is.null(adata)){
                selectInput('obsm_select', 'Embedded Objects', choices = c('-'))
            }else{
                menulist = adata$obsm_keys()
                strlist = c()
                for (m in menulist){
                    for (d in 2:dim(adata$obsm[m])[2]){
                        if (d>3){
                            break
                        }
                        strlist = c(strlist,toupper(paste0(m,' ',d,'D')))
                    }
                }
                strlist = str_replace(strlist,'X_','')
                selectInput('obsm_select', 'Embedded Objects', strlist)
            }
        })

        output$obs_menu = renderUI({
            if (is.null(adata)){
                selectInput('obs_select', 'Observations', choices = c('-'))
            }else{
                menulist = adata$obs_keys()
                selectInput('obs_select', 'Observations', menulist)
            }
        })

        output$recal_menu = renderUI({
            if (is.null(adata)){
                selectInput('obsm_recal_select', 'Process', choices = c('-'))
            }else{
                menulist = adata$obsm_keys()
                strlist = c()
                for (m in menulist){
                    for (d in 2:dim(adata$obsm[m])[2]){
                        if (d>3){
                            break
                        }
                        strlist = c(strlist,toupper(paste0(m,' ',d,'D')))
                    }
                }
                strlist = str_replace(strlist,'X_','')
                strlist = setdiff(all_obsm,strlist)
                selectInput('obsm_recal_select', 'To calculate', strlist)
            }
        })

        observeEvent(input$obsm_recal_select, {
            select_recal = str_split(input$obsm_recal_select,' ',simplify = TRUE)[1]
            if (select_recal != '-'){
                cal_obsm(adata,select_recal)

                menulist = adata$obsm_keys()
                strlist = c()
                for (m in menulist){
                    for (d in 2:dim(adata$obsm[m])[2]){
                        if (d>3){
                            break
                        }
                        strlist = c(strlist,toupper(paste0(m,' ',d,'D')))
                    }
                }
                strlist = str_replace(strlist,'X_','')
                diffstrlist = setdiff(all_obsm,strlist)

                updateSelectInput(inputId = "obsm_select", choices = strlist, selected = input$obsm_recal_select)
                updateSelectInput(inputId = "obsm_recal_select", choices = diffstrlist, selected = diffstrlist[1])

            }
        })

        observeEvent(input$file1, {
            filename1 <<- input$file1$datapath
            ext <- tools::file_ext(filename1)

            req(filename1)
            validate(need(ext == "h5ad", "Please upload a h5ad file"))

            adata <<- get_anndata(filename1)

            menulist_obs = adata$obs_keys()
            menulist = adata$obsm_keys()
            strlist = c()
            for (m in menulist){
                for (d in 2:dim(adata$obsm[m])[2]){
                    if (d>3){
                        break
                    }
                    strlist = c(strlist,toupper(paste0(m,' ',d,'D')))
                }
            }
            strlist = str_replace(strlist,'X_','')
            diffstrlist = setdiff(all_obsm,strlist)

            updateSelectInput(inputId = "obsm_select", choices = strlist, selected = strlist[1])
            updateSelectInput(inputId = 'obs_select', choices = menulist_obs, selected = menulist_obs[1])
            updateSelectInput(inputId = "obsm_recal_select", choices = diffstrlist, selected = diffstrlist[1])
        })

        output$mainInterative <- renderPlotly({

            if (input$obsm_select != '-' & input$obsm_select %in% all_obsm){
                select_obsm = str_split(tolower(input$obsm_select),' ',simplify = TRUE)
                obsm_idx = paste0('X_',select_obsm[1])
                dimension = select_obsm[2]

                if (dimension == "2d"){
                    p <- plot_ly(x = adata$obsm[obsm_idx][,1], y = adata$obsm[obsm_idx][,2], key = adata$obs_names$to_list(), color = adata$obs[,input$obs_select]) %>%
                        toWebGL()
                }else{
                    p <- plot_ly(x = adata$obsm[obsm_idx][,1], y = adata$obsm[obsm_idx][,2], z = adata$obsm[obsm_idx][,3], key = adata$obs_names$to_list(), color = adata$obs[,input$obs_select]) %>%
                        add_markers(customdata = row.names(adata$obsm[obsm_idx])) %>%
                        toWebGL()
                }
            }
        })

        output$saveAllData <-  downloadHandler(
            filename = function() {
                ext <- tools::file_ext(filename1)
                paste0("download_all_data.", ext)
            },
            content = function(file) {
                save_all_h5ad(adata,file)
            }
        )

        output$saveSelectedData <-  downloadHandler(
            filename = function() {
                ext <- tools::file_ext(filename1)
                paste0("download_selected_data.", ext)
            },
            content = function(file) {
                selected_data <<- event_data("plotly_selected")
                print(selected_data)
                save_sub_h5ad(adata, file, selected_data$key)
            }
        )
    }

    # Run the application
    shinyApp(ui = ui, server = server)

}
