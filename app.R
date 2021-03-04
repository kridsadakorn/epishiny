#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# https://docs.rstudio.com/tutorials/user/using-python-with-rstudio-and-reticulate/
# https://rstudio.github.io/shinydashboard/
#
# create new virtual environment
# virtualenv .venv
# activate virtual environment
# source .venv/bin/activate


library(shiny)
library(reticulate)
library(stringr)
library(plotly)

source_python("interface_anndata.py")
options(shiny.maxRequestSize=100*1024^2)

# Define UI for application that draws a histogram
ui <- bootstrapPage(

    # Application title
    titlePanel("Display options"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose h5ad file", accept = ".h5ad"),
            uiOutput('obsm_menu'),
            uiOutput('obs_menu'),
            uiOutput('recal_menu')
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("mainPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    adata = get_anndata("/Users/kris/Desktop/Test/processed_windows_CG_luo_et_al_nov2020_paper_resubmission.h5ad")
    #adata = NULL
    
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
            selectInput('obsm_recal_select', 'To calculate', choices = c('-'))
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
            all_obsm = c('-','PCA 2D','PCA 3D','UMAP 2D','UMAP 3D','TSNE 2D','DRAW_GRAPH_FA 2D')
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
            all_obsm = c('-','PCA 2D','PCA 3D','UMAP 2D','UMAP 3D','TSNE 2D','DRAW_GRAPH_FA 2D')
            diffstrlist = setdiff(all_obsm,strlist)
            
            updateSelectInput(inputId = "obsm_select", choices = strlist, selected = input$obsm_recal_select)
            updateSelectInput(inputId = "obsm_recal_select", choices = diffstrlist, selected = diffstrlist[1])
            
        }
    })
    
    observeEvent(input$file1, {
        ext <- tools::file_ext(input$file1$datapath)
        
        req(input$file1$datapath)
        validate(need(ext == "h5ad", "Please upload a h5ad file"))
        
        adata <- get_anndata(input$file1$datapath)
        
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
        all_obsm = c('-','PCA 2D','PCA 3D','UMAP 2D','UMAP 3D','TSNE 2D','DRAW_GRAPH_FA 2D')
        diffstrlist = setdiff(all_obsm,strlist)
        
        updateSelectInput(inputId = "obsm_select", choices = strlist, selected = strlist[1])
        updateSelectInput(inputId = 'obs_select', choices = menulist_obs, selected = menulist_obs[1])
        updateSelectInput(inputId = "obsm_recal_select", choices = diffstrlist, selected = diffstrlist[1])
    })
    
    output$mainPlot <- renderPlotly({
        if (input$obsm_select != '-'){
            select_obsm = str_split(tolower(input$obsm_select),' ',simplify = TRUE)
            obsm_idx = paste0('X_',select_obsm[1])
            dimension = select_obsm[2]
            
            if (dimension == "2d"){
                plot_ly(x = adata$obsm[obsm_idx][,1], y = adata$obsm[obsm_idx][,2], color = adata$obs[,input$obs_select])
            }else{
                plot_ly(x = adata$obsm[obsm_idx][,1], y = adata$obsm[obsm_idx][,2], z = adata$obsm[obsm_idx][,3], color = adata$obs[,input$obs_select])
            }
        }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
