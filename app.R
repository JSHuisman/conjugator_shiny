###########################################################
# Shiny App for Conjugation Rate Estimation
# Author: Jana Huisman
###########################################################

library(shiny)
library(deSolve)
library(DT) # >= 0.6 needed

library(tidyverse)
library(reshape2)
library(ggrepel)
library(gridExtra) # to arrange multiple ggplots
library(cowplot)

#library(rmarkdown)
library(shinythemes)
library(shinyBS)
library(shinycssloaders)


library(conjugator)
paper_git_dir = "../../conjugator_paper/"

source(paste0(paper_git_dir, "R/model_functions.R"))
source(paste0(paper_git_dir, "R/plotting_functions.R"))
source(paste0(paper_git_dir, "R/utils.R"))
source(paste0(paper_git_dir, "R/shiny_utils.R"))

###########################################################
example_data <- read.table("data/example.csv",
                           header = TRUE,
                           sep = ",", stringsAsFactors = FALSE,
                           colClasses=c("character", rep("double",9)))

colnames_meas_var <- c("psi.D", "psi.R", "psi.T", "D.0",
                       "R.0", "T.0", "D.t", "R.t", "T.t", "t")

complete_dataframe <- data.frame()
complete_dataframe["ID"] <- as.character()
for (var in colnames_meas_var) complete_dataframe[[var]] <- as.numeric()

##################################
calculate_C0 <- function(input, tstat = NULL, Dfinal = NULL){
  
  if(!is.null(Dfinal)){
    tstat = log(Dfinal/(10**as.numeric(input$logD0)))/input$psiD
  }
  
  C0 = 1*(10**as.numeric(input$logD0) + 10**as.numeric(input$logR0) + 
            ifelse(10**as.numeric(input$nonzero_logT0), input$logT0, 0))*
    (exp(max(input$psiD, input$psiR, input$psiT, input$psi)*tstat) - 1)
  
  return(C0)
}

input_user_data <- function(file_path){
  df <- try(read.csv(file_path, 
                     header = TRUE, sep = ",",
                     stringsAsFactors = FALSE), silent = TRUE)
  # some error handling
  if (class(df) == "try-error"){
    showNotification(paste0("This input file can not be read as a ",
                            "csv file. Please try with a different file."), 
                     action = NULL, duration = 5, closeButton = TRUE,
                     id = "csv_input_error", type = "error",
                     session = getDefaultReactiveDomain())
  } 
  req(class(df) != "try-error")
  df <- df %>% 
    mutate_if(is.integer, as.double) %>%
    mutate_at("ID", as.character)
  
  return(df)
}

###########################################################
# expSlider javascript function
JS.expify <-
  "
// function to exponentiate a sliderInput
function expSlider (sliderId, sci = false) {
  $('#'+sliderId).data('ionRangeSlider').update({
  'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
  })
}"

# call expSlider for each relevant sliderInput
JS.onload <-
  "
// execute upon document loading
$(document).ready(function() {
// wait a few ms to allow other scripts to execute
setTimeout(function() {
// include call for each slider
expSlider('log_gammaD', sci = true)
expSlider('log_gammaT', sci = true)
expSlider('log_gamma', sci = true)
expSlider('logT0', sci = true)
expSlider('logR0', sci = true)
expSlider('logD0', sci = true)
expSlider('logC0', sci = true)
}, 5);
// make sure layout is updated when slider is
$('#C0_widget_setval').click(function() {
setTimeout(function() {
// include call for each slider
expSlider('logC0', sci = true)
}, 100);
})
})
"

###########################################################
options(shiny.sanitize.errors = TRUE)

ui <- fluidPage(theme = shinytheme("flatly"),
                withMathJax(),
                tags$head(tags$script(HTML(JS.expify))),
                tags$head(tags$script(HTML(JS.onload))),
                tags$head(tags$style(".rightAlign{float:right;}")),
                navbarPage("Conjugation Rate Estimator",
                           ##################################
                           # Application home ----
                           tabPanel("Home",
                                    h1("Conjugation Rate Estimator"),
                                    br(),
                                    withMathJax(includeMarkdown("app_home.md")),
                                    br(),
                                    br(),
                                    fluidRow(
                                      column(width = 3, plotOutput("ethlogo", inline = FALSE)),
                                      column(width = 3, plotOutput("snflogo", inline = FALSE))
                                    )
                                    
                                    
                           ),
                           ##################################
                           # Data analysis ----
                           tabPanel("Analyse experimental data",
                                    #h1("Analyse experimental data"),
                                    withMathJax(includeMarkdown("analyse_data_intro.md")),
                                    # Input experimental data ----
                                    tabsetPanel(type = "tabs",
                                                tabPanel("DRT",
                                                         sidebarLayout(
                                                           sidebarPanel(
                                                             div(id = "info_DRT_data", class = 'rightAlign', icon(name = "info-circle")),
                                                             h3("Upload Data"),
                                                             # Input: Select a file 
                                                             fileInput("userdataDRT", "Choose CSV File",
                                                                       multiple = FALSE,
                                                                       accept = c("text/csv",
                                                                                  "text/comma-separated-values,text/plain",
                                                                                  ".csv"),
                                                                       placeholder = "No file selected"),
                                                             selectInput("id_cols", "Identifying columns (to be included in output)",
                                                                         multiple = TRUE,
                                                                         choices = colnames(example_data), 
                                                                         selected = "ID"),
                                                             p('If previously uploaded data is not removed, 
                                                                the app will attempt a full join of the old and the new data.'),
                                                             bsButton("trash_dataDRT", label = "Remove DRT data", icon = icon("trash")),
                                                             h3("Example DRT Data"),
                                                             checkboxInput("exampledata", "Show Example", value = TRUE)
                                                           ),
                                                           mainPanel(
                                                             DTOutput("editdataDRT")
                                                           )) # Sidebarlayout 
                                                ), #tabpanel
                                                tabPanel("TRT",
                                                         sidebarLayout(
                                                           sidebarPanel(
                                                             h3("Upload Data"),
                                                             # Input: Select a file 
                                                             fileInput("userdataTRT", "Choose CSV File",
                                                                       multiple = FALSE,
                                                                       accept = c("text/csv",
                                                                                  "text/comma-separated-values,text/plain",
                                                                                  ".csv"),
                                                                       placeholder = "No file selected"),
                                                             bsButton("trash_dataTRT", label = "Remove TRT data", icon = icon("trash"))
                                                           ),
                                                           mainPanel(
                                                             DTOutput("editdataTRT")
                                                           )) # Sidebarlayout 
                                                ) #tabPanel
                                    ), #tabsetPanel
                                    # Displaying outputs ----
                                    #h1("Analyses"),
                                    tabsetPanel(type = "tabs",
                                                tabPanel("Conjugation rates",
                                                         h3("Estimated conjugation rates"),
                                                         fluidRow(
                                                           column(width = 4,
                                                                  wellPanel(
                                                                    div(id = "info_conj_estimates", class = 'rightAlign', icon(name = "info-circle")),
                                                                    checkboxGroupInput("selected_estimates", 
                                                                                       "Select estimates to compute", 
                                                                                       selected = c("SM", "ASM"),
                                                                                       choiceNames = list("SM", 
                                                                                                          "ASM",
                                                                                                          "T/D", "T/DR",
                                                                                                          "T/(R+T)", "Dionisio: T/sqrt(DR)", "Gama: log(T/sqrt(DR))"),
                                                                                       choiceValues = list("SM", "ASM", "TD", "T_DR",
                                                                                                           "T_RT", "Dionisio", "Gama")),
                                                                    conditionalPanel(condition = "input.userdataTRT != null",
                                                                                     checkboxGroupInput("selected_estimatesTRT", 
                                                                                                        "Select estimates to compute", 
                                                                                                        selected = c("gamma.T"),
                                                                                                        choiceNames = list("\\(\\gamma_T\\)"),
                                                                                                        #"Minimal critical time",
                                                                                                        #"Time at which stationary phase is reached"), 
                                                                                                        choiceValues = list("gamma.T")),
                                                                                     #,  "tcrit", "tstat"
                                                                                     br()
                                                                    ),
                                                                    conditionalPanel(condition = "output.conj_estimates != null",
                                                                                     downloadButton("downloadData", "Download"),
                                                                                     bsTooltip("downloadData",
                                                                                               'Download a csv file with the estimated conjugation rates.',
                                                                                               placement = "right", trigger = "hover")
                                                                    ) #conditionalpanel
                                                                  )), #wellpanel + column
                                                           column(width = 8, 
                                                                  withMathJax(dataTableOutput("conj_estimates")))
                                                         ) # fluidrow
                                                ), #tabpanel
                                                # tabPanel("Boxplot",
                                                #          plotOutput("boxplot")),
                                                tabPanel("Critical times",
                                                         h3("Estimated critical times"),
                                                         fluidRow(
                                                           column(width = 4,
                                                                  wellPanel(
                                                                    HTML("The critical times indicate when the ",
                                                                             "assumptions underlying the (A)SM break down:<br>",
                                                                             "- tcrit1 describes when conjugation from transconjugants outweighs conjugation from donors.<br>",
                                                                             "- tcrit2, tcrit3 describe when the recipient population is depleted by ",
                                                                             "conjugation from donors or transconjugants respectively.<br><br>",
                                                                             "Estimates obtained after the first critical time is reached, ",
                                                                             "will most likely deviate from the true value of the conjugation rate."),
                                                                    br(),
                                                                    downloadButton("download_crittime", "Download critical time estimates"),
                                                                    bsTooltip("download_crittime",
                                                                              'Download a csv file with the estimated critical times.',
                                                                              placement = "right", trigger = "hover")
                                                                  )), #wellpanel + column 4
                                                           column(width = 8, 
                                                                  #div(id = "info_crittime_estimates", class = 'leftAlign', icon(name = "info-circle")),
                                                                  conditionalPanel(condition = "input.userdataTRT == null",
                                                                                   p("To estimate all critical times, ",
                                                                                     "please upload data for the TRT experiment.")),
                                                                  withMathJax(tableOutput("time_estimates"))
                                                           ) #column width 8
                                                         ), # fluidrow
                                                         br(),
                                                         p("Prediction of the minimal critical time",
                                                           " based only on the data of the DRT experiment."),
                                                         plotOutput("tcrit_plot"),
                                                         br()
                                                ) # tabpanel
                                    ) # tabsetpanel (estimates)
                           ), #tabpanel experimental data
                           ##################################
                           # Simulation panel ----
                           tabPanel("Simulate population dynamics",
                                    fluidRow(
                                      column(width = 12, withMathJax(includeMarkdown("simulation_plot_intro.md")))
                                    ),
                                    fluidRow(
                                      column(width = 6, withMathJax(includeMarkdown("simulation_plot_explanation.md"))),
                                      ## Model inputs ----
                                      column(width = 6, wellPanel(
                                        div(id = "info_model", class = 'rightAlign', icon(name = "info-circle")),
                                        selectInput("model", "Model to use",
                                                    choices = list("Simonsen Model" = "model_SM", 
                                                                   "Simonsen Extended Model" = "model_ESM",
                                                                   "Approximate Extended Simonsen Model" = "model_ASM"), 
                                                    selected = "model_ESM"),
                                        sliderInput("tend", label = "Measurement time point (hr)",
                                                    min = 0, max = 48, value = 24, step = 0.2),
                                        div(id = "info_crittimes", class = 'rightAlign', icon(name = "info-circle")),
                                        checkboxInput("t.crit", "Display critical times", value = TRUE),
                                        downloadButton("saveplot", "Plots"),
                                        downloadButton("save_output", "Simulation Output\n and Parameters")
                                      )) # column/wellpanel
                                    ),
                                    fluidRow(
                                      column(width = 4, h3("Population dynamics", align = "center"), 
                                             withSpinner(plotOutput("interactive_dyn_plot"), 
                                                         size = .5, type = 3, 
                                                         color = '#2c3e50', color.background = 'white')),
                                      column(width = 4, h3("Estimated growth rate", align = "center"),
                                             withSpinner(plotOutput("interactive_growth_plot"), size = .5,
                                                         type = 3,
                                                         color = '#2c3e50', color.background = 'white')),
                                      column(width = 4, h3("Estimated conjugation rate", align = "center"),
                                             withSpinner(plotOutput("interactive_conj_plot"),
                                                         size = .5, type = 3,
                                                         color = '#2c3e50', color.background = 'white'))
                                    ),
                                    
                                    fluidRow(
                                      ## Population size inputs ----
                                      column(width = 3, wellPanel(
                                        div(id = "info_init_pop", class = 'rightAlign', icon(name = "info-circle")),
                                        sliderInput("logD0", label = "Initial donor population size \\(D_0\\)",
                                                    min = 4, max = 8, value = 6, step = 0.5),
                                        sliderInput("logR0", label = "Initial recipient population size \\(R_0\\)",
                                                    min = 4, max = 8, value = 6, step = 0.5),
                                        checkboxInput("nonzero_logT0", "Nonzero initial transconjugant population", value = FALSE),
                                        conditionalPanel(condition = "input.nonzero_logT0",
                                                         sliderInput("logT0", label = "Initial transconjugant population size \\(T_0\\)",
                                                                     min = 0, max = 8, value = 0, step = 0.5)
                                        )
                                      )), # column + wellpanel
                                      
                                      ## Growth rate inputs ----
                                      column(width = 3, wellPanel(
                                        div(id = "info_growthrates", class = 'rightAlign', icon(name = "info-circle")),
                                        conditionalPanel(condition = "input.model != 'model_SM'",
                                                         sliderInput("psiD", label = "Donor growth rate \\(\\psi_D\\)",
                                                                     min = 0, max = 2, value = 1, step = 0.1),
                                                         sliderInput("psiR", label = "Recipient growth rate \\(\\psi_R\\)",
                                                                     min = 0, max = 2, value = 1, step = 0.1),
                                                         sliderInput("psiT", label = "Transconjugant growth rate \\(\\psi_T\\)",
                                                                     min = 0, max = 2, value = 1, step = 0.1)
                                        ),
                                        conditionalPanel(condition = "input.model == 'model_SM'",
                                                         sliderInput("psi", label = "Growth rate \\(\\psi_{max}\\)",
                                                                     min = 0.1, max = 2, value = 1, step = 0.1)
                                        )
                                      )),
                                      ## Conjugation rate inputs ----
                                      column(width = 3, wellPanel(
                                        div(id = "info_conjrates", class = 'rightAlign', icon(name = "info-circle")),
                                        conditionalPanel(condition = "input.model != 'model_SM'",
                                                         sliderInput("log_gammaD", 
                                                                     label = "Donor conjugation rate \\(\\gamma_D\\)",
                                                                     min = -14, max = -5, value = -11, step = 0.5),
                                                         sliderInput("log_gammaT", 
                                                                     label = "Transconjugant conjugation rate \\(\\gamma_T\\)",
                                                                     min = -14, max = -5, value = -11, step = 0.5)
                                        ),
                                        conditionalPanel(condition = "input.model == 'model_SM'",
                                                         sliderInput("log_gamma", 
                                                                     label = "Conjugation rate \\(\\gamma_{max}\\)",
                                                                     min = -14, max = -5, value = -11, step = 0.5)          
                                        )
                                      )),
                                      ## Resource inputs ----
                                      
                                      conditionalPanel(condition = "input.model != 'model_ASM'",
                                                       column(width = 3, wellPanel(
                                                         div(id = "info_init_resource", class = 'rightAlign', icon(name = "info-circle")),
                                                         sliderInput("logC0", label = "Initial resource concentration \\(C_0\\)",
                                                                     min = 5, max = 20, value = 12, step = 0.5, width = "100%"),
                                                         br(),
                                                         actionButton("show_C0_widget", "Need help to select \\(C_0\\)?"),
                                                         ## C0 widget ----
                                                         conditionalPanel(condition = "input.show_C0_widget % 2 == 1",
                                                                          hr(),
                                                                          div(id = "info_C0_widget", class = 'rightAlign', icon(name = "info-circle")),
                                                                          selectInput("C0_widget_var", "Desired variable to fix:",
                                                                                      c("Time to stationary phase" = "tstat",
                                                                                        "Final Donor population size" = "Dfinal"),
                                                                                      selected = 'Dfinal'),
                                                                          numericInput("C0_widget_val", "Fix variable to: ", '10E7'),
                                                                          htmlOutput("C0_widget_out"),
                                                                          br(),
                                                                          br(),
                                                                          actionButton("C0_widget_setval", "Set C0 slider to this value")
                                                                          
                                                         ) # conditional panel
                                                       )) # outer C0 column and wellpanel
                                      ) # conditional panel
                                      
                                    ), # fluidrow
                                    
                                    
                                    
                           ), # tabpanel
                           ## Navbar page options ----
                           id = "navbar", selected = "Home", collapsible = TRUE
                ) #navbarpage
) #fluidpage

#############################################################################################
# Define server logic 
server <- function(input, output, session) {
  
  ##################################
  # input data ----
  
  # like variables but notify downstream dependencies
  # data is displayed DRT, user_dataDRT is saved user upload
  dataDRT <- reactiveVal(complete_dataframe)
  user_dataDRT <- reactiveVal(complete_dataframe)
  user_dataTRT <- reactiveVal(complete_dataframe)
  
  # read new user data
  new_user_inputDRT <- eventReactive(input$userdataDRT$datapath, {
    return(input_user_data(input$userdataDRT$datapath))
  })
  new_user_inputTRT <- eventReactive(input$userdataTRT$datapath, {
    return(input_user_data(input$userdataTRT$datapath))
  })
  
  # update displayed data (once user uploads data)
  observeEvent(input$userdataDRT$datapath, {
    joint_df <- try(join_and_overwrite_df(user_dataDRT(), new_user_inputDRT(), id_cols = input$id_cols))
    if (class(joint_df) == "try-error"){
      showNotification(paste0("This input file can not be appended ",
                              "to the existing data. Please click ",
                              "\"Remove Data\" and try again.\n",
                              "This error may be related to non-numeric data",
                              " in a measurement column."), 
                       action = NULL, duration = 5, closeButton = TRUE,
                       id = "append_input_error", type = "error",
                       session = getDefaultReactiveDomain())
    } else {
      user_dataDRT(joint_df)
      dataDRT(user_dataDRT())
      updateSelectInput(session, "id_cols", choices = colnames(dataDRT()), selected = colnames(dataDRT())[1])
    }
  })
  # same for TRT
  observeEvent(input$userdataTRT$datapath, {
    joint_df <- try(join_and_overwrite_df(user_dataTRT(), new_user_inputTRT(), id_cols = input$id_cols))
    if (class(joint_df) == "try-error"){
      showNotification(paste0("This input file can not be appended ",
                              "to the existing data. Please click ",
                              "\"Remove Data\" and try again.\n",
                              "This error may be related to non-numeric data",
                              " in a measurement column."), 
                       action = NULL, duration = 5, closeButton = TRUE,
                       id = "append_input_error", type = "error",
                       session = getDefaultReactiveDomain())
    } else {
      user_dataTRT(joint_df)
    }
  })
  
  # uncheck exampledata if new user data is uploaded
  # This is only for DRT!
  observeEvent(input$userdataDRT$datapath, {
    req(input$exampledata)
    updateCheckboxInput(session, "exampledata", value = FALSE)
  }, ignoreInit = TRUE)
  
  # remove data
  observeEvent(input$trash_dataDRT, {
    user_dataDRT(complete_dataframe)
    dataDRT(user_dataDRT())
    updateSelectInput(session, "id_cols", choices = colnames(dataDRT()), selected = "ID")
  })
  observeEvent(input$trash_dataTRT, {
    user_dataTRT(complete_dataframe)
  })
  
  # uncheck exampledata if "remove data" is clicked
  observeEvent(input$trash_dataDRT, {
    req(input$exampledata)
    updateCheckboxInput(session, "exampledata", value = FALSE)
  }, ignoreInit = TRUE)
  
  ## show results on example data
  # defaulting back to userdata is why we can't turn this
  # into an ActionButton!
  observeEvent(input$exampledata, {
    if (input$exampledata){
      # example data selected, show example data
      dataDRT(example_data)
    } else if (!is.null(input$userdataDRT)){
      # example data unselected, user data loaded, show user data
      dataDRT(user_dataDRT())
    } else {
      dataDRT(complete_dataframe)
    }
    updateSelectInput(session, "id_cols", choices = colnames(dataDRT()), selected = colnames(dataDRT())[1])
  }) 
  
  # Allow editable datatables
  observeEvent(input$editdataDRT_cell_edit, {
    dataDRT(editData(dataDRT(), input$editdataDRT_cell_edit, 'editdataDRT'))
  })
  observeEvent(input$editdataTRT_cell_edit, {
    user_dataTRT(editData(user_dataTRT(), input$editdataTRT_cell_edit, 'editdataTRT'))
  })
  
  ##################################
  # display input data ----
  output$editdataDRT <- renderDT({dataDRT()}, 
                                 selection = 'none', 
                                 options  = list(dom = 't', pageLength = -1, ordering = FALSE), 
                                 rownames = FALSE,
                                 server = TRUE, editable = 'cell')
  
  output$editdataTRT <- renderDT({user_dataTRT()}, 
                                 selection = 'none', 
                                 options  = list(dom = 't', pageLength = -1, ordering = FALSE), 
                                 rownames = FALSE,
                                 server = TRUE, editable = 'cell')
  
  ##################################
  # output estimates crit time ----
  
  data_time_crit <- reactive({
    # Do not compute if data() is empty
    req(dim(dataDRT())[1]>1)
    
    time_crit <- suppressWarnings(estimate_crit_time(dataDRT(), user_dataTRT(), tol_factor = 10, id_cols = input$id_cols))
    
    return(time_crit)
  })
  
  ## table ##
  output$time_estimates <- renderTable({
    req(data_time_crit(), cancelOutput = TRUE)
    data_time_crit()
  })
  
  ## plot ##
  output$tcrit_plot <- renderPlot({
    # Do not compute if data() is empty
    req(dim(dataDRT())[1]>1)
    req(length(input$id_cols)>0)
    
    plot_crittimes_sweep(dataDRT(), id_cols = input$id_cols)
  })
  
  ##################################
  # output estimates conjugation rates ----
  data_conj_rates <- reactive({
    # Do not compute if dataDRT() is empty
    req(dim(dataDRT())[1]>1)
    if (length(input$id_cols)<1){
      ids = "ID"
    } else {ids = input$id_cols}
    
    long_conj_rates <- estimate_conj_rate(dataDRT(), method = input$selected_estimates, id_cols = input$id_cols)
    conj_rates <- pivot_wider(long_conj_rates, id_cols = all_of(ids), 
                              names_from = 'method', values_from = 'estimate')
    
    return(conj_rates)
  })
  
  output$conj_estimates <- renderDataTable({
    DT::datatable(data_conj_rates(),
                  rownames = FALSE,
                  options = list(dom = 't', pageLength = -1, ordering = FALSE,
                                 rowCallback = JS(
                                   "function(row, data) {
                                   for (i = 1; i < data.length; i++) {
                                   $('td:eq('+i+')', row).html(parseFloat(data[i]).toExponential(1));
                                   }}"
                                 ) 
                  )
    )
  })
  
  ##################################
  # boxplot ----
  # output$boxplot <- renderPlot({
  #   req(all(c('ASM', 'SM') %in% colnames(data_conj_rates())))
  #   boxplot_scan(data_conj_rates(), 'ID', 'ASM', 'SM',
  #                logplot = FALSE, add_means = FALSE, plotfile = NULL)
  # })
  
  ##################################
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = "conjugation_rates.csv",
    content = function(file) {
      write.csv(data_conj_rates(), file, row.names = FALSE)
    }
  )
  
  output$download_crittime <- downloadHandler(
    filename = "crit_time_estimates.csv",
    content = function(file) {
      write.csv(data_time_crit(), file, row.names = FALSE)
    }
  )
  
  ####################################################################
  ## Linked sliders simonsen model ----
  observeEvent(c(input$model), {
    req((input$model == 'model_SM'))
    
    # Control the value of psiD, psiR, and psiT
    updateSliderInput(session, "psiD", value = input$psi)
    updateSliderInput(session, "psiR", value = input$psi)
    updateSliderInput(session, "psiT", value = input$psi)
  })
  
  ###########################
  ## C0 widget ----
  suggested_C0 <- reactive({
    if (input$C0_widget_var == 'tstat'){
      suggested_C0 <- calculate_C0(input, tstat = input$C0_widget_val, Dfinal = NULL)
    } else {
      suggested_C0  <- calculate_C0(input, tstat = NULL, Dfinal = input$C0_widget_val)
    }
    return(suggested_C0)
  })
  
  output$C0_widget_out <- reactive({ 
    HTML('Suggested C0: 10<sup>', round(log10(suggested_C0()), digits = 1),'</sup>' )
  })
  
  observeEvent(input$C0_widget_setval, {
    # Control the value of C0
    updateSliderInput(session, "logC0", value = log10(suggested_C0()))
    
    JS("$('#'+'logC0').data('ionRangeSlider').update({
          'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
        })"
    )
  })
  
  ###########################
  # interactive simulation plot ----
  
  ## Inputs ##
  vars <- reactive({
    v = c(r = 10^as.numeric(input$logR0),
          t = ifelse(input$nonzero_logT0, 10^as.numeric(input$logT0), 0),
          d = 10^as.numeric(input$logD0))
    
    if (input$model == "model_ASM"){
      return(v)
    } else {
      return(c(v, c=10^as.numeric(input$logC0)))
    }
  })
  
  parms <- reactive({
    growth_rates = c(psi.R = as.numeric(input$psiR),
                     psi.T = as.numeric(input$psiT),
                     psi.D = as.numeric(input$psiD))
    conj_rates = c(gamma.t = 10^(as.numeric(input$log_gammaT)),
                   gamma.d = 10^(as.numeric(input$log_gammaD)))
    
    if (input$model == "model_ASM"){
      # no resource dynamics
      return(c(growth_rates, conj_rates))
    } else if (input$model == "model_SM"){
      # linked parameters
      growth_rates = c(psi.R=as.numeric(input$psi),
                       psi.T=as.numeric(input$psi),
                       psi.D=as.numeric(input$psi))
      conj_rates = c(gamma.t=10^(as.numeric(input$log_gamma)),
                     gamma.d=10^(as.numeric(input$log_gamma)))
      return(c(growth_rates, conj_rates, q=1E2))
    } else {
      # simonsen extended model
      return(c(growth_rates, conj_rates, q=1E2))
    }
  })
  
  ## Results ##
  out <- reactive({
    t.vec <- seq(0, input$tend, length = 100)
    out <- integrate_model(vars(), t.vec, input$model, parms())
    return(out)
  })
  
  crit_time <- reactive({
    if(input$t.crit) {
      crit_times <- estimate_crittime_from_sim(parms(), vars(), tol_factor = 10)
    } else {
      crit_times = NULL
    }
    return(crit_times)
  })
  
  ## Plot Outputs ##
  dyn_plot <- reactive({plot_dynamics_ggplot(out(), crit_time() )})
  output$interactive_dyn_plot <- renderPlot({ dyn_plot() })
  
  growth_plot <- reactive({ plot_growth_rate_ggplot(out(), parms(), crit_time() ) })
  output$interactive_growth_plot <- renderPlot({ growth_plot() })
  
  conj_plot <- reactive({ plot_conj_rate_ggplot(out(), parms(), crit_time() ) })
  output$interactive_conj_plot <- renderPlot({ conj_plot() }) 
  
  
  ###########################
  # Downloadable plots of simulation ----
  output$saveplot <- downloadHandler(
    filename = function() {
      paste('popdynamic_simulation', ".pdf", sep = "")
    },
    content = function(file) {
      all_p_row <- plot_grid( dyn_plot(), growth_plot(), conj_plot(),
                              nrow = 1, ncol = 3,
                              label_size = 15) 
      
      title <- ggdraw() +
        draw_label(
          "Produced using Conjugator. Huisman et al. 2020",
          fontface = 'bold',
          x = 0,
          hjust = 0
        ) +
        theme(
          # add margin on the left of the drawing canvas,
          # so title is aligned with left edge of first plot
          plot.margin = margin(0, 0, 0, 7)
        )
      
      all_p <- plot_grid(
        all_p_row, title, 
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(1, 0.1)
      )   
      
      ggsave(filename = file,
             plot = all_p, width = 14, height = 7)
    }
  )
  
  ###########################
  # Download csv of simulation output ----
  output$save_output <- downloadHandler(
    filename = "popdynamic_simulation.csv",
    content = function(file) {
      # all initial parameters; final output;
      # and measurement timepoint + model
      current_vars <- vars()
      if(length(current_vars)== 4){
        names(current_vars) <- c("R0", "T0", "D0", "C0")
      } else {
        names(current_vars) <- c("R0", "T0", "D0")
      }
      
      sim_params <- c(current_vars, 
                      parms()[c("psi.R", "psi.T", "psi.D",
                                "gamma.t","gamma.d")],
                      out()[dim(out())[1], c("r", "t", "d")],
                      input$tend, input$model)
      
      
      
      write.csv(sim_params, file, row.names = FALSE)
    }
  )
  
  ######################################################
  # Logos ----
  output$ethlogo <- renderImage({list(src = "www/logoethz.png", width = '75%')}, deleteFile = FALSE)
  output$snflogo <- renderImage({list(src = "www/SNF_logo_v2.png", width = '75%')}, deleteFile = FALSE)
  
  ###########################
  ## Info on data analysis ----
  addPopover(session, "info_DRT_data", "Upload DRT Data", 
             paste0("Please upload a csv file, following the example columns shown.<br>",
             "The required columns vary between methods. ",
             "D.t, R.t, T.t are final population sizes;<br>",
             "D.0, R.0 are initial pop sizes;<br>",
             "psi.D, psi.R, psi.T are growth rates;<br>",
             "t is the time of measurement.<br>",
             "Optional columns: T.0, psi.max"), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_conj_estimates", "Select conjugation rate estimates", 
             paste0("The SM method requires a mating population growth rate psi.max, and D.t, R.t, T.t, D.0, R.0.<br>",
                    "The ASM method requires psi.D, psi.R, psi.T, D.t, R.t, T.t, D.0, R.0.<br>",
                    "The other methods require final population sizes D.t, R.t, and/or T.t."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  # addPopover(session, "info_crittime_estimates", "Critical time estimates", 
  #            paste0("The critical times indicate when the ",
  #                   "assumptions underlying the (extended) Simonsen model",
  #                   " are expected to break down.<br>",
  #                   "tcrit1 describes when conjugation from transconjugants outweighs conjugation from donors.<br>",
  #                   "tcrit2, tcrit3 describe when the recipient population is depleted by ",
  #                   "conjugation from donors or transconjugants respectively.<br>",
  #                   "Estimates obtained after the first critical time is reached, ",
  #                   "will most likely deviate from the true value of the conjugation rate."), 
  #            placement = "top",
  #            trigger = "hover", options = list(container = "body"))
  
  ## Info on simulations ----
  addPopover(session, "info_init_pop", "Initial conditions", 
             paste0("These parameters govern the initial population sizes",
                    " in CFU/mL.<br>"), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_init_resource", "Initial resources", 
             paste0("These parameters govern the initial resources",
                    " in ug/mL (where we assume the amount of resource needed to",
                    " produce one new cell is 1 ug/CFU).<br>",
                    "If the ASM is selected, the initial",
                    " resources are assumed to be infinite.<br>",
                    "If the resource concentrations is not known, or one ",
                    "expects stationary phase not to be reached, we suggest ",
                    "to use the ASM or to set the maximal value ",
                    "for C0 when using the Simonsen/extended Simonsen method."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_C0_widget", "Initial resources", 
             # paste0("It is often easier to measure the final density of donors D",
             #        " or specify the time at which the stationary growth phase starts,",
             #        " than to know the resource concentration needed to reach that density",
             #        " or timing. Thus we offer a rough calculation to translate between these ",
             #        "measures. The final density of recipients R or transconjugants T can not",
             #        " be selected because this changes as a function of the conjugation process."), 
             paste0("Fill in the desired time when stationary phase is ",
                    "reached (e.g. 10 hrs) or the desired population size ",
                    "of donors at the measurement time (e.g. 10E7 CFU/mL), ",
                    "and we will suggest a C0 value to set for the simulations."),
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_growthrates", "Growth rates", 
             paste0("These parameters govern the growth rates",
                    " of the bacterial populations. Rates are measured",
                    " per hour.<br>",
                    "If the Simonsen model is selected,",
                    " growth rates are assumed equal for all populations.<br>",
                    "In the ESM and ASM, ",
                    "growth rates can be set independently for all populations."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_conjrates", "Conjugation rates", 
             paste0("These parameters govern the conjugation rates.<br>",
                    "If the Simonsen model is selected, ",
                    "conjugation rates are assumed equal for all populations.<br>",
                    "In the ESM and ASM, ",
                    "conjugation rates can be set independently for all populations."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_model", "Simulation model", 
             paste0("Three population dynamic models can be selected:<br>",
                    "the Simonsen model, which assumes all growth and conjugation",
                    " rates are the same,<br>",
                    "the extended Simonsen model, which relaxes these assumptions,<br>",
                    "and the ASM, which simplifies the extended Simonsen ",
                    "model to a case with infinite resources."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
  addPopover(session, "info_crittimes", "Critical times", 
             paste0("The critical times indicate when the ",
                    "assumptions underlying the (extended) Simonsen model",
                    " are expected to break down.<br>",
                    "Estimates obtained after the first critical time is reached, ",
                    "will most likely deviate from the true value of the conjugation rate."), 
             placement = "top",
             trigger = "hover", options = list(container = "body"))
  
}

###########################################################
# Run the application 
shinyApp(ui = ui, server = server)
