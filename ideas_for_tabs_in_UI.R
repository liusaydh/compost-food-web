## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- message=FALSE-------------------------------------------------------------------------------
require(shiny)
require(shinyWidgets)

## ---- message = FALSE-----------------------------------------------------------------------------
require(deSolve)

## --state variebles, model-state is passes as a vector------------------------------------------------------
Pads.tot <- 0
x.fast <- 0.5
state <- c(Fe = 10.3, 
           S = 0, 
           Sads = 0, 
           P = 0, 
           Pads.fast = Pads.tot*x.fast, 
           Pads.slow = Pads.tot*(1-x.fast)
) # both in [mmol/L], Fe=conc.in the reactor, 
#P=initial P, Pads=adsorbed P porpotional to Fe loading
#!Important! the order of these state variables is important, i.e. time-derivatives
#must follow the same order(and maybe the other parts)
#default.parms is set as a list for the reason to be updated in shiny
default.parms <- list(
  D         = 1/25, # [/min]     Q = 1 mL/min, V = 25 mL
  m         = 0.5, #1/2,  # [-]
  n         = 0.8, #2/3,  # [-]
  y         = 1, # y=exponent S^(y)
  kFe       = 10^(-2.76),    # [mmol Fe min-1]
  kPads      = 1e-4,
  kPads2kdes = 0.33, # kdes = kPads*kPads2kdes
  kSads      = 1e-10,
  kSads2kdes = 1, #?
  kS.fast   = sqrt(10^(1.3)),
  kS.slow   = sqrt(10^(1.3)),
  Sin       = 1.8*1,           # [mmol S], inflow Stot concentration, should be the max(breakthrough curve)
  P.ini     = state[["P"]],
  Fe.ini    = state[["Fe"]],   #initial mineral loading
  x.fast    = x.fast,
  Pads.ini  = state[["Pads.fast"]] + state[["Pads.slow"]]
)

FeSPmodel <- function(t, state, pars) {
  with (as.list(c(state, pars)),{
    
    # to avoid that the power of a tiny negative number gives NaN
    Fe <- max(0, Fe)
    S <- max(0, S)
    P <- max(0, P)
    Sads <- max(0,Sads)
    Pads.slow <- max(0, Pads.slow)
    Pads.fast <- max(0, Pads.fast)
    
    # rate expressions [mol/m3/d]
    R    <- kFe * S^m * Fe.ini^(1-n) * Fe^n
    Rads <- kPads * P  * Fe.ini^(1-n) * Fe^n
    kdes <- kPads/kPads2kdes
    RSads <- kSads * S  * Fe.ini^(1-n) * Fe^n
    kSdes <- kSads/kSads2kdes
    RSdes <- kSdes * Sads
    Rdes.fast <- kdes * Pads.fast
    Rdes.slow <- kdes * Pads.slow
    RHS.fast <- kS.fast^y * Pads.fast * S^y
    RHS.slow <- kS.slow^y * Pads.slow * S^y
    
    #x.fast <- ifelse(Pads.fast+Pads.slow>0, Pads.fast/(Pads.fast + Pads.slow), 0.5)  
    #?cause trouble when x.fast != 0.5 and returns x.fast = 0.5, bug condition unclear(? only shows up when x.fast = 0.5?)
    
    # Time-derivatives:
    dFe.dt   <- -R
    dS.dt    <- -3/2*R + D*Sin - D*S - RSads + RSdes
    dSads.dt <- RSads - RSdes
    dP.dt    <-     -Rads + Rdes.fast + Rdes.slow + RHS.fast + RHS.slow - D*P
    dPads.fast.dt <- Rads*x.fast      - Rdes.fast - RHS.fast
    dPads.slow.dt <- Rads*(1-x.fast)  - Rdes.slow - RHS.slow
    
    
    list(c(dFe.dt, dS.dt, dSads.dt, dP.dt, dPads.fast.dt, dPads.slow.dt), # vector with derivatives 
         # other output: process rates
         R = R,
         Rads = Rads,
         Rdes.slow = Rdes.slow,
         Rdes.fast = Rdes.fast,
         RHS.slow = RHS.slow,
         RHS.fast = RHS.fast,
         x.fast = x.fast,
         RSads = RSads,
         RSdes = RSdes
    )
  })
}


## ----data import--------------------------------------------------------------------------------
outtimes <- seq(from = 0, to = 450, length.out = 451)

#initial conditions:
Default  <- ode(y=state, parms=default.parms, func=FeSPmodel, times=outtimes, method = "vode")

## specify experimental data

# base-directory with the experimental data:
# "~/Desktop/R_ptrap/" for Mingkai
# "~/3_github/P_trap/src/FeSPads/data/"
bdir <- "src/FeSPads/data/"

# input file names,include all data sets we have
data_files <- c("pH7373R1.csv",
                "pH7373R2.csv",
                "pH7373R3.csv",
                "pH7373R4.csv",
                "pH7373R5.csv",
                "pH7383R1.csv",
                "pH7383R2.csv",
                "pH7383R3.csv",
                "pH7383R4.csv",
                "pH7383R5.csv",
                "pH8373R1.csv",
                "pH8373R2.csv",
                "pH8373R3.csv",
                "pH8373R4.csv",
                "pH8373R5.csv",
                "pH8383R1.csv",
                "pH8383R2.csv",
                "pH8383R3.csv",
                "pH8383R4.csv",
                "pH8383R5.csv",
                "buffer7373.csv",
                "buffer7383.csv",
                "buffer8373.csv",
                "buffer8383.csv")
#fun for importing data as dataframe

importdata <- function(i){
  bdir <- "src/FeSPads/data/"
  data <- read.csv(file = paste(bdir, data_files[i], sep=""))
  FeSPads.data <- data.frame(data)
  return(FeSPads.data)
}

#creat 24 datasets as import exp data.
for(i in 1:length(data_files)){
  assign(paste("FeSPads.data",i,sep = ""),importdata(i = i))
}

## ---UI-----------------------------------------------------------------------------
UI.FeSP <- shinyUI(fluidPage(fluid=TRUE,      # Define UI (user interface)
                             # Application title
                             headerPanel(strong("FePS model")),
                             sidebarLayout(
                               sidebarPanel(fluid = T,width = 8,
                                            fluidRow(fluid = T,
                                                     column(4, offset = 0,
                                                            uiOutput("tabpanelUI")),
                                                     column(4, offset = 0,
                                                            prettyCheckboxGroup(inputId = "plot_check",choiceValues = c(1:24),#inline = T,
                                                                                choiceNames = data_files,
                                                                                selected = "1",
                                                                                label=strong("check to plot these exp datasets"))),
                                                     column(4, offset = 0,
                                                            prettyCheckbox(inputId="appendix",value=FALSE,shape = "curve",outline = FALSE,
                                                                           animation = "tada",label=strong("Show equations")),
                                                            prettyCheckbox(inputId = "show_parm",value = FALSE,shape = "curve",outline = FALSE,
                                                                           animation = "tada",label =strong("show fitted parms.")),
                                                            prettyCheckbox(inputId = "show_plot",value = TRUE,shape = "curve",outline = FALSE,
                                                                           animation = "tada",label =strong("show plot")),
                                                            
                                                     ),
                                                     column(4, offset = 0,
                                                            prettyCheckboxGroup(inputId="fit_select",shape = "round",outline = F,fill = T,
                                                                                choices = c("fit1","fit2","fit3","fit4","fit5"),
                                                                                label = strong("select for fitting"),selected =c("fit1","fit2","fit3","fit4","fit5"))
                                                     ),
                                                     column(4,offset =0,
                                                            prettyCheckboxGroup(inputId="plot_elements",shape = "curve",outline = F,fill = T,
                                                                                choices = colnames(Default),
                                                                                label = strong("select plot elements"),
                                                                                selected =c("Fe", "S", "P", "Pads.fast",  "R", "RSdes", "RSads","Sads")))
                                            )),
                               mainPanel(
                                 actionBttn(inputId = "save_fit_parm1",size = "sm",style = "unite",no_outline = T,
                                            label =strong("save Fit1")),
                                 actionBttn(inputId = "save_fit_parm2",size = "sm",style = "unite",no_outline = T,
                                            label =strong("save Fit2")),
                                 actionBttn(inputId = "save_fit_parm3",size = "sm",style = "unite",no_outline = T,
                                            label =strong("save Fit3")),
                                 actionBttn(inputId = "save_fit_parm4",size = "sm",style = "unite",no_outline = T,
                                            label =strong("save Fit4")),
                                 actionBttn(inputId = "save_fit_parm5",size = "sm",style = "unite",no_outline = T,
                                            label =strong("save Fit5")),
                               ) 
                               
                             ),
                             HTML('<hr style="color: purple;">'),    #add a line to split the two blocks
                             mainPanel(width = 12, 
                                       conditionalPanel(condition = "input.show_plot",
                                                        plotOutput("PlotFeSP"),#plot the output
                                       ),
                                       verbatimTextOutput("verb"),   #show helptext, with the checkboxInput
                                       conditionalPanel(condition = "input.show_parm",
                                                        dataTableOutput("table_parmfit"),
                                       ),
                                       
                             ),
)
)

## ---server----------------------------------------------------------------------------------------
Server.FeSP <- shinyServer(function(input, output, session) {
  # 
  
  
  #renderUI, put UI back to server, we could also do it in UI part
  CreatUI <- renderUI({
    mainPanel(width = 12,
              tabsetPanel(id = "tabselected", 
                          tabPanel(strong("ini. conditions"), value = 1, helpText("initial conditions...and settings")),
                          tabPanel(strong("sulfidation parm"), value =2, helpText("sulfidation paramaeters...")),
                          tabPanel(strong("phosphate parm"),value = 3, helpText("phosphate release parameters..."))
              ),
              conditionalPanel(condition = "input.tabselected == 1",
                               column(4, offset = 2,
                                      numericInput(width = '100%', inputId="tau",
                                                   label = "t (min), V/Q",
                                                   min = 15, max = 35, step = 1, value = 1/default.parms$D),
                                      numericInput(width = '100%', inputId="Sin",
                                                   label = "Sin",
                                                   min = 0, max = 2, step = 0.1, value = default.parms$Sin),
                                      numericInput(width = '100%', inputId="Fe.ini",
                                                   label = "Fe.ini",
                                                   min = 0, max = 20, step = 0.1, value = default.parms$Fe.ini)
                               ),
                               column(4, offset = 2,
                                      pickerInput(label = "click4default fit",multiple = F,
                                                  inputId = "expdataID",
                                                  choices = as.list(data_files),
                                                  selected = "buffer8383.csv")
                               )
              ),
              conditionalPanel(condition = "input.tabselected == 2",
                               fluidRow(
                                 column(4, offset = 2,
                                        numericInput(width = '100%', inputId="n", 
                                                     label = "Exponent n",
                                                     min = 0.1, max = 2, step = 0.05, value = default.parms$n),
                                        numericInput(width = '100%', inputId="m", 
                                                     label = "Exponent m",
                                                     min = 0.3, max = 2, step = 0.05, value = default.parms$m),
                                        numericInput(width = '100%', inputId="kFe", 
                                                     label = "log(kFe)", 
                                                     min = -4, max = 1, step = 0.02, value = log10(default.parms$kFe)),
                                 ),
                                 column(4, offset = 2,
                                        numericInput(width = '100%', inputId="kSads2kdes", 
                                                     label = "kSads/kdes",
                                                     min = 0.0001, max = 5, step = 0.01, value = default.parms$kSads2kdes),
                                        numericInput(width = '100%', inputId="kSads", 
                                                     label = "log(kSads)",
                                                     min = -16, max = 2, step = 0.05, value = log10(default.parms$kSads)),
                                        numericInput(width = '100%', inputId="y", 
                                                     label = "exponent y",
                                                     min = 0.3, max = 4, step = 0.05, value = default.parms$y)
                                 ),
                               ),
              ),
              conditionalPanel(condition = "input.tabselected == 3",
                               fluidRow(
                                 column(4, offset = 1,
                                        numericInput(width = '100%', inputId="kPads", 
                                                     label = "log(kPads)", 
                                                     min = -5, max = -1, step = 0.02, value = log10(default.parms$kPads)),
                                        numericInput(width = '100%', inputId="kPads2kdes", 
                                                     label = "kPads/kdes",
                                                     min = 0.2, max = 1, step = 0.01, value = default.parms$kPads2kdes)
                                 ),
                                 column(4, offset = 2,
                                        numericInput(width = '100%', inputId="kS.fast", 
                                                     label = "log(kS.fast)", 
                                                     min = -4, max = 4, step = 0.02, value = log10(default.parms$kS.fast)),
                                        numericInput(width = '100%', inputId="x.fast", 
                                                     label = "x.fast",
                                                     min = 0.0, max = 1, step = 0.05, value = default.parms$x.fast),
                                        numericInput(width = '100%', inputId="kS.slow", 
                                                     label = "log(kS.slow)", 
                                                     min = -4, max = 4, step = 0.02, value = log10(default.parms$kS.slow)),
                                 ),
                                 column(4,offset = 2,
                                        numericInput(width = '100%', inputId="P.ini", 
                                                     label = "P.ini",
                                                     min = 0, max = 0.03, step = 0.001, value = default.parms$P.ini),
                                        numericInput(width = '100%', inputId="Pads.ini", 
                                                     label = "Pads.ini",
                                                     min = 0, max = 0.6, step = 0.01, value = default.parms$Pads.ini)
                                 )),
              ),
    )
  }
  
  
  )
  
  output$tabpanelUI <-renderUI({  
    CreatUI})
  
  
  # Get the model parameters, as defined in the UI, this is a function,return list(parms)
  getparms_interface <- reactive( {#update parms from the interface
    parms        <- default.parms 
    parms$D     <- 1/input$tau
    parms$n     <- input$n
    parms$m     <- input$m
    parms$Fe.ini   <- input$Fe.ini
    parms$x.fast     <- input$x.fast
    parms$y     <- input$y
    parms$kFe   <- 10^input$kFe
    parms$kPads   <- 10^input$kPads
    parms$kPads2kdes <- input$kPads2kdes
    parms$kS.fast   <- 10^input$kS.fast
    parms$kS.slow   <- 10^input$kS.slow
    parms$Sin <- input$Sin
    parms$P.ini <- input$P.ini
    parms$Pads.ini <- input$Pads.ini
    parms$kSads  <- 10^input$kSads
    parms$kSads2kdes  <- input$kSads2kdes
    #parms$Fe <- input$Fe
    parms
  })
  
  
  
  
  #set conditions for parm.1-5,state.1-5 
  
  parm_reactive <- reactiveValues(parm_fit1 = default.parms,
                                  parm_fit2 = default.parms,
                                  parm_fit3 = default.parms,
                                  parm_fit4 = default.parms,
                                  parm_fit5 = default.parms)
  
  
  
  observeEvent(input$save_fit_parm1,{
    updatePickerInput(session, inputId = "expdataID",selected = "none")#fake none, return nothing
    parm.IM.1 <- getparms_interface()
    parm_reactive$parm_fit1 <- parm.IM.1
  })
  observeEvent(input$save_fit_parm2,{
    updatePickerInput(session, inputId = "expdataID",selected = "none")
    parm.IM.2 <- getparms_interface()
    parm_reactive$parm_fit2 <- parm.IM.2
  })
  observeEvent(input$save_fit_parm3,{
    updatePickerInput(session, inputId = "expdataID",selected = "none")
    parm.IM.3 <- getparms_interface()
    parm_reactive$parm_fit3 <- parm.IM.3
  })
  observeEvent(input$save_fit_parm4,{
    updatePickerInput(session, inputId = "expdataID",selected = "none")
    parm.IM.4 <- getparms_interface()
    parm_reactive$parm_fit4 <- parm.IM.4
  })
  observeEvent(input$save_fit_parm5,{
    updatePickerInput(session, inputId = "expdataID",selected = "none")
    parm.IM.5 <- getparms_interface()
    parm_reactive$parm_fit5 <- parm.IM.5
  })
  
  state_reactive_1 <- reactive({
    state_fit1 <- state
    state_fit1[["P"]] <- parm_reactive$parm_fit1$P.ini
    state_fit1[["Pads.fast"]] <- parm_reactive$parm_fit1$Pads.ini *
      parm_reactive$parm_fit1$x.fast
    state_fit1[["Pads.slow"]] <- parm_reactive$parm_fit1$Pads.ini * 
      (1-parm_reactive$parm_fit1$x.fast)
    state_fit1[["Fe"]] <- parm_reactive$parm_fit1$Fe.ini
    state_fit1
  })
  state_reactive_2 <- reactive({
    state_fit2 <- state
    state_fit2[["P"]] <- parm_reactive$parm_fit2$P.ini
    state_fit2[["Pads.fast"]] <- parm_reactive$parm_fit2$Pads.ini *
      parm_reactive$parm_fit2$x.fast
    state_fit2[["Pads.slow"]] <- parm_reactive$parm_fit2$Pads.ini * 
      (1-parm_reactive$parm_fit2$x.fast)
    state_fit2[["Fe"]] <- parm_reactive$parm_fit2$Fe.ini
    state_fit2
  })
  state_reactive_3 <- reactive({
    state_fit3 <- state
    state_fit3[["P"]] <- parm_reactive$parm_fit3$P.ini
    state_fit3[["Pads.fast"]] <- parm_reactive$parm_fit3$Pads.ini *
      parm_reactive$parm_fit3$x.fast
    state_fit3[["Pads.slow"]] <- parm_reactive$parm_fit3$Pads.ini * 
      (1-parm_reactive$parm_fit3$x.fast)
    state_fit3[["Fe"]] <- parm_reactive$parm_fit3$Fe.ini
    state_fit3
  })
  state_reactive_4 <- reactive({
    state_fit4 <- state
    state_fit4[["P"]] <- parm_reactive$parm_fit4$P.ini
    state_fit4[["Pads.fast"]] <- parm_reactive$parm_fit4$Pads.ini *
      parm_reactive$parm_fit4$x.fast
    state_fit4[["Pads.slow"]] <- parm_reactive$parm_fit4$Pads.ini * 
      (1-parm_reactive$parm_fit4$x.fast)
    state_fit4[["Fe"]] <- parm_reactive$parm_fit4$Fe.ini
    state_fit4
  })
  state_reactive_5 <- reactive({
    state_fit5 <- state
    state_fit5[["P"]] <- parm_reactive$parm_fit5$P.ini
    state_fit5[["Pads.fast"]] <- parm_reactive$parm_fit5$Pads.ini *
      parm_reactive$parm_fit5$x.fast
    state_fit5[["Pads.slow"]] <- parm_reactive$parm_fit5$Pads.ini * 
      (1-parm_reactive$parm_fit5$x.fast)
    state_fit5[["Fe"]] <- parm_reactive$parm_fit5$Fe.ini
    state_fit5
  })
  
  
  # -------------------
  # the 'Plot' tab
  # ------------------- 
  
  output$PlotFeSP <- renderPlot({     # will be visible in the main panel
    
    parm.1 <- parm_reactive$parm_fit1
    state.1 <- state_reactive_1()
    parm.2 <- parm_reactive$parm_fit2
    state.2 <- state_reactive_2()
    parm.3 <- parm_reactive$parm_fit3
    state.3 <- state_reactive_3()
    parm.4 <- parm_reactive$parm_fit4
    state.4 <- state_reactive_4()
    parm.5 <- parm_reactive$parm_fit5
    state.5 <- state_reactive_5()
    
    
    parm_table <- rbind.data.frame("parm.1",parm.1,
                                   "parm.2",parm.2,
                                   "parm.3",parm.3,
                                   "parm.4",parm.4,
                                   "parm.5",parm.5,make.row.names = T)
    
    observeEvent(input$show_parm,{
      output$table_parmfit <- renderDataTable(parm_table)
    })
    
    
    
    # evaluate the model
    out1   <- ode(y=state.1, parms=parm.1, func=FeSPmodel, times=outtimes, method="vode")
    out2   <- ode(y=state.2, parms=parm.2, func=FeSPmodel, times=outtimes, method="vode")
    out3   <- ode(y=state.3, parms=parm.2, func=FeSPmodel, times=outtimes, method="vode")
    out4   <- ode(y=state.4, parms=parm.3, func=FeSPmodel, times=outtimes, method="vode")
    out5   <- ode(y=state.5, parms=parm.4, func=FeSPmodel, times=outtimes, method="vode")
    # outplot <- cbind(out1,out2,out3,out4,out5)
    # parm_table <- cbind.data.frame(out1,out2,out3,out4,out5)
    
    # # choose input experimental data
    data_plotall <- list(
      FeSPads.data1,
      FeSPads.data2,
      FeSPads.data3,
      FeSPads.data4,
      FeSPads.data5,
      FeSPads.data6,
      FeSPads.data7,
      FeSPads.data8,
      FeSPads.data9,
      FeSPads.data10,
      FeSPads.data11,
      FeSPads.data12,
      FeSPads.data13,
      FeSPads.data14,
      FeSPads.data15,
      FeSPads.data16,
      FeSPads.data17,
      FeSPads.data18,
      FeSPads.data19,
      FeSPads.data20,
      FeSPads.data21,
      FeSPads.data22,
      FeSPads.data23,
      FeSPads.data24)
    
    # ----pickerinput for return saved initial fit parameters----------------------------
    
    observeEvent(input$expdataID, {
      if (input$expdataID == "pH7373R1.csv"){
        updateNumericInput(session, "tau",     value = 25)
        updateNumericInput(session, "n",       value = 0.8)
        updateNumericInput(session, "m",       value = 0.5)
        updateNumericInput(session, "Fe.ini",  value = 0)
        updateNumericInput(session, "x.fast",       value = 0.5)
        updateNumericInput(session, "y",       value = 1)
        updateNumericInput(session, "kFe",     value = -1.58)
        updateNumericInput(session, "kPads",    value = -4)
        updateNumericInput(session, "kPads2kdes", value = 0.33)
        updateNumericInput(session, "kS.fast",      value = 0.65)
        updateNumericInput(session, "kS.slow",      value = 0.65)
        updateNumericInput(session, "Sin",     value = 1.6)
        updateNumericInput(session, "P.ini",   value = 0)
        updateNumericInput(session, "Pads.ini", value = 0)
        updateNumericInput(session, "kSads", value = -10)
        updateNumericInput(session, "kSads2kdes", value = 1)
      }
      if (input$expdataID == "pH7373R2.csv"){
        updateNumericInput(session, "tau",     value = 25)
        updateNumericInput(session, "n",       value = 0.8)
        updateNumericInput(session, "m",       value = 0.5)
        updateNumericInput(session, "Fe.ini",  value = 12)
        updateNumericInput(session, "x.fast",       value = 0.5)
        updateNumericInput(session, "y",       value = 1)
        updateNumericInput(session, "kFe",     value = -1.564)
        updateNumericInput(session, "kPads",    value = -4)
        updateNumericInput(session, "kPads2kdes", value = 0.33)
        updateNumericInput(session, "kS.fast",      value = 0.65)
        updateNumericInput(session, "kS.slow",      value = 0.65)
        updateNumericInput(session, "Sin",     value = 1.6)
        updateNumericInput(session, "P.ini",   value = 0)
        updateNumericInput(session, "Pads.ini", value = 0)
        updateNumericInput(session, "kSads", value = -10)
        updateNumericInput(session, "kSads2kdes", value = 1)
      }
      if (input$expdataID == "pH7373R3.csv"){
        updateNumericInput(session, "tau",     value = 25)
        updateNumericInput(session, "n",       value = 0.8)
        updateNumericInput(session, "m",       value = 0.5)
        updateNumericInput(session, "Fe.ini",  value = 12.4)
        updateNumericInput(session, "x.fast",       value = 0.5)
        updateNumericInput(session, "y",       value = 1)
        updateNumericInput(session, "kFe",     value = -1.74)
        updateNumericInput(session, "kPads",    value = -4)
        updateNumericInput(session, "kPads2kdes", value = 0.33)
        updateNumericInput(session, "kS.fast",      value = 0.65)
        updateNumericInput(session, "kS.slow",      value = 0.65)
        updateNumericInput(session, "Sin",     value = 1.6)
        updateNumericInput(session, "P.ini",   value = 0)
        updateNumericInput(session, "Pads.ini", value = 0)
        updateNumericInput(session, "kSads", value = -10)
        updateNumericInput(session, "kSads2kdes", value = 1)
      }
      if (input$expdataID == "pH7373R4.csv"){
        updateNumericInput(session, "tau",     value = 25)
        updateNumericInput(session, "n",       value = 0.8)
        updateNumericInput(session, "m",       value = 0.5)
        updateNumericInput(session, "Fe.ini",  value = 9.4)
        updateNumericInput(session, "x.fast",       value = 0.23)
        updateNumericInput(session, "y",       value = 1)
        updateNumericInput(session, "kFe",     value = -1.62)
        updateNumericInput(session, "kPads",    value = -4)
        updateNumericInput(session, "kPads2kdes", value = 0.33)
        updateNumericInput(session, "kS.fast",      value = -0.02)
        updateNumericInput(session, "kS.slow",      value = -0.94)
        updateNumericInput(session, "Sin",     value = 1.6)
        updateNumericInput(session, "P.ini",   value = 0.01)
        updateNumericInput(session, "Pads.ini", value = 0.09)
        updateNumericInput(session, "kSads", value = -10)
        updateNumericInput(session, "kSads2kdes", value = 1)
      }
      if (input$expdataID == "pH7373R5.csv"){
        updateNumericInput(session, "tau",     value = 25)
        updateNumericInput(session, "n",       value = 0.8)
        updateNumericInput(session, "m",       value = 0.5)
        updateNumericInput(session, "Fe.ini",  value = 11.6)
        updateNumericInput(session, "x.fast",       value = 0.3)
        updateNumericInput(session, "y",       value = 1)
        updateNumericInput(session, "kFe",     value = -1.62)
        updateNumericInput(session, "kPads",    value = -4)
        updateNumericInput(session, "kPads2kdes", value = 0.33)
        updateNumericInput(session, "kS.fast",      value = -0.14)
        updateNumericInput(session, "kS.slow",      value = -0.94)
        updateNumericInput(session, "Sin",     value = 1.6)
        updateNumericInput(session, "P.ini",   value = 0.01)
        updateNumericInput(session, "Pads.ini", value = 0.11)
        updateNumericInput(session, "kSads", value = -10)
        updateNumericInput(session, "kSads2kdes", value = 1)
      }
      
      
    })
    # -----------------------------------------------------
    # function blocked, as this requires larger memory
    # plotsize <- reactiveValues(a = 2,
    #                            b =4)
    # 
    # observeEvent(input$plot_elements, {
    #   if (length(input$plot_elements) <= 4){
    #     plotsize$a = 2 
    #     plotsize$b =2
    #   }
    #   else if(length(input$plot_elements)<= 6) {
    #     plotsize$a = 2
    #     plotsize$b =3
    #   }
    #   else if(length(input$plot_elements)<= 8){
    #     plotsize$a = 2
    #     plotsize$b = 4
    #   }
    #   else if(length(input$plot_elements) == 9){
    #     plotsize$a = 3
    #     plotsize$b = 3
    #   }
    #   else{plotsize$a = 4
    #   plotsize$b = 4}
    # })
    # 
    
    # --------------------------------------------------------------------
    # display model results + exp. data
    plot (out1,out2,out3,out4,out5,
          lwd = 2, las = 1, lty = 1, mfrow=c(2,4),
          obs = data_plotall[as.numeric(input$plot_check)],
          which=input$plot_elements,
          cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
    legend("topright", legend = c("Fit1", "Fit2","Fit3","Fit4","Fit5"),
           cex = 1.5, col = 1:2, lty = 1, bty = "n")  
  })                             # end ouput$plot
  
  
  
  output$verb <- renderText(
    if(input$appendix == TRUE){
      paste("
            Rate expressions:
            R <- kFe * S^m * Fe.ini^(1-n) * Fe^n; 
            Rads <- kPads * P  * Fe.ini^(1-n) * Fe^n;
            kdes <- kPads/kPads2kdes;
            RSads <- kSads * S  * Fe.ini^(1-n) * Fe^n;
            kSdes <- kSads/kSads2kdes;RSdes <- kSdes * Sads;
            Rdes.fast <- kdes * Pads.fast;Rdes.slow <- kdes * Pads.slow;
            RHS.fast <- kS.fast^y * Pads.fast * S^y;
            RHS.slow <- kS.slow^y * Pads.slow * S^y,
            Mass balances:
            dFe.dt <- -R;dS.dt <- -3/2*R + D*Sin - D*S - RSads + RSdes;
            dSads.dt <- RSads - RSdes;
            dP.dt <- -Rads + Rdes.fast + Rdes.slow + RHS.fast + RHS.slow - D*P;
            dPads.fast.dt <- Rads*x.fast - Rdes.fast - RHS.fast;
            dPads.slow.dt <- Rads*(1-x.fast)  - Rdes.slow - RHS.slow")
    }
    
  )    # end of the display of the equations
})     # end of the definition of shinyServer

## To run as ShinyApp:
#shinyApp(ui = UI.FeSP, server = Server.FeSP)


