#install.packages("shiny")
#install.packages("ggplot2")
library(shiny)
library(ggplot2)

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Estimation of SP Protective efficacy"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      tags$h3("Inputs:"),
      textInput("txt1", "Incidence (per person per year):", "10"),
      textInput("txt2", "Frequency of dhps _AKAA:", "0.1"),
      textInput("txt3", "Frequency of dhps _GKAA:", "0.5"),
      textInput("txt4", "Frequency of dhps _GEAA:", "0.3"),
      textInput("txt5", "Frequency of dhps _GEGA:", "0.05"),
      textInput("txt6", "Frequency of other dhps combinations:", "0.05"),
      h3("Note: above frequencies should add up to 1"),
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
       # Output: Histogram ----
      plotOutput(outputId = "SP_PE_plot"),
      # 
      h3("30-day Protective efficacy:"),
      
      verbatimTextOutput("protective_efficacy"),
      tags$head(tags$style(HTML("
                            #protective_efficacy {
                              font-size: 20px;
                            }
                            "))),
      
      h3("Median duration of protection:"),
      verbatimTextOutput("days_med"),
      tags$head(tags$style(HTML("
                            #days_med {
                              font-size: 20px;
                            }
                            ")))
      
      
    )
  )
)

server <- function(input, output) {
    
  output$SP_PE_plot<- renderPlot({
    
  inc_py<- as.numeric(input$txt1)
  
  freq_trip<-as.numeric(input$txt2)
  freq_quadr<-as.numeric(input$txt3)     
  freq_quint<-as.numeric(input$txt4)  
  freq_sext<- as.numeric(input$txt5)   
  freq_other<- as.numeric(input$txt6)
 
  
  reported_dur<-30
  
  inc_pd<-inc_py/365
  
  ############################################# to be fixed
  
  lambda_SP_trip<-59.41
  w_SP_trip<- 8.44
  lambda_SP_quint<-18.43
  w_SP_quint<-2.56
  lambda_SP_quadr<-32.96
  w_SP_quadr<- 4.91
  lambda_SP_sext<-12.86
  w_SP_sext<-3.61
  lambda_SP_other<-23
  w_SP_other<-4.5
  
  
  ###lambda= 30/ (gamma(1+(1/5)))
  ## time step for deterministic model (half a day but could use smaller)
  dt=0.1  ## 
  
  ## create drug protection weibull curve
  time<-seq(from=0,to=63,by=dt)      #### Ran for 60 days, Protective Efficacy will depend on which time interval is being reported.
  
  
  PEdata<-data.frame(time)
  
  PEdata$SP_PE<-NA
  
  p_protect_SP_trip<- exp(-(time/lambda_SP_trip)^w_SP_trip)
  p_protect_SP_quadr<- exp(-(time/lambda_SP_quadr)^w_SP_quadr)
  p_protect_SP_quint<- exp(-(time/lambda_SP_quint)^w_SP_quint)
  p_protect_SP_sext<- exp(-(time/lambda_SP_sext)^w_SP_sext)
  p_protect_SP_other<- exp(-(time/lambda_SP_other)^w_SP_other)
  
  prob_inf<-1-exp(-inc_pd*dt)  ## prob of infection at each time step
  
  ## simulate reinfection in a control vs two chemoprevention groups in a simple difference equation set up.
  
  control<-treated_SP<- treated_SP_I_trip<- treated_SP_I_quadr<-treated_SP_I_quint<-treated_SP_I_sext<-treated_SP_I_other<-vector(length=length(time))   
  
  ##start with everyone being uninfected (at risk) at the first time point
  
  treated_SP_I_trip[1]<-0
  treated_SP_I_quadr[1]<-0
  treated_SP_I_sext[1]<-0
  treated_SP_I_quint[1]<-0
  treated_SP_I_other[1]<-0
  
  treated_SP[1]<-1
  control[1]<-1  
  PEdata$SP_PE[1]<-1
  
  
  for(i in 2:length(time)) {
    control[i]<-control[i-1] - prob_inf*control[i-1]  
    treated_SP_I_quadr[i]<-treated_SP_I_trip[i-1]+ (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  ## proportion of new infections with S in treated group
    treated_SP_I_quadr[i]<-treated_SP_I_quadr[i-1]+ (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  ## proportion of new infections with S in treated group
    treated_SP_I_quint[i]<-treated_SP_I_quint[i-1]+ (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))   ## proportion of new infections with S in treated group
    treated_SP_I_sext[i]<-treated_SP_I_sext[i-1]+ (prob_inf*freq_sext*treated_SP[i-1]*(1-p_protect_SP_sext[i])) 
    treated_SP_I_other[i]<-treated_SP_I_other[i-1]+ (prob_inf*freq_other*treated_SP[i-1]*(1-p_protect_SP_other[i])) 
    PEdata$SP_PE[i]<- (freq_trip*(p_protect_SP_trip[i])) + 
      (freq_quadr*(p_protect_SP_quadr[i]))+
      (freq_quint*(p_protect_SP_quint[i]))+
      (freq_sext*(p_protect_SP_sext[i]))+
      (freq_other*(p_protect_SP_other[i]))
    
    treated_SP[i]<-treated_SP[i-1] -  
      (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  -  
      (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  -  
      (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))  - 
      (prob_inf*freq_sext *treated_SP[i-1]*(1-p_protect_SP_sext[i]))   -
      (prob_inf*freq_other *treated_SP[i-1]*(1-p_protect_SP_other[i]))   
  }
  
  
  ### calculate protective efficacy (will depend on what duration you report for - e.g 30 days)
  
  protective_efficacy<- round((1-((1-treated_SP[length(treated_SP[1:(1+(reported_dur/dt))])])/(1-control[length(control[1:(1+(reported_dur/dt))])]))), digits=5)
  
 days_med<-PEdata$time[which(abs(PEdata$SP_PE-0.5)==min(abs(PEdata$SP_PE-0.5)))]  ### this line will tell you at which day protective efficacy is 50% (ie. median duration of protection)
  
  
 
  
 if (round((freq_trip + freq_quadr +freq_quint+freq_sext+freq_other), digits = 5)==1 ) {
   ggplot(PEdata) + theme(axis.text = element_text(size = 20)) +
      geom_line(aes(x=time, y=SP_PE), linewidth=2, colour="purple") +
      theme_bw() +ylab("Protective Efficacy") +xlab("Time since SP dose") +theme(
        text = element_text(size = 20)
      )

  ##})
  ##print(c("Number of days to reach 50% Protective efficacy:",  days_med))
} 
  })
  
  output$protective_efficacy <- renderText({
    inc_py<- as.numeric(input$txt1)
    
    freq_trip<-as.numeric(input$txt2)
    freq_quadr<-as.numeric(input$txt3)     
    freq_quint<-as.numeric(input$txt4)  
    freq_sext<- as.numeric(input$txt5)   
    freq_other<- as.numeric(input$txt6)
    
    
    reported_dur<-30
    
    inc_pd<-inc_py/365
    
    ############################################# to be fixed
    
    lambda_SP_trip<-59.57659
    w_SP_trip<- 8.435971
    lambda_SP_quint<-18.55328
    w_SP_quint<-2.4840752
    lambda_SP_quadr<-33.05391
    w_SP_quadr<- 4.862126
    lambda_SP_sext<- 12.81186
    w_SP_sext<-3.691953
    lambda_SP_other<-23
    w_SP_other<-4.5
    
    ###lambda= 30/ (gamma(1+(1/5)))
    ## time step for deterministic model (half a day but could use smaller)
    dt=0.1  ## 
    
    ## create drug protection weibull curve
    time<-seq(from=0,to=63,by=dt)      #### Ran for 60 days, Protective Efficacy will depend on which time interval is being reported.
    
    
    PEdata<-data.frame(time)
    
    PEdata$SP_PE<-NA
    
    p_protect_SP_trip<- exp(-(time/lambda_SP_trip)^w_SP_trip)
    p_protect_SP_quadr<- exp(-(time/lambda_SP_quadr)^w_SP_quadr)
    p_protect_SP_quint<- exp(-(time/lambda_SP_quint)^w_SP_quint)
    p_protect_SP_sext<- exp(-(time/lambda_SP_sext)^w_SP_sext)
    p_protect_SP_other<- exp(-(time/lambda_SP_other)^w_SP_other)
    
    prob_inf<-1-exp(-inc_pd*dt)  ## prob of infection at each time step
    
    ## simulate reinfection in a control vs two chemoprevention groups in a simple difference equation set up.
    
    control<-treated_SP<- treated_SP_I_trip<- treated_SP_I_quadr<-treated_SP_I_quint<-treated_SP_I_sext<-treated_SP_I_other<-vector(length=length(time))   
    
    ##start with everyone being uninfected (at risk) at the first time point
    
    treated_SP_I_trip[1]<-0
    treated_SP_I_quadr[1]<-0
    treated_SP_I_sext[1]<-0
    treated_SP_I_quint[1]<-0
    treated_SP_I_other[1]<-0
    
    treated_SP[1]<-1
    control[1]<-1  
    PEdata$SP_PE[1]<-1
    
    
    for(i in 2:length(time)) {
      control[i]<-control[i-1] - prob_inf*control[i-1]  
      treated_SP_I_quadr[i]<-treated_SP_I_trip[i-1]+ (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  ## proportion of new infections with S in treated group
      treated_SP_I_quadr[i]<-treated_SP_I_quadr[i-1]+ (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  ## proportion of new infections with S in treated group
      treated_SP_I_quint[i]<-treated_SP_I_quint[i-1]+ (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))   ## proportion of new infections with S in treated group
      treated_SP_I_sext[i]<-treated_SP_I_sext[i-1]+ (prob_inf*freq_sext*treated_SP[i-1]*(1-p_protect_SP_sext[i])) 
      treated_SP_I_other[i]<-treated_SP_I_other[i-1]+ (prob_inf*freq_other*treated_SP[i-1]*(1-p_protect_SP_other[i])) 
      PEdata$SP_PE[i]<- (freq_trip*(p_protect_SP_trip[i])) + 
        (freq_quadr*(p_protect_SP_quadr[i]))+
        (freq_quint*(p_protect_SP_quint[i]))+
        (freq_sext*(p_protect_SP_sext[i]))+
        (freq_other*(p_protect_SP_other[i]))
      
      treated_SP[i]<-treated_SP[i-1] -  
        (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  -  
        (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  -  
        (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))  - 
        (prob_inf*freq_sext *treated_SP[i-1]*(1-p_protect_SP_sext[i]))   -
        (prob_inf*freq_other *treated_SP[i-1]*(1-p_protect_SP_other[i]))   
    }
    
    
    ### calculate protective efficacy (will depend on what duration you report for - e.g 30 days)
    
    protective_efficacy<- round((1-((1-treated_SP[length(treated_SP[1:(1+(reported_dur/dt))])])/(1-control[length(control[1:(1+(reported_dur/dt))])]))), digits=5)
    
    days_med<-PEdata$time[which(abs(PEdata$SP_PE-0.5)==min(abs(PEdata$SP_PE-0.5)))]  ### this line will tell you at which day protective efficacy is 50% (ie. median duration of protection)
    
    
    if (round((freq_trip + freq_quadr +freq_quint+freq_sext+freq_other), digits = 5)==1 ) {
      paste( protective_efficacy )
    } else {paste("Frequencies should add up to 1")}
  })
  
  output$days_med <- renderText({
    inc_py<- as.numeric(input$txt1)
    
    freq_trip<-as.numeric(input$txt2)
    freq_quadr<-as.numeric(input$txt3)     
    freq_quint<-as.numeric(input$txt4)  
    freq_sext<- as.numeric(input$txt5)   
    freq_other<- as.numeric(input$txt6)
    
    
    reported_dur<-30
    
    inc_pd<-inc_py/365
    
    ############################################# to be fixed
    
    lambda_SP_trip<-59.57659
    w_SP_trip<- 8.435971
    lambda_SP_quint<-18.55328
    w_SP_quint<-2.4840752
    lambda_SP_quadr<-33.05391
    w_SP_quadr<- 4.862126
    lambda_SP_sext<- 12.81186
    w_SP_sext<-3.691953
    lambda_SP_other<-23
    w_SP_other<-4.5
    
    
    ###lambda= 30/ (gamma(1+(1/5)))
    ## time step for deterministic model (half a day but could use smaller)
    dt=0.1  ## 
    
    ## create drug protection weibull curve
    time<-seq(from=0,to=63,by=dt)      #### Ran for 60 days, Protective Efficacy will depend on which time interval is being reported.
    
    
    PEdata<-data.frame(time)
    
    PEdata$SP_PE<-NA
    
    p_protect_SP_trip<- exp(-(time/lambda_SP_trip)^w_SP_trip)
    p_protect_SP_quadr<- exp(-(time/lambda_SP_quadr)^w_SP_quadr)
    p_protect_SP_quint<- exp(-(time/lambda_SP_quint)^w_SP_quint)
    p_protect_SP_sext<- exp(-(time/lambda_SP_sext)^w_SP_sext)
    p_protect_SP_other<- exp(-(time/lambda_SP_other)^w_SP_other)
    
    prob_inf<-1-exp(-inc_pd*dt)  ## prob of infection at each time step
    
    ## simulate reinfection in a control vs two chemoprevention groups in a simple difference equation set up.
    
    control<-treated_SP<- treated_SP_I_trip<- treated_SP_I_quadr<-treated_SP_I_quint<-treated_SP_I_sext<-treated_SP_I_other<-vector(length=length(time))   
    
    ##start with everyone being uninfected (at risk) at the first time point
    
    treated_SP_I_trip[1]<-0
    treated_SP_I_quadr[1]<-0
    treated_SP_I_sext[1]<-0
    treated_SP_I_quint[1]<-0
    treated_SP_I_other[1]<-0
    
    treated_SP[1]<-1
    control[1]<-1  
    PEdata$SP_PE[1]<-1
    
    
    for(i in 2:length(time)) {
      control[i]<-control[i-1] - prob_inf*control[i-1]  
      treated_SP_I_quadr[i]<-treated_SP_I_trip[i-1]+ (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  ## proportion of new infections with S in treated group
      treated_SP_I_quadr[i]<-treated_SP_I_quadr[i-1]+ (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  ## proportion of new infections with S in treated group
      treated_SP_I_quint[i]<-treated_SP_I_quint[i-1]+ (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))   ## proportion of new infections with S in treated group
      treated_SP_I_sext[i]<-treated_SP_I_sext[i-1]+ (prob_inf*freq_sext*treated_SP[i-1]*(1-p_protect_SP_sext[i])) 
      treated_SP_I_other[i]<-treated_SP_I_other[i-1]+ (prob_inf*freq_other*treated_SP[i-1]*(1-p_protect_SP_other[i])) 
      PEdata$SP_PE[i]<- (freq_trip*(p_protect_SP_trip[i])) + 
        (freq_quadr*(p_protect_SP_quadr[i]))+
        (freq_quint*(p_protect_SP_quint[i]))+
        (freq_sext*(p_protect_SP_sext[i]))+
        (freq_other*(p_protect_SP_other[i]))
      
      treated_SP[i]<-treated_SP[i-1] -  
        (prob_inf*freq_trip*treated_SP[i-1]*(1-p_protect_SP_trip[i]))  -  
        (prob_inf*freq_quadr*treated_SP[i-1]*(1-p_protect_SP_quadr[i]))  -  
        (prob_inf*freq_quint*treated_SP[i-1]*(1-p_protect_SP_quint[i]))  - 
        (prob_inf*freq_sext *treated_SP[i-1]*(1-p_protect_SP_sext[i]))   -
        (prob_inf*freq_other *treated_SP[i-1]*(1-p_protect_SP_other[i]))   
    }
    
    
    ### calculate protective efficacy (will depend on what duration you report for - e.g 30 days)
    
    protective_efficacy<- round((1-((1-treated_SP[length(treated_SP[1:(1+(reported_dur/dt))])])/(1-control[length(control[1:(1+(reported_dur/dt))])]))), digits=5)
    days_med<-PEdata$time[which(abs(PEdata$SP_PE-0.5)==min(abs(PEdata$SP_PE-0.5)))]  ### this line will tell you at which day protective efficacy is 50% (ie. median duration of protection)
    
    
    if (round((freq_trip + freq_quadr +freq_quint+freq_sext+freq_other), digits = 5)==1 ) {
      
    paste(days_med )
    } else {paste("Frequencies should add up to 1")}
  })
 
  
   
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
