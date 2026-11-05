library(shiny)
library(bslib)
library(rhandsontable)
library(DT)
library(dplyr)

#############################################################################
# Define functions for calculating sample size
#############################################################################

calculate_sample_size_stage1 <- function(design,
                                         planned_aim,
                                         stdeffsize,
                                         alpha_test,
                                         target_power,
                                         response_probability,
                                         has_baseline,
                                         rho,
                                         prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1 - rho^2)*(1/(1 - prob_missing))
  }else{
    rho <- NA
    q <- 4
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1/(1 - prob_missing))
  }
  
  value_sample_size <- ceiling(value_sample_size)
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = value_sample_size,
                                     Power = target_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_sample_size = value_sample_size)
  
  return(list_output)
}

calculate_sample_size_design1_stage2nonresp <- function(design,
                                                        planned_aim,
                                                        stdeffsize,
                                                        alpha_test,
                                                        target_power,
                                                        response_probability,
                                                        has_baseline,
                                                        rho,
                                                        prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1 - rho^2)*(1/(1 - prob_missing))
  }else{
    rho <- NA
    q <- 4/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1/(1 - prob_missing))
  }
  
  value_sample_size <- ceiling(value_sample_size)
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = value_sample_size,
                                     Power = target_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_sample_size = value_sample_size)
  
  return(list_output)
}

calculate_sample_size_design1_stage2resp <- function(design,
                                                     planned_aim,
                                                     stdeffsize,
                                                     alpha_test,
                                                     target_power,
                                                     response_probability,
                                                     has_baseline,
                                                     rho,
                                                     prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/response_probability
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1 - rho^2)*(1/(1 - prob_missing))
  }else{
    rho <- NA
    q <- 4/response_probability
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1/(1 - prob_missing))
  }
  
  value_sample_size <- ceiling(value_sample_size)
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = value_sample_size,
                                     Power = target_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_sample_size = value_sample_size)
  
  return(list_output)
}

calculate_sample_size_design2_stage2nonresp <- function(design,
                                                        planned_aim,
                                                        stdeffsize,
                                                        alpha_test,
                                                        target_power,
                                                        response_probability,
                                                        has_baseline,
                                                        rho,
                                                        prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1 - rho^2)*(1/(1 - prob_missing))
  }else{
    rho <- NA
    q <- 4/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1/(1 - prob_missing))
  }
  
  value_sample_size <- ceiling(value_sample_size)
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = value_sample_size,
                                     Power = target_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_sample_size = value_sample_size)
  
  return(list_output)
}

calculate_sample_size_design3_stage2nonresp <- function(design,
                                                        planned_aim,
                                                        stdeffsize,
                                                        alpha_test,
                                                        target_power,
                                                        response_probability,
                                                        has_baseline,
                                                        rho,
                                                        prob_missing)
{
  if(has_baseline == "yes"){
    q <- 8/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1 - rho^2)*(1/(1 - prob_missing))
  }else{
    rho <- NA
    q <- 8/(1 - response_probability)
    const <- (qnorm(1 - alpha_test/2) + qnorm(target_power))/stdeffsize
    value_sample_size <- (const^2)*q*(1/(1 - prob_missing))
  }
  
  value_sample_size <- ceiling(value_sample_size)
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = value_sample_size,
                                     Power = target_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_sample_size = value_sample_size)
  
  return(list_output)
}

#############################################################################
# Define functions for calculating power
#############################################################################

calculate_power_stage1 <- function(design,
                                   planned_aim,
                                   stdeffsize,
                                   alpha_test,
                                   target_sample_size,
                                   response_probability,
                                   has_baseline,
                                   rho,
                                   prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q * (1 - rho^2))) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }else{
    rho <- NA
    q <- 4
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q)) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = target_sample_size,
                                     Power = value_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_power = value_power)
  
  return(list_output)
}

calculate_power_design1_stage2nonresp <- function(design,
                                                  planned_aim,
                                                  stdeffsize,
                                                  alpha_test,
                                                  target_sample_size,
                                                  response_probability,
                                                  has_baseline,
                                                  rho,
                                                  prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q * (1 - rho^2))) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }else{
    rho <- NA
    q <- 4/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q)) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = target_sample_size,
                                     Power = value_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_power = value_power)
  
  return(list_output)
}

calculate_power_design1_stage2resp <- function(design,
                                               planned_aim,
                                               stdeffsize,
                                               alpha_test,
                                               target_sample_size,
                                               response_probability,
                                               has_baseline,
                                               rho,
                                               prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/response_probability
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q * (1 - rho^2))) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }else{
    rho <- NA
    q <- 4/response_probability
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q)) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = target_sample_size,
                                     Power = value_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_power = value_power)
  
  return(list_output)
}

calculate_power_design2_stage2nonresp <- function(design,
                                                  planned_aim,
                                                  stdeffsize,
                                                  alpha_test,
                                                  target_sample_size,
                                                  response_probability,
                                                  has_baseline,
                                                  rho,
                                                  prob_missing)
{
  if(has_baseline == "yes"){
    q <- 4/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q * (1 - rho^2))) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }else{
    rho <- NA
    q <- 4/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q)) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = target_sample_size,
                                     Power = value_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_power = value_power)
  
  return(list_output)
}

calculate_power_design3_stage2nonresp <- function(design,
                                                  planned_aim,
                                                  stdeffsize,
                                                  alpha_test,
                                                  target_sample_size,
                                                  response_probability,
                                                  has_baseline,
                                                  rho,
                                                  prob_missing)
{
  if(has_baseline == "yes"){
    q <- 8/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q * (1 - rho^2))) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }else{
    rho <- NA
    q <- 8/(1 - response_probability)
    const <- stdeffsize*sqrt((target_sample_size*(1/(1 - prob_missing))) / (q)) - qnorm(1 - alpha_test/2)
    value_power <- pnorm(const)
  }
  
  # We save inputs and output to facilitate checking whether the reactivity of
  # the app is functioning as expected
  dat_input_and_output <- data.frame(SMARTDesign = design,
                                     Aim = planned_aim,
                                     SampleSize = target_sample_size,
                                     Power = value_power,
                                     CohensD = stdeffsize,
                                     TypeOneErrorRate = alpha_test,
                                     ResponseProbability = response_probability,
                                     WhetherHasBaseline = has_baseline,
                                     PrePostCorr = rho,
                                     ProbabilityMissing = prob_missing)
  
  list_output <- list(dat_input_and_output = dat_input_and_output, value_power = value_power)
  
  return(list_output)
}

#############################################################################
# Define user interface
#############################################################################

app_theme <- bs_theme(version = 5)

app_theme <- bs_add_rules(
  app_theme,
  "
  .card .card-header {
    background-color: #e2e8f0;
    color: #1e293b;
    font-size: 1.2rem;
    font-weight: 600;
  }

  .card .card-body {
    background-color: #f8fafc;
    color: #1e293b;
    font-size: 0.90rem;
  }
  
  /* Colors for one specific card */
  #result_card .card-header {
    background-color: #00274c;
    color: #ffffff;
  }
  
  #result_card .card-body {
  background-color: #fff8dc;
  color: #212529;
  }
  
  .field-note {
  color: #52647A;
  font-size: 0.9rem;
  line-height: 1.4;
  }
  "
)

input_cards <- list(
  card(card_header(tagList(icon(lib="font-awesome", name="dice"), "Type of SMART and Aim")), 
       h6("Type of SMART"),
       radioButtons(
         inputId = "design",
         label = NULL,
         choices = c("Rerandomize everybody (Design I)" = "design1", 
                     "Rerandomize non-responders only (Prototypical; Design II)" = "design2",
                     "Rerandomize arm 1 non-responders only (Design III)" = "design3"),
         selected = "design2",
         width = '100%'),
       
       h6("Primary Research Aim"),
       tabs_design_aims <- tabsetPanel(
         id = "tabs_design_aims",
         type = "hidden",
         tabPanel("design1",
                  radioButtons(
                    inputId = "design1aims",
                    label = NULL,
                    choices = c("Compare first-stage intervention options" = "stage1", 
                                "Compare second-stage intervention options among responders" = "stage2resp",
                                "Compare second-stage intervention options among non-responders" = "stage2nonresp"),
                    selected = "stage1",
                    width = '100%'
                  )
         ),
         tabPanel("design2",
                  radioButtons(
                    inputId = "design2aims",
                    label = NULL,
                    choices = c("Compare first-stage intervention options" = "stage1", 
                                "Compare second-stage intervention options among non-responders" = "stage2nonresp"),
                    selected = "stage1",
                    width = '100%'
                  )
         ),
         tabPanel("design3",
                  radioButtons(
                    inputId = "design3aims",
                    label = NULL,
                    choices = c("Compare first-stage intervention options" = "stage1", 
                                "Compare second-stage intervention options among arm 1 non-responders" = "stage2nonresp"),
                    selected = "stage1",
                    width = '100%'
                  )
         )
         )
       ), ##### End of card
  
  card(card_header("What You Want to Calculate"),
       card_body(fillable = FALSE,
                 radioButtons(
                   inputId = "need",
                   label = NULL,
                   choices = c("Find sample size for a given power" = "sampsize", 
                               "Find power for a given sample size" = "power"),
                   selected = "sampsize"),
                 
                 tabs_need <- tabsetPanel(
                   id = "tabs_need",
                   type = "hidden",
                   tabPanel("sampsize",
                            h6("Target power"),
                            rHandsontableOutput("target_power"),
                            tags$span(class = "field-note", "Corresponds to 80% target power")
                   ),
                   tabPanel("power",
                            h6("Target total enrolled participants"),
                            rHandsontableOutput("target_sample_size")
                   )))
       ), #### End of card
  
  card(card_header("Standardized effect size (Cohen's d)"), 
       card_body(fillable = FALSE,
                 rHandsontableOutput("stdeffsize"),
                 tags$span(class = "field-note", "Corresponds to the anticipated standardized difference in end-of-study means between the two options under comparison"))
       ), #### End of card
  
  card(
    card_header("Type I error rate"),
    card_body(fillable = FALSE,
              rHandsontableOutput("alpha_test"),
              tags$span(class = "field-note", textOutput("alpha_test_note")))
  ), #### End of card
  
  card(card_header("Proportion Responders"), 
       card_body(fillable = FALSE,
                 rHandsontableOutput("response_probability"),
                 tags$span(class = "field-note", textOutput("response_probability_note")))
  ), #### End of card
  
  card(card_header("Attrition Proportion"), 
       card_body(fillable = FALSE,
                 rHandsontableOutput("prob_missing"),
                 tags$span(class = "field-note", textOutput("prob_missing_note")))
  ), #### End of card
  
  card(card_header("Baseline Measurement"),
       radioButtons(
         inputId = "has_baseline",
         label = NULL,
         choices = c("None" = "no",
                     "As covariate" = "yes"),
         selected = "yes"
       ),
       
       tabs_has_baseline <- tabsetPanel(
         id = "tabs_has_baseline",
         type = "hidden",
         tabPanel("yes",
                  fluidRow(
                    column(width=12,
                           h6("Correlation between baseline and end-of-study measurements"),
                           rHandsontableOutput("rho"),
                           tags$span(class = "field-note", "Corresponds to a Pearson correlation; useful for potentially improving statistical efficiency by adjusting for a baseline measurement"))
                  )),
         tabPanel("no",
                  fluidRow(
                    column(width=12,
                           h4(""))
                  )
         )
       )
       )  #### End of card
)

output_cards <- list(
  card(card_header(tagList(icon(lib="font-awesome", name="firstdraft"), "SMART Schematic (Scroll Down to View)")), 
       full_screen = TRUE, 
       
       tabs_smart_diagram <- tabsetPanel(
         id = "tabs_smart_diagram",
         type = "hidden",
         tabPanel("design1",
                  tags$span(class = "field-note", "In a Design I SMART, responders and non-responders are re-randomized, leading to eight experimental conditions C-J in the figure below. In this app, we assume that at each stage of the SMART, randomized participants have equal probability of being assigned to each intervention option."),
                  card_body(
                    tags$img(
                      src = file.path("images", "design1.png"),
                      alt = "Design I SMART",
                      style = "
                              width: 100%;
                              height: 100%;
                              min-height: 0;
                              object-fit: contain;
                              "
                      )
                    )
                  #tags$img(src = file.path("images", "design1.png"), width = 400, height = 352, units = "inches", alt = "Design I SMART")
                  #tags$img(src = file.path("images", "design1.png"), width = 612, height = 540, units = "inches", alt = "Design I SMART")
         ),
         tabPanel("design2",
                  tags$span(class = "field-note", "In a Design II (Prototypical) SMART, non-responders are re-randomized, leading to six experimental conditions C-H in the figure below. In this app, we assume that at each stage of the SMART, randomized participants have equal probability of being assigned to each intervention option."),
                  card_body(
                    tags$img(
                      src = file.path("images", "design2.png"),
                      alt = "Design II SMART",
                      style = "
                              width: 100%;
                              height: 100%;
                              min-height: 0;
                              object-fit: contain;
                              "
                    )
                  )
                  #tags$img(src = file.path("images", "design2.png"), width = 400, height = 328, units = "px", alt = "Design II (Prototypical) SMART")
                  #tags$img(src = file.path("images", "design2.png"), width = 620, height = 514, units = "px", alt = "Design II (Prototypical) SMART")
         ),
         tabPanel("design3",
                  tags$span(class = "field-note", "In a Design III SMART, only arm 1 non-responders are re-randomized, leading to five experimental conditions C-G in the figure below. In this app, we assume that at each stage of the SMART, randomized participants have equal probability of being assigned to each intervention option."),
                  card_body(
                    tags$img(
                      src = file.path("images", "design3.png"),
                      alt = "Design III SMART",
                      style = "
                              width: 100%;
                              height: 100%;
                              min-height: 0;
                              object-fit: contain;
                              "
                    )
                  )
                  #tags$img(src = file.path("images", "design3.png"), width = 400, height = 300, units = "inches", alt = "Design III SMART")
                  #tags$img(src = file.path("images", "design3.png"), width = 680, height = 514, units = "inches", alt = "Design III SMART")
         )
       )
       ), #### End of card
  
  card(card_header(tagList(icon(lib="font-awesome", name="trophy"), "Current Estimate")), 
       id = "result_card",
       tabs_what_output <- tabsetPanel(
         id = "tabs_what_output",
         type = "hidden",
         tabPanel("power",
                  h6("Power that can be attained is estimated at"),
                  htmlOutput("value_power"),
                  h6("based on the selected design and current assumptions."),
                  tags$span(
                    class = "field-note",
                    textOutput("power_assumptions")
                  )
         ),
         tabPanel("sampsize",
                  h6("The required total enrolled participants is estimated at"),
                  htmlOutput("value_sample_size"),
                  h6("based on the selected design and current assumptions."),
                  tags$span(
                    class = "field-note",
                    textOutput("sampsize_assumptions")
                  )
         )
       )
       ),

  
  card(card_header(tagList(icon(lib="font-awesome", name="clock-rotate-left"), "History")),
       card_body(
         fillable = FALSE,
         tags$span(class = "field-note", "Hover over the table and slide it to the right to view all columns. Enter any text or numeric value in the search box to only view the rows containing those values."),
         DTOutput("history_table"),
         downloadButton("download_history", "Download Calculation History (.csv)")
       )
       ),
  
  card(card_header(tagList(icon(lib="font-awesome", name="circle-info"), "Acknowledgement and Related Resources")),
       p(strong("About this R Shiny app")),
       p(strong("By Jamie R. T. Yap and John J. Dziak")),
       p(class = "text-muted",
         "(c) 2026, d3Center, University of Michigan"),
       p(strong("Funding: This work was supported by NIH awards P50 DA 054039 and NIH R01 DA 039901 from the National Institute on Drug Abuse.")),
       
       # target = "_blank": Forces the link to open in a new browser tab so users don't accidentally navigate away from our app
       p("For a version for pilot studies see",
         tags$a("here", href = "https://d3center.shinyapps.io/SMART_Pilot_Sample_Size/", target = "_blank"),
         "(by Daniel Almirall and team)"),
       p("For binary outcomes or for more options see",
         tags$a("here", href = "https://umich-biostat-smartsize.share.connect.posit.cloud/", target = "_blank"),
         "(by Kelley Kidwell and team)"),
       p("For tradeoffs see",
         tags$a("here", href = "https://d3center.shinyapps.io/smart-tradeoffs-01/", target = "_blank"),
         "and",
         tags$a("here", href = "https://d3center.shinyapps.io/smart-tradeoffs-02/", target = "_blank"),
         "(by Jamie Yap and team)"),
       p("For other randomized trial planning resources see",
         tags$a("here", href = "https://d3c.isr.umich.edu/code-library/", target = "_blank"), "(d3Center code library)")
       )
)

ui <- page_fluid(
  titlePanel("Simplified Sample Size and Power Calculator for SMARTs"),
  theme = app_theme,
  
  tags$p(
    "Adjust settings in any order to explore sample size and power. ",
    "The current estimate updates automatically.",
    class = "workflow-note"
  ),
  
  layout_columns(col_widths = c(9,3),
                 card(card_header(tagList(icon(lib="font-awesome", name="wrench"), "Settings")),
                      layout_columns(col_widths = c(4,4,4),
                        layout_columns(col_widths = c(12,12), input_cards[[1]], input_cards[[2]]),
                        layout_columns(col_widths = c(12,12,12), input_cards[[3]], input_cards[[4]], input_cards[[5]]),
                        layout_columns(col_widths = c(12,12), input_cards[[7]], input_cards[[6]])
                      )
                      ),
                 layout_columns(output_cards[[2]])
                 ),
  layout_columns(output_cards[[1]], output_cards[[3]]),
  layout_columns(output_cards[[4]])
)

#############################################################################
# Define server
#############################################################################

server <- function(input, output, session) {
  #############################################################################
  # Hide or reveal research aims depending on type of SMART
  #############################################################################
  observeEvent(input$design, {
    # The list of aims displayed depend on the type of SMART
    updateTabsetPanel(inputId = "tabs_design_aims", selected = input$design)
  })
  
  observeEvent(input$need, {
    # The user is asked to provide target power if they say that want to plan 
    # sample size. On the other hand, the user is asked to provide target
    # sample size if they say that they want to plan power.
    # Note that the user is always asked to type-1 error rate.
    updateTabsetPanel(inputId = "tabs_need", selected = input$need)
  })
  
  observeEvent(input$has_baseline, {
    # Whether the user is asked to provide the correlation between the baseline
    # measurement and the measurement at follow-up depends on whether they 
    # reported whether a baseline covariate will be accounted for in sample
    # size/power planning
    updateTabsetPanel(inputId = "tabs_has_baseline", selected = input$has_baseline)
  })
  
  observeEvent(input$need, {
    # The right panel will show minimum attainable power if the user says they 
    # want to plan power. On the other hand, the right panel will show minimum
    # required sample size if the user says they want to plan sample size
    updateTabsetPanel(inputId = "tabs_what_output", selected = input$need)
  })
  
  observeEvent(input$design, {
    # The diagram displayed depend on the type of SMART
    updateTabsetPanel(inputId = "tabs_smart_diagram", selected = input$design)
  })
  
  #############################################################################
  # Create hands on table for standardized effect size (Cohen's d)
  #############################################################################
  output$stdeffsize <- renderRHandsontable({
    rhandsontable(data = data.frame(stdeffsize = 0.50), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0001,          # Minimum allowed value (inclusive)
        max = 0.9999,          # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for Type-1 Error Rate
  #############################################################################
  output$alpha_test <- renderRHandsontable({
    rhandsontable(data = data.frame(alpha_test = 0.05), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0001,          # Minimum allowed value (inclusive)
        max = 0.50,            # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for target power
  #############################################################################
  output$target_power <- renderRHandsontable({
    rhandsontable(data = data.frame(target_power = 0.80), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0001,          # Minimum allowed value (inclusive)
        max = 0.9999,          # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for target sample size
  #############################################################################
  output$target_sample_size <- renderRHandsontable({
    rhandsontable(data = data.frame(target_sample_size = as.integer(200)), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      # Note that this code will still allow non-integer (i.e., decimal) values
      # The answer is not "wrong" per se and the app will not grind to a halt
      # if a non-integer value is provided. Rather than include a validation
      # check for integer/non-integer values, we simply advise the user to
      # provide an integer valued target sample size.
      hot_validate_numeric(
        cols = 1,
        min = 8,          # Minimum allowed value (inclusive)
        max = 1000,       # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for response probability
  #############################################################################
  output$response_probability <- renderRHandsontable({
    rhandsontable(data = data.frame(response_probability = 0.40), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0001,          # Minimum allowed value (inclusive)
        max = 0.9999,          # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for proportion of participants with missing outcome 
  # data
  #############################################################################
  output$prob_missing <- renderRHandsontable({
    rhandsontable(data = data.frame(prob_missing = 0.0000), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0000,          # Minimum allowed value (inclusive)
        max = 1.0000,          # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # Create hands on table for correlation between baseline measurement and
  # end-of-study measurement
  #############################################################################
  output$rho <- renderRHandsontable({
    rhandsontable(data = data.frame(rho = 0.30), 
                  stretchH = "all", height = 50,
                  width = 150,
                  useTypes = TRUE, 
                  digits = 3, 
                  colHeaders = FALSE,
                  rowHeaders = FALSE) %>%
      hot_validate_numeric(
        cols = 1,
        min = 0.0001,          # Minimum allowed value (inclusive)
        max = 0.9999,          # Maximum allowed value (inclusive)
        allowInvalid = FALSE
      )
  })
  
  #############################################################################
  # We lay out the logic for calculating minimum required sample size
  # that will be shown as an output in the right panel
  #############################################################################
  
  react_calculate_sample_size <- reactive({
    # Ensure inputs are not NULL before processing (important for startup)
    req(input$design, 
        input$design1aims, input$design2aims, input$design3aims,
        input$stdeffsize, input$alpha_test, input$target_power, input$response_probability, input$has_baseline)
    
    # Convert handsontable input data (JSON format) back to an R data frame
    dat_stdeffsize <- hot_to_r(input$stdeffsize)
    dat_alpha_test <- hot_to_r(input$alpha_test)
    dat_target_power <- hot_to_r(input$target_power)
    dat_response_probability <- hot_to_r(input$response_probability)
    dat_rho <- hot_to_r(input$rho)
    dat_prob_missing <- hot_to_r(input$prob_missing)
    
    if(input$design == "design1"){
      if(input$design1aims == "stage1"){
        # Sample size for comparing first-stage options
        list_output <- calculate_sample_size_stage1(design = input$design,
                                                    planned_aim = input$design1aims,
                                                    has_baseline = input$has_baseline, 
                                                    stdeffsize = dat_stdeffsize$stdeffsize,
                                                    alpha_test = dat_alpha_test$alpha_test,
                                                    target_power = dat_target_power$target_power,
                                                    response_probability = dat_response_probability$response_probability,
                                                    rho = dat_rho$rho,
                                                    prob_missing = dat_prob_missing$prob_missing)
      }else if(input$design1aims == "stage2nonresp"){
        # Sample size for comparing second-stage options among non-responders
        list_output <- calculate_sample_size_design1_stage2nonresp(design = input$design,
                                                                   planned_aim = input$design1aims,
                                                                   has_baseline = input$has_baseline, 
                                                                   stdeffsize = dat_stdeffsize$stdeffsize,
                                                                   alpha_test = dat_alpha_test$alpha_test,
                                                                   target_power = dat_target_power$target_power,
                                                                   response_probability = dat_response_probability$response_probability,
                                                                   rho = dat_rho$rho,
                                                                   prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Sample size for comparing second-stage options among responders
        list_output <- calculate_sample_size_design1_stage2resp(design = input$design,
                                                                planned_aim = input$design1aims,
                                                                has_baseline = input$has_baseline, 
                                                                stdeffsize = dat_stdeffsize$stdeffsize,
                                                                alpha_test = dat_alpha_test$alpha_test,
                                                                target_power = dat_target_power$target_power,
                                                                response_probability = dat_response_probability$response_probability,
                                                                rho = dat_rho$rho,
                                                                prob_missing = dat_prob_missing$prob_missing)
      }
    }else if(input$design == "design2"){
      if(input$design2aims == "stage1"){
        # Sample size for comparing first-stage options
        list_output <- calculate_sample_size_stage1(design = input$design,
                                                    planned_aim = input$design2aims,
                                                    has_baseline = input$has_baseline, 
                                                    stdeffsize = dat_stdeffsize$stdeffsize,
                                                    alpha_test = dat_alpha_test$alpha_test,
                                                    target_power = dat_target_power$target_power,
                                                    response_probability = dat_response_probability$response_probability,
                                                    rho = dat_rho$rho,
                                                    prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Sample size for comparing second-stage options among non-responders
        list_output <- calculate_sample_size_design2_stage2nonresp(design = input$design,
                                                                   planned_aim = input$design2aims,
                                                                   has_baseline = input$has_baseline, 
                                                                   stdeffsize = dat_stdeffsize$stdeffsize,
                                                                   alpha_test = dat_alpha_test$alpha_test,
                                                                   target_power = dat_target_power$target_power,
                                                                   response_probability = dat_response_probability$response_probability,
                                                                   rho = dat_rho$rho,
                                                                   prob_missing = dat_prob_missing$prob_missing)
      }
    }else{
      if(input$design3aims == "stage1"){
        # Sample size for comparing first-stage options
        list_output <- calculate_sample_size_stage1(design = input$design,
                                                    planned_aim = input$design3aims,
                                                    has_baseline = input$has_baseline, 
                                                    stdeffsize = dat_stdeffsize$stdeffsize,
                                                    alpha_test = dat_alpha_test$alpha_test,
                                                    target_power = dat_target_power$target_power,
                                                    response_probability = dat_response_probability$response_probability,
                                                    rho = dat_rho$rho,
                                                    prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Sample size for comparing second-stage options among arm 1 non-responders
        list_output <- calculate_sample_size_design3_stage2nonresp(design = input$design,
                                                                   planned_aim = input$design3aims,
                                                                   has_baseline = input$has_baseline,
                                                                   stdeffsize = dat_stdeffsize$stdeffsize,
                                                                   alpha_test = dat_alpha_test$alpha_test,
                                                                   target_power = dat_target_power$target_power,
                                                                   response_probability = dat_response_probability$response_probability,
                                                                   rho = dat_rho$rho,
                                                                   prob_missing = dat_prob_missing$prob_missing)
      }
    }
    
    return(list_output)
  })
  
  #############################################################################
  # This is the resulting sample size given all the current inputs
  #############################################################################
  output$value_sample_size <- renderUI({
    list_output <- react_calculate_sample_size()
    tags$p(sprintf("%.0f", list_output$value_sample_size), style = "font-size: 70px; font-weight: 1000;")
  })
  
  output$sampsize_assumptions <- renderText({
    list_output <- react_calculate_sample_size()
    
    smart_design_text <- case_when(
      list_output$dat_input_and_output$SMARTDesign == "design1" ~ "Design I",
      list_output$dat_input_and_output$SMARTDesign == "design2" ~ "Design II",
      list_output$dat_input_and_output$SMARTDesign == "design3" ~ "Design III",
      .default = NULL
    )
    
    aim_text <- case_when(
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage2resp" ~ "Compare second-stage intervention options among responders",
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among non-responders",
      list_output$dat_input_and_output$SMARTDesign == "design2" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design2" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among non-responders",
      list_output$dat_input_and_output$SMARTDesign == "design3" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design3" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among arm 1 non-responders",
      .default = NULL
    )
    
    paste(paste("SMART Design: ", smart_design_text, sep = ""),
          paste("Aim: ", aim_text, sep = ""),
          paste("Target Power: ", list_output$dat_input_and_output$Power, sep = ""),
          paste("Effect Size (Cohen's d): ", list_output$dat_input_and_output$CohensD, sep = ""),
          paste("Type I Error Rate: ", round(list_output$dat_input_and_output$TypeOneErrorRate*100, 0), "%", sep = ""),
          paste("Proportion Responders: ", round(list_output$dat_input_and_output$ResponseProbability*100, 0), "%", sep = ""),
          paste("Attrition Proportion: ", round(list_output$dat_input_and_output$ProbabilityMissing*100, 0), "%", sep = ""),
          paste("Correlation between Baseline and End-of-Study Measurements: ", ifelse(is.na(list_output$dat_input_and_output$PrePostCorr), "Not Applicable", list_output$dat_input_and_output$PrePostCorr), sep = ""),
          sep = ", ")
  })
  
  #############################################################################
  # We lay out the logic for calculating minimum attainable power 
  # that will be shown as an output in the right panel
  #############################################################################
  
  react_calculate_power <- reactive({
    # Ensure inputs are not NULL before processing (important for startup)
    req(input$design, 
        input$design1aims, input$design2aims, input$design3aims,
        input$stdeffsize, input$alpha_test, input$target_sample_size, input$response_probability, input$has_baseline)
    
    # Convert handsontable input data (JSON format) back to an R data frame
    dat_stdeffsize <- hot_to_r(input$stdeffsize)
    dat_alpha_test <- hot_to_r(input$alpha_test)
    dat_target_sample_size <- hot_to_r(input$target_sample_size)
    dat_response_probability <- hot_to_r(input$response_probability)
    dat_rho <- hot_to_r(input$rho)
    dat_prob_missing <- hot_to_r(input$prob_missing)
    
    if(input$design == "design1"){
      if(input$design1aims == "stage1"){
        # Power for comparing first-stage options
        list_output <- calculate_power_stage1(design = input$design,
                                              planned_aim = input$design1aims,
                                              has_baseline = input$has_baseline, 
                                              stdeffsize = dat_stdeffsize$stdeffsize,
                                              alpha_test = dat_alpha_test$alpha_test,
                                              target_sample_size = dat_target_sample_size$target_sample_size,
                                              response_probability = dat_response_probability$response_probability,
                                              rho = dat_rho$rho,
                                              prob_missing = dat_prob_missing$prob_missing)
      }else if(input$design1aims == "stage2nonresp"){
        # Power for comparing second-stage options among non-responders
        list_output <- calculate_power_design1_stage2nonresp(design = input$design,
                                                             planned_aim = input$design1aims,
                                                             has_baseline = input$has_baseline, 
                                                             stdeffsize = dat_stdeffsize$stdeffsize,
                                                             alpha_test = dat_alpha_test$alpha_test,
                                                             target_sample_size = dat_target_sample_size$target_sample_size,
                                                             response_probability = dat_response_probability$response_probability,
                                                             rho = dat_rho$rho,
                                                             prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Power for comparing second-stage options among responders
        list_output <- calculate_power_design1_stage2resp(design = input$design,
                                                          planned_aim = input$design1aims,
                                                          has_baseline = input$has_baseline, 
                                                          stdeffsize = dat_stdeffsize$stdeffsize,
                                                          alpha_test = dat_alpha_test$alpha_test,
                                                          target_sample_size = dat_target_sample_size$target_sample_size,
                                                          response_probability = dat_response_probability$response_probability,
                                                          rho = dat_rho$rho,
                                                          prob_missing = dat_prob_missing$prob_missing)
      }
    }else if(input$design == "design2"){
      if(input$design2aims == "stage1"){
        # Power for comparing first-stage options
        list_output <- calculate_power_stage1(design = input$design,
                                              planned_aim = input$design2aims,
                                              has_baseline = input$has_baseline, 
                                              stdeffsize = dat_stdeffsize$stdeffsize,
                                              alpha_test = dat_alpha_test$alpha_test,
                                              target_sample_size = dat_target_sample_size$target_sample_size,
                                              response_probability = dat_response_probability$response_probability,
                                              rho = dat_rho$rho,
                                              prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Power for comparing second-stage options among non-responders
        list_output <- calculate_power_design2_stage2nonresp(design = input$design,
                                                             planned_aim = input$design2aims,
                                                             has_baseline = input$has_baseline, 
                                                             stdeffsize = dat_stdeffsize$stdeffsize,
                                                             alpha_test = dat_alpha_test$alpha_test,
                                                             target_sample_size = dat_target_sample_size$target_sample_size,
                                                             response_probability = dat_response_probability$response_probability,
                                                             rho = dat_rho$rho,
                                                             prob_missing = dat_prob_missing$prob_missing)
      }
    }else{
      if(input$design3aims == "stage1"){
        # Power for comparing first-stage options
        list_output <- calculate_power_stage1(design = input$design,
                                              planned_aim = input$design3aims,
                                              has_baseline = input$has_baseline, 
                                              stdeffsize = dat_stdeffsize$stdeffsize,
                                              alpha_test = dat_alpha_test$alpha_test,
                                              target_sample_size = dat_target_sample_size$target_sample_size,
                                              response_probability = dat_response_probability$response_probability,
                                              rho = dat_rho$rho,
                                              prob_missing = dat_prob_missing$prob_missing)
      }else{
        # Power for comparing second-stage options among arm 1 non-responders
        list_output <- calculate_power_design3_stage2nonresp(design = input$design,
                                                             planned_aim = input$design3aims,
                                                             has_baseline = input$has_baseline, 
                                                             stdeffsize = dat_stdeffsize$stdeffsize,
                                                             alpha_test = dat_alpha_test$alpha_test,
                                                             target_sample_size = dat_target_sample_size$target_sample_size,
                                                             response_probability = dat_response_probability$response_probability,
                                                             rho = dat_rho$rho,
                                                             prob_missing = dat_prob_missing$prob_missing)
      }
    }
    
    return(list_output)
  })
  
  #############################################################################
  # This is the resulting power given all current inputs
  #############################################################################
  output$value_power <- renderUI({
    list_output <- react_calculate_power()
    tags$p(sprintf("%.2f", list_output$value_power), style = "font-size: 70px; font-weight: 1000;")
  })
  
  output$power_assumptions <- renderText({
    list_output <- react_calculate_power()
    
    smart_design_text <- case_when(
      list_output$dat_input_and_output$SMARTDesign == "design1" ~ "Design I",
      list_output$dat_input_and_output$SMARTDesign == "design2" ~ "Design II",
      list_output$dat_input_and_output$SMARTDesign == "design3" ~ "Design III",
      .default = NULL
    )
    
    aim_text <- case_when(
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage2resp" ~ "Compare second-stage intervention options among responders",
      list_output$dat_input_and_output$SMARTDesign == "design1" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among non-responders",
      list_output$dat_input_and_output$SMARTDesign == "design2" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design2" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among non-responders",
      list_output$dat_input_and_output$SMARTDesign == "design3" & list_output$dat_input_and_output$Aim == "stage1" ~ "Compare first-stage intervention options",
      list_output$dat_input_and_output$SMARTDesign == "design3" & list_output$dat_input_and_output$Aim == "stage2nonresp" ~ "Compare second-stage intervention options among arm 1 non-responders",
      .default = NULL
    )
    
    paste(paste("SMART Design: ", smart_design_text, sep = ""),
          paste("Aim: ", aim_text, sep = ""),
          paste("Target Sample Size: ", list_output$dat_input_and_output$SampleSize, sep = ""),
          paste("Effect Size (Cohen's d): ", list_output$dat_input_and_output$CohensD, sep = ""),
          paste("Type I Error Rate: ", round(list_output$dat_input_and_output$TypeOneErrorRate*100, 0), "%", sep = ""),
          paste("Proportion Responders: ", round(list_output$dat_input_and_output$ResponseProbability*100, 0), "%", sep = ""),
          paste("Attrition Proportion: ", round(list_output$dat_input_and_output$ProbabilityMissing*100, 0), "%", sep = ""),
          paste("Correlation between Baseline and End-of-Study Measurements: ", ifelse(is.na(list_output$dat_input_and_output$PrePostCorr), "Not Applicable", list_output$dat_input_and_output$PrePostCorr), sep = ""),
          sep = ", ")
  })
  
  #############################################################################
  # This is for the dynamically updated brief explanation of the inputs
  #############################################################################
  output$alpha_test_note <- renderText({
    
    if(input$need == "sampsize"){
      list_output <- react_calculate_sample_size()
      paste("Corresponds to a ", round(list_output$dat_input_and_output$TypeOneErrorRate*100, 0), "% significance level",sep = "")
    }else{
      list_output <- react_calculate_power()
      paste("Corresponds to a ", round(list_output$dat_input_and_output$TypeOneErrorRate*100, 0), "% significance level",sep = "")
    }
  })
  
  output$response_probability_note <- renderText({
    
    if(input$need == "sampsize"){
      list_output <- react_calculate_sample_size()
      paste("Corresponds to ", round(list_output$dat_input_and_output$ResponseProbability*100, 0), "% of total enrolled participants expected to respond",sep = "")
    }else{
      list_output <- react_calculate_power()
      paste("Corresponds to ", round(list_output$dat_input_and_output$ResponseProbability*100, 0), "% of total enrolled participants expected to respond",sep = "")
    }
  })
  
  output$prob_missing_note <- renderText({
    
    if(input$need == "sampsize"){
      list_output <- react_calculate_sample_size()
      paste("Corresponds to ", round(list_output$dat_input_and_output$ProbabilityMissing*100, 0), "% of total enrolled participants expected to lack end-of-study outcome",sep = "")
    }else{
      list_output <- react_calculate_power()
      paste("Corresponds to ", round(list_output$dat_input_and_output$ProbabilityMissing*100, 0), "% of total enrolled participants expected to lack end-of-study outcome",sep = "")
    }
  })
  
  #############################################################################
  # We will show calculation history
  #############################################################################
  history_data <- reactiveVal(NULL)
  
  observeEvent(react_calculate_sample_size()$dat_input_and_output, {
    req(input$need)
    if(input$need == "sampsize"){
      updated_log <- rbind(history_data(), react_calculate_sample_size()$dat_input_and_output)
      updated_log <- as.data.frame(updated_log)
      row.names(updated_log) <- 1:nrow(updated_log)
      history_data(updated_log)
    }
  })
  
  observeEvent(react_calculate_power()$dat_input_and_output, {
    req(input$need)
    if(input$need == "power"){
      updated_log <- rbind(history_data(), react_calculate_power()$dat_input_and_output)
      updated_log <- as.data.frame(updated_log)
      row.names(updated_log) <- 1:nrow(updated_log)
      history_data(updated_log)
    }
  })
  
  observeEvent(input$need, {
    req(input$need)
    if(input$need == "sampsize"){
      updated_log <- rbind(history_data(), react_calculate_sample_size()$dat_input_and_output)
      updated_log <- as.data.frame(updated_log)
      row.names(updated_log) <- 1:nrow(updated_log)
      history_data(updated_log)
    }else{
      updated_log <- rbind(history_data(), react_calculate_power()$dat_input_and_output)
      updated_log <- as.data.frame(updated_log)
      row.names(updated_log) <- 1:nrow(updated_log)
      history_data(updated_log)
    }
  })
  
  output$history_table <- renderDT({
    
    display_data <- history_data()
    display_data <- display_data %>%
      as.data.frame(.) %>%
      mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))
    
    datatable(
      display_data,
      rownames  = TRUE,
      extensions = "FixedColumns",
      options = list(
        scrollX = TRUE,
        #dom = "t", # commenting this argument out allows results to be displayed as pages
        fixedColumns = list(leftColumns = 3),
        lengthChange = FALSE,
        pageLength = 20
      )
    )
  })
  
  output$download_history <- downloadHandler(
    filename = function() {
      # Returns a string specifying the default file name
      paste0("calculation-history-", Sys.time(), ".csv")
    },
    content = function(file) {
      # req() ensures the app does not crash if they download an empty history
      req(history_data()) 
      
      display_data <- history_data()
      display_data <- display_data %>%
        as.data.frame(.) %>%
        mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))
      
      # Pull the current state of the history data frame and write it to 'file'
      write.csv(display_data, file, row.names = TRUE)
    }
  )
}

shinyApp(ui, server)

