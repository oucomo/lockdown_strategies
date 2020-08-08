## source: https://github.com/biouea/covidhm

## set-up
devtools::install_github("biouea/covidhm", dependencies = TRUE)

## Run a single scenario and plot a network

library(covidhm)

#Load association matrices
load(file.path("..", "..", "resources", "covidhm-master", "data-raw/am_list.RData"))

#First item in the list is data across all days
m <- am_list[[1]]

#Plot network
plot_network(
  am = m,
  day = 20,
  num.initial.cases = 1,
  prop.asym = 0.4,
  delay_shape =  1,
  delay_scale = 1.4,
  prop.ascertain = 0.8,
  presymrate = 0.4,
  R = 6.5,
  outside = 0.001,
  sensitivity = "high",
  testing = "none",
  s = 333,
  isolation = FALSE,
  secondary = FALSE,
  tracing = FALSE,
  quarantine = FALSE)

net1 <- format_network(m)
case_data <- outbreak_model2(net = net1,
                            num.initial.cases = 1,
                            prop.ascertain = 0.8,
                            cap_max_days = 20,
                            R = 6.5, presymrate = 0.4,
                            delay_shape = 1,
                            delay_scale = 1.4, prop.asym = 0.4,
                            quarantine = FALSE, isolation = FALSE,
                            tracing = FALSE, secondary = FALSE,
                            outside = 0.001, sensitivity = "high",
                            testing = "none", cap_max_tests = NULL,
                            weekly = FALSE, s = NULL)


# functions ---------------------------------------------------------------

inf_prob
function (day = NULL, inc_samp = NULL, theta = NULL, R = NULL,
          contactrate = NULL, infasym = NULL)
{
  presym <- rbernoulli(length(inc_samp), p = theta)
  presym_inds <- NA
  postsym_inds <- NA
  if (sum(presym) > 0) {
    presym_inds <- sn::dsn(x = day[presym], xi = inc_samp[presym],
                           omega = (inc_samp[presym]/max(inc_samp[presym])) *
                             3, alpha = -Inf)
  }
  if (sum(!presym) > 0) {
    postsym_inds <- sn::dsn(x = day[!presym], xi = inc_samp[!presym],
                            omega = 2, alpha = Inf)
  }
  out <- rep(NA, length(inc_samp))
  out[presym] <- presym_inds
  out[!presym] <- postsym_inds
  out2 <- 1 - exp(1)^(-(out * R * contactrate * infasym))
  return(out2)
}

outbreak_setup <- function(net, num.initial.cases, incfn, delayfn, prop.asym, isolation) {

  # Set up table of population
  popsize <- length(unique(c(net$Var1,net$Var2)))
  case_data <- tibble(exposure = NA, # Exposure time of 0 for all initial cases
                      asym = purrr::rbernoulli(popsize, p = prop.asym),
                      caseid = unique(c(net$Var1,net$Var2)), # set case id
                      infector = NA,
                      onset = NA,
                      isolated_time = Inf,
                      quarantine_time = Inf,
                      test_time = Inf,
                      release_time = NA,
                      recovery_time = NA,
                      status = "S",
                      isolated = FALSE,
                      quarantined = FALSE)

  # Set up initial cases
  initial_cases <- sample(1:popsize,num.initial.cases)
  case_data$exposure[initial_cases] <- 0
  case_data$onset[initial_cases] <- incfn(num.initial.cases)
  case_data$recovery_time[initial_cases] <- case_data$onset[initial_cases] + 7
  case_data$status[initial_cases] <- "I"

  if(isolation){
    #Isolation times for symptomatic cases: onset + delay
    sym_cases <- initial_cases[!case_data$asym[initial_cases]]
    case_data$isolated_time[sym_cases] <- case_data$onset[sym_cases] +
      delayfn(length(sym_cases))
    case_data$release_time[sym_cases] <- case_data$isolated_time[sym_cases] + 14
  }


  # return
  return(case_data)
}

outbreak_model2 <- function(net = haslemere,
                           num.initial.cases,
                           prop.ascertain,
                           cap_max_days,
                           R , presymrate, delay_shape,
                           delay_scale, prop.asym,
                           quarantine, isolation,
                           tracing, secondary,
                           outside, sensitivity = "high",
                           testing = "none", cap_max_tests = NULL,
                           weekly = TRUE, s = NULL) {

  browser()

  # Set up functions to sample from distributions
  # incubation period sampling function
  incfn <- dist_setup(dist_shape = 2.322737,
                      dist_scale = 6.492272)
  # incfn <- dist_setup(dist_shape = 3.303525,dist_scale = 6.68849) # incubation function for ECDC run
  # onset to isolation delay sampling function
  delayfn <- dist_setup(delay_shape,
                        delay_scale)


  # Set initial values for loop indices
  total.cases <- num.initial.cases
  latest.onset <- 0
  extinct <- FALSE
  popsize <- nrow(net)
  cday <- 1
  daily_isolated <- 0 #none isolated on day 0
  daily_quarantined <- 0
  daily_tests <- 0




  # Initial setup
  if(exists("s")){set.seed(s)}
  case_data <- outbreak_setup(net = net,
                              num.initial.cases = num.initial.cases,
                              incfn = incfn,
                              prop.asym = prop.asym,
                              delayfn = delayfn,
                              isolation = isolation)



  # Model loop
  if(exists("s")){set.seed(s)}
  while (cday < cap_max_days & total.cases < popsize & !extinct) {
    case_data <- outbreak_step(case_data = case_data,
                               day = cday,
                               net = net,
                               incfn = incfn,
                               delayfn = delayfn,
                               prop.ascertain = prop.ascertain,
                               R = R,
                               presymrate = presymrate,
                               quarantine = quarantine,
                               isolation = isolation,
                               tracing = tracing,
                               secondary = secondary,
                               prop.asym = prop.asym,
                               outside = outside,
                               sensitivity = sensitivity,
                               testing = testing,
                               cap_max_tests = cap_max_tests)


    total.cases <- sum(!is.na(case_data$exposure))
    extinct <- all(case_data$isolated,na.rm = TRUE)
    daily_isolated <- c(daily_isolated,sum(case_data$isolated))
    daily_quarantined <- c(daily_quarantined,sum(case_data$quarantined))
    daily_tests <- c(daily_tests,sum(floor(case_data$test_time) == cday))
    cday <- cday + 1

  }

  # Prepare output, group into weeks
  weekly_isolation <- c()
  weekly_quarantine <- c()
  weekly_tests <- c()

  for(i in seq(1,cap_max_days+1,7)){
    weekly_isolation <- c(weekly_isolation,
                          mean(daily_isolated[i:(i+6)],na.rm = TRUE))

    weekly_quarantine <- c(weekly_quarantine,
                           mean(daily_quarantined[i:(i+6)],na.rm = TRUE))

    weekly_tests <- c(weekly_tests,
                      sum(daily_tests[i:(i+6)],na.rm = TRUE))
  }

  weekly_cases <- tibble(week = unique(floor((1:cap_max_days)/7)),
                         weekly_isolation,
                         weekly_tests,
                         weekly_quarantine) %>%
    left_join(case_data %>%
                dplyr::mutate(week = floor(onset/7)) %>%
                dplyr::group_by(week) %>%
                dplyr::summarise(weekly_cases = n()),
              by = "week") %>%
    mutate(weekly_cases = tidyr::replace_na(weekly_cases,0),
           weekly_isolation = tidyr::replace_na(weekly_isolation,0))


  # order and sum up, cut at max_week and add effective R0
  weekly_cases %<>%
    dplyr::arrange(week) %>%
    dplyr::mutate(cumcases = cumsum(weekly_cases),
                  cumiso = cumsum(weekly_isolation))

  # return
  if(weekly) {
    return(weekly_cases)
  } else {
    return(case_data)
  }
}

plot_network <- function(am,
                         am.layout = NULL,
                         day,
                         num.initial.cases = NULL,
                         prop.ascertain = NULL,
                         R = NULL, presymrate = NULL, delay_shape = NULL,
                         delay_scale = NULL, prop.asym = NULL,
                         quarantine = NULL, isolation = NULL,
                         tracing = NULL, secondary = NULL,
                         outside = NULL, sensitivity = NULL,
                         testing = NULL, cap_max_tests = NULL, s = NULL) {

  net1 <- format_network(am)
  case_data <- outbreak_model(net = net1,
                              num.initial.cases = num.initial.cases,
                              prop.ascertain = prop.ascertain,
                              cap_max_days = day,
                              R = R, presymrate = presymrate,
                              delay_shape = delay_shape,
                              delay_scale = delay_scale, prop.asym = prop.asym,
                              quarantine = quarantine, isolation = isolation,
                              tracing = tracing, secondary = secondary,
                              outside = outside, sensitivity = sensitivity,
                              testing = testing, cap_max_tests = cap_max_tests,
                              weekly = FALSE, s = s)
  set.seed(s)
  draw.contagion(am = am,
                 am.layout = am.layout,
                 use.df = case_data,
                 day = day)

}


