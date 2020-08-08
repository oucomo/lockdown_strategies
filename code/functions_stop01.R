library(tidyverse)
library(tidyr)


# misc --------------------------------------------------------------------

## copy from covidhm
dist_setup <- function(dist_shape = NULL, dist_scale = NULL) {
  out <- purrr::partial(rweibull,
                        shape = dist_shape,
                        scale = dist_scale)
  return(out)
}


# social contact ----------------------------------------------------------

## use Haslemere network data
# load(file.path("..", "..", "resources", "covidhm-master", "data-raw/am_list.RData"))
# m <- am_list[[1]]
# diag(m) <- NA
# rownames(m) <- c(1:nrow(m))
# colnames(m) <- c(1:nrow(m))
# contactnet2 <- as_tibble(reshape2::melt(m)) %>%
#   filter(!is.na(value))

cluster_sim <- function(n, cmax, cmean, seed = NULL) {
  ## calculate distribution of cluster (binomial distribution)
  cp <- dbinom(x = 0:cmax, size = cmax, prob = cmean/cmax)
  ### combine cluster with size 0 & 1
  cp_01 <- cp[1] + cp[2]
  cp <- c(cp_01, cp[-c(1:2)])

  ## calculate number of cluster
  cn <- floor(n/(c(1:cmax) %*% cp))

  ## simulate clusters (multinomial distribution)
  if (!is.null(seed)) {set.seed(seed)}
  cluster <- rmultinom(n = 1, size = cn, prob = cp)
  csize <- rep(1:cmax, time = cluster)

  ## make sure the total subjects is fixed
  if (sum(csize) > n) {
    l <- max(which(cumsum(csize) <= n))
    csize2 <- c(csize[1:l], rep(1, n - sum(csize[1:l])))
  } else {
    csize2 <- c(csize, rep(1, n - sum(csize)))
  }

  ## output
  return(data.frame(case_id = 1:n,
                    cluster_id = rep(1:length(csize2), time = csize2)))
}

pop_sim <- function(n,
                    hh_max = 8,
                    hh_mean = 4,
                    ws_max = 100,
                    ws_mean = 50,
                    seed = NULL){
  ## simulate household cluster
  hh_cluster <- cluster_sim(n = n, cmax = hh_max, cmean = hh_mean, seed = seed)
  names(hh_cluster)[2] <- "hh_id"

  ## simulate work/school cluster
  ws_cluster <- cluster_sim(n = n, cmax = ws_max, cmean = ws_mean, seed = seed)
  names(ws_cluster)[2] <- "ws_id"

  ## output
  return(merge(hh_cluster, ws_cluster, by = "case_id"))
}

contact_sim <- function(pop,
                        assume = cbind(hh = c(0.1, 0.2, 0.4, 0.3),
                                       ws = c(0.2, 0.4, 0.3, 0.1),
                                       ot = c(0.65, 0.2, 0.1, 0.05)),
                        seed = NULL) {
  assume <- cbind(assume, hw = apply(assume[, 1:2], 1, mean))
  pop1 <- pop; names(pop1) <- c("case_1", "hh_1", "ws_1")
  pop2 <- pop; names(pop2) <- c("case_2", "hh_2", "ws_2")
  tmp <- expand.grid(case_1 = pop$case_id,
                     case_2 = pop$case_id)
  if (!is.null(seed)) {set.seed(seed)}
  tmp1 <- filter(merge(merge(tmp, pop1, by = "case_1", all.x = TRUE), pop2, by = "case_2", all.x = TRUE), case_1 != case_2) %>%
    filter(case_1 > case_2) %>%
    mutate(value = ifelse((hh_1 == hh_2) & (ws_1 == ws_2), apply(rmultinom(n = n(), size = 1, prob = assume[, 4]), 2, function(x) which(x == 1)),
                          ifelse(hh_1 == hh_2, apply(rmultinom(n = n(), size = 1, prob = assume[, 1]), 2, function(x) which(x == 1)),
                                 ifelse(ws_1 == ws_2, apply(rmultinom(n = n(), size = 1, prob = assume[, 2]), 2, function(x) which(x == 1)),
                                        apply(rmultinom(n = n(), size = 1, prob = assume[, 3]), 2, function(x) which(x == 1)))))) %>%
    mutate(value = value - 1) %>%
    select(case_1, case_2, hh_1, ws_1, hh_2, ws_2, value)
  tmp2 <- tmp1; names(tmp2) <- c("case_2", "case_1", "hh_2", "ws_2", "hh_1", "ws_1", "value")
  return(rbind(tmp1, tmp2[, c("case_1", "case_2", "hh_1", "ws_1", "hh_2", "ws_2", "value")]) %>% arrange(case_1, case_2))
}

contact_step <- function(contact_setup, day, p = 0.9) {
  q <- (1-p)/3
  out <- cbind(day = day, contact_setup) %>%
    mutate(value2 = ifelse(value == 0, 0,
                           ifelse(value == 1, sample(x = 0:3, size = n(), prob = c(q, p, q, q), replace = TRUE),
                                  ifelse(value == 2, sample(x = 0:3, size = n(), prob = c(q, q, p, q), replace = TRUE),
                                         sample(x = 0:3, size = n(), prob = c(q, q, q, p), replace = TRUE))))) %>%
    select(day, case_1, case_2, value2)
  names(out)[4] <- "value"
  return(out)
}

# epidemic model ----------------------------------------------------------

## setup
outbreak_setup0 <- function(net, n_initial, incfn, p_asym, d_recovery = 14) {
  ## net: contact network
  ## n_initial: number of initial cases
  ## incfn: function of incubation time period
  ## p_asym: probability of asymptomatic
  ## d_recovery: duration from onset to recovery

  # Set up table of population
  popsize <- length(unique(net$case_1))
  case_data <- tibble(
    # day of exposure: 0 for all initial cases, NA for uninfected cases
    time_exposure = NA,
    # asymptomatic status: 0/1
    asym = purrr::rbernoulli(popsize, p = p_asym),
    # id of all cases in the contact network
    caseid = unique(net$case_1),
    #
    infector = NA,
    # day of disease onset
    time_onset = NA,
    # day of isolation: Inf for cases without isolation
    time_isolate = Inf,
    # day of quarantine: Inf for cases without quarantine
    time_quarantine = Inf,
    # day of release: 14 days after isolation/quarantine
    time_release = NA,
    # day of recovery: 14 days after onset
    time_recovery = NA,
    # infection status: S = 0/I = 1/R = 2
    status = 0,
    # indication of being isolated
    isolated = 0,
    # indication of being quarantine
    quarantined = 0,
    # indication of being tested
    tested = 0)

  # Set up initial cases
  initial_cases <- sample(1:popsize, n_initial)
  case_data$time_exposure[initial_cases] <- 0
  case_data$time_onset[initial_cases] <- incfn(n_initial)
  case_data$time_recovery[initial_cases] <- case_data$time_onset[initial_cases] + d_recovery
  case_data$status[initial_cases] <- 1

  # return
  return(case_data)
}

## Samples the serial interval for given incubation period samples - gives an infection prob for a given day
inf_prob <- function(day = NULL, inc_samp = NULL, theta = NULL, R = NULL, contactrate = NULL, infasym = NULL) {
  ## day: day of simulation
  ## inc_samp: vector of samples from the incubation period distribution
  ## theta: probability of presymptomatic transmission
  ## R: scaling factor
  ## contactrate: ?
  ## infasym: vector of weights based on whether inds are asymptomatic

  presym <- rbernoulli(length(inc_samp), p = theta)

  presym_inds <- NA
  postsym_inds <- NA

  if(sum(presym) > 0)
  {
    presym_inds <- sn::dsn(x = day[presym],
                           xi = inc_samp[presym],
                           omega = (inc_samp[presym]/max(inc_samp[presym]))*3,
                           alpha = -Inf)
  }
  if(sum(!presym) > 0)
  {
    postsym_inds <- sn::dsn(x = day[!presym],
                            xi = inc_samp[!presym],
                            omega = 2,
                            alpha = Inf)
  }

  out <- rep(NA, length(inc_samp))
  out[presym] <- presym_inds
  out[!presym] <- postsym_inds

  out2 <- 1 - exp(1)^(-(out*R*contactrate*infasym))

  return(out2)
}

## step
outbreak_step0 <- function(day, case_data, net = haslemere,
                          p_asym, incfn, p_presym, p_asymtrans, R,
                          p_outside, d_recovery = 14) {
  ## day: day of investigation
  ## case_data: case_data at the begining of the day
  ## net: contact network
  ## p_asym: probability of asymptomatic
  ## incfn: function of incubation time period
  ## p_presym: probability of presymptomatic transmission,
  ## p_asymtrans: Infectiousness of asymptomatic individuals (relative to symptomatic cases)
  ## R: ??,
  ## p_outside: probability of infection from outside,
  ## d_recovery: duration from onset to recovery
  ## d_release: duration from isolate/quarantine to release

  # rename some variables ---------------------------------------------------

  newnet <- net[, c("case_1", "case_2", "value")]
  colnames(newnet) <- c("caseid", "contact", "rate")

  # update recovery status -----------------------------------
  # assign recovered status to recovered individuals
  recovered <- which(day > case_data$time_recovery)
  case_data$status[recovered] <- 2

  # add infections from outside ---------------------------------------------

  if(p_outside > 0) {
    potential_new_inf<- which(case_data$status == 0 & !case_data$isolated & !case_data$quarantined)
    new_infections <- potential_new_inf[rbernoulli(length(potential_new_inf), p = p_outside)]

    case_data$time_exposure[new_infections] <- day
    case_data$time_onset[new_infections] <- day + incfn(length(new_infections))
    case_data$time_recovery[new_infections] <- case_data$time_onset[new_infections] + d_recovery
    case_data$status[new_infections] <- 1
  }

  # pull out all infectious inds
  infectors <- dplyr::filter(case_data, status == 1, !isolated)

  # new potential cases  ----------------------------

  # get contacts of infectious inds from network
  new_cases <- dplyr::filter(newnet, caseid %in% infectors$caseid, rate > 0)

  new_inf_rows <- match(new_cases$contact, case_data$caseid)

  # only keep susceptible contacts
  new_cases <- new_cases[case_data$status[new_inf_rows] == 0 & !case_data$isolated[new_inf_rows] & !case_data$quarantined[new_inf_rows],]

  # generate new infections -------------------------------------------------

  #browser()

  if(nrow(new_cases) > 0) {

    # filter based on probability that each contact is infected
    infector_rows <- match(new_cases$caseid, case_data$caseid)

    asymrate <- ifelse(case_data$asym[infector_rows], p_asymtrans, 1)

    infected <- rbernoulli(nrow(new_cases),
                           p = inf_prob(day = day - case_data$time_exposure[infector_rows],
                                        inc_samp = case_data$time_onset[infector_rows],
                                        contactrate = new_cases$rate,
                                        theta = p_presym,
                                        infasym = asymrate,
                                        R = R))

    # each contact can only be infected once
    new_cases <- new_cases[infected,]

    if (any(infected)) {
      new_cases <- new_cases %>%
        group_by(contact) %>%
        sample_n(1)
    }
  }

  # compile data for all new infections -------------------------------------

  # compile a data frame for all new cases, new_cases is the amount of people that each infector has infected

  if (nrow(new_cases) > 0){

    prob_samples <- match(new_cases$contact, case_data$caseid)
    case_data$infector[prob_samples] <- case_data$caseid[match(new_cases$caseid, case_data$caseid)]
    case_data$time_exposure[prob_samples] <- day
    case_data$status[prob_samples] <- 1
    case_data$time_onset[prob_samples] <- day + incfn(length(prob_samples))
    case_data$time_recovery[prob_samples] <- case_data$time_onset[prob_samples] + d_recovery
  }
  return(case_data)
}


# outbreak_step <- function(day, case_data, net = haslemere,
#                           p_asym, incfn, delayfn,
#                           prop.ascertain, p_presym, p_asymtrans, R,
#                           p_quarantine, p_isolate, p_trace,
#                           secondary, p_outside,
#                           sensitivity = "high",
#                           d_recovery = 14, d_release = 14) {
#   ## day: day of investigation
#   ## case_data: case_data at the begining of the day
#   ## net: contact network
#   ## p_asym: probability of asymptomatic
#   ## incfn: function of incubation time period
#   ## delayfn: function of time period from onset to isolation
#   ## prop.ascertain: ???
#   ## p_presym: probability of presymptomatic transmission,
#   ## p_asymtrans: Infectiousness of asymptomatic individuals (relative to symptomatic cases)
#   ## R: ??,
#   ## p_quarantine: probability of being quarantined,
#   ## p_isolate: probability of being isolated,
#   ## p_trace: probability of being traced,
#   ## secondary: secondary contact tracing?,
#   ## p_outside: probability of infection from outside,
#   ## sensitivity = "high"???
#   ## d_recovery: duration from onset to recovery
#   ## d_release: duration from isolate/quarantine to release
#
#   # rename some variables ---------------------------------------------------
#
#   newnet <- net
#   colnames(newnet) <- c("caseid", "contact", "rate")
#
#   # update isolation, quarantine and infection status -----------------------------------
#   if (p_isolate > 0) {
#     new_isolations <- which(day > case_data$time_isolate)
#
#     if (p_quarantine > 0){
#       new_quarantines <- which(day > case_data$time_quarantine)
#     }
#
#     new_releases <- which(day > case_data$time_release)
#
#     case_data$isolated[new_isolations] <- TRUE
#     if (p_quarantine > 0){case_data$quarantined[new_quarantines] <- TRUE}
#     case_data$isolated[new_releases] <- FALSE
#     case_data$quarantined[new_releases] <- FALSE
#
#     # reset isolation, release for individuals who didn't undergo full isolation or quarantine
#     early_releases <- c(new_releases[case_data$release_time[new_releases] < (case_data$time_isolate[new_releases] + d_release)],
#                         new_releases[case_data$release_time[new_releases] < (case_data$time_quarantine[new_releases] + d_release)])
#     case_data$time_isolate[early_releases] <- Inf
#     case_data$time_release[early_releases] <- Inf
#     case_data$time_quarantine[early_releases] <- Inf
#   }
#
#   # assign recovered status to recovered individuals
#   recovered <- which(day > case_data$time_recovery)
#   case_data$status[recovered] <- "R"
#
#   # add infections from outside ---------------------------------------------
#
#   if(p_outside > 0) {
#     potential_new_inf<- which(case_data$status == "S" & !case_data$isolated & !case_data$quarantined)
#     new_infections <- potential_new_inf[rbernoulli(length(potential_new_inf), p = p_outside)]
#
#     case_data$time_exposure[new_infections] <- day
#     case_data$time_onset[new_infections] <- day + incfn(length(new_infections))
#     case_data$time_recovery[new_infections] <- case_data$time_onset[new_infections] + d_recovery
#     case_data$status[new_infections] <- "I"
#
#     # isolation times for symptomatic new infections
#     if (p_isolate > 0)
#     {
#       sym_cases <- new_infections[!case_data$asym[new_infections]]
#       isolate_cases <- sym_cases[rbernoulli(length(sym_cases), p = p_isolate)]
#       case_data$time_isolate[isolate_cases] <- case_data$time_onset[isolate_cases] + delayfn(length(isolate_cases))
#       case_data$time_release[isolate_cases] <- case_data$time_isolate[isolate_cases] + d_release
#     }
#   }
#
#   # pull out infectious inds who are not isolated
#   infectors <- dplyr::filter(case_data, status == "I", !isolated)
#
#   # new potential cases  ----------------------------
#
#   # get contacts of infectious inds from network
#   new_cases <- dplyr::filter(newnet, caseid %in% infectors$caseid, rate > 0)
#
#   new_inf_rows <- match(new_cases$contact, case_data$caseid)
#
#   # only keep susceptible contacts who are not isolated
#   new_cases <- new_cases[case_data$status[new_inf_rows] == "S" & !case_data$isolated[new_inf_rows] & !case_data$quarantined[new_inf_rows],]
#
#   # generate new infections -------------------------------------------------
#
#   if(nrow(new_cases) > 0) {
#
#     # filter based on probability that each contact is infected
#     infector_rows <- match(new_cases$caseid, case_data$caseid)
#
#     asymrate <- ifelse(case_data$asym[infector_rows], p_asymtrans, 1)
#
#     infected <- rbernoulli(nrow(new_cases),
#                            p = inf_prob(day = day - case_data$time_exposure[infector_rows],
#                                         inc_samp = case_data$time_onset[infector_rows],
#                                         contactrate = new_cases$rate,
#                                         theta = p_presym,
#                                         infasym = asymrate,
#                                         R = R))
#
#     # each contact can only be infected once
#     new_cases <- new_cases[infected,] %>%
#       group_by(contact) %>%
#       sample_n(1)
#   }
#
#   # compile data for all new infections -------------------------------------
#
#   # compile a data frame for all new cases, new_cases is the amount of people that each infector has infected
#
#   if (nrow(new_cases) > 0){
#
#     prob_samples <- match(new_cases$contact, case_data$caseid)
#     case_data$infector[prob_samples] <- case_data$caseid[match(new_cases$caseid, case_data$caseid)]
#     case_data$time_exposure[prob_samples] <- day
#     case_data$status[prob_samples] <- "I"
#     case_data$time_onset[prob_samples] <- day + incfn(length(prob_samples))
#     case_data$time_recovery[prob_samples] <- case_data$time_onset[prob_samples] + d_recovery
#
#     # isolation times for symptomatic cases
#     if (p_isolate > 0) {
#       sym_cases <- prob_samples[!case_data$asym[prob_samples]]
#       isolate_cases <- sym_cases[rbernoulli(length(sym_cases), p = p_isolate)]
#       case_data$time_isolate[isolate_cases] <- case_data$time_onset[isolate_cases] + delayfn(length(isolate_cases))
#       case_data$time_release[isolate_cases] <- case_data$time_isolated[isolate_cases] + d_release
#     }
#   }
#
#   # contact tracing ---------------------------------------------------------
#
#   # empty vector of contacts
#   traced_contacts <- c()
#
#
#   if (p_trace > 0) {
#     # get contacts of symptomatic infectors who onset yesterday
#     new_ill <- filter(case_data, status == "I", !asym, time_onset < day, time_onset > (day - 1))
#
#     if (nrow(new_ill) > 0){
#       # get contacts of infectious inds from network
#       case_contacts <- filter(newnet,
#                               caseid %in% new_ill$caseid,
#                               rate > ifelse(sensitivity == "high", 0, 1))
#
#       new_contact_rows <- match(case_contacts$contact, case_data$caseid)
#
#       # only keep contacts who are not isolated and not recovered
#       traced_contacts <- case_contacts[case_data$status[new_contact_rows] != "R" &
#                                          !case_data$isolated[new_contact_rows] &
#                                          !case_data$quarantined[new_contact_rows],]
#       traced_contacts <- traced_contacts$contact[rbernoulli(length(traced_contacts$contact), p = p_trace)]
#     }
#   }
#
#   # secondary contact tracing ----------------
#
#   # estimate whether secondary contacts are traced
#   if (secondary) {
#     if(nrow(new_ill) > 0) {
#       if(nrow(case_contacts) > 0) {
#         sec_contacts <- filter(newnet,
#                                caseid %in% case_contacts$contact,
#                                rate > ifelse(sensitivity == "high", 0, 1))
#
#         new_contact_rows <- match(sec_contacts$contact, case_data$caseid)
#
#
#         # only keep secondary contacts who are not isolated and not recovered
#         traced_sec_contacts <- sec_contacts[case_data$status[new_contact_rows] != "R" &
#                                               !case_data$isolated[new_contact_rows] &
#                                               !case_data$quarantined[new_contact_rows],]
#
#         # filter based on ascertainment rate
#         traced_sec_contacts <- traced_sec_contacts$contact[rbernoulli(length(traced_sec_contacts$contact), p = p_trace)]
#
#         # add to list of traced contacts
#         traced_contacts <- unique(c(traced_contacts, traced_sec_contacts))
#       }
#     }
#   }
#
#   # set quarantine times based on tracing -----------------------------------
#
#   if((p_isolate > 0) & (p_trace > 0) & (length(traced_contacts) > 0)) {
#
#     if (p_quarantine > 0) {
#       # if you are recovered and asymptomatic and traced, you isolate
#       recovered_traced <- which(case_data$caseid %in% traced_contacts &
#                                   case_data$status == "R" &
#                                   case_data$asym)
#       case_data$time_quarantine[recovered_traced] <- day + delayfn(length(recovered_traced))
#
#       # if you are susceptible you quarantine on being traced
#       susceptible_traced <- which(case_data$caseid %in% traced_contacts & case_data$status == "S")
#
#       case_data$time_quarantine[susceptible_traced] <- day + delayfn(length(susceptible_traced))
#
#       # if you are infectious you isolate if the delay is shorter than your current iso time
#       infectious_traced <- which(case_data$caseid %in% traced_contacts & case_data$status == "I")
#       new_iso_time <- day + delayfn(length(infectious_traced))
#       case_data$time_quarantine[infectious_traced] <-  day + delayfn(length(infectious_traced))
#
#     }
#
#     case_data$time_release <- ifelse(case_data$time_isolated < case_data$time_quarantine,
#                                      case_data$time_isolated + d_release,
#                                      case_data$time_quarantine + d_release)
#
#   }
#
#   return(case_data)
# }

## overall
outbreak_model <- function(pop,
                           net,
                           n_initial,
                           p_asym,
                           d_recovery,
                           p_presym,
                           p_asymtrans,
                           R,
                           p_outside,
                           p_contact,
                           cap_max_days = 90,
                           seed = NULL,
                           outdir = "simulation",
                           control = c("none", "full_lockdown"),
                           control_t = NULL,
                           control_assume = NULL) {
  # control
  if (length(control) > 1){control <- "none"}

  # create output directory
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  dt <- gsub(pattern = " |-|:|+", replacement = "", x = Sys.time())
  filename <- paste0(paste("case_data", dt, sep = "_"), ".Rdata")
  ptm <- proc.time()
  cat("start:", as.character(Sys.time()), "\n")

  # Set up functions to sample from distributions
  # incubation period sampling function
  incfn <- dist_setup(dist_shape = 2.322737,
                      dist_scale = 6.492272)
  # incfn <- dist_setup(dist_shape = 3.303525,dist_scale = 6.68849) # incubation function for ECDC run

  # Set initial values for loop indices
  popsize <- length(unique(c(net$case_1, net$case_2)))
  n_susceptible <- popsize - n_initial
  cday <- 1

  # save contact data
  save(net, file = file.path(outdir, paste0(paste("contact", dt, sep = "_"), ".Rdata")))
  contact <- net

  # Initial setup
  if (!is.null(seed)){set.seed(seed)}
  case_data <- outbreak_setup0(net = net,
                               n_initial = n_initial,
                               incfn = incfn,
                               p_asym = p_asym,
                               d_recovery = d_recovery)
  out <- array(dim = c(nrow(case_data), ncol(case_data), cap_max_days + 1),
               dimnames = list(rownames(case_data), colnames(case_data), 0:cap_max_days))
  out[, , 1] <- as.matrix(case_data)
  save(out, file = file.path(outdir, filename))

  #browser()
  # Model loop
  if (!is.null(seed)){set.seed(seed)}
  while (cday < cap_max_days & n_susceptible > 0) {
    cat("\r", cday)

    if (control != "none") {
      if (cday == control_t) {
        net <- contact_sim(pop = pop,
                           assume = control_assume)
      }
    }

    ## update contact data
    newnet <- contact_step(contact_setup = net, day = cday, p = p_contact)
    contact[, paste0("D", cday)] <- newnet$value
    save(contact, file = file.path(outdir, paste0(paste("contact", dt, sep = "_"), ".Rdata")))

    case_data <- outbreak_step0(day = cday,
                                case_data = case_data,
                                net = newnet,
                                incfn = incfn,
                                p_presym = p_presym,
                                p_asymtrans = p_asymtrans,
                                R = R,
                                p_outside = p_outside,
                                d_recovery = d_recovery)
    out[, , cday + 1] <- as.matrix(case_data)
    n_susceptible <- sum(case_data$status == 0)
    cday <- cday + 1
    save(out, file = file.path(outdir, filename))
  }
  ptm <- proc.time() - ptm
  cat("\n finish:", as.character(Sys.time()), "\n elapse: ", ptm[["elapsed"]]/60, " min \n")

  # return
  return(out)
}


# sum(test[, 10, 92] == 0)

# control strategies ------------------------------------------------------


# summary -----------------------------------------------------------------

outbreak_summary <- function(file_case_data){
  ## load case data
  load(file_case_data)
  ## daily summary
  daily_sum <- as.data.frame(do.call(rbind, lapply(1:dim(out)[3], function(i){
    cbind(day = i - 1,
          infection = sum(out[, 10, i] == 1, na.rm = TRUE),
          susceptible = sum(out[, 10, i] == 0, na.rm = TRUE),
          recovery = sum(out[, 10, i] == 2, na.rm = TRUE),
          isolated = sum(out[, 11, i], na.rm = TRUE),
          quarantined = sum(out[, 11, i], na.rm = TRUE))
  }))) %>%
    mutate(t = day,
           infection_new = pmax(diff(c(0, infection)), 0),
           infection_cum = cumsum(infection_new),
           recovery_new = pmax(diff(c(0, recovery)), 0),
           recovery_cum = cumsum(recovery_new),
           isolated_new = pmax(diff(c(0, isolated)), 0),
           isolated_cum = cumsum(isolated_new),
           quarantined_new = pmax(diff(c(0, quarantined)), 0),
           quarantined_cum = cumsum(quarantined_new))

  ## weekly summary
  weekly_sum <- daily_sum %>%
    mutate(week = ceiling(day/7)) %>%
    group_by(week) %>%
    summarise(infection = mean(infection),
              infection_new = mean(infection_new),
              infection_cum = mean(infection_cum),
              susceptible = mean(susceptible),
              recovery = mean(recovery),
              recovery_new = mean(recovery_new),
              recovery_cum = mean(recovery_cum),
              isolated = mean(isolated),
              isolated_new = mean(isolated_new),
              isolated_cum = mean(isolated_cum),
              quarantined = mean(quarantined),
              quarantined_new = mean(quarantined_new),
              quarantined_cum = mean(quarantined_cum),
              .groups = "drop") %>%
    mutate(t = week)
  ## return
  return(list(daily = daily_sum,
       weekly = weekly_sum))
}

scenario_summary <- function(dir_case_data) {
  tmp <- list.files(path = dir_case_data, pattern = "case_data")
  daily <- weekly <- vector("list", length = length(tmp))
  for (i in c(1:length(tmp))) {
    tmpi <- outbreak_summary(file_case_data = file.path(dir_case_data, tmp[i]))
    daily[[i]] <- cbind(sim = i, tmpi$daily)
    weekly[[i]] <- cbind(sim = i, tmpi$weekly)
  }
  return(list(daily = do.call(rbind, daily),
              weekly = do.call(rbind, weekly)))
}

outbreak_plot <- function(summary_data, y = c("Incidence of infection"), x = c("Week", "Day"), period = c(0, Inf), type = c("line", "bar")) {
  require(ggplot2)
  if (length(y) > 1) y <- y[1]
  if (length(x) > 1) x <- x[1]
  if (length(type) > 1) type <- type[1]

  ## choose x
  tmp <- switch(x,
                "Day" = summary_data$daily,
                "Week" = summary_data$weekly)

  ## choose y
  tmp$y <- switch(y,
                  "Incidence of infection" = tmp$infection_new)

  ## choose label
  xlabel <- x
  ylabel <- y

  p <- ggplot(data = filter(tmp, t >= period[1], t <= period[2]), aes(x = t, y = y)) +
    theme_bw() +
    xlab(xlabel) + ylab(ylabel)

  if (type == "line"){p <- p + geom_line()}
  if (type == "bar"){p <- p + geom_bar(stat = "identity")}

  p
}

# outbreak_plot(summary_data = tmp, y = c("Incidence of infection"), x = "Week", period = c(0, 10))
# outbreak_plot(summary_data = tmp, y = c("Incidence of infection"), x = "Week", period = c(0, 20), type = "bar")

# cost-effectiveness ------------------------------------------------------


