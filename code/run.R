library(tidyverse)
library(tidyr)

## simulate population
pop <- pop_sim(n = 500, hh_max = 8, hh_mean = 4, ws_max = 100, ws_mean = 50, seed = 1)

## simulate contact
contactnet <- contact_sim(pop = pop,
                          assume = cbind(hh = c(0.1, 0.2, 0.4, 0.3),
                                         ws = c(0.69, 0.2, 0.1, 0.01),
                                         ot = c(0.889, 0.1, 0.01, 0.001)),
                          seed = 1)

## setup
ncore <- 4
n_sim <- 1000
net <- contactnet
n_initial <- 1
p_asym <- 0.4
d_recovery <- 14
p_presym <- 0.4
p_asymtrans <- 0.5
R <- c(2, 4, 6.5, 8)
p_outside <- 0.001
p_contact <- 0.9
cap_max_days <- 180
outdir <- "simulation"
control_t <- c(1, 5, 10)
full_lockdown_assume <- list(cbind(hh = c(0.01, 0.02, 0.5, 0.47),
                              ws = c(0.9889, 0.01, 0.001, 0.0001),
                              ot = c(0.99889, 0.001, 0.0001, 0.00001)),
                             cbind(hh = c(0.01, 0.02, 0.5, 0.47),
                                   ws = c(1, 0, 0, 0),
                                   ot = c(1, 0, 0, 0)))

## no control
for (i in c(1:length(R))) {
  assign(paste("nocontrol", i, sep = "_"),
         purrr::map(.x = 1:n_sim,
                    ~ outbreak_model(
                      pop = pop,
                      net = net,
                      n_initial = n_initial,
                      p_asym = p_asym,
                      d_recovery = d_recovery,
                      p_presym = p_presym,
                      p_asymtrans = p_asymtrans,
                      R = R[i],
                      p_outside = p_outside,
                      p_contact = p_contact,
                      cap_max_days = cap_max_days,
                      outdir = file.path(outdir, "nocontrol", paste("scenario", i, sep = "_")),
                      control = "none"
                    )))
  assign(paste("res_nocontrol", i, sep = "_"),
         scenario_summary(dir_case_data = file.path(outdir, "nocontrol", paste("scenario", i, sep = "_"))))
  save(list = paste("res_nocontrol", i, sep = "_"), file = file.path(outdir, "nocontrol", paste0("res_nocontrol_", i, ".Rdata")))
}

## full lockdown
tmp <- expand.grid(t = 1:length(control_t), assume = 1:length(full_lockdown_assume))

for (i in c(1:nrow(tmp))) {
  assign(paste("full_lockdown", i, sep = "_"),
         purrr::map(.x = 1:n_sim,
                    ~ outbreak_model(
                      pop = pop,
                      net = net,
                      n_initial = n_initial,
                      p_asym = p_asym,
                      d_recovery = d_recovery,
                      p_presym = p_presym,
                      p_asymtrans = p_asymtrans,
                      R = R[3],
                      p_outside = p_outside,
                      p_contact = p_contact,
                      cap_max_days = cap_max_days,
                      outdir = file.path(outdir, "full_lockdown", paste("scenario", i, sep = "_")),
                      control = "full_lockdown",
                      control_t = control_t[tmp$t[i]],
                      control_assume = full_lockdown_assume[[tmp$assume[i]]]
                    )))
  assign(paste("res_full_lockdown", i, sep = "_"),
         scenario_summary(dir_case_data = file.path(outdir, "full_lockdown", paste("scenario", i, sep = "_"))))
  save(list = paste("res_full_lockdown", i, sep = "_"), file = file.path(outdir, "full_lockdown", paste0("res_full_lockdown_", i, ".Rdata")))
}


## summary
outbreak_plot(summary_data = tmp, y = c("Incidence of infection"), x = "Day", period = c(0, 10))
load(file.path(outdir, "full_lockdown", paste("scenario", 1, sep = "_"), "contact_20200802214700.Rdata"))
