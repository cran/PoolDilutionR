## ---- include = FALSE---------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PoolDilutionR)

OneSample_dat <- subset(Morris2023, subset = id == "71")
print(OneSample_dat)

## ----Solve for P and k--------------------------------------------------------
results <- pdr_optimize(
  time = OneSample_dat$time_days, # time as a vector
  m = OneSample_dat$cal12CH4ml + OneSample_dat$cal13CH4ml, # total pool size
  n = OneSample_dat$cal13CH4ml, # pool size of heavy isotopologue
  P = 0.1, # inital production rate for optim()
  pool = "CH4", # indicates use of default fractionation constants for methane
  m_prec = 0.001, # instrument precision for total pool size, as standard deviation
  ap_prec = 0.01, # instrument precision for atom percent, as standard deviation
)

## ----print(results)-----------------------------------------------------------
print(results)

## ----OneSample Results Dataframe----------------------------------------------
pdr_optimize_df(
  time = OneSample_dat$time_days,
  m = OneSample_dat$cal12CH4ml + OneSample_dat$cal13CH4ml,
  n = OneSample_dat$cal13CH4ml,
  P = 0.1,
  pool = "CH4",
  m_prec = 0.001,
  ap_prec = 0.01
)

## ----Multi-sample Example Input, fig.height = 5, fig.width = 8, fig.align='center'----
# create long data for plotting
library(tidyr)
long_Morris2023 <- pivot_longer(Morris2023, cols = cal12CH4ml:AP_obs)

library(ggplot2)
ggplot(
  data = long_Morris2023,
  aes(time_days, value, color = id)
) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y")

## ----lapply(PoolDilutionR), message = FALSE-----------------------------------
samples <- split(Morris2023, Morris2023$id)

all_results <- lapply(samples, FUN = function(samples) {
  pdr_optimize_df(
    time = samples$time_days,
    m = samples$cal12CH4ml + samples$cal13CH4ml,
    n = samples$cal13CH4ml,
    P = 0.1,
    pool = "CH4",
    m_prec = 0.001,
    ap_prec = 0.01,
    quiet = TRUE,
    include_progress = FALSE
  )
})
all_results <- do.call(rbind, all_results)

## ----Six samples with all preds, echo=FALSE, message = FALSE, results='hide'----
pk_results <- list()
incdat_out <- list()
all_predictions <- list()

library(tibble)

for (i in unique(Morris2023$id)) {
  message("------------------- ", i)
  # Isolate this sample's data
  dat <- Morris2023[Morris2023$id == i, ]


  result <- pdr_optimize(
    time = dat$time_days,
    m = dat$cal12CH4ml + dat$cal13CH4ml,
    n = dat$cal13CH4ml,
    P = 0.1,
    pool = "CH4",
    m_prec = 0.001,
    ap_prec = 0.01,
    quiet = TRUE,
    include_progress = TRUE
  )

  # Save progress details separately so they don't print below
  progress_detail <- result$progress
  result$progress <- NULL

  P <- result$par["P"]
  id <- dat$id[1]
  pk_results[[i]] <- tibble(
    id = id,
    P = P,
    k = result$par["k"],
    k0 = result$initial_par["k"],
    convergence = result$convergence,
    message = result$message
  )

  # Predict based on the optimized parameters
  pred <- pdr_predict(
    time = dat$time_days,
    m0 = dat$cal12CH4ml[1] + dat$cal13CH4ml[1],
    n0 = dat$cal13CH4ml[1],
    P = P,
    k = result$par["k"],
    pool = "CH4"
  )
  dat <- cbind(dat, pred)

  # Predict based on ALL the models that were tried

  x <- split(progress_detail, seq_len(nrow(progress_detail)))
  all_preds <- lapply(x, FUN = function(x) {
    y1 <- data.frame(
      P = x$P[1],
      k = x$k[1],
      time = seq(min(dat$time_days), max(dat$time_days), length.out = 20)
    )
    y2 <- pdr_predict(
      time = y1$time,
      m0 = dat$cal12CH4ml[1] + dat$cal13CH4ml[1],
      n0 = dat$cal13CH4ml[1],
      P = x$P[1],
      k = x$k[1],
      pool = "CH4"
    )
    cbind(y1, y2)
  })
  all_predictions[[i]] <- do.call(rbind, all_preds)
  all_predictions[[i]]$id <- id
  # Calculate implied consumption (ml) based on predictions
  # Equation 4: dm/dt = P - C, so C = P - dm/dt
  total_methane <- dat$cal12CH4ml + dat$cal13CH4ml
  change_methane <- c(0, diff(total_methane))
  change_time <- c(0, diff(dat$time_days))
  dat$Pt <- P * change_time # P is ml/day
  # amount of methane produced at time (t) of this incubation, a volume in mL
  dat$Ct <- dat$Pt - change_methane

  incdat_out[[i]] <- dat
}

pk_results <- do.call(rbind, pk_results)
incdat_out <- do.call(rbind, incdat_out)
all_predictions <- do.call(rbind, all_predictions)

## ----Future Note, include = FALSE, echo = FALSE-------------------------------
# For consideration, include link here to additional supplementary materials that include code for making the plots below

## ----Multisample Fit APE, echo = FALSE, fig.height = 5.5, fig.width = 8.25----
# ----- Plot AP results -----
ggplot(incdat_out, aes(time_days, AP_obs, color = id)) +
  geom_point(aes(shape = ""), size = 4) +
  geom_line(
    data = all_predictions,
    aes(time, AP_pred, group = paste(id, P, k)), color = "grey", linetype = 2
  ) +
  geom_line(aes(y = AP_pred, group = id, linetype = ""),
    size = 1.5
  ) +
  scale_linetype_manual(
    name = "Prediction",
    values = "dotted"
  ) +
  scale_shape_manual(
    name = "Observations",
    values = 20
  ) +
  scale_color_discrete(guide = "none") +
  facet_wrap(~id, scales="free_y") +
  xlab("\n Timestep \n") +
  ylab("\n (13C-CH4/Total CH4) x 100 \n") +
  ggtitle("\n Atom% 13C \n") +
  theme(legend.position = "bottom")

## ----Multisample Fit Total Pool, echo = FALSE, fig.height = 5.5, fig.width = 8.25----
ggplot(incdat_out, aes(time_days, cal12CH4ml + cal13CH4ml, color = id)) +
  geom_point(aes(shape = ""), size = 4) +
  geom_line(
    data = all_predictions,
    aes(time, mt, group = paste(id, P, k)), color = "grey", linetype = 2
  ) +
  geom_line(aes(y = mt, group = id, linetype = ""),
    size = 1.5
  ) +
  scale_linetype_manual(
    name = "Prediction",
    values = "dotted"
  ) +
  scale_shape_manual(
    name = "Observations",
    values = 20
  ) +
  scale_color_discrete(guide = "none") +
  facet_wrap(~id, scales="free_y") +
  xlab("\n Timestep \n") +
  ylab("\n Volume (mL) \n") +
  ggtitle("\n Total Methane \n") +
  theme(legend.position = "bottom")

## ----Optimize frac rates, fig.height = 4, fig.width = 6-----------------------
# In this example, we assume that P and k are known, but that fractionation rates are unknown.

Sample64_dat <- subset(Morris2023, subset = id == "64")

new_fit <- pdr_optimize_df(
  time = Sample64_dat$time_days,
  m = Sample64_dat$cal12CH4ml + Sample64_dat$cal13CH4ml,
  n = Sample64_dat$cal13CH4ml,
  P = 8,
  k = 1,
  params_to_optimize = c("frac_P", "frac_k"),
  m_prec = 0.001,
  ap_prec = 0.01
)

newFitFracP <- new_fit[new_fit$par == "frac_P", ]$value
newFitFrack <- new_fit[new_fit$par == "frac_k", ]$value

# dataframe with new fit (optimizing P_frac and k_frac)
y_new <- pdr_predict(
  time = Sample64_dat$time_days,
  m0 = Sample64_dat$cal12CH4ml[1] + Sample64_dat$cal13CH4ml[1],
  n0 = Sample64_dat$cal13CH4ml[1],
  P = 8,
  k = 1,
  frac_P = newFitFracP,
  frac_k = newFitFrack
)


dat_new <- cbind(Sample64_dat, y_new)
dat_new$oldAPpred <- incdat_out[incdat_out$id == "64", ]$AP_pred


# graph new fit
ggplot(data = dat_new) +
  geom_point(aes(time_days, AP_pred, shape = "New Prediction", color = "New Prediction"),
    size = 3, fill = "#619CFF"
  ) +
  geom_point(aes(time_days, AP_obs, shape = "Observations", color = "Observations"),
    size = 3
  ) +
  geom_point(aes(time_days, oldAPpred, shape = "Old Prediction", color = "Old Prediction"),
    size = 3
  ) +
  scale_shape_manual(
    name = "Sample 64",
    breaks = c("New Prediction", "Observations", "Old Prediction"),
    values = c("New Prediction" = 21, "Observations" = 16, "Old Prediction" = 1)
  ) +
  scale_color_manual(
    name = "Sample 64",
    breaks = c("New Prediction", "Observations", "Old Prediction"),
    values = c("New Prediction" = "black", "Observations" = "#619CFF", "Old Prediction" = "#619CFF")
  ) +
  ylab("Atom%13C\n")

print(new_fit)

