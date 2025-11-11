---
title: "R-Code Walker et al 2025 - Companion associations"
author: "Jett Walker"
date: "2025-11-05"
output: html_document
---

# Purpose
#   - Fit and select a hurdle NB model for zero-inflated count data.
#   - Pre-screen with pscl::hurdle, refit with glmmTMB, rank by AICc/BIC.
#   - Report diagnostics (ROC/AUC, 10-fold CV, calibration), and effects/partials.

# Reproduction tips
#   - Use R >= 4.3. Consider `renv::init(); renv::snapshot()` for OA reproducibility.
#   - Output CSVs/figures land in ./outputs

# Inputs
#   - No_Outliers.xlsx   (expects columns: Comp_Abund, Species, windU, windV, SST,
#                        logDepth, logSalinity, logPP, logPPsquared, logDFS, logDTP, Oxy, logCHl, etc.)

---
#Packages
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr); library(tidyr); library(tibble); library(purrr)
    library(glmmTMB)
    library(AICcmodavg)
    library(glue)
    library(furrr); library(future)
    library(R.utils)
    library(pscl)
    library(readr)
    library(caret)  
    library(pROC)
    library(ggplot2); library(patchwork)
    library(scales)
  })

set.seed(123)
options(future.rng.onMisuse = "ignore")



##########################################################################################################################################
# ------------------------------- CONFIG ---------------------------------------
##########################################################################################################################################

DATA_PATH   <- "No_Outliers.xlsx"
OUT_DIR     <- "outputs"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Candidate predictors for best-subsets search
predictors <- c("logDepth","windU","windV","SST","logSalinity","logPP","logDFS","logDTP","logPPsquared")
cand_zero  <- predictors
cand_count <- predictors

# Subset limits + ranking
max_k_zero     <- 9
max_k_count    <- 9
keep_intercept <- TRUE
rank_by        <- "AICc"

# FAST pre-screen with pscl::hurdle
FAST_MODE            <- TRUE
TOP_M                <- 600
DELTA_CUT            <- 8
PRE_SCREEN_SUBSAMPLE <- 1.0

# glmmTMB refit controls + parallel
TIMEOUT_SEC <- 45
glmm_ctrl   <- glmmTMBControl(optCtrl = list(iter.max = 200, eval.max = 200))
workers     <- max(1, future::availableCores() - 1)
future::plan(multisession, workers = workers)

# Diagnostics toggles
MAKE_ROC            <- TRUE
MAKE_CV10           <- TRUE
MAKE_FOREST_IQR     <- TRUE
MAKE_PARTIALS_DELTA <- TRUE

# ---- Policy variant toggles --------------------------------------------------

RUN_BASELINE_BESTSUBSETS <- TRUE          # main best-subsets search (no DV*)
RUN_VARIANT_DVMPA        <- TRUE          # add DVMPA to both parts
RUN_VARIANT_DVHP_DVPP    <- TRUE          # include DVHP + DVPP
CV_ON_DVHP_DVPP          <- TRUE          # do 10-fold CV + calibration for DVHP+DVPP

# lock the figures to known “final” terms,
# set them here (used by ROC/partials/forest when baseline runs):
final_count_vars <- c("windU","logSalinity","logCHl","logDFS")
final_zero_vars  <- c("windU","Oxy","logSalinity","logCHl","logDFS","logDTP")

# ------------------------------ UTILITIES -------------------------------------

all_subsets <- function(vars, K, keep1 = TRUE) {
  out <- list()
  if (K >= 1) for (m in 1:K) out <- c(out, combn(vars, m, simplify = FALSE))
  if (keep1) out <- c(list(character(0)), out) # "" == intercept-only
  out
}
n_subsets <- function(p, k, keep1 = TRUE) (if (keep1) 1L else 0L) + sum(choose(p, 1:k))

build_formulas <- function(cvars, zvars) {
  crhs <- if (length(cvars)==0) "1" else paste(cvars, collapse=" + ")
  zrhs <- if (length(zvars)==0) "1" else paste(zvars, collapse=" + ")
  list(
    f_cond = as.formula(glue("Comp_Abund ~ {crhs}")),
    f_zi   = as.formula(glue("~ {zrhs}")),
    f_pscl = as.formula(glue("Comp_Abund ~ {crhs} | {zrhs}"))
  )
}

fit_one_glmm <- function(cvars, zvars, data, timeout_sec = TIMEOUT_SEC, rank_by = c("AICc","BIC")) {
  rank_by <- match.arg(rank_by)
  f <- build_formulas(cvars, zvars)
  fit <- try(
    R.utils::withTimeout(
      glmmTMB(f$f_cond, ziformula = f$f_zi,
              family = truncated_nbinom2(), data = data, control = glmm_ctrl),
      timeout = timeout_sec, onTimeout = "silent"
    ),
    silent = TRUE
  )
  if (!inherits(fit, "glmmTMB")) return(NULL)
  metric <- if (rank_by == "AICc") AICcmodavg::AICc(fit) else BIC(fit)
  tibble(
    count_k = length(cvars),
    zero_k  = length(zvars),
    count_predictors = paste(cvars, collapse = if (length(cvars)) " + " else ""),
    zero_predictors  = paste(zvars,  collapse = if (length(zvars))  " + " else ""),
    metric = metric,
    .fit = list(fit)
  )
}

nice_axis_label <- function(v) {
  s <- v
  s <- gsub("_", " ", s)
  s <- sub("^log", "log ", s)
  s <- sub("^SST$", "Sea Surface Temperature (°C)", s)
  s <- sub("^logDFS$", "Distance to Shore (log)", s)
  s <- sub("^logDTP$", "Distance to Port (log)", s)
  s <- sub("^logPP$", "Primary production (log)", s)
  s <- sub("^windU$", "Wind U (m s\u207b\u00b9)", s)
  s <- sub("^windV$", "Wind V (m s\u207b\u00b9)", s)
  s <- sub("^logSalinity$", "Salinity (0.001 psu)", s)
  s <- sub("^Oxy$", "Dissolved oxygen", s)
  s
}

print_stats_list <- function(s) {
  stopifnot(is.list(s))
  cat(sprintf(
    "N = %d\nlogLik = %.3f\nAIC = %.3f\nAICc = %.3f\nBIC = %.3f\nMcFadden R^2 = %.4f\n",
    s$N, s$logLik, s$AIC, s$AICc, s$BIC, s$McFadden_R2
  ))
  invisible(s)
}

# ------------------------------ LOAD DATA -------------------------------------

raw <- read_excel(DATA_PATH)
mydata <- raw %>%
  mutate(
    Species    = factor(Species),
    Comp_Abund = as.integer(Comp_Abund)
  )
message("Rows: ", nrow(mydata))
message("Proportion zeros: ", round(mean(mydata$Comp_Abund == 0, na.rm = TRUE), 3))



##########################################################################################################################################
# ============================ BASELINE: BEST-SUBSETS ==========================
##########################################################################################################################################

if (RUN_BASELINE_BESTSUBSETS) {
  count_sets <- all_subsets(cand_count, max_k_count, keep1 = keep_intercept)
  zero_sets  <- all_subsets(cand_zero,  max_k_zero,  keep1 = keep_intercept)
  
  S_count <- n_subsets(length(cand_count), max_k_count, keep1 = keep_intercept)
  S_zero  <- n_subsets(length(cand_zero),  max_k_zero,  keep1 = keep_intercept)
  cat("Planned combos: ", S_count, " x ", S_zero, " = ", S_count * S_zero, "\n", sep = "")
  
  grid_all <- tidyr::expand_grid(i = seq_along(count_sets),
                                 j = seq_along(zero_sets))
  
  if (FAST_MODE) {
    message("FAST pre-screen (pscl::hurdle) ...")
    mydata_quick <- if (PRE_SCREEN_SUBSAMPLE < 1) {
      ix <- sample.int(nrow(mydata), size = floor(PRE_SCREEN_SUBSAMPLE * nrow(mydata)))
      mydata[ix, , drop = FALSE]
    } else mydata
    
    lst <- furrr::future_pmap(
      list(grid_all$i, grid_all$j),
      function(i, j) {
        f <- build_formulas(count_sets[[i]], zero_sets[[j]])
        fit <- try(pscl::hurdle(f$f_pscl, data = mydata_quick,
                                dist = "negbin", zero.dist = "binomial", link = "logit"),
                   silent = TRUE)
        if (!inherits(fit, "hurdle")) return(NULL)
        tibble(i = i, j = j, quick_metric = AICcmodavg::AICc(fit))
      },
      .options  = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
    
    quick_rows <- purrr::compact(lst) %>% bind_rows()
    if (!nrow(quick_rows)) stop("Pre-screen failed for all combos.")
    quick_rows <- quick_rows %>% arrange(quick_metric) %>%
      mutate(delta = quick_metric - min(quick_metric))
    
    keep_cut <- dplyr::filter(quick_rows, delta <= DELTA_CUT)
    grid <- dplyr::bind_rows(
      keep_cut[, c("i","j")],
      dplyr::slice_head(quick_rows[, c("i","j")], n = TOP_M)
    ) %>% distinct()
    
    cat("Refitting ", nrow(grid), " combos with glmmTMB (ΔAICc≤", DELTA_CUT, " and/or top ", TOP_M, ").\n", sep = "")
  } else {
    grid <- grid_all
  }
  
  rows <- furrr::future_pmap(
    list(grid$i, grid$j),
    function(i, j) fit_one_glmm(count_sets[[i]], zero_sets[[j]], data = mydata, rank_by = rank_by),
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  ) %>% purrr::compact() %>% dplyr::bind_rows()
  
  if (!nrow(rows)) stop("No baseline models fit (failed or timed out).")
  
  best_table <- rows %>%
    arrange(metric) %>%
    mutate(delta = metric - dplyr::first(metric),
           weight = exp(-0.5 * delta) / sum(exp(-0.5 * delta)))
  
  readr::write_csv(
    best_table %>% transmute(count_k, zero_k, count_predictors, zero_predictors,
                             !!rank_by := metric, delta, weight),
    file.path(OUT_DIR, "baseline_best_subsets_ranked.csv")
  )
  
  disp <- best_table %>%
    transmute(
      count_k, zero_k,
      count_predictors = ifelse(count_predictors == "", "(1)", count_predictors),
      zero_predictors  = ifelse(zero_predictors  == "", "(1)", zero_predictors),
      metric = round(metric, 3), delta = round(delta, 3), weight = round(weight, 3)
    )
  nm <- names(disp); nm[nm=="metric"] <- rank_by; names(disp) <- nm
  cat("\n--- Baseline best-subsets (", rank_by, "), top 20 ---\n", sep = "")
  print(utils::head(disp, 20), n = 20)
  
  # Winner refit cleanly
  f_count <- as.formula(paste("Comp_Abund ~", ifelse(best_table$count_predictors[1]=="","1",best_table$count_predictors[1])))
  f_zero  <- as.formula(paste("~",           ifelse(best_table$zero_predictors[1]=="","1",best_table$zero_predictors[1])))
  best_mod <- glmmTMB(f_count, ziformula = f_zero, family = truncated_nbinom2(), data = mydata)
  saveRDS(best_mod, file.path(OUT_DIR, "baseline_best_model_glmmTMB.rds"))
}

# ============================ VARIANT B: DVHP + DVPP ==========================

if (RUN_VARIANT_DVHP_DVPP && all(c("DVHP","DVPP") %in% names(raw))) {
  dat_dvhp_pp <- raw %>% 
    transmute(
      Comp_Abund = as.integer(Comp_Abund),
      windU, logDepth, windV, logPP, logSalinity, logCHl, logDFS,
      SST, logDTP, Oxy, DVHP, DVPP
    ) %>%
    mutate(
      DVHP = factor(DVHP, levels = c(0,1), labels = c("Other","IUCN2")),
      DVPP = factor(DVPP, levels = c(0,1), labels = c("Other","IUCN45"))
    ) %>%
    drop_na()

  # Base fit (earlier structure)
  count_vars_b <- c("windU","logSalinity","logPP","DVHP","DVPP")
  zero_vars_b  <- c("windU","SST","logSalinity","logPP","logDFS","DVPP")
  fm_b <- as.formula(paste("Comp_Abund ~", paste(count_vars_b, collapse = " + "),
                           "|", paste(zero_vars_b, collapse = " + ")))
  mod_b <- pscl::hurdle(fm_b, data = dat_dvhp_pp, dist = "negbin", zero.dist = "binomial", link = "logit")
  null_b <- pscl::hurdle(Comp_Abund ~ 1 | 1, data = dat_dvhp_pp,
                         dist = "negbin", zero.dist = "binomial", link = "logit")
  stats_b <- list(
    N = nrow(dat_dvhp_pp),
    logLik = as.numeric(logLik(mod_b)),
    AIC = AIC(mod_b),
    AICc = AICcmodavg::AICc(mod_b),
    BIC = BIC(mod_b),
    McFadden_R2 = 1 - (as.numeric(logLik(mod_b)) / as.numeric(logLik(null_b)))
  )
  cat("\n--- Variant B: DVHP + DVPP (pscl::hurdle) ---\n")
  print_stats_list(stats_b)
  capture.output(summary(mod_b), file = file.path(OUT_DIR, "variant_DVHP_DVPP_summary.txt"))
  
  # Optional: glmmTMB CV on a slightly expanded zero set
  if (CV_ON_DVHP_DVPP) {
    count_cv <- c("windU","logSalinity","logPP","logDFS","DVHP","DVPP")
    zero_cv  <- c("windU","SST","logSalinity","logPP","logDFS","DVHP","DVPP")
  }
    
    dat_cv <- raw %>%
      mutate(
        .row = dplyr::row_number(),
        DVHP = factor(DVHP, levels = c(0,1), labels = c("Other","IUCN2")),
        DVPP = factor(DVPP, levels = c(0,1), labels = c("Other","IUCN45"))
      ) %>%
      dplyr::select(.row, Comp_Abund, dplyr::all_of(unique(c(count_cv, zero_cv)))) %>%
      tidyr::drop_na() %>%
      as.data.frame()
    
    form_cond <- as.formula(paste("Comp_Abund ~", paste(count_cv, collapse = " + ")))
    form_zi   <- as.formula(paste("~",               paste(zero_cv,  collapse = " + ")))
    set.seed(123); K <- 10
    folds <- caret::createFolds(dat_cv$Comp_Abund, k = K, list = TRUE, returnTrain = FALSE)
    
    per_fold <- vector("list", K); oof_pred <- vector("list", K)
    for (k in seq_len(K)) {
      test_idx  <- folds[[k]]
      train_dat <- dat_cv[-test_idx, , drop = FALSE]
      test_dat  <- dat_cv[ test_idx, , drop = FALSE]
      fit_k <- glmmTMB(form_cond, ziformula = form_zi,
                       family = truncated_nbinom2(), data = train_dat)
      pred_EY_k   <- as.numeric(predict(fit_k, newdata = test_dat, type = "response"))
      p_nonzero_k <- as.numeric(1 - predict(fit_k, newdata = test_dat, type = "zprob"))
      obs_k       <- as.numeric(test_dat$Comp_Abund)
      pres_k      <- as.integer(obs_k > 0)
      per_fold[[k]] <- tibble(fold = k,
                              RMSE = caret::RMSE(pred_EY_k, obs_k),
                              MAE  = caret::MAE (pred_EY_k, obs_k),
                              AUC  = tryCatch(as.numeric(pROC::auc(pres_k, p_nonzero_k)),
                                              error = function(e) NA_real_),
                              n = nrow(test_dat))
      oof_pred[[k]] <- tibble(fold = k, .row = test_dat$.row,
                              pred_EY = pred_EY_k, p_nonzero = p_nonzero_k, obs = obs_k)
    }
    per_fold <- bind_rows(per_fold); oof <- bind_rows(oof_pred)
    overall <- tibble(
      RMSE = caret::RMSE(oof$pred_EY, oof$obs),
      MAE  = caret::MAE (oof$pred_EY, oof$obs),
      AUC  = as.numeric(pROC::auc(as.integer(oof$obs > 0), oof$p_nonzero))
    )
    write_csv(per_fold, file.path(OUT_DIR, "variant_DVHP_DVPP_cv10_per_fold.csv"))
    write_csv(overall,   file.path(OUT_DIR, "variant_DVHP_DVPP_cv10_overall.csv"))
    print(overall)
    
    # binned calibration
    make_ntile_equal <- function(x, n_bins = 10L, tie_break = NULL) {
      if (is.null(tie_break)) tie_break <- seq_along(x)
      ord <- order(x, tie_break, na.last = NA); n <- length(ord)
      base <- floor(n / n_bins); r <- n - base * n_bins
      sizes <- rep(base, n_bins); if (r > 0) sizes[1:r] <- sizes[1:r] + 1
      bins_sorted <- rep(seq_len(n_bins), times = sizes)
      out <- integer(n); out[ord] <- bins_sorted; out
    }
    oof <- oof %>% mutate(bin = make_ntile_equal(pred_EY, n_bins = 10L, tie_break = .row))
    bins <- oof %>%
      group_by(bin) %>% summarise(
        n = n(),
        pred_mean = mean(pred_EY),
        obs_mean  = mean(obs),
        obs_sd    = sd(obs),
        .groups   = "drop"
      ) %>%
      mutate(obs_se = obs_sd / sqrt(n),
             obs_lo = obs_mean - 1.96 * obs_se,
             obs_hi = obs_mean + 1.96 * obs_se)
    
}


##########################################################################################################################################
# ============================== 10-fold CV with equal-frequency bins ==============================
##########################################################################################################################################

if (MAKE_CV10 && all(c("DVHP","DVPP") %in% names(raw))) {
  # Model terms
  count_vars <- c("windU","logSalinity","logPP","logDFS","DVHP","DVPP")
  zero_vars  <- c("windU","SST","logSalinity","logPP","logDFS","DVHP","DVPP")
  
  # Data for CV (keep original row id for tie-breaking)
  dat <- raw %>%
    dplyr::mutate(
      .row = dplyr::row_number(),
      DVHP = factor(DVHP, levels = c(0,1), labels = c("Other","IUCN2")),
      DVPP = factor(DVPP, levels = c(0,1), labels = c("Other","IUCN45"))
    ) %>%
    dplyr::select(.row, Comp_Abund, dplyr::all_of(unique(c(count_vars, zero_vars)))) %>%
    tidyr::drop_na() %>%
    as.data.frame()
  
  form_cond <- as.formula(paste("Comp_Abund ~", paste(count_vars, collapse = " + ")))
  form_zi   <- as.formula(paste("~",               paste(zero_vars,  collapse = " + ")))
  
  # 10 folds, stratified by presence/absence
  set.seed(123); K <- 10
  pres  <- as.integer(dat$Comp_Abund > 0)
  folds <- caret::createFolds(factor(pres), k = K, list = TRUE, returnTrain = FALSE)
  
  per_fold <- vector("list", K)
  oof_pred <- vector("list", K)
  
  for (k in seq_len(K)) {
    test_idx  <- folds[[k]]
    train_dat <- dat[-test_idx, , drop = FALSE]
    test_dat  <- dat[ test_idx, , drop = FALSE]
    
    fit_k <- glmmTMB(form_cond, ziformula = form_zi,
                     family = truncated_nbinom2(), data = train_dat)
    
    pred_EY   <- as.numeric(predict(fit_k, newdata = test_dat, type = "response"))  # E[Y]
    p_nonzero <- as.numeric(1 - predict(fit_k, newdata = test_dat, type = "zprob")) # Pr(Y>0)
    obs       <- as.numeric(test_dat$Comp_Abund)
    pres_k    <- as.integer(obs > 0)
    
    per_fold[[k]] <- tibble::tibble(
      fold = k,
      RMSE = caret::RMSE(pred_EY, obs),
      MAE  = caret::MAE (pred_EY, obs),
      AUC  = tryCatch(as.numeric(pROC::auc(pres_k, p_nonzero)), error = function(e) NA_real_),
      n    = nrow(test_dat)
    )
    oof_pred[[k]] <- tibble::tibble(fold = k, .row = test_dat$.row,
                                    pred_EY = pred_EY, p_nonzero = p_nonzero, obs = obs)
  }
  
  per_fold <- dplyr::bind_rows(per_fold)
  oof      <- dplyr::bind_rows(oof_pred)
  
  overall <- tibble::tibble(
    RMSE = caret::RMSE(oof$pred_EY, oof$obs),
    MAE  = caret::MAE (oof$pred_EY, oof$obs),
    AUC  = as.numeric(pROC::auc(as.integer(oof$obs > 0), oof$p_nonzero))
  )
  readr::write_csv(per_fold, file.path(OUT_DIR, "variant_DVHP_DVPP_cv10_per_fold.csv"))
  readr::write_csv(overall,   file.path(OUT_DIR, "variant_DVHP_DVPP_cv10_overall.csv"))
  print(overall)
  
  # Equal-frequency bins (for pooled OOF predictions)
  make_ntile_equal <- function(x, n_bins = 10L, tie_break = NULL) {
    if (is.null(tie_break)) tie_break <- seq_along(x)
    ord <- order(x, tie_break, na.last = NA); n <- length(ord)
    base <- floor(n / n_bins); r <- n - base * n_bins
    sizes <- rep(base, n_bins); if (r > 0) sizes[1:r] <- sizes[1:r] + 1
    bins_sorted <- rep(seq_len(n_bins), times = sizes)
    out <- integer(n); out[ord] <- bins_sorted; out
  }
  
  oof <- dplyr::mutate(oof, bin = make_ntile_equal(pred_EY, n_bins = 10L, tie_break = .row))
  
  # Binned calibration summary (+ 95% CI on observed mean)
  bins <- oof %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      n         = dplyr::n(),
      pred_mean = mean(pred_EY),
      obs_mean  = mean(obs),
      obs_sd    = sd(obs),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      obs_se = obs_sd / sqrt(n),
      obs_lo = obs_mean - 1.96 * obs_se,
      obs_hi = obs_mean + 1.96 * obs_se,
      ratio  = obs_mean / pmax(pred_mean, .Machine$double.eps)
    )
  
  readr::write_csv(bins, file.path(OUT_DIR, "calibration_bins.csv"))
  
  # Calibration plot
  rng <- range(c(bins$pred_mean, bins$obs_mean), na.rm = TRUE)
  p_cal <- ggplot2::ggplot(bins, ggplot2::aes(x = pred_mean, y = obs_mean)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = obs_lo, ymax = obs_hi), width = 0) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(ggplot2::aes(label = bin), nudge_x = -0.1, vjust = -1, size = 3) +
    ggplot2::coord_equal(xlim = rng, ylim = rng) +
    ggplot2::labs(title = "10-fold CV binned calibration (OOF predictions)",
                  x = "Mean predicted E[Y] (bin)",
                  y = "Mean observed (bin)") +
    ggplot2::theme_minimal(base_size = 12)
  
  ggplot2::ggsave(filename = file.path(OUT_DIR, "calibration_plot.png"),
                  plot = p_cal, width = 6.5, height = 5, dpi = 300)
}


##########################################################################################################################################
# ============================ PARTIAL PLOTS (concise) ============================
##########################################################################################################################################


# Use the Variant B hurdle fit for partial plots
mod <- mod_b
dat <- dat_dvhp_pp

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(cowplot)
})

# ---------- helpers ----------
mode_val <- function(x) if (is.factor(x) || is.character(x)) {
  ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]
} else stats::median(x, na.rm = TRUE)

make_grid <- function(var, data, n = 200, trim = c(0.02, 0.98)) {
  grid <- data[0, , drop = FALSE]; grid[1, ] <- lapply(data, mode_val)
  x <- data[[var]]
  if (is.numeric(x)) {
    lo <- quantile(x, trim[1], na.rm = TRUE); hi <- quantile(x, trim[2], na.rm = TRUE)
    grid <- grid[rep(1, n), , drop = FALSE]; grid[[var]] <- seq(lo, hi, length.out = n)
    grid$type__ <- "num"
  } else {
    levs <- if (is.factor(x)) levels(x) else unique(x)
    grid <- grid[rep(1, length(levs)), , drop = FALSE]
    grid[[var]] <- if (is.factor(x)) factor(levs, levels = levs) else levs
    grid$type__ <- "fac"
  }
  # carry factor levels from data
  for (nm in names(data)) if (is.factor(data[[nm]])) grid[[nm]] <- factor(grid[[nm]], levels(data[[nm]]))
  grid
}

# unified accessors for pscl::hurdle and glmmTMB
get_terms  <- function(m, part = c("count","zero")) {
  part <- match.arg(part)
  if (inherits(m,"hurdle")) {
    stats::delete.response(stats::terms(m, model = ifelse(part=="count","count","zero")))
  } else if (inherits(m,"glmmTMB")) {
    f <- stats::formula(m)[[ifelse(part=="count","cond","zi")]]
    stats::delete.response(stats::terms(f, data = model.frame(m)))
  } else stop("Unsupported model class")
}
get_fixef  <- function(m, part = c("count","zero")) {
  part <- match.arg(part)
  if (inherits(m,"hurdle")) {
    stats::coef(m, model = ifelse(part=="count","count","zero"))
  } else { # glmmTMB
    if (part=="count") glmmTMB::fixef(m)$cond else glmmTMB::fixef(m)$zi
  }
}
get_vcov   <- function(m, part = c("count","zero")) {
  part <- match.arg(part); bnm <- names(get_fixef(m, part))
  if (inherits(m,"hurdle")) {
    V <- stats::vcov(m, model = ifelse(part=="count","count","zero"))
  } else {
    V <- try(stats::vcov(m, component = ifelse(part=="count","cond","zi")), silent = TRUE)
    if (inherits(V,"try-error") || is.list(V)) V <- if (part=="count") stats::vcov(m)$cond else stats::vcov(m)$zi
  }
  as.matrix(V)[bnm, bnm, drop = FALSE]
}
get_frame  <- function(m, dat, part=c("count","zero")) {
  trm <- get_terms(m, part); stats::model.frame(trm, data = dat, na.action = stats::na.pass)
}

# build one partial (returns on response scale: EY for count; Pr(Y=0) for zero)
build_partial <- function(var, m, dat, part = c("count","zero"), n = 200) {
  part <- match.arg(part)
  trm <- get_terms(m, part); b <- get_fixef(m, part); V <- get_vcov(m, part)
  ref <- get_frame(m, dat, part); stopifnot(var %in% names(ref))
  grid <- make_grid(var, ref, n = n)
  
  X <- stats::model.matrix(trm, data = grid)
  X <- X[, intersect(colnames(X), names(b)), drop = FALSE]
  miss <- setdiff(names(b), colnames(X))
  if (length(miss)) X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames = list(NULL, miss)))
  X <- X[, names(b), drop = FALSE]
  
  eta <- as.numeric(X %*% b)
  # count part may have offset; hold at typical value
  if (part=="count") {
    off <- try(stats::model.offset(ref), silent=TRUE)
    off_val <- if (!inherits(off,"try-error") && !is.null(off)) stats::median(off, na.rm=TRUE) else 0
    eta <- eta + off_val
  }
  se <- sqrt(rowSums((X %*% V) * X)); z <- qnorm(0.975)
  
  if (part=="count") {
    tibble(var=var, x=grid[[var]], type=unique(grid$type__), fit=exp(eta),
           lwr=exp(eta - z*se), upr=exp(eta + z*se))
  } else {
    tibble(var=var, x=grid[[var]], type=unique(grid$type__), fit=plogis(eta),
           lwr=plogis(eta - z*se), upr=plogis(eta + z*se))
  }
}

plot_partial <- function(df_one, xlab, ylab) {
  if (unique(df_one$type)=="num") {
    ggplot(df_one, aes(x, fit)) +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.35) +
      geom_line(linewidth=.9) +
      labs(x=xlab, y=ylab) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text=element_text(colour="black"))
  } else {
    ggplot(df_one, aes(as.factor(x), fit)) +
      geom_pointrange(aes(ymin=lwr, ymax=upr)) +
      labs(x=xlab, y=ylab) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text=element_text(colour="black"))
  }
}

# 2-level bar for binary policy vars (count→EY, zero→Pr(Y=0))
binary_bar <- function(var, m, dat, part=c("count","zero"), labels=c("outside","inside"),
                       xlab, ylab) {
  part <- match.arg(part)
  trm <- get_terms(m, part); b <- get_fixef(m, part); V <- get_vcov(m, part)
  ref <- get_frame(m, dat, part); base <- ref[0,]; base[1,] <- lapply(ref, mode_val)
  
  # two rows: baseline vs nonref (or 0/1)
  bnames <- names(b); grid <- base[rep(1,2), , drop=FALSE]
  if (is.factor(ref[[var]])) {
    nonref <- sub(paste0("^", var), "", grep(paste0("^",var), bnames, value=TRUE)[1])
    grid[[var]] <- factor(c(levels(ref[[var]])[1], nonref),
                          levels = c(levels(ref[[var]])[1], nonref))
  } else grid[[var]] <- c(0,1)
  
  X <- stats::model.matrix(trm, data=grid)
  X <- X[, intersect(colnames(X), names(b)), drop=FALSE]
  miss <- setdiff(names(b), colnames(X))
  if (length(miss)) X <- cbind(X, matrix(0, nrow(X), length(miss), dimnames=list(NULL, miss)))
  X <- X[, names(b), drop=FALSE]
  
  eta <- as.numeric(X %*% b)
  if (part=="count") {
    off <- try(stats::model.offset(ref), silent=TRUE)
    eta <- eta + if (!inherits(off,"try-error") && !is.null(off)) stats::median(off, na.rm=TRUE) else 0
    tr <- function(v) exp(v)
  } else tr <- function(v) plogis(v)
  
  se <- sqrt(rowSums((X %*% V) * X)); z <- qnorm(0.975)
  out <- tibble(group=factor(labels, levels=labels),
                fit=tr(eta), lwr=tr(eta - z*se), upr=tr(eta + z*se))
  
  ggplot(out, aes(group, fit, fill = group)) +
    geom_col(width=.7, alpha=.85, show.legend = FALSE) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=.15) +
    labs(x=xlab, y=ylab) +
    theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(),
          axis.text=element_text(colour="black"))
}

# ---------- ORCHESTRATORS ----------
# labels
lab_count <- c(windU="Wind (U component)", logSalinity="Salinity (log)",
               logPP="Primary production (log)")
lab_zero  <- c(windU="Wind U", SST="Sea surface temperature (°C)",
               logSalinity="Salinity (log)", logDFS="Distance from shore (log)")

# which variables to plot
count_want <- c("windU","logSalinity","logPP")
zero_want  <- c("windU","SST","logSalinity","logDFS")

# build partials
count_df <- map_dfr(count_want, ~build_partial(.x, m = mod, dat = dat, part = "count"))
zero_df  <- map_dfr(zero_want,  ~build_partial(.x, m = mod, dat = dat, part = "zero"))

count_plots <- imap(split(count_df, count_df$var),
                    ~plot_partial(.x, xlab = if (!is.na(lab_count[.y])) lab_count[.y] else .y, ylab = "Exp. Count"))
zero_plots  <- imap(split(zero_df, zero_df$var),
                    ~plot_partial(.x, xlab = if (!is.na(lab_zero[.y]))  lab_zero[.y]  else .y, ylab = "Pr(Count = 0)"))

# binary bars
p_DVHP <- binary_bar("DVHP", m = mod, dat = dat, part = "count",
                     labels = c("outside","inside"),
                     xlab = "High protection", ylab = "Exp. Count")
p_DVPP <- binary_bar("DVPP", m = mod, dat = dat, part = "zero",
                     labels = c("0","1"), xlab = "DVPP", ylab = "Pr(Count = 0)")

# stitch to match layouts
p_count <- plot_grid(count_plots$windU, count_plots$logSalinity, count_plots$logPP,
                     nrow = 1, labels = c("a)","b)","c)")) +
  theme(plot.margin = margin(5,5,5,5))
p_count <- plot_grid(p_count, p_DVHP, ncol = 1, rel_heights = c(1, .85), labels = c("", "d)"))

p_zero  <- plot_grid(zero_plots$windU, zero_plots$SST, zero_plots$logSalinity,
                     nrow = 1, labels = c("a)","b)","c)"))
p_zero  <- plot_grid(p_zero, plot_grid(zero_plots$logDFS, p_DVPP, NULL, nrow = 1,
                                       rel_widths = c(1,1,1), labels = c("d)","e)","")),
                     ncol = 1, rel_heights = c(1, .85))

# save
ggsave(file.path(OUT_DIR, "count_partials.png"), p_count, width = 10, height = 7, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "zero_partials.png"),  p_zero,  width = 10, height = 7, dpi = 300, bg = "white")

             
