# Setup ####

# load libraries
library(tidyverse)
library(fs)
library(MplusAutomation)
library(here)
library(corx)
library(lavaan)
library(semTools)


## Functions ####

modify_mplus <- function(x, name, old = NULL, new = NULL, int = FALSE) {
  if (name %in% c("TITLE", "VARIABLE", "DEFINE", "MODEL", "ANALYSIS", "MODELCONSTRAINT")) {
    if (!is.null(x)) {
      if (!is.null(old) & !is.null(new)) {
        x <- stringr::str_replace_all(x, old, new)
      }

      if (int) {
        x <- stringr::str_replace_all(x, ";{0,1}!@{1}([A-Z]| )", "\\1")
      }
    }
  }
  x
}

calculate_mplus <- function(x, y, dir) {
  if (length(class(x)) == 1) {
    class(x) <- c("mplusObject", "list")
  }

  MplusAutomation::mplusModeler(
    object = x,
    dataout = "../data/data.dat",
    modelout = here::here(paste0("mplus/", dir, "/", y, ".inp")),
    writeData = "never",
    check = FALSE,
    run = TRUE,
    hashfilename = FALSE
  )
}

add_constraint <- function(x, y, z) {
  x$MODELTEST <- NULL
  constraint <- which(y == z$Model)
  if (!purrr::is_empty(constraint)) {
    x$MODELCONSTRAINT <- paste0(x$MODELCONSTRAINT, "\n ", z$Constraint[constraint])
  }
  return(x)
}


remove_quote <- function(x, y) {
  if (y %in% c("VARIABLE", "MODEL", "ANALYSIS")) {
    stringr::str_remove_all(x, ";{0,1}!@")
  } else {
    x
  }
}

stars_fun <- function(x) {
  dplyr::case_when(
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    .default = ""
  )
}


diff_fun <- function(m0, m1, corrected = TRUE, as_vector = FALSE) {
  
  if (purrr::pluck_exists(m0$results)) {
    m0 <- m0$results$summaries
  } else {
    m0 <- m0$summaries
  }
  
  if (purrr::pluck_exists(m1$results)) {
    m1 <- m1$results$summaries
  } else {
    m1 <- m1$summaries
  }
  
  df <- abs(m0$Parameters - m1$Parameters)
  
  lambda <- -2 * (m0$LL - m1$LL)
  
  if (corrected) {
    cd <- (m0$Parameters * m0$LLCorrectionFactor - m1$Parameters * m1$LLCorrectionFactor)
    cd <- cd / (m0$Parameters - m1$Parameters)
    lambda <- lambda / cd
  }
  
  lrtP <- pchisq(lambda, df = df, lower.tail = FALSE)
  
  if (as_vector) {
    c(df = df, ChiSqDiff = lambda, pval = lrtP)
  } else {
    paste0(
      "$\\chi^2$(", df, ") = ", papaja::apa_num(lambda), ", p ",
      papaja::apa_p(lrtP, add_equals = TRUE)
    )
  }
}

create_mi_syntax_fun <- function(mod, var, pat) {
  mod$VARIABLE <- stringr::str_remove_all(mod$VARIABLE, var)
  mod$MODEL <- stringr::str_remove_all(mod$MODEL, var)
  mod$MODELCONSTRAINT <- stringr::str_remove_all(mod$MODELCONSTRAINT, var)
  mod$MODELCONSTRAINT <- stringr::str_remove_all(mod$MODELCONSTRAINT, pat)
  return(mod)
}

## data ####
# load data
dat <- readr::read_csv("dat_mplus.csv")

# recode gender
dat <- dat |>
  dplyr::mutate(
    sex = dplyr::case_match(
      sex,
      "M" ~ 0,
      "F" ~ 1,
      .default = NA
    )
  )

# create mplus dir
fs::dir_create(here::here("mplus"))
fs::dir_create(here::here("mplus/data"))
fs::dir_create(here::here("mplus/cfa"))
fs::dir_create(here::here("mplus/sem"))

# Mplus Analysis
# Write data
MplusAutomation::prepareMplusData(
  df = dat,
  filename = here::here("mplus/data", "data.dat"),
  writeData = "always",
  interactive = FALSE,
  inpfile = FALSE,
  quiet = TRUE
)

# Mplus Analysis ####

## CFA: Single Time Point #####

### Subject Interest ####
# 4f model
sg_cfa_iv_4f <- MplusAutomation::mplusObject(
  TITLE = "Intrinsic Value - 4F-CFA model",
  VARIABLE = "
  IDVARIABLE = ID;
  USEVAR =
    IV_CH_1 - IV_CH_5
    IV_GE_1 - IV_GE_5
    IV_MA_1 - IV_MA_5
    IV_PH_1 - IV_PH_5;",
  ANALYSIS = "
  ESTIMATOR = MLR;",
  MODEL = "
  IV_CH BY IV_CH_1*
    IV_CH_2 - IV_CH_5;

  IV_GE BY IV_GE_1*
    IV_GE_2 - IV_GE_5;

  IV_MA BY IV_MA_1*
    IV_MA_2 - IV_MA_5;

  IV_PH BY IV_PH_1*
    IV_PH_2 - IV_PH_5;

  IV_CH-IV_PH@1;",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)

fit_sg_cfa_iv_4f <- MplusAutomation::mplusModeler(
  object = sg_cfa_iv_4f,
  dataout = "../data/data.dat",
  modelout = here::here("mplus/cfa", "sg_cfa_iv_4f.inp"),
  writeData = "never",
  check = FALSE,
  run = TRUE,
  hashfilename = FALSE
)

# bifactor model
sg_cfa_iv_bf <- MplusAutomation::mplusObject(
  TITLE = "Intrinsic Value - BF-CFA model",
  VARIABLE = "
  IDVARIABLE = ID;
  USEVAR =
    IV_CH_1 - IV_CH_5
    IV_GE_1 - IV_GE_5
    IV_MA_1 - IV_MA_5
    IV_PH_1 - IV_PH_5;",
  ANALYSIS = "
  ESTIMATOR = MLR;
  MODEL = NOCOVARIANCES;
  STARTS = 20;",
  MODEL = "
  IV_CH BY IV_CH_1*
    IV_CH_2 - IV_CH_5;

  IV_GE BY IV_GE_1*
    IV_GE_2 - IV_GE_5;

  IV_MA BY IV_MA_1*
    IV_MA_2 - IV_MA_5;

  IV_PH BY IV_PH_1*
    IV_PH_2 - IV_PH_5;

  IV_CH-IV_PH@1;

  IV BY IV_CH_1*
    IV_CH_2 - IV_PH_5;

  IV@1;
  ",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)

fit_sg_cfa_iv_bf <- MplusAutomation::mplusModeler(
  object = sg_cfa_iv_bf,
  dataout = "../data/data.dat",
  modelout = here::here("mplus/cfa", "sg_cfa_iv_bf.inp"),
  writeData = "never",
  check = FALSE,
  run = TRUE,
  hashfilename = FALSE
)

### Conscientiousness ####
sg_cfa_bfc_1f <- MplusAutomation::mplusObject(
  TITLE = "Big five Consciousness - 1F-CFA model",
  VARIABLE = "
  IDVARIABLE = ID;
  USEVAR =
    BF_C_1 - BF_C_4;",
  ANALYSIS = "
  ESTIMATOR = MLR;",
  MODEL = "
  BF_C BY BF_C_1*
    BF_C_2 - BF_C_4;

  BF_C@1;",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)

fit_sg_cfa_bfc_1f <- MplusAutomation::mplusModeler(
  object = sg_cfa_bfc_1f,
  dataout = "../data/data.dat",
  modelout = here::here("mplus/cfa", "sg_cfa_bfc_1f.inp"),
  writeData = "never",
  check = FALSE,
  run = TRUE,
  hashfilename = FALSE
)


mods_sg_cfa <- MplusAutomation::readModels(
  "mplus/cfa",
  filefilter = "sg_cfa_(iv|bfc)",
  what = c("summaries", "warn_err")
)

mods_sum_sg_cfa <- purrr::map(
  mods_sg_cfa, function(x) data.frame(x$summaries)
) |>
  dplyr::bind_rows()

mods_sum_sg_cfa |>
  dplyr::select(Filename, CFI, TLI, RMSEA_Estimate, SRMR) |>
  knitr::kable()

rm(list = ls(envir = .GlobalEnv)[grepl("sg_cfa", ls(envir = .GlobalEnv))])


## CFA: Multiple Time Points - Measurement Invariance #####

sg_cfa_mitime <- MplusAutomation::mplusObject(
  TITLE = "CFA FOR MI Testing over TIME",
  VARIABLE = "
  IDVARIABLE = ID;
  CLUSTER = CLASSID;
  USEVAR =
    !@VMAA_MSI_2 - AA_MSI_6
    !@VMBA_MSI_2 - BA_MSI_6
    !@VEAB_EF_1 - AB_EF_3
    !@VEBB_EF_1 - BB_EF_3
    ;",
  ANALYSIS = "
  ESTIMATOR = MLR;
  TYPE = COMPLEX;",
  MODEL = "
  !@VMAA_MSI BY AA_MSI_2* (lamsi1)
  !@VM  AA_MSI_4 - AA_MSI_6 (lamsi2 lamsi3);

  !@VM[AA_MSI_2 - AA_MSI_6] (iamsi1 - iamsi3);
  !@VMAA_MSI_2 - AA_MSI_6 (eamsi1 - eamsi3);

  !@VM[AA_MSI@0];
  !!@VMAA_MSI@1;

  !@VMBA_MSI BY BA_MSI_2* (lbmsi1)
  !@VM  BA_MSI_4 - BA_MSI_6 (lbmsi2 lbmsi3);

  !@VM[BA_MSI_2 - BA_MSI_6] (ibmsi1 - ibmsi3);
  !@VMBA_MSI_2 - BA_MSI_6 (ebmsi1 - ebmsi3);

  !@VM[BA_MSI] (mbamsi);
  
  !@VMAA_MSI_2 - AA_MSI_6 PWITH BA_MSI_2 - BA_MSI_6;

  ! Effort
  !@VEAB_EF BY AB_EF_1* (laef1)
  !@VE  AB_EF_2 - AB_EF_3 (laef2 laef3);

  !@VE[AB_EF_1 - AB_EF_3] (iaef1 - iaef3);
  !@VEAB_EF_1 - AB_EF_3 (eaef1 - eaef3);

  !@VE[AB_EF@0];
  !!@VEAB_EF@1;

  !@VEBB_EF BY BB_EF_1* (lbef1)
  !@VE  BB_EF_2 - BB_EF_3 (lbef2 lbef3);

  !@VE[BB_EF_1 - BB_EF_3] (ibef1 - ibef3);
  !@VEBB_EF_1 - BB_EF_3 (ebef1 - ebef3);

  !@VE[BB_EF] (mbbef);

  !@VEAB_EF_1 - AB_EF_3 PWITH BB_EF_1 - BB_EF_3;
  ",
  MODELCONSTRAINT = "
  ! factor loadings;
  !@VE!laef1 = 3 - (laef2 + laef3);
  !@VE!!@Clbef1 = 3 - (lbef2 + lbef3);
  !@VM!lamsi1 = 3 - (lamsi2 + lamsi3);
  !@VM!!@Clbmsi1 = 3 - (lbmsi2 + lbmsi3);
  
  !@VE0 = 1 - laef1;
  !@VM0 = 1 - lamsi2;
  !@VE!@C0 = 1 - lbef1;
  !@VM!@C0 = 1 - lbmsi2;
  !@VM!@CM0 = 1 - mbamsi;
  !@VE!@CM0 = 1 - mbbef;

  ! intercepts;
  !@VE!iaef1 = 0 - (iaef2 + iaef3);
  !@VE!!@CMibef1 = 0 - (ibef2 + ibef3);
  !@VM!iamsi1 = 0 - (iamsi2 + iamsi3);
  !@VM!!@CMibmsi1 = 0 - (ibmsi2 + ibmsi3);

  ! time invariance
  !@VE!@M0 = laef1 - lbef1;
  !@VE!@M0 = laef2 - lbef2;
  !@VE!@M0 = laef3 - lbef3;
  !@VE!@S0 = iaef1 - ibef1;
  !@VE!@S0 = iaef2 - ibef2;
  !@VE!@S0 = iaef3 - ibef3;
  !@VE!@E0 = eaef1 - ebef1;
  !@VE!@E0 = eaef2 - ebef2;
  !@VE!@E0 = eaef3 - ebef3;
  !@VM!@M0 = lamsi1 - lbmsi1;
  !@VM!@M0 = lamsi2 - lbmsi2;
  !@VM!@M0 = lamsi3 - lbmsi3;
  !@VM!@S0 = iamsi1 - ibmsi1;
  !@VM!@S0 = iamsi2 - ibmsi2;
  !@VM!@S0 = iamsi3 - ibmsi3;
  !@VM!@E0 = eamsi1 - ebmsi1;
  !@VM!@E0 = eamsi2 - ebmsi2;
  !@VM!@E0 = eamsi3 - ebmsi3;
  ",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)


sg_cfa_mitime_msi_list <- list(
  sg_cfa_mitime_config_msi = create_mi_syntax_fun(sg_cfa_mitime, "!@VM", "!@CM{0,1}"),
  sg_cfa_mitime_metric_msi = create_mi_syntax_fun(sg_cfa_mitime, "!@VM", "!@C*M"),
  sg_cfa_mitime_scalar_msi = create_mi_syntax_fun(sg_cfa_mitime, "!@VM", "!@(M|S)"),
  sg_cfa_mitime_strict_msi = create_mi_syntax_fun(sg_cfa_mitime, "!@VM", "!@(M|S|E)"),
  sg_cfa_mitime_config_eff = create_mi_syntax_fun(sg_cfa_mitime, "!@VE", "!@CM{0,1}"),
  sg_cfa_mitime_metric_eff = create_mi_syntax_fun(sg_cfa_mitime, "!@VE", "!@C*M"),
  sg_cfa_mitime_scalar_eff = create_mi_syntax_fun(sg_cfa_mitime, "!@VE", "!@(M|S)"),
  sg_cfa_mitime_strict_eff = create_mi_syntax_fun(sg_cfa_mitime, "!@VE", "!@(M|S|E)"),
  sg_cfa_mitime_config_ful = create_mi_syntax_fun(sg_cfa_mitime, "!@VE|!@VM", "!@CM{0,1}"),
  sg_cfa_mitime_metric_ful = create_mi_syntax_fun(sg_cfa_mitime, "!@VE|!@VM", "!@C*M"),
  sg_cfa_mitime_scalar_ful = create_mi_syntax_fun(sg_cfa_mitime, "!@VE|!@VM", "!@(M|S)"),
  sg_cfa_mitime_strict_ful = create_mi_syntax_fun(sg_cfa_mitime, "!@VE|!@VM", "!@(M|S|E)")
)


sg_cfa_mitime_fit <- purrr::map2(
  sg_cfa_mitime_msi_list,
  names(sg_cfa_mitime_msi_list),
  calculate_mplus,
  dir = "cfa"
)


sg_cfa_mitime_sums <- purrr::map(sg_cfa_mitime_fit, ~ .x$results$summaries) |>
  dplyr::bind_rows(.id = "Model") |>
  dplyr::mutate(
    IV = toupper(stringr::str_remove(Model, "^sg_cfa_mitime_[a-z]*_")),
    Model = paste(
      toupper(stringr::str_remove_all(Model, "^sg_cfa_mitime_|_[a-z]{3}(_[a-z]{2}){0,1}$")),
      "MODEL"
    )
  ) |>
  dplyr::arrange(IV, Model) |>
  dplyr::mutate(
    .by = c(IV),
    dRMSEA = c(NA, abs(diff(RMSEA_Estimate))),
    dCFI = c(NA, abs(diff(CFI))),
    dSRMR = c(NA, abs(diff(SRMR)))
  ) |>
  dplyr::rename(RMSEA = RMSEA_Estimate) |>
  dplyr::select(
    IV,
    Model,
    CFI, TLI, RMSEA, SRMR,
    dCFI, dRMSEA, dSRMR
  )

sg_cfa_mitime_sums |>
  knitr::kable()

rm(list = ls(envir = .GlobalEnv)[grepl("sg_cfa", ls(envir = .GlobalEnv))])

## SEM MODELS ####

# - M1: Main Effects of DSI and CON on Effort: 
#       DSI -> EFF 
#       CON -> EFF
# - M2: Interaction of DSI and CON on Effort
#       DSI -> EFF 
#       CON -> EFF 
#       CON x DSI -> EFF
# - M3: Main Effects of MSI and CON on Effort (check)
#       MSI -> EFF 
#       CON -> EFF
# - M4: Interaction of MSI and CON on Effort (check)
#       MSI -> EFF 
#       CON -> EFF 
#       CON x MSI -> EFF
# - M5: Mediation: Main Effects of DSI, MSI and CON on Effort (check)
#       DSI -> EFF 
#       MSI -> EFF 
#       CON -> EFF
#       DSI -> MSI
# - M6: Interaction of DSI and CON and Main Effects of MSI on Effort (easy)
#       DSI -> EFF 
#       CON -> EFF 
#       CON x DSI -> EFF
#       MSI -> EFF
#       DSI -> MSI
# - M7: Interaction of MSI and CON and Main Effects of DSI on Effort (check)
#       MSI -> EFF 
#       CON -> EFF 
#       CON x MSI -> EFF
#       DSI -> EFF
#       DSI -> MSI
# - M8: Interaction of MSI/DSI and CON on Effort (Mediation) (check)
#       MSI -> EFF 
#       CON -> EFF 
#       CON x MSI -> EFF
#       DSI -> EFF
#       DSI -> MSI
#       CON x DSI -> EFF

### M1: Main Effects of DSI and CON on Effort ####
sg_sem_m1_iv <- MplusAutomation::mplusObject(
  TITLE = "Models for IV
    - M1: Main Effects of DSI and CON on EffortH
    - M2: Interaction of DSI and CON on Effort",
  VARIABLE = "
  IDVARIABLE = ID;
  CLUSTER = CLASSID;
  USEVAR =
    BF_C_1 - BF_C_4
    AB_EF_1 - AB_EF_3
    BB_EF_1 - BB_EF_3
    IV_GE_1 - IV_GE_5
    GRADE SCHOOL
    KFT_WLE SEX
    SEQ KEY;",
  DEFINE = "
  STANDARDIZE
    !AB_ACH BB_ACH
    GRADE
    KFT_WLE;",
  ANALYSIS = "
  PROCESSORS = 16;
  TYPE = COMPLEX;!@ RANDOM;
  !@ALGORITHM = INTEGRATION;
  !@MITERATIONS = 5000;
  ESTIMATOR = MLR;
  MODEL = NOCOVARIANCES;
  ITERATIONS = 200000;",
  MODEL = "
  ! Conscientiousness
  CON BY BF_C_1*
    BF_C_2 - BF_C_4;

  CON@1;

  ! Interest
  IV_GE BY IV_GE_1*
    IV_GE_2 - IV_GE_5;

  IV_GE@1;

  ! Effort
  AB_EF BY AB_EF_1* (laef1)
    AB_EF_2 - AB_EF_3 (laef2 laef3);

  [AB_EF_1 - AB_EF_3] (iaef1 - iaef3);
  AB_EF_1 - AB_EF_3 (eaef1 - eaef3);

  [AB_EF@0];
  AB_EF@1;

  BB_EF BY BB_EF_1* (lbef1)
    BB_EF_2 - BB_EF_3 (lbef2 lbef3);

  [BB_EF_1 - BB_EF_3] (ibef1 - ibef3);
  BB_EF_1 - BB_EF_3 (ebef1 - ebef3);

  [BB_EF];
  BB_EF;

  AB_EF_1 - AB_EF_3 PWITH BB_EF_1 - BB_EF_3;
  
  ! Interaction
  !@CONXINT | CON XWITH IV_GE;

  ! Correlations
  CON WITH IV_GE;
  
  !Regression
  !@AB_EF ON CONXINT (raexi);
  AB_EF ON IV_GE CON (raeg raec);

  !@BB_EF ON CONXINT (rbexi);
  BB_EF ON IV_GE CON (rbeg rbec);
  
  BB_EF WITH AB_EF;

  ! Covariates
  AB_EF ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (racov1-racov6);
  BB_EF ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (rbcov1-rbcov6);
  CON ON GRADE SCHOOL KFT_WLE SEX SEQ KEY;
  IV_GE ON GRADE SCHOOL KFT_WLE SEX SEQ KEY;
  !GRADE KFT_WLE WITH GRADE KFT_WLE;
  [KEY SEX];
  ",
  MODELCONSTRAINT = "
  ! factor loadings;
  ! time invariance
  0 = laef1 - lbef1;
  0 = laef2 - lbef2;
  0 = laef3 - lbef3;
  0 = iaef1 - ibef1;
  0 = iaef2 - ibef2;
  0 = iaef3 - ibef3;
  0 = eaef1 - ebef1;
  0 = eaef2 - ebef2;
  0 = eaef3 - ebef3;
  ",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)


eq_test <- c("rbeg", "rbec", paste0("rbcov", 1:5))
eq_test <- c(NA, paste0("0 = ", eq_test, " - ", stringr::str_replace(eq_test, "^([a-z])b", "\\1a"), ";"))


fit_sg_sem_m1_iv_eq_test <- purrr::map2(eq_test, 1:length(eq_test), \(eq, no, model) {
  if (!is.na(eq)) {
    model$MODELTEST <- eq
  } else {
    model$MODELTEST <- NULL
  }
  MplusAutomation::mplusModeler(
    object = model,
    dataout = "../data/data.dat",
    modelout = here::here(paste0("mplus/sem/sg_sem_m1_iv_eq_test_t", no - 1, ".inp")),
    writeData = "never",
    check = FALSE,
    run = TRUE,
    hashfilename = FALSE
  )
}, model = sg_sem_m1_iv)

fit_sg_sem_m1_iv_eq_test_sum <- purrr::map(fit_sg_sem_m1_iv_eq_test, \(x) {
  if (purrr::is_null(x$MODELTEST)) {
    x$MODELTEST <- NA
  }
  data.frame(x$results$summaries, TEST = x$MODELTEST)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Model = stringr::str_remove_all(Filename, "_t1{0,1}[0-9]\\.out$")) |>
  dplyr::select(Filename, TEST, CFI, TLI, RMSEA_Estimate, SRMR, WaldChiSq_PValue)

fit_sg_sem_m1_iv_eq_test_sum |>
  dplyr::filter(WaldChiSq_PValue < 0.05) |>
  knitr::kable()


sg_sem_m1_iv_eq_test_constraints <- fit_sg_sem_m1_iv_eq_test_sum |>
  dplyr::filter(WaldChiSq_PValue > 0.05) |>
  dplyr::summarise(Constraint = paste0(TEST, collapse = "\n"))


if (sg_sem_m1_iv_eq_test_constraints$Constraint != "") {
  sg_sem_m1_iv$MODELCONSTRAINT <- paste0(
    sg_sem_m1_iv$MODELCONSTRAINT, 
    sg_sem_m1_iv_eq_test_constraints$Constraint, 
    collapse = "\n")
}

fit_sg_sem_m1_iv <- calculate_mplus(sg_sem_m1_iv, "sg_sem_m1_iv", dir = "sem")


### M2: Interaction of MSI and CON on Effort ####

sg_sem_m2_iv_int <- purrr::map2(sg_sem_m1_iv, names(sg_sem_m1_iv), modify_mplus, int = TRUE)


# Test Interaction to be equal across time
sg_sem_m2_iv_int_rsex <- sg_sem_m2_iv_int
sg_sem_m2_iv_int_rsex$MODELTEST <- "0 = rbexi - raexi;"
sg_sem_m2_iv_int_rsex$MODELINDIRECT <- NULL

fit_sg_sem_m2_iv_int_rsex <- calculate_mplus(sg_sem_m2_iv_int_rsex, "sg_sem_m2_iv_int_rsex", "sem")


# Constraints
sg_sem_eq_test_sum <- list(
  data.frame(
    fit_sg_sem_m2_iv_int_rsex$results$summaries,
    Test = fit_sg_sem_m2_iv_int_rsex$MODELTEST)) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Model = stringr::str_remove_all(Filename, "_rb[ae]x\\.out$")) |>
  dplyr::select(Filename, Model, Test, WaldChiSq_PValue) |>
  dplyr::filter(WaldChiSq_PValue > 0.05)

if (any(sg_sem_eq_test_sum$WaldChiSq_PValue > 0.05)) {
  sg_sem_m2_iv_int$MODELCONSTRAINT <- paste0(
    sg_sem_m2_iv_int$MODELCONSTRAINT,
    paste0(sg_sem_eq_test_sum$Test, collapse = "\n"),
    collapse = "\n")
}

fit_sg_sem_m2_iv_int <- calculate_mplus(sg_sem_m2_iv_int, "sg_sem_m2_iv_int", dir = "sem")

diff_fun(fit_sg_sem_m1_iv, fit_sg_sem_m2_iv_int_rsex, corrected = FALSE) # no
diff_fun(fit_sg_sem_m1_iv, fit_sg_sem_m2_iv_int, corrected = FALSE) # no

##########################################################################
### M3: Main Effects of MSI and CON on Effort ####
sg_sem_m3_msi <- MplusAutomation::mplusObject(
  TITLE = "Models for MSI
    - M0: Main Effects of MSI and CON on Effort/ACH
    - M1: Interaction of MSI and CON on Effort/ACH
    - m0: Main Effects of MSI and CON and IV on Effort/ACH (Mediation)
    - M3: Interaction of MSI and CON  and IV on Effort (Mediation)
    - M4: Interaction of MSI/IV and CON on Effort (Mediation)",
  VARIABLE = "
  IDVARIABLE = ID;
  CLUSTER = CLASSID;
  USEVAR =
    AA_MSI_2 - AA_MSI_6
    BA_MSI_2 - BA_MSI_6
    BF_C_1 - BF_C_4
    AB_EF_1 - AB_EF_3
    BB_EF_1 - BB_EF_3
    !AB_ACH BB_ACH
    !@@{IV}_1 - {IV}_5
    GRADE SCHOOL
    KFT_WLE SEX
    SEQ KEY;",
  DEFINE = "
  STANDARDIZE
    !AB_ACH BB_ACH
    GRADE
    KFT_WLE;",
  ANALYSIS = "
  PROCESSORS = 16;
  TYPE = COMPLEX;!@ RANDOM;
  !@ALGORITHM = INTEGRATION;
  !@MITERATIONS = 5000;
  !@!@!@@INTEGRATION = MONTECARLO;
  ESTIMATOR = MLR;
  MODEL = NOCOVARIANCES;
  ITERATIONS = 200000;",
  MODEL = "
  ! Conscientiousness
  CON BY BF_C_1*
    BF_C_2 - BF_C_4;

  CON@1;

  ! Interest
  AA_MSI BY AA_MSI_2* (lamsi1)
    AA_MSI_4 - AA_MSI_6 (lamsi2 lamsi3);

  [AA_MSI_2 - AA_MSI_6] (iamsi1 - iamsi3);
  AA_MSI_2 - AA_MSI_6 (eamsi1 - eamsi3);

  [AA_MSI@0];
  AA_MSI@1 (vamsi);

  BA_MSI BY BA_MSI_2* (lbmsi1)
    BA_MSI_4 - BA_MSI_6 (lbmsi2 lbmsi3);

  [BA_MSI_2 - BA_MSI_6] (ibmsi1 - ibmsi3);
  BA_MSI_2 - BA_MSI_6 (ebmsi1 - ebmsi3);

  [BA_MSI];
  BA_MSI (vbmsi);

  AA_MSI_2 - AA_MSI_6 PWITH BA_MSI_2 - BA_MSI_6;

  ! Interest
  !@@{IV} BY {IV}_1*
    !@@{IV}_2 - {IV}_5;

  !@@{IV}@1 (viv);

  ! Effort
  AB_EF BY AB_EF_1* (laef1)
    AB_EF_2 - AB_EF_3 (laef2 laef3);

  [AB_EF_1 - AB_EF_3] (iaef1 - iaef3);
  AB_EF_1 - AB_EF_3 (eaef1 - eaef3);

  [AB_EF@0];
  AB_EF@1 (vaef);

  BB_EF BY BB_EF_1* (lbef1)
    BB_EF_2 - BB_EF_3 (lbef2 lbef3);

  [BB_EF_1 - BB_EF_3] (ibef1 - ibef3);
  BB_EF_1 - BB_EF_3 (ebef1 - ebef3);

  [BB_EF];
  BB_EF (vbef);

  AB_EF_1 - AB_EF_3 PWITH BB_EF_1 - BB_EF_3;
  
  ! Achievement
  !AB_A BY AB_ACH@1;
  !AB_A@1;
  !AB_ACH@0;

  !BB_A BY BB_ACH@1;
  !BB_A@1;
  !BB_ACH@0;

  ! Interaction
  !@ACONXINT | CON XWITH AA_MSI;

  ! Interaction
  !@BCONXINT | CON XWITH BA_MSI;
  
  ! Interaction
  !@!@!@@CONXINT | CON XWITH {IV};

  ! Correlations
  CON WITH AA_MSI BA_MSI (cacov cbcov);
  !@@CON WITH {IV};
  
  !AB_A WITH AB_EF (aacov);
  !BB_A WITH BB_EF (abcov);

  !Regression
  AB_EF ON AA_MSI CON (raei raec);
  !AB_A ON AA_MSI CON (raai raac);
  !@AB_EF ON ACONXINT (raex);
  !!@AB_A ON ACONXINT (raax);
  !@!@!@@AB_EF ON CONXINT (raexi);
  !!@!@!@@AB_A ON CONXINT (raexi);
  !@@AB_EF ON {IV} (raeg);
  !!@@AB_A ON {IV} (raag);

  BB_EF ON BA_MSI CON (rbei rbec);
  !BB_A ON BA_MSI CON (rbai rbac);
  !@BB_EF ON BCONXINT (rbex);
  !!@BB_A ON BCONXINT (rbax);
  !@!@!@@BB_EF ON CONXINT (rbexi);
  !!@!@!@@BB_A ON CONXINT (rbexi);
  !@@BB_EF ON {IV} (rbeg);
  !!@@BB_A ON {IV} (rbag);
  
  !BB_EF BB_A ON AB_EF AB_A;
  BB_EF WITH AB_EF;
  BA_MSI WITH AA_MSI;
  BB_EF WITH AA_MSI;

  !@@AA_MSI BA_MSI ON {IV} (rami rbmi);

  ! Covariates
  AB_EF AA_MSI ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (racov1-racov12);
  BB_EF BA_MSI ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (rbcov1-rbcov12);
  CON ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (rccov1-rccov6);
  !@@{IV} ON GRADE SCHOOL KFT_WLE SEX SEQ KEY (ricov1-ricov6);
  GRADE KFT_WLE WITH GRADE KFT_WLE;
  !GRADE SCHOOL KFT_WLE SEX SEQ KEY (vcov1-vcov6);
  [KEY SEX];
  ",
  MODELCONSTRAINT = "
  ! factor loadings;
  !laef1 = 3 - (laef2 + laef3);
  !lamsi1 = 3 - (lamsi2 + lamsi3);

  ! intercepts;
  !iaef1 = 0 - (iaef2 + iaef3);
  !iamsi1 = 0 - (iamsi2 + iamsi3);

  ! time invariance
  0 = laef1 - lbef1;
  0 = laef2 - lbef2;
  0 = laef3 - lbef3;
  0 = iaef1 - ibef1;
  0 = iaef2 - ibef2;
  0 = iaef3 - ibef3;
  0 = eaef1 - ebef1;
  0 = eaef2 - ebef2;
  0 = eaef3 - ebef3;
  0 = lamsi1 - lbmsi1;
  0 = lamsi2 - lbmsi2;
  0 = lamsi3 - lbmsi3;
  0 = iamsi1 - ibmsi1;
  0 = iamsi2 - ibmsi2;
  0 = iamsi3 - ibmsi3;
  0 = eamsi1 - ebmsi1;
  0 = eamsi2 - ebmsi2;
  0 = eamsi3 - ebmsi3;
  ",
  OUTPUT = "STANDARDIZED;",
  rdata = dat,
  usevariables = names(dat)
)


sg_sem_m3_msi_eq_test <- sg_sem_m3_msi
sg_sem_m3_msi_eq_test <- purrr::map(sg_sem_m3_msi_eq_test, \(x) {
  if (is.vector(x) && length(x) == 1) {
    if (is.character(x)) {
    x <- stringr::str_glue(x, IV = "IV_GE")
    }
  }
  x
})
class(sg_sem_m3_msi_eq_test) <- c("mplusObject", "list")

eq_test <- c("rbei", "rbec", paste0("rbcov", 1:12), "cbcov")
eq_test <- c(NA, paste0("0 = ", eq_test, " - ", stringr::str_replace(eq_test, "^([a-z])b", "\\1a"), ";"))


fit_sg_sem_m3_msi_eq_test <- purrr::map2(eq_test, 1:length(eq_test), \(eq, no, model) {
  if (!is.na(eq)) {
    model$MODELTEST <- eq
  } else {
    model$MODELTEST <- NULL
  }
  MplusAutomation::mplusModeler(
    object = model,
    dataout = "../data/data.dat",
    modelout = here::here(paste0("mplus/sem/sg_sem_m3_msi_eq_test_t", no - 1, ".inp")),
    writeData = "never",
    check = FALSE,
    run = TRUE,
    hashfilename = FALSE
  )
}, model = sg_sem_m3_msi_eq_test)

fit_sg_sem_m3_msi_eq_test_sum <- purrr::map(fit_sg_sem_m3_msi_eq_test, \(x) {
  if (purrr::is_null(x$MODELTEST)) {
    x$MODELTEST <- NA
  }
  data.frame(x$results$summaries, TEST = x$MODELTEST)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Model = stringr::str_remove_all(Filename, "_t1{0,1}[0-9]\\.out$")) |>
  dplyr::select(Filename, TEST, CFI, TLI, RMSEA_Estimate, SRMR, WaldChiSq_PValue)

fit_sg_sem_m3_msi_eq_test_sum |>
  dplyr::filter(WaldChiSq_PValue < 0.05) |>
  knitr::kable()


sg_sem_m3_msi_eq_test_constraints <- fit_sg_sem_m3_msi_eq_test_sum |>
  dplyr::filter(WaldChiSq_PValue > 0.05) |>
  dplyr::summarise(Constraint = paste0(TEST, collapse = "\n"))


if (sg_sem_m3_msi_eq_test_constraints$Constraint != "") {
  sg_sem_m3_msi$MODELCONSTRAINT <- paste0(
    sg_sem_m3_msi$MODELCONSTRAINT, 
    sg_sem_m3_msi_eq_test_constraints$Constraint, 
    collapse = "\n")
}

fit_sg_sem_m3_msi <- calculate_mplus(sg_sem_m3_msi, "sg_sem_m3_msi", dir = "sem")


### M4: Interaction of MSI and CON on Effort ####

sg_sem_m4_msi_int <- purrr::map2(sg_sem_m3_msi, names(sg_sem_m3_msi), modify_mplus, int = TRUE)


# Test Interaction to be equal across time
sg_sem_m4_msi_int_rsex <- sg_sem_m4_msi_int
sg_sem_m4_msi_int_rsex$MODELTEST <- "0 = rbex - raex;"
sg_sem_m4_msi_int_rsex$MODELINDIRECT <- NULL

fit_sg_sem_m4_msi_int_rsex <- calculate_mplus(sg_sem_m4_msi_int_rsex, "sg_sem_m4_msi_int_rsex", "sem")


# Constraints
sg_sem_eq_test_sum <- list(
  data.frame(
    fit_sg_sem_m4_msi_int_rsex$results$summaries,
    Test = fit_sg_sem_m4_msi_int_rsex$MODELTEST)) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Model = stringr::str_remove_all(Filename, "_rb[ae]x\\.out$")) |>
  dplyr::select(Filename, Model, Test, WaldChiSq_PValue) |>
  dplyr::filter(WaldChiSq_PValue > 0.05)

if (any(sg_sem_eq_test_sum$WaldChiSq_PValue > 0.05)) {
  sg_sem_m4_msi_int$MODELCONSTRAINT <- paste0(
    sg_sem_m4_msi_int$MODELCONSTRAINT, "\n",
    paste0(sg_sem_eq_test_sum$Test, collapse = "\n"),
    collapse = "\n")
}

fit_sg_sem_m4_msi_int <- calculate_mplus(sg_sem_m4_msi_int, "sg_sem_m4_msi_int", dir = "sem")



### M5: Main Effects of MSI and CON on Effort and IV (Mediation) ####

sg_sem_m5_msi_iv <- sg_sem_m3_msi %>%
  purrr::map2(.x = ., .y = names(.), ~ modify_mplus(.x, name = .y, old = "!@@", new = ""))

sg_sem_m5_msi_iv$VARIABLE <- stringr::str_glue(sg_sem_m5_msi_iv$VARIABLE, IV = "IV_GE")
sg_sem_m5_msi_iv$MODEL <- stringr::str_glue(sg_sem_m5_msi_iv$MODEL, IV = "IV_GE")

sg_sem_m5_msi_iv_rami <- sg_sem_m5_msi_iv
sg_sem_m5_msi_iv_rami$MODELTEST <- "0 = rami - rbmi;"
fit_sg_sem_m5_msi_iv_rami <- calculate_mplus(sg_sem_m5_msi_iv_rami, "sg_sem_m5_msi_iv_rami", dir = "sem")
fit_sg_sem_m5_msi_iv_rami$results$summaries$WaldChiSq_PValue

sg_sem_m5_msi_iv_raeg <- sg_sem_m5_msi_iv
sg_sem_m5_msi_iv_raeg$MODELTEST <- "0 = raeg - rbeg;"
fit_sg_sem_m5_msi_iv_raeg <- calculate_mplus(sg_sem_m5_msi_iv_raeg, "sg_sem_m5_msi_iv_raeg", dir = "sem")
fit_sg_sem_m5_msi_iv_raeg$results$summaries$WaldChiSq_PValue


# indirect path
const_ind_std <- "
NEW (iv_var con_var amsi_var bmsi_var 
     aef_var bef_var a_in_std a_to_std 
     b_in_std b_to_std);

iv_var =  ricov1**2*vcov1 + 
          ricov2**2*vcov2 +
          ricov3**2*vcov3 +
          ricov4**2*vcov4 +
          ricov5**2*vcov5 + 
          ricov6**2*vcov6 +
          1;


con_var = rccov1**2*vcov1 + 
          rccov2**2*vcov2 +
          rccov3**2*vcov3 +
          rccov4**2*vcov4 +
          rccov5**2*vcov5 + 
          rccov6**2*vcov6 +
          1;

amsi_var = racov6**2*vcov1 + 
           racov7**2*vcov2 +
           racov8**2*vcov3 +
           racov9**2*vcov4 +
           racov10**2*vcov5 + 
           rami**2*iv_var +
           1;

bmsi_var = rbcov6**2*vcov1 + 
           rbcov7**2*vcov2 +
           rbcov8**2*vcov3 +
           rbcov9**2*vcov4 +
           rbcov10**2*vcov5 + 
           rbmi**2*iv_var +
           vbmsi;

         
aef_var = racov1**2*vcov1 +
          racov2**2*vcov2 +
          racov3**2*vcov3 +
          racov4**2*vcov4 +
          racov5**2*vcov5 +
          raei**2*amsi_var +
          raec**2*con_var + 
          raeg**2*iv_var +
          1;    


bef_var = rbcov1**2*vcov1 +
          rbcov2**2*vcov2 +
          rbcov3**2*vcov3 +
          rbcov4**2*vcov4 +
          rbcov5**2*vcov5 +
          rbei**2*bmsi_var +
          rbec**2*con_var + 
          rbeg**2*iv_var +
          vbef;      


a_in_std = rami * raei * 1/sqrt(aef_var);
b_in_std = rbmi * rbei * 1/sqrt(bef_var);
a_to_std = a_in_std + raeg * 1/sqrt(aef_var);
b_to_std = b_in_std + rbeg * 1/sqrt(bef_var);
"

const_ind <- "
NEW (a_in a_to b_in b_to);

a_in = rami * raei;
b_in = rbmi * rbei;
a_to = a_in + raeg;
b_to = b_in + rbeg;
"

sg_sem_m5_msi_iv$MODELINDIRECT <- "AB_EF ind IV_GE;\nBB_EF ind IV_GE"
sg_sem_m5_msi_iv$MODELCONSTRAINT <- paste0(sg_sem_m5_msi_iv$MODELCONSTRAINT, "\n0 = raeg - rbeg;")
fit_sg_sem_m5_msi_iv <- calculate_mplus(sg_sem_m5_msi_iv, "sg_sem_m5_msi_iv", dir = "sem")


fit_sg_sem_m5_msi_iv$results$indirect$stdyx.standardized$overall |>
  knitr::kable()


### M6: Interaction of DSI and CON and Main Effects of MSI on Effort ####

sg_sem_m6_msi_iv_int <- sg_sem_m5_msi_iv
sg_sem_m6_msi_iv_int$MODELINDIRECT <- NULL
sg_sem_m6_msi_iv_int$MODEL <- stringr::str_remove_all(sg_sem_m6_msi_iv_int$MODEL, "!@!@")
sg_sem_m6_msi_iv_int$ANALYSIS <- stringr::str_remove_all(sg_sem_m6_msi_iv_int$ANALYSIS, "(;| )!@")

sg_sem_m6_msi_iv_int_raexi <- sg_sem_m6_msi_iv_int
sg_sem_m6_msi_iv_int_raexi$MODELTEST <- "0 = raexi - rbexi;"
fit_sg_sem_m6_msi_iv_int_raexi <- calculate_mplus(sg_sem_m6_msi_iv_int_raexi, "sg_sem_m6_msi_iv_int_raexi", dir = "sem")
fit_sg_sem_m6_msi_iv_int_raexi$results$summaries$WaldChiSq_PValue

sg_sem_m6_msi_iv_int$MODELCONSTRAINT <- paste0(sg_sem_m6_msi_iv_int$MODELCONSTRAINT, "\n0 = raexi - rbexi;")
fit_sg_sem_m6_msi_iv_int <- calculate_mplus(sg_sem_m6_msi_iv_int, "sg_sem_m6_msi_iv_int", dir = "sem")

### M7: Interaction of MSI and CON on Effort and IV (Mediation) ####

sg_sem_m7_msi_int_iv <- sg_sem_m5_msi_iv %>%
  purrr::map2(.x = ., .y = names(.), ~ modify_mplus(.x, name = .y, old = ";{0,1}!@( {0,1})([A-Z])", new = "\\1\\2"))

sg_sem_m7_msi_int_iv$MODELINDIRECT <- NULL
sg_sem_m7_msi_int_iv$MODELCONSTRAINT <- paste0(sg_sem_m7_msi_int_iv$MODELCONSTRAINT, "\n0 = rbex - raex;")
sg_sem_m7_msi_int_iv_std <- sg_sem_m7_msi_int_iv
sg_sem_m7_msi_int_iv$MODELCONSTRAINT <- paste0(sg_sem_m7_msi_int_iv$MODELCONSTRAINT, const_ind)
sg_sem_m7_msi_int_iv_std$MODELCONSTRAINT <- paste0(sg_sem_m7_msi_int_iv_std$MODELCONSTRAINT, const_ind_std)
sg_sem_m7_msi_int_iv_std$MODEL <- stringr::str_replace(
  sg_sem_m7_msi_int_iv_std$MODEL, 
  "!GRADE SCHOOL KFT_WLE SEX SEQ KEY", 
  "GRADE SCHOOL KFT_WLE SEX SEQ KEY"
)
fit_sg_sem_m7_msi_int_iv <- calculate_mplus(sg_sem_m7_msi_int_iv, "sg_sem_m7_msi_int_iv", dir = "sem")
fit_sg_sem_m7_msi_int_iv_std <- calculate_mplus(sg_sem_m7_msi_int_iv_std, "sg_sem_m7_msi_int_iv_std", dir = "sem")

### M8: Interaction of MSI/IV and CON on Effort (Mediation) ####
sg_sem_m8_msi_iv_int2 <- sg_sem_m7_msi_int_iv %>%
  purrr::map2(.x = ., .y = names(.), ~ modify_mplus(.x, name = .y, int = TRUE))

cat(sg_sem_m8_msi_iv_int2$MODEL)
sg_sem_m8_msi_iv_int2$MODELCONSTRAINT <- paste0(sg_sem_m8_msi_iv_int2$MODELCONSTRAINT, "\n0 = rbexi - raexi;")
fit_sg_sem_m8_msi_iv_int2 <- calculate_mplus(sg_sem_m8_msi_iv_int2, "sg_sem_m8_msi_iv_int2", dir = "sem")


### Model Comparisons ####

# fit_m2 <- MplusAutomation::readModels("mplus/sem/sg_sem_msi_m2_iv.out")
# fit_m3 <- MplusAutomation::readModels("mplus/sem/sg_sem_m7_msi_int_iv.out")
# fit_m8 <- MplusAutomation::readModels("mplus/sem/sg_sem_m8_msi_iv_int2.out")
# diff_fun(fit_m2, fit_m3, corrected = FALSE)
# diff_fun(fit_m3, fit_m4, corrected = FALSE)
# Test

fit_m1 <- MplusAutomation::readModels("mplus/sem/sg_sem_m1_iv.out")
fit_m2 <- MplusAutomation::readModels("mplus/sem/sg_sem_m2_iv_int.out")
fit_m3 <- MplusAutomation::readModels("mplus/sem/sg_sem_m3_msi.out")
fit_m4 <- MplusAutomation::readModels("mplus/sem/sg_sem_m4_msi_int.out")
fit_m5 <- MplusAutomation::readModels("mplus/sem/sg_sem_m5_msi_iv.out")
fit_m6 <- MplusAutomation::readModels("mplus/sem/sg_sem_m6_msi_iv_int.out")
fit_m7 <- MplusAutomation::readModels("mplus/sem/sg_sem_m7_msi_int_iv.out")
fit_m8 <- MplusAutomation::readModels("mplus/sem/sg_sem_m8_msi_iv_int2.out")

diff_fun(fit_m1, fit_m2, corrected = FALSE)
diff_fun(fit_m3, fit_m4, corrected = FALSE)
diff_fun(fit_m5, fit_m6, corrected = FALSE)
diff_fun(fit_m5, fit_m7, corrected = FALSE)
#diff_fun(fit_m6, fit_m7, corrected = FALSE)
diff_fun(fit_m7, fit_m8, corrected = FALSE)
diff_fun(fit_m6, fit_m8, corrected = FALSE)
