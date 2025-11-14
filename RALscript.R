#2021 POPPK Model 
###############################
model <- 
  "
$PARAM
F = 1             // Bioavailability fixed in 1
KA = 0.741        // Absorption rate (h⁻¹)
MAT = 0.336       // Mean absortion time fasting (h)
thetaMAT_FED = 1.6     // Change in MAT for fed status (160%)
V1F = 44.3         // Apparent central Vd (L; V1/F)
CLF = 55.8         // apparent Clearance (L/h; CL/F)
thetaCLF_ATV = -0.17    // Change in CL for atazanavir
QF = 5.68          // apparent intercompartimental Clearance  (L/h; Q/F)
V2F = 92.8         // Apparent peripheric Vd (L; V2/F)
thetaFLOWFAT = -0.459  // Change in F for lowfat diet
thetaFPREG = -0.487    // Change in F for pregnancy
thetaFEFV = -0.167     // Change in F for efavirenz
OCC = 1

$PARAM @annotated @covariates
WT : 70 :          // weight (kg)
PREG : 0  :        // pregnancy (0 = N, 1 = Y)
ATV : 0  :         // atazanavir (0 = n, 1 = y)
EFV : 0   :        // efavirenz (0 = n, 1 = y)
FED : 0  :         // fed status (0 = fast, 1 = fed)
LOWFAT : 0 :       // low fat diet (0 = n, 1 = y)

$OMEGA @block @name IIV @labels eta_V1 eta_CL eta_Q eta_V2
0.3959594                                  // IIV_V1 (69.7% CV) 
0        0.07862262                         // IIV_CL (28.6% CV) 
0        0.03243243 0.4129206              // IIV_Q (71.5% CV) 
0        0          0.3484313       0.8446246     // IIV_V2 (115.2% CV)

$OMEGA @name IOV @labels eta_F eta_MAT
0.8138774    // IOV_F (112.1% CV)
1.089916    // IOV_MAT (140.5% CV)

$SIGMA @name ERR @labels eps_prop_leq3h eps_prop_gt3h
0.01 // 0.1733018   // proportional Error t ≤ 3 h (43.5% CV)
0.01 // 0.08075015    // proportional Error t > 3 h (29.0% CV)

$CMT GUT CENTRAL PERIPH TRANS1 TRANS2 TRANS3

$MAIN
double CLi = CLF * pow((WT / 70), 0.75) * (1 + ATV * thetaCLF_ATV) * exp(eta_CL); 
double V1i = V1F * (WT / 70) * exp(eta_V1); 
double Qi = QF * pow((WT / 70), 0.75) * exp(eta_Q);
double V2i = V2F * (WT / 70) * exp(eta_V2);
double MATi = MAT * (1 + FED * thetaMAT_FED) * exp(eta_MAT);
double KTRi = 3 / MATi;
double Fi = (1 * (1 + LOWFAT * thetaFLOWFAT) * (1 + PREG * thetaFPREG) * (1 + EFV * thetaFEFV)) * exp(eta_F);
F_GUT = Fi;

$ODE
dxdt_GUT = -KTRi * GUT;
dxdt_TRANS1 = KTRi * GUT - KTRi * TRANS1;
dxdt_TRANS2 = KTRi * TRANS1 - KTRi * TRANS2;
dxdt_TRANS3 = KTRi * TRANS2 - KA * TRANS3;
dxdt_CENTRAL = KA * TRANS3 - (CLi / V1i) * CENTRAL - (Qi / V1i) * CENTRAL + (Qi / V2i) * PERIPH;
dxdt_PERIPH = (Qi / V1i) * CENTRAL - (Qi / V2i) * PERIPH;

$TABLE
double CP = CENTRAL / V1i; // Concentração plasmática;
double TIME_REL = fmod(TIME, 12); //  12 hour cycle
double DV;

if (TIME_REL <= 3) {
    DV = CP * (1 + eps_prop_leq3h); //  proportional Error  t ≤ 3 hours
} else {
    DV = CP * (1 + eps_prop_gt3h); // proportional Error t > 3 hours
}


$CAPTURE CP DV CLi V1i KTRi Qi V2i Fi MATi;

"
mod <- mcode("final_model", model)
# WT distribution
set.seed(15611)
WT = tibble(ID = 1:4500) %>% 
  mutate(WT = rtruncnorm(n(), a=45, b=120, mean=70, sd=10))

# PREG (0 = no, 1 = yes)
set.seed(4231)
PREG = tibble(ID = 1:4500) %>% 
  mutate(PREG = rbinom(n(), 1, 0.4))  # 15% das mulheres são grávidas

# Atazanavir (ATV)
set.seed(7984)
ATV = tibble(ID = 1:4500) %>% 
  mutate(ATV = rbinom(n(), 1, 0.10))  # 10% usam Atazanavir

# Efavirenz (EFV)
set.seed(54612)
EFV = tibble(ID = 1:4500) %>% 
  mutate(EFV = rbinom(n(), 1, 0.10))  # 10% usam Efavirenz

# --- FED ---
set.seed(4658)
fed_data = tibble(ID = 1:4500) %>%
  mutate(FED = rbinom(n(), 1, 0.50)) # A coluna interna ainda se chama FED

# --- LOWFAT ---
set.seed(9108)
lowfat_data <- fed_data %>% # Usar o tibble renomeado
  mutate(LOWFAT = ifelse(FED == 1, rbinom(n(), 1, 0.50), 0)) %>%
  select(ID, LOWFAT)

#Dose simulation- 400 mg BID e 1200 mg QD

set.seed(4040)
dose_events <- ev(ID = 1:4500, amt = 400, ii = 12, addl = 10)
dose <- as_tibble(dose_events) %>%
  select(ID, time, amt, ii, addl, evid, cmt) %>%
  mutate(dose_group = "400mg BID") %>%
  uncount(addl + 1, .id = "OCC") %>%
  mutate(time = time + (OCC - 1) * ii) %>%  
  select(-ii, -addl) %>%  
  arrange(ID, time)

pac <- left_join(WT, PREG, by = "ID") %>%
  left_join(ATV, by = "ID") %>%
  left_join(EFV, by = "ID") %>%
  left_join(fed_data, by = "ID") %>%  
  left_join(lowfat_data, by = "ID") %>%
  left_join(dose, by = "ID")

#model simulation
set.seed(2)

out <- mod %>%
  data_set(pac) %>%
  Req(CP, DV, CLi, V2i, KTRi, Qi, V1i, Fi, MATi) %>%
  mrgsim(end = 120, delta = 0.5, carry_out = "OCC")

simul_ral <- as_tibble(out) %>%
  left_join(pac %>%
              select(ID, OCC, dose_group, PREG, ATV, EFV, FED, WT, LOWFAT) %>% distinct(ID, OCC, .keep_all = TRUE),
            by = c("ID", "OCC"))

simul_ral2 <- simul_ral %>% filter(time>107)
aucs <- simul_ral2 %>%
  group_by(ID) %>%
  mutate(AUCt = trapz(time, DV))
dataconc <- aucs %>% mutate(Group = "Bukkems et al. (2021)") %>% select(ID, time,Group, DV, AUCt)

simul_ral3 <- aucs %>% filter(time==107.5 | time==109 |time==110 | time==111|time==112 | time==108.5)
pivot_conc <- simul_ral3 %>% 
  select(ID, time, DV) %>% 
  pivot_wider(id_cols = ID, names_from = time, names_prefix = "conc_", values_from = DV)

base_prediction_auc <- simul_ral3 %>% 
  filter(time==107.5) %>% 
  select( -time, -CP) %>% left_join(pivot_conc)

#Remove <1 >99 percentiles
centile_1_auc <- quantile(base_prediction_auc$AUCt, 0.01)
centile_99_auc <- quantile(base_prediction_auc$AUCt, 0.99)

dataconc <- dataconc %>% group_by(ID) %>% 
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc) %>% ungroup()
base_prediction_auc1 <- base_prediction_auc %>%
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc)

# data <- fread("ralsim2021.csv")
# fwrite(base_prediction_auc1, file = "ralsim2021.csv")

#2012 POPPK Model

model <- '
$PARAM
CL_base = 60.2,    // CL/F  (L/h)
V1_base = 223,     // V1/F  (L)
V2_base = 113,    // V2/F (L)
Q_base = 8.5,           // Q/F (L/h)
KA_hiv_pos = 0.21, // absorption rate (h⁻¹) for HIV+
F1_hiv_pos = 0.75, // Bioavailability for HIV+
KA_hiv_neg = 0.65, // absorption rate (h⁻¹) for HIV-
F1_hiv_neg = 1.0,  // Bioavailability for HIV-
theta_female = 0.55, // female sex effect
theta_bilirubin = 0.36, // bilirubin effect
theta_race = 0.59 // race effect
theta_atv = 0.39 // atv effect



$PARAM @annotated @covariates
sex : 1 :   0 = female 1=male
bilirubin : 12 : median umol/liter 
race : 1 : 0=caucasian 1=not caucasian
hiv: 1 : 0=hiv neg 1 = hiv pos
atv : 1 : 0= not atv 1= atv

$OMEGA @name IIV @labels eta_V1 eta_KA eta_F1
0.5675137  // Variance eta_V1 (87.4% CV) log((0.874^2)+1)
0.6351811  // Variance eta_KA (94.2% CV)
0.5605805  // Variance eta_F1 (86.7% CV)

$SIGMA 
0.01 // 0.3074847   // residual variance (60% CV)

$PK
// variables

// --- HIV dependant ---
double KA_actual_base; // base KA
double F1_actual_base; // base F1

if (hiv == 0) { // HIV Ne
  KA_actual_base = KA_hiv_neg; // 0.65
  F1_actual_base = F1_hiv_neg; // 1.0
} else { // HIV Positivo (hiv == 1)
  KA_actual_base = KA_hiv_pos; // 0.21
  F1_actual_base = F1_hiv_pos; // 0.75
}

double KA = KA_actual_base * exp(eta_KA);
double F1 = ((F1_actual_base * (1 + (theta_bilirubin * (bilirubin - 12)/bilirubin))) * (1.55 - theta_female * sex)* (1 + theta_atv * atv)) * exp(eta_F1);
double V1 = (V1_base * (0.41 + theta_race * race)) * exp(eta_V1);
double CL = CL_base ;
double V2 = V2_base ;
double Q = Q_base ;
F_GUT = F1;

$CMT GUT CENT PERIPH

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (CL/V1)*CENT + (Q/V2)*PERIPH - (Q/V1)*CENT;
dxdt_PERIPH = (Q/V1)*CENT - (Q/V2)*PERIPH;

$TABLE
capture CP = (CENT / V1)*(1 + EPS(1));

$CAPTURE CP CL V1 KA F1 eta_V1 eta_KA eta_F1
'

# Compile
mod <- mcode("final_model", model)

i=1
#Sex
set.seed(str_c(12,i))
sex = tibble(ID = 1:4500) %>% 
  mutate(sex = rbinom(n(),1,0.786))

#Race
set.seed(str_c(12,3))
race = tibble(ID = 1:4500) %>% 
  mutate(race = rbinom(n(),1,0.09))
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5377560/

#HIV
set.seed(15744)
hiv = tibble(ID = 1:4500) %>% 
  mutate(hiv = rbinom(n(),1,0.85))

#ATV
set.seed(4442)
atv <- hiv %>%
  mutate(atv = ifelse(hiv == 1, rbinom(n(), 1, 0.09), 0))
set.seed(str_c(12,i))

#Bilirrubin
bilirubin = tibble(ID = 1:4000) %>% 
  mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=12, sd=4)) %>% 
  rbind(tibble(ID = 4001:4450) %>% 
          mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=30, sd=8))) %>% 
  rbind (tibble(ID = 4451:4500) %>% 
           mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=70, sd=12)))

#DOSE 400mg BID in both studies 800 daily
dose <- as_tibble(ev(ID = 1:4500, amt = 400, ii = 12, addl = 10)) %>%
  mutate(dose_group = "400mg BID")
pac <- left_join(sex, race) %>% left_join(bilirubin) %>% left_join(hiv) %>% left_join(atv) %>%  left_join(dose)  

#model simulation
out <- mod %>%
  data_set(pac) %>%
  Req(CP, CL, F1) %>%
  mrgsim(end = 120, delta = 0.5)

simul_ral <- as_tibble(out) %>% 
  left_join(pac %>% 
              select(ID, dose_group, bilirubin, sex, race, hiv,atv), by = c("ID" = "ID"))

simul_ral <- simul_ral %>% 
  mutate(dose = parse_number(dose_group),
         auc = dose/CL)
simul_ral2 <- simul_ral %>% filter(time>107)
aucs <- simul_ral2 %>%
  group_by(ID) %>%
  mutate(AUCt = trapz(time, CP))
dfconc <- aucs %>% mutate(Group = "Arab-Alameddine et al. (2012)") %>% rename(DV = CP) %>% select(ID, time,Group, DV, AUCt)

simul_ral3 <- aucs %>% filter(time==107.5 | time==109 |time==110 | time==111|time==112 | time==108.5)
pivot_conc <- simul_ral3 %>% 
  select(ID, time, CP) %>% 
  pivot_wider(id_cols = ID, names_from = time, names_prefix = "conc_", values_from = CP)

base_prediction_auc <- simul_ral3 %>% 
  filter(time==107.5) %>% 
  select( -time, -CP) %>% left_join(pivot_conc)

#Remove <1 >99 percentiles
centile_1_auc <- quantile(base_prediction_auc$AUCt, 0.01)
centile_99_auc <- quantile(base_prediction_auc$AUCt, 0.99)
dfconc <- dfconc %>% group_by(ID) %>% 
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc) %>% ungroup()
base_prediction_auc1 <- base_prediction_auc %>%
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc)
#   fwrite(base_prediction_auc1, file = "ralsim2012.csv")

# Parallel
n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
data <- as.data.frame(fread("ralsim2021.csv"))

data <- data %>% select(-ID,-dose_group) %>% 
  rename(C0 = conc_107.5, C1 = conc_109, C2 = conc_110, C3 = conc_111,
         C4= conc_112, C05= conc_108.5) 

# Metrics
calculate_metrics <- function(pred_df) {
  if (nrow(pred_df) == 0) return(NULL)
  
  pred_df %>%
    mutate(
      biais_brut = .pred - AUCt,
      bias_rel = (.pred - AUCt) / AUCt,
      bias_rel_square = bias_rel^2,
      biais_brut_sqr = biais_brut^2
    ) %>%
    summarise(
      biais_brut = mean(biais_brut),
      rmse = sqrt(mean(biais_brut_sqr)),
      biais_rel = mean(bias_rel),
      relative_rmse = sqrt(mean(bias_rel_square)),
      biais_out_20percent = mean((bias_rel) > 0.2),
      nb_out_20percent = sum((bias_rel) > 0.2),
      n = n()
    ) %>%
    mutate(
      R2 = summary(lm(AUCt ~ .pred, data = pred_df))$r.squared
    )
}

#2023 data arrangement
df <- as.data.frame(fread("moreira2023.csv"))
df <- df %>% filter(CMT == 2)
df$DV2 <- df$DV*444.4163/1000
df <- df %>% select(ID, TTIME, DV2, tratamento, IDD)
df <- df %>%
  group_by(IDD) %>%
  mutate(AUCt = trapz(TTIME, DV2))
df2conc <- df %>% mutate(Group = "Moreira et al. (2023)") %>% select(-ID) %>%
  rename(DV = DV2, ID = IDD, time = TTIME) %>% select(ID, time,Group, DV, AUCt)
df <- df %>% mutate(TTIME = case_when(
  TTIME == 0 ~ "C0",
  TTIME == 0.5 ~ "C05",
  TTIME == 1 ~ "C1",
  TTIME == 2 ~ "C2",
  TTIME == 3 ~ "C3",
  TTIME == 4 ~ "C4",
  TTIME == 6 ~ "C6",
  TTIME == 8 ~ "C8",
  TTIME == 12 ~ "C12",
  TRUE ~ as.character(TTIME)
))

pivot_df <- df %>% 
  select(IDD, AUCt, TTIME, DV2) %>% 
  pivot_wider(id_cols = IDD, names_from = TTIME, values_from = DV2)

df2 <- df %>% 
  filter(TTIME=="C0") %>% 
  select( -TTIME, -DV2) %>% left_join(pivot_df)
df2 <- as.data.frame(df2)

df <- as.data.frame(fread("ralsim20121710sigma001.csv"))
df <- df %>% rename(C0 = conc_107.5, C1 = conc_109, C2 = conc_110, C3 = conc_111,
                    C4= conc_112, C05= conc_108.5) 


### Combination function
test_combinations <- function(data, time_predictors, fixed_predictors = NULL, response,
                              external_dfs = list(), combo_sizes = c(2, 3),
                              include_fixed_predictors = TRUE,
                              include_ratio_predictors = FALSE, 
                              include_diff_predictors = FALSE,
                              algorithms = c("rf", "xgb", "svm", "glmnet"),
                              return_fits = FALSE,
                              return_workflows = FALSE) {
  
  # Lists
  all_results <- list()
  all_best_params <- list()
  all_tabs <- list()
  all_fits <- list()
  all_workflows <- list()
  
  # Control
  ctrl <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
  race_ctrl <- control_race(save_pred = TRUE)
  
  # Algorithms
  valid_algs <- c("rf", "xgb", "svm", "glmnet")
  algorithms <- intersect(algorithms, valid_algs)
  
  # Derived Predictors
  create_derived_predictors <- function(df, predictors, 
                                        include_ratio, include_diff) {
    new_df <- df
    ratio_names <- c()
    diff_names <- c()
    
    if (include_ratio && length(predictors) >= 2) {
      ratios <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_div_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := ifelse(!!sym(pair[1]) == 0, NA, !!sym(pair[2]) / !!sym(pair[1])))
          name
        })
      ratio_names <- unlist(ratios)
    }
    
    if (include_diff && length(predictors) >= 2) {
      diffs <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_sub_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := !!sym(pair[2]) - !!sym(pair[1]))
          name
        })
      diff_names <- unlist(diffs)
    }
    
    list(data = new_df, ratio_names = ratio_names, diff_names = diff_names)
  }
  
  # Algs specs
  model_specs <- list(
    rf = rand_forest(
      mtry = tune(),
      min_n = tune(),
      trees = tune()
    ) %>% 
      set_engine("ranger", nthread = 10) %>% 
      set_mode("regression"),
    
    xgb = boost_tree(
      mode = "regression",
      mtry = tune(),
      trees = tune(),
      min_n = tune(),
      sample_size = tune(),
      tree_depth = tune(),
      learn_rate = tune()
    ) %>% 
      set_engine("xgboost", nthread = 10),
    
    svm = svm_linear(
      mode = "regression",
      cost = tune(),
      margin = tune()
    ) %>% 
      set_engine("kernlab", nthread = 10),
    
    glmnet = linear_reg(
      penalty = tune(),
      mixture = tune()
    ) %>% 
      set_engine("glmnet", nthread = 10)
  )
  
  # Combination process
  for (n in combo_sizes) {
    size_results <- list()
    size_best_params <- list()
    size_tabs <- list()
    size_fits <- list()
    size_workflows <- list()
    
    # Time combinations
    time_combos <- combn(time_predictors, n, simplify = FALSE)
    
    for (combo in time_combos) {
      # Comb names
      combo_name <- paste(combo, collapse = "_")
      suffix <- ""
      
      # Suffix
      if (include_fixed_predictors && length(fixed_predictors) > 0) {
        suffix <- paste0(suffix, "_fixed(", paste(fixed_predictors, collapse = "+"), ")")
      }
      
      derived <- create_derived_predictors(
        data, combo, include_ratio_predictors, include_diff_predictors
      )
      main_data <- derived$data
      all_predictors <- c(combo, 
                          if(include_fixed_predictors) fixed_predictors,
                          derived$ratio_names, 
                          derived$diff_names)
      
      # Derived predictor suffix
      if (include_ratio_predictors && length(derived$ratio_names) > 0) {
        suffix <- paste0(suffix, "_ratio")
      }
      if (include_diff_predictors && length(derived$diff_names) > 0) {
        suffix <- paste0(suffix, "_diff")
      }
      
      combo_base_name <- paste0(combo_name, suffix)
      cat("\n--- Processando combinação:", combo_base_name, "---\n")
      
      # Data
      df_subset <- main_data %>% 
        select(all_of(c(all_predictors, response))) %>%
        na.omit()
      
      # Train/Test split
      set.seed(123)
      split <- initial_split(df_subset, strata = !!sym(response), prop = 0.75)
      train <- training(split)
      test <- testing(split)
      
      # Recipe
      rec <- recipe(as.formula(paste(response, "~ .")), data = train) %>% 
      step_normalize(all_numeric_predictors()) %>%
        step_zv(all_numeric_predictors())
      
      # Alg selection
      for (alg in algorithms) {
        alg_name <- paste0(combo_base_name,alg)
        cat("  Algoritmo:", alg, "\n")
        
        tryCatch({
          # WF
          wf <- workflow() %>% 
            add_recipe(rec) %>% 
            add_model(model_specs[[alg]])
          
          # CV10F
          set.seed(1234)
          folds <- vfold_cv(train, strata = !!sym(response))
          
          # Tuning
          tune_res <- tune_race_anova(
            wf, 
            resamples = folds, 
            grid = 60,  
            metrics = metric_set(rmse),
            control = race_ctrl
          )
          
          # Best hyperparameters
          best_params <- select_best(tune_res, metric = "rmse")
          final_wf <- finalize_workflow(wf, best_params)
          
          # Fit
          final_fit <- fit(final_wf, train)
          
          # Save WF 
          if (return_fits) {
            size_fits[[alg_name]] <- final_fit
          }
          if (return_workflows) {
            size_workflows[[alg_name]] <- final_wf
          }
          
          # CV Metrics
          cv_res <- fit_resamples(final_wf, folds, control = ctrl)
          cv_preds <- collect_predictions(cv_res)
          cv_metrics <- calculate_metrics(cv_preds) %>% 
            mutate(val = "CV10F", Combinacao = alg_name, Algoritmo = alg)
          
          # Test Metrics
          test_preds <- predict(final_fit, test) %>% bind_cols(test)
          test_metrics <- calculate_metrics(test_preds) %>% 
            mutate(val = "Test", Combinacao = alg_name, Algoritmo = alg)
          
          # External Metrics
          external_metrics <- list()
          for (df_name in names(external_dfs)) {
            df <- external_dfs[[df_name]]
            
            # Derived predictors for external data
            derived_ext <- create_derived_predictors(
              df, combo, include_ratio_predictors, include_diff_predictors
            )
            df_ext <- derived_ext$data
            
            # Verify
            missing_vars <- setdiff(all_predictors, colnames(df_ext))
            if (length(missing_vars) > 0) {
              warning(paste("Variáveis faltantes em", df_name, ":", paste(missing_vars, collapse = ", ")))
              next
            }
            
            preds <- predict(final_fit, df_ext) %>% 
              bind_cols(df_ext %>% select(all_of(response)))
            metrics <- calculate_metrics(preds) %>% 
              mutate(val = df_name, Combinacao = alg_name, Algoritmo = alg)
            
            external_metrics[[df_name]] <- metrics
          }
          
          # Save results
          size_results[[alg_name]] <- cv_metrics$rmse
          size_best_params[[alg_name]] <- best_params
          all_metrics <- c(list(cv_metrics, test_metrics), external_metrics)
          size_tabs[[alg_name]] <- bind_rows(all_metrics)
          
        }, error = function(e) {
          message("Erro no algoritmo ", alg, " para combinação ", combo_base_name, ": ", e$message)
        })
      }
    }
    
    all_results[[paste0("size_", n)]] <- size_results
    all_best_params[[paste0("size_", n)]] <- size_best_params
    all_tabs[[paste0("size_", n)]] <- bind_rows(size_tabs)
    if (return_fits) {
      all_fits[[paste0("size_", n)]] <- size_fits
    }
    if (return_workflows) {
      all_workflows[[paste0("size_", n)]] <- size_workflows
    }
  }
  
  if (return_fits && return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits,
      workflows = all_workflows
    ))
  } else if (return_fits) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits
    ))
  } else if (return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      workflows = all_workflows
    ))
  } else {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs
    ))
  }
}


# Run function
 results <- test_combinations(
  data = data,
  time_predictors = c("C05","C2","C4"),
  fixed_predictors = c(),
  response = "AUCt",
  external_dfs = list(
    "Fernanda" = df2
  ),
  combo_sizes = c(3),
  include_fixed_predictors = F,
  include_ratio_predictors = F,
  include_diff_predictors = F,
  return_fits = T,
  return_workflows = T,
  algorithms = c("xgb"))
final_results <- bind_rows(results$all_tabs)
#fwrite(final_results, file = "models.csv")

final_data <- data %>% select(C05, C2, C4, AUCt)

# Training/testing split with stratification
set.seed(123)
final_split <- initial_split(final_data, strata = AUCt, prop = 0.75)
final_train <- training(final_split)
final_test <- testing(final_split)

set.seed(123)

ral_split <- initial_split(data, 
                           strata = AUCt, 
                           prop= 0.75)
ral_ml_train  <- training(ral_split)
ral_ml_test  <- testing(ral_split)

# Create preprocessing recipe
final_rec <- recipe(AUCt ~ ., data = final_train) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_numeric_predictors())

# Importance plot
ral_ml_train <- ral_ml_train %>% select(C05, C2, C4, AUCt, PREG, ATV, FED, LOWFAT, EFV, WT)
final_rec2 <- recipe(AUCt ~ ., data = ral_ml_train)
ral_ml_rec_prep <-  prep(final_rec2)

final_xgb_spec <- boost_tree(mode = "regression",
                             mtry = tune(),
                             trees = tune(),
                             min_n = tune(),
                             sample_size = tune(),
                             tree_depth = tune(),
                             learn_rate = tune()) %>%
  set_engine("xgboost", nthread = 10)

final_xgb_spec <- finalize_model(final_xgb_spec,results[["best_params_list"]][["size_3"]][["C05_C2_C4xgb"]])

final_xgb_spec %>% set_engine("xgboost", importance = "permutation") %>%
  fit(AUCt ~ ., data = juice(ral_ml_rec_prep)) %>%
  vip::vip(geom = "col") + theme_bw() +  theme(axis.title.x = element_text(size = 20),
                                                axis.text.x = element_text(size = 17),
                                                axis.text.y = element_text(size = 17),
                                                axis.title.y = element_text(size = 20),
                                                legend.text = element_text(size = 17),
                                                plot.title = element_text(size = 16))

final_wf <- results$workflows$size_3$C05_C2_C4xgb

#Save model
final_wf_xgb <- results$fits$size_3$C05_C2_C4xgb
#saveRDS(final_wf_xgb, "model_xgboost2021.rds")

explainer_external <- explain_tidymodels(
  model = final_wf_xgb, 
  data = final_train, 
  y = final_train$AUCt,
  label = "xgb")
#saveRDS(explainer_external, file ="explainer_externalxgb2021.rds" )

# Predictions on sets
# final_predictions <- final_wf_xgb %>% predict(final_train) %>% bind_cols(final_train)

#resample
set.seed(1234)
folds <- vfold_cv(final_train, strata = AUCt)

#10 fold CV
xgb_rs<- fit_resamples(object = final_wf, 
                       resamples = folds, 
                       control = control_resamples(verbose=T, save_pred = T))

xgb_pred_obs <- xgb_rs%>% collect_predictions()
CV10F <- calculate_metrics(xgb_pred_obs)
CV10F$val <- "CV10F"

Train_pred_obs <- xgb_pred_obs %>% 
  left_join(ral_ml_train, by="AUCt") 
Train_pred_obs$PREG <- as.factor(Train_pred_obs$PREG)
Train_pred_obs <- Train_pred_obs %>% mutate(PREG = case_when(
  PREG == 0 ~ "Not Pregnant",
  PREG == 1 ~ "Pregnant"))

Figure_2a <- Train_pred_obs %>%
  ggplot(mapping = aes(x = .pred, y = AUCt, colour = PREG)) + 
  geom_point() +
  geom_smooth(method=lm, color="black") + 
  labs(x="Reference AUC (mg*h/L)", y= "Predicted raltegravir AUC (mg*h/L)",
       color="", title = "A") + 
  theme_bw()+
  theme(axis.title.x = element_text(size = 23),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 23),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 17),
        # Adicionar bordas para melhor visualização
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        # Mover legenda para dentro do gráfico (canto superior direito)
        legend.position = c(0.85, 0.265),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "transparent")
  )
Figure_2a

Figure_2c <- xgb_pred_obs%>%
  ggplot(mapping = aes(x = AUCt, y = AUCt - .pred)) +
  geom_point() + geom_smooth(method=lm) +
  labs(x = "Reference raltegravir AUC (mg*h/L)", y = "Reference - predicted AUC", 
       title="(C)") + 
  theme_bw()

Figure_2c

final_res <- final_wf %>%
  last_fit(final_split,type = "conf_int")

#Performances biais imprecision test
final_res_predictions <- final_res %>% collect_predictions()
ab <- calculate_metrics(final_res_predictions)
ab$val <- "Test"

Test_pred_obs <- final_res_predictions %>% left_join(ral_ml_test, by="AUCt")

Test_pred_obs$PREG <- as.factor(Test_pred_obs$PREG)
Test_pred_obs <- Test_pred_obs %>% mutate(PREG = case_when(
  PREG == 0 ~ "Not Pregnant",
  PREG == 1 ~ "Pregnant"))
Figure_2b <- Test_pred_obs %>%
  ggplot(mapping = aes(x = AUCt, y = .pred, colour = PREG)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (mg*h/L)", y = "Predicted raltegravir AUC (mg*h/L)", 
       color ="", title ="B") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 23),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 23),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 17),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85, 0.265),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "transparent")
  )
Figure_2b

Figure_2d <- final_res_predictions %>%
  ggplot(mapping = aes(x = AUCt, y = AUCt - .pred)) +
  geom_point() +
  geom_smooth(method=lm)+
  labs(x = "Reference raltegravir AUC (mg*h/L)", y = "Test predicted AUC", title="(D)") + 
  theme_bw()
Figure_2d

#2023 External validation organizing, metrics and plots
predictions <- predict(final_wf_xgb, new_data = df2)
df2$AUC_pred <- predictions$.pred
df2$.pred <- predictions$.pred

df2 <- df2 %>% mutate(tratamento = case_when(
  tratamento == 0 ~ "Third trimester",
  tratamento == 1 ~ "Postpartum"))
df2 %>%
  ggplot(mapping = aes(x = AUCt, y = AUC_pred, colour = tratamento)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (mg*h/L)", y = "Predicted raltegravir AUC (mg*h/L)", 
       color="", title = "B") + 
  theme_bw()+
  theme(axis.title.x = element_text(size = 23),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 23),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 17),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85, 0.265),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "transparent")
  ) 

#Plots
df2 %>%
  ggplot(mapping = aes(x = AUCt, y = AUCt - AUC_pred)) +
  geom_point() +
  geom_smooth(method=lm)+
  labs(x = "Reference raltegravir AUC (mg*h/L)", y = "Pred AUC XGB 2 conc 2021 Fernanda", title="(D)") + 
  theme_bw()

ae2 <- calculate_metrics(df2)
ae2$val <- "Fernanda"

#2012 External validation organizing, metrics and plots
predictions <- predict(final_wf_xgb, new_data = df)
df$AUC_pred <- predictions$.pred
df$.pred <- predictions$.pred

be2 <- calculate_metrics(df)
be2$val <- "2012"

df$sex <- as.factor(df$sex)
df <- df %>% mutate(sex = case_when(
  sex == 0 ~ "Female",
  sex == 1 ~ "Male"))

#Plots
df %>%
  ggplot(mapping = aes(x = AUCt, y = AUC_pred, colour = sex)) +
  geom_point() +
  geom_smooth(method=lm,color = "black")+
  labs(x = "Reference AUC (mg*h/L)", y = "Predicted raltegravir AUC (mg*h/L)", 
       color ="", title ="A") +  
  theme_bw() +
  theme(axis.title.x = element_text(size = 23),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 23),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 17),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = c(0.85, 0.265),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "transparent")
  )

df %>%
  ggplot(mapping = aes(x = AUCt, y = AUCt - AUC_pred)) +
  geom_point() +
  geom_smooth(method=lm)+
  labs(x = "Reference raltegravir AUC (mg*h/L)", y = "Pred AUC GLM 3 conc 2021 Fernanda", title="(D)") + 
  theme_bw()

#Final metrics
tab <- rbind(CV10F,ab,ae2,be2)

#Correlation pairs plot
datagraf <- data %>% select(LOWFAT,PREG,ATV,EFV,FED,WT,AUCt,C05,C2,C4)
datagraf <- datagraf %>%
  mutate(across(c(LOWFAT,PREG, ATV, EFV, FED), as.factor))

#table
ral_ml_train  <- training(ral_split)

table1 <- CreateTableOne(vars = c("C0", "C05","C1","C2","C3","C4","AUCt", "WT","PREG", "ATV","EFV","FED","LOWFAT"), factorVars = c("PREG", "ATV","EFV","FED","LOWFAT") ,
                          data = ral_ml_train)

tableOne2b<-print(table1, nonnormal = c("C0", "C05","C1","C2","C3","C4","AUCt", "WT"), printToggle=F, minMax=F)

kableone(tableOne2b)

ral_ml_train  <- training(ral_split)

table1 <- CreateTableOne(factorVars = c("atv", "sex") ,
                         data = df)

tableOne2b<-print(table1, nonnormal = c("C0", "C05","C1","C2","C3","C4","AUCt", "WT"), printToggle=F, minMax=F)
tab1 <- tab %>% select(val,rmse,R2,biais_rel,relative_rmse, biais_out_20percent,nb_out_20percent)
kableone(tableOne2b)
# Id numeric and binary variables
bin_vars <- c("PREG", "ATV", "EFV", "FED", "LOWFAT")  #  binary
num_vars <- c("WT", "AUCt", "C05", "C2", "C4")  # numeric

# ggpairs customization
custom_plot <- function(data, mapping, ...) {
  x_name <- as.character(mapping$x)[2]
  y_name <- as.character(mapping$y)[2]
  
  if ((x_name %in% bin_vars & y_name %in% num_vars) || (x_name %in% num_vars & y_name %in% bin_vars)) {
    ggplot(data, mapping) + geom_boxplot()
  } else if (x_name %in% bin_vars & y_name %in% bin_vars) {
    data %>%
      count(!!sym(x_name), !!sym(y_name)) %>%
      ggplot(aes(x = !!sym(x_name), y = n, fill = as.factor(!!sym(y_name)))) +
      geom_col(position = "dodge") +
      scale_fill_brewer(palette = "Set1") +
      theme_bw() +
      labs(fill = y_name)
  } else {
    ggplot(data, mapping) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = FALSE)
  }
}

# Plot
p1 <- ggpairs(datagraf, 
              columns = c(bin_vars, num_vars),
              lower = list(continuous = wrap("points"), combo = wrap(custom_plot), discrete = wrap(custom_plot)),
              upper = list(continuous = wrap("cor")))

print(p1)

##Correlation Ctrough AUC

lm <- lm((df2$AUCt) ~ (df2$C0))
lm <- lm((df$AUCt) ~ (df$C0))
lm <- lm((data$AUCt) ~ (data$C0))
summary(lm)

ggplot(data, aes(x = `C0`, y = `AUCt`)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", col = "black", se=T, linewidth = 1.5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "A   R² = 0.36, P < 0.01" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_bw(base_size = 22)

ggplot(df, aes(x = `C0`, y = `AUCt`)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", col = "black", se=T, linewidth = 1.5) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "B   R² = 0.55, P < 0.01" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_bw(base_size = 22)

ggplot(df2, aes(x = `C0`, y = `AUCt`)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", col = "black", se=T, linewidth = 1.5) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,21), breaks = c(0,5,10,15,20)) +
  labs(title = "C   R² = -0.04, P > 0.05" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_bw(base_size = 22)

# Data distribution
statcal <- function(df, col_test, col_group, method = "median") {
  df_results <- df %>%
    group_by(across(all_of(col_group))) %>%
    summarise(across(all_of(col_test), 
                     list(estatistica = ~ {
                       if (method == "median") {
                         sprintf("%.3f (%.3f - %.3f)", median(.x, na.rm = TRUE),
                                 quantile(.x, probs = 0.25, na.rm = TRUE),
                                 quantile(.x, probs = 0.75, na.rm = TRUE))
                       } else if (method == "geom") {
                         gmean <- Gmean(.x, na.rm = TRUE, conf.level = 0.95)
                         sprintf("%.3f (%.3f - %.3f)", gmean[1], gmean[2], gmean[3])
                       } else if (method == "mean") {
                         sprintf("%.3f ± %.3f", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE))
                       } else {
                         NA  
                       }
                     }),
                     .names = "{.col}")) %>%
    ungroup()
  
  return(df_results)
}

df$group <- 1
df2$group <- 1
set.seed(123)

ral_split <- initial_split(data, 
                           strata = AUCt, 
                           prop= 0.75)
ral_ml_train  <- training(ral_split)
ral_ml_train$group <- 1
ral_ml_test$group <- 1

# Results

resultsstat2012 <- statcal(df, col_test = c("C0", "C05", "C1", "C2", "C3", "C4", "AUCt"),
                                                col_group = c("group"), method = "median") %>% t()
resultsstattest <- statcal(ral_ml_test, col_test = c("C0","C05","C1", "C2","C3","C4", "AUCt"),
                                                col_group = c("group"), method = "median") %>% t()
resultsstattrain <- statcal(ral_ml_train, col_test = c("C0","C05","C1", "C2","C3","C4", "AUCt"),
                           col_group = c("group"), method = "median") %>% t()
resultsstat2023 <- statcal(df2, col_test = c("C0", "C05", "C1", "C2", "C3", "C4", "AUCt"),
                           col_group = c("group"), method = "median") %>% t()
resusltsstatss <- rbind(resultsstattrain,resultsstattest,resultsstat2023,resultsstat2012)

# Feature selection

data <- as.data.frame(fread("ralsim2021FINAL.csv"))
response <- "AUCt"
data <- data %>% select(-ID,-dose_group) %>% 
  mutate(Diff_C1C0 = conc_61 - conc_59.5, ratioc0c1 = conc_61/conc_59.5) %>% 
  rename(C0 = conc_59.5, C1 = conc_61, C2 = conc_62, C3 = conc_63,
         C4= conc_64, C05= conc_60.5)
concs <- data %>% select(C0,C05,C1,C2,C3,C4)
covs <- data %>%  select(LOWFAT, FED, PREG, EFV, ATV, WT)

#Ratios and subtractions calculations

combos <- combn(concs, 2, simplify = FALSE)
col_names <- names(concs)
name_combos <- combn(col_names, 2, simplify = FALSE)
ratio_list <- map(name_combos, function(pair) {
  col_name1 <- pair[1] 
  col_name2 <- pair[2]
  
  
  vector1 <- concs[[col_name1]]
  vector2 <- concs[[col_name2]]
  
  ratio_result <- vector2 / vector1
  
  return(ratio_result)
})


ratio_names <- map_chr(name_combos, function(pair) {
  col_name1 <- pair[1]
  col_name2 <- pair[2]
  paste0(col_name2, "/", col_name1)
})

names(ratio_list) <- ratio_names
ratio_df <- as_tibble(ratio_list)

sub_list <- map(name_combos, function(pair) {
  col_name1 <- pair[1] 
  col_name2 <- pair[2]
  
  
  vector1 <- concs[[col_name1]]
  vector2 <- concs[[col_name2]]
  
  sub_result <- vector2 - vector1
  
  return(sub_result)
})

sub_names <- map_chr(name_combos, function(pair) {
  col_name1 <- pair[1]
  col_name2 <- pair[2]
  paste0(col_name2, "-", col_name1)
})

names(sub_list) <- sub_names
sub_df <- as_tibble(sub_list)

AUCt <- data$AUCt

final_data <- cbind(ratio_df, concs, AUCt,covs, sub_df)

# Random Forest

rf_model <- rand_forest(mode = "regression", trees = 500) %>%
  set_engine("ranger", importance = "permutation")

rf_fit <- rf_model  %>% 
  fit(AUCt ~ ., data = final_data)
vip(rf_fit)

# GLMNet

glmnet_model <- linear_reg(mode = "regression", mixture = 1, penalty = 1) %>%
  set_engine("glmnet")
glm_net <- glmnet_model  %>% 
  fit(AUCt ~ ., data = final_data)
vip(glm_net)

# Boruta

boruta_model <- Boruta(AUCt ~ ., data = final_data, doTrace = 2)
stats_boruta <- attStats(boruta_model)

# Cattaneo external validation

df2$.pred <- 2.082*df2$C0 + 0.821*df2$C1 + 1.238*df2$C2 - 0.210*df2$C3 + 4.280*df2$C4 + 0.7831
data$.pred <- 2.082*data$C0 + 0.821*data$C1 + 1.238*data$C2 - 0.210*data$C3 + 4.280*data$C4 + 0.7831
df$.pred <- 2.082*df$C0 + 0.821*df$C1 + 1.238*df$C2 - 0.210*df$C3 + 4.280*df$C4 + 0.7831
aa1 <- calculate_metrics(data)
aa2 <-calculate_metrics(df)
aa3 <-calculate_metrics(df2)
aa4 <- rbind(aa1,aa2,aa3)

data$.pred <-0.876*data$C1 + 1.201*data$C2 + 3.994*data$C4 + 1.8242 
df$.pred<- 0.876*df$C1 + 1.201*df$C2 + 3.994*df$C4 + 1.8242 
df2$.pred<- 0.876*df2$C1 + 1.201*df2$C2 + 3.994*df2$C4 + 1.8242 
aa1 <- calculate_metrics(data)
aa2 <-calculate_metrics(df)
aa3 <-calculate_metrics(df2)
aa4 <- rbind(aa1,aa2,aa3)

data$.pred <- 2.289*data$C0 + 1.5153*data$C2 + 4.026*data$C4 + 1.5881
df$.pred<- 2.289*df$C0 + 1.5153*df$C2 + 4.026*df$C4 + 1.5881
df2$.pred<- 2.289*df2$C0 + 1.5153*df2$C2 + 4.026*df2$C4 + 1.5881
aa1 <- calculate_metrics(data)
aa2 <-calculate_metrics(df)
aa3 <-calculate_metrics(df2)
aa4 <- rbind(aa1,aa2,aa3)

#Typical concentration-time curves for each dataset
concss <- rbind(dataconc,dfconc,df2conc) 
concss <- concss %>% ungroup %>%  mutate(time = case_when(
  time == 107.5 ~ 0,
  time == 108.5 ~ 0.5,
  time == 109 ~ 1,
  time == 110 ~ 2,
  time == 111 ~ 3,
  time == 112 ~ 4,
  time == 114 ~ 6,
  time == 116 ~ 8,
  time == 120 ~ 12,
  TRUE ~ time)) %>% filter(time %in% c(0,0.5,1,2,3,4,6,8,12))
  concs <- concss %>% 
    group_by(ID, time, Group) %>% 
    slice_max(DV, n = 1, with_ties = TRUE) %>% ungroup() %>% pivot_wider(id_cols = c("ID","Group"), names_from = time, values_from = DV)
  
  concs <- concss %>% 
    group_by(ID, time, Group) %>% 
    slice_max(DV, n = 1, with_ties = TRUE) %>% ungroup() 
    
  concs1 <- concs %>% 
    group_by(time, Group) %>% 
    summarise(
      median_DV = median(DV, na.rm = TRUE),
      Q1 = quantile(DV, 0.25, na.rm = TRUE),  
      Q3 = quantile(DV, 0.75, na.rm = TRUE),  
      IQR = IQR(DV, na.rm = TRUE),            
      n = n(),                                 
      .groups = 'drop'  
    )
  
  concs1 %>% 
    ggplot(mapping = aes(x = time, y = median_DV, colour = Group)) + 
    geom_point(position = position_dodge(width = 0.2), size = 3) +
    geom_line(position = position_dodge(width = 0.2), linewidth = 1.2) +
    geom_errorbar(
      aes(ymin = Q1, ymax = Q3), 
      position = position_dodge(width = 0.2), 
      width = 0.3,
      size = 0.8
    ) +
    scale_y_log10(breaks = c(0.01, 0.1, 1), limits = c(0.008,2.5)) +
    labs(
      x = "Time (h)", 
      y = "Plasma concentration (mg/L)",
      color = ""
    ) + 
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 17),
      axis.text.y = element_text(size = 17),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 17),
      plot.title = element_text(size = 16),
        panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.grid.minor = element_line(color = "gray95", size = 0.1),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = c(0.6, 0.4),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "transparent")
      )
