#2021 POPPK Model 

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

$OMEGA @name IIV @labels eta_V1 eta_CL eta_Q eta_V2
0.3959594    // IIV_V1 (69.7% CV) 
0.07862262    // IIV_CL (28.6% CV)  
0.4129206    // IIV_Q (71.5% CV) 
0.8446246    // IIV_V2 (115.2% CV)

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
double KTRi = 4 / MATi;
double Fi = (1 * (1 + LOWFAT * thetaFLOWFAT) * (1 + PREG * thetaFPREG) * (1 + EFV * thetaFEFV)) * exp(eta_F);
F_GUT = Fi;


$ODE
dxdt_GUT = -KA * GUT;
dxdt_TRANS1 = KA * GUT - KTRi * TRANS1;
dxdt_TRANS2 = KTRi * TRANS1 - KTRi * TRANS2;
dxdt_TRANS3 = KTRi * TRANS2 - KTRi * TRANS3;
dxdt_CENTRAL = KTRi * TRANS3 - (CLi / V1i) * CENTRAL - (Qi / V1i) * CENTRAL + (Qi / V2i) * PERIPH;
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

# --- FED (Renomeado para 'fed_data') ---
set.seed(4658)
fed_data = tibble(ID = 1:4500) %>%
  mutate(FED = rbinom(n(), 1, 0.50)) # A coluna interna ainda se chama FED

# --- LOWFAT (Renomeado para 'lowfat_data' e usando 'fed_data') ---
set.seed(9108)
lowfat_data <- fed_data %>% # Usar o tibble renomeado
  mutate(LOWFAT = ifelse(FED == 1, rbinom(n(), 1, 0.50), 0)) %>%
  select(ID, LOWFAT)

#Dose simulation- 400 mg BID e 1200 mg QD

set.seed(4040)
dose_events <- ev(ID = 1:4500, amt = 400, ii = 12, addl = 6)
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
  mrgsim(end = 72, delta = 0.5, carry_out = "OCC")

simul_ral <- as_tibble(out) %>%
  left_join(pac %>%
              select(ID, OCC, dose_group, PREG, ATV, EFV, FED, WT, LOWFAT) %>% distinct(ID, OCC, .keep_all = TRUE),
            by = c("ID", "OCC"))

simul_ral2 <- simul_ral %>% filter(time>59)
aucs <- simul_ral2 %>%
  group_by(ID) %>%
  mutate(AUCt = trapz(time, DV))

simul_ral3 <- aucs %>% filter(time==59.5 | time==61 |time==62 | time==63|time==64 | time==60.5)
pivot_conc <- simul_ral3 %>% 
  select(ID, time, DV) %>% 
  pivot_wider(id_cols = ID, names_from = time, names_prefix = "conc_", values_from = DV)

base_prediction_auc <- simul_ral3 %>% 
  filter(time==59.5) %>% 
  select( -time, -CP) %>% left_join(pivot_conc)

#Remove <1 >99 percentiles
centile_1_auc <- quantile(base_prediction_auc$AUCt, 0.01)
centile_99_auc <- quantile(base_prediction_auc$AUCt, 0.99)

base_prediction_auc1 <- base_prediction_auc %>%
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc)

data22 <- fread("ralsim2021.csv")
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
0.01 //0.3074847   // residual variance (60% CV)

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
set.seed(str_c(12,i))
sex = tibble(ID = 1:4500) %>% 
  mutate(sex = rbinom(n(),1,0.786))
set.seed(str_c(12,3))
race = tibble(ID = 1:4500) %>% 
  mutate(race = rbinom(n(),1,0.09))
#https://pmc.ncbi.nlm.nih.gov/articles/PMC5377560/
set.seed(15744)
hiv = tibble(ID = 1:4500) %>% 
  mutate(hiv = rbinom(n(),1,0.85))
set.seed(4442)
atv <- hiv %>%
  mutate(atv = ifelse(hiv == 1, rbinom(n(), 1, 0.09), 0))
set.seed(str_c(12,i))
bilirubin = tibble(ID = 1:4000) %>% 
  mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=12, sd=4)) %>% 
  rbind(tibble(ID = 4001:4450) %>% 
          mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=30, sd=8))) %>% 
  rbind (tibble(ID = 4451:4500) %>% 
           mutate(bilirubin = rtruncnorm(n(),a=5, b=91, mean=70, sd=12)))
#DOSE 400mg BID in both studies 800 daily
dose <- as_tibble(ev(ID = 1:4500, amt = 400, ii = 12, addl = 5)) %>%
  mutate(dose_group = "400mg BID")
pac <- left_join(sex, race) %>% left_join(bilirubin) %>% left_join(hiv) %>% left_join(atv) %>%  left_join(dose)  



#model simulation
out <- mod %>%
  data_set(pac) %>%
  Req(CP, CL, F1) %>%
  mrgsim(end = 72, delta = 0.5)

simul_ral <- as_tibble(out) %>% 
  left_join(pac %>% 
              select(ID, dose_group, bilirubin, sex, race, hiv,atv), by = c("ID" = "ID"))

simul_ral <- simul_ral %>% 
  mutate(dose = parse_number(dose_group),
         auc = dose/CL)
simul_ral2 <- simul_ral %>% filter(time>59)
aucs <- simul_ral2 %>%
  group_by(ID) %>%
  mutate(AUCt = trapz(time, CP))

simul_ral3 <- aucs %>% filter(time==59.5 | time==61 |time==62 | time==63|time==64 | time==60.5)
pivot_conc <- simul_ral3 %>% 
  select(ID, time, CP) %>% 
  pivot_wider(id_cols = ID, names_from = time, names_prefix = "conc_", values_from = CP)

base_prediction_auc <- simul_ral3 %>% 
  filter(time==59.5) %>% 
  select( -time, -CP) %>% left_join(pivot_conc)

#Remove <1 >99 percentiles
centile_1_auc <- quantile(base_prediction_auc$AUCt, 0.01)
centile_99_auc <- quantile(base_prediction_auc$AUCt, 0.99)

base_prediction_auc1 <- base_prediction_auc %>%
  filter(AUCt >= centile_1_auc, AUCt <= centile_99_auc)

fwrite(base_prediction_auc1, file = "ralsim2012.csv")

# Parallel
n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
data <- as.data.frame(fread("ralsim2021.csv"))
all_predictors <- c( "C0","C05", "C1", "C2", "C3", "C4")
response <- "AUCt"

data <- data %>% select(-ID,-dose_group) %>% 
  mutate(Diff_C1C0 = conc_61 - conc_59.5) %>%
  rename(C0 = conc_59.5, C1 = conc_61, C2 = conc_62, C3 = conc_63,
         C4= conc_64, C05= conc_60.5) 

# Arrange dataframes
df <- as.data.frame(fread("moreira2023.csv"))
df <- df %>% filter(CMT == 2)
df$DV2 <- df$DV*444.4163/1000
df <- df %>% select(ID, TTIME, DV2, tratamento, IDD)
df <- df %>%
  group_by(IDD) %>%
  mutate(AUCt = trapz(TTIME, DV2))
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

df <- as.data.frame(fread("ralsim2012.csv"))
df <- df %>% mutate(Diff_C1C0 = conc_61 - conc_59.5) %>%
  rename(C0 = conc_59.5, C1 = conc_61, C2 = conc_62, C3 = conc_63,
         C4= conc_64, C05= conc_60.5)


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
#saveRDS(final_wf_xgb, "model_xgboost2021FIM.rds")

explainer_external <- explain_tidymodels(
  model = final_wf_xgb, 
  data = final_train, 
  y = final_train$AUCt,
  label = "xgb")
#saveRDS(explainer_external, file ="explainer_externalxgb2021FIM.rds" )

# Predictions on sets
# final_predictions <- final_wf_xgb %>% predict(final_train) %>% bind_cols(final_train)

#resample
set.seed(456)
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
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16))
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
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16))
Figure_2b

Figure_2d <- final_res_predictions %>%
  ggplot(mapping = aes(x = AUCt, y = AUCt - .pred)) +
  geom_point() +
  geom_smooth(method=lm)+
  labs(x = "Reference raltegravir AUC (mg*h/L)", y = "Test predicted AUC", title="(D)") + 
  theme_bw()
Figure_2d

#2023 External validation organizing, metrics and plots


final_wf_xgb <- readRDS("modelo_xgboost20213concsemiovCERTO.rds")
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
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16)) 


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
  theme(axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17),
        plot.title = element_text(size = 16))

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


# Id numeric and binary variables
bin_vars <- c("PREG", "ATV", "EFV", "FED", "LOWFAT")  #  bináry
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
      theme_minimal() +
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
  geom_point() +
  geom_smooth(method = "lm", col = "black", se=T) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "A   R² = 0.39, P < 0.01" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_classic(base_size = 22)

ggplot(df, aes(x = `C0`, y = `AUCt`)) +
  geom_point() +
  geom_smooth(method = "lm", col = "black", se=T) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "B   R² = 0.54, P < 0.01" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_classic(base_size = 22)

ggplot(df2, aes(x = `C0`, y = `AUCt`)) +
  geom_point() +
  geom_smooth(method = "lm", col = "black", se=T) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "C   R² = -0.03, P > 0.05" ,
       x = expression(C["trough"] * " " * (mg/L)),
       y = expression(AUC["0-12"] * " " * (mg.h/L)
       )) +
  theme_classic(base_size = 22)

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
                           col_group = c("group"), method = "median")
resultsstattest <- statcal(ral_ml_test, col_test = c("C0","C05","C1", "C2","C3","C4", "AUCt"),
                           col_group = c("group"), method = "median")
resultsstattrain <- statcal(ral_ml_train, col_test = c("C0","C05","C1", "C2","C3","C4", "AUCt"),
                            col_group = c("group"), method = "median")
resultsstat2023 <- statcal(df2, col_test = c("C0", "C05", "C1", "C2", "C3", "C4", "AUCt"),
                           col_group = c("group"), method = "median")
resusltsstatss <- rbind(resultsstattrain,resultsstattest,resultsstat2023,resultsstat2012)

# Feature selection

data <- as.data.frame(fread("ralsim2021.csv"))
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