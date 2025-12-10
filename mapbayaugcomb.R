if(!require(furrr)) install.packages("furrr")
library(furrr)  
if(!require(future)) install.packages("future")
library(future)  

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
        biais_rel = mean(bias_rel),
        relative_rmse = sqrt(mean(bias_rel_square)),
        rmse = sqrt(mean(biais_brut_sqr)),
        biais_out_20percent = mean((bias_rel) > 0.2),
        nb_out_20percent = sum((bias_rel) > 0.2),
        n = n()
      ) %>%
      mutate(
        R2 = summary(lm(AUCt ~ .pred, data = pred_df))$r.squared
      )
  }

  test_purrr <- function(data, timestart, timeend, combo_sizes, time_predictors, returnAUC = FALSE) {
    dataAUCt <- data %>% select(ID, AUCt)
    data <- data %>% select(-AUCt)
    
    # Parallel
    n_cores <- availableCores() - 2
    plan(multisession, workers = n_cores)
    on.exit(plan(sequential))
    
    # Combinations
    all_combos <- map(combo_sizes, ~combn(time_predictors, .x, simplify = FALSE)) %>% flatten()
    names(all_combos) <- map_chr(all_combos, ~paste(.x, collapse = "_"))
    
    # Parallel combination processing
    process_combo <- function(combo, combo_name) {
      tryCatch({
        # Filter
        datatime <- data %>% 
          filter(time %in% combo | EVID == 1) %>%
          distinct(ID, time, .keep_all = TRUE)
        
        
        unique_ids <- unique(datatime$ID)
        
        # Parallel IDs
        pred_df <- future_map_dfr(unique_ids, ~{
          tryCatch({
            dataID <- datatime %>% filter(ID == .x)
            estim <- mapbayest(mod, data = dataID)
            aug <- mapbayr::augment(estim, start = timestart, end = timeend, delta = 0.5)
            
            auc_pred <- aug$aug_tab %>%
              filter(type == "IPRED", between(time, timestart, timeend)) %>%
              summarise(auc = trapz(time, value)) %>% pull(auc)
            
            tibble(ID = .x, .pred = auc_pred)
          }, error = function(e) { 
            tibble(ID = .x, .pred = NA_real_) 
          })
        }, .options = furrr_options(seed = TRUE)) %>% 
          left_join(dataAUCt, by = "ID") %>% 
          filter(!is.na(.pred), !is.na(AUCt)) %>%
          mutate(combo = combo_name, n_points = length(combo)) %>% 
          distinct()
        
        # Metrics
        metrics <- calculate_metrics(pred_df) %>%
          mutate(combo = combo_name, 
                 n_points = length(combo),
                 n_success = nrow(pred_df),
                 n_total = length(unique_ids),
                 success_rate = n_success / n_total)
        
        list(metrics = metrics, predictions = pred_df)
        
      }, error = function(e) {
        list(metrics = tibble(combo = combo_name, n_points = length(combo), error = "Erro"))
      })
    }
    
    # Execute parallel combinations
    results <- future_map2(
      all_combos, 
      names(all_combos), 
      process_combo,
      .options = furrr_options(seed = TRUE, scheduling = 2)
    )
    
    # Predictions and metrics
    metrics_df <- map_dfr(results, ~.x$metrics)
    
    if (returnAUC) {
      predictions_list <- map(results, ~.x$predictions)
      names(predictions_list) <- names(all_combos)
      list(metrics = metrics_df, predictions = predictions_list)
    } else {
      list(metrics = metrics_df)
    }
  }
  
  #result
  results <- test_purrr(
    data = df1,
    timestart = 107.5,
    timeend = 120,
    combo_sizes = 3,
    time_predictors = c(108,108.5,109,110,111,112),
    returnAUC = T
  ) 

  df2est <- df2est %>% filter(ID == 1)
estim <- mapbayest(mod, data = df2est)

aug <- mapbayr::augment(estim, start = 107.5 ,end = 120, delta = 0.5) 
aug$aug_tab
individual_esitmate <- aug$aug_tab %>%
  filter(type == "IPRED",
         between(time, 107.5, 120)) %>% 
  select(ID, time, DV = value)  %>% 
  group_by(ID) %>% 
  mutate(.pred = trapz(time, DV)) %>% 
  slice_min(DV) %>% 
  select(.pred) %>%
  distinct(.pred, .keep_all = TRUE)
