## Purpose: Run models for the botany data
## Project: FireBioOptima

botmod_run <- function(decay_rate,
                       buff_width,
                       tws, # trait weights in percent, must be in order of severity, fri, patch size
                       n.fold = NULL,
                       save_model = T, # logical, do we save the full model or just the summary output
                       save_waic = F, # logical, do we save the waic output
                       return = F #logical, do we return the model object
) 
{
  library(tidyverse)
  library(brms)
  library(tidybayes)
  
  ## read in y and covariate data
  d <- read_csv('Data/plant_richness.csv') %>% 
    mutate(point = tolower(plot),
           point = case_when(point == "ynp_1" ~ 'ynp_01',
                             point == 'ynp_2' ~ 'ynp_02',
                             point == 'ynp_3' ~ 'ynp_03',
                             point == 'ynp_4' ~ 'ynp_04',
                             point == 'ynp_7' ~ 'ynp_07',
                             TRUE ~ point))

  covs <- read_rds(paste0('Data/ucb_covs', buff_width, '.rds'))
  
  ## Pull out severity data
  sev <- filter(covs$sev, decay == decay_rate, win == 27, stat == "mean") %>% 
    pivot_wider(id_cols = c(point), names_from = c(ref_yr), names_prefix = "sev_", values_from = cbi) %>% 
    ## make sure all sampled single digits have zeros before
    mutate(point = case_when(point == "ynp_2" ~ "ynp_02",
                             point == "ynp_1" ~ "ynp_01",
                             point == "ynp_4" ~ "ynp_04",
                             point == "ynp_3" ~ "ynp_03",
                             point == "ynp_6" ~ "ynp_06",
                             point == "ynp_5" ~ "ynp_05",
                             point == "ynp_7" ~ "ynp_07",
                             point == "ynp_8" ~ "ynp_08",
                             point == "ynp_9" ~ "ynp_09",
                             point == "aaw_7" ~ "aaw_07",
                             point == "aaw_6" ~ "aaw_06",
                             point == "aaw_1" ~ "aaw_01",
                             point == "aaw_5" ~ "aaw_05",
                             point == "aaw_9" ~ "aaw_09",
                             TRUE ~ point)) %>% 
    filter(point %in% d$point) 
  
  ## pull average TSLF
  ## pull average frequency
  fri <- filter(covs$freq, buffer == buff_width, win == 27, 
                decay == decay_rate, stat == "mean") %>% 
    pivot_wider(id_cols = c(point), names_from = c(metric, ref_yr), values_from = value) %>% 
    mutate(point = case_when(point == "ynp_2" ~ "ynp_02",
                             point == "ynp_1" ~ "ynp_01",
                             point == "ynp_4" ~ "ynp_04",
                             point == "ynp_3" ~ "ynp_03",
                             point == "ynp_6" ~ "ynp_06",
                             point == "ynp_5" ~ "ynp_05",
                             point == "ynp_7" ~ "ynp_07",
                             point == "ynp_8" ~ "ynp_08",
                             point == "ynp_9" ~ "ynp_09",
                             point == "aaw_7" ~ "aaw_07",
                             point == "aaw_6" ~ "aaw_06",
                             point == "aaw_1" ~ "aaw_01",
                             point == "aaw_5" ~ "aaw_05",
                             point == "aaw_9" ~ "aaw_09",
                             TRUE ~ point)) %>% 
    filter(point %in% d$point) 
  
  ## pull out pyrodiversity
  pd <- filter(covs$pd, radius == buff_width, decay == decay_rate, wind == 27, 
               sev_w == tws[1], fri_w == tws[2], pat_w == tws[3]) %>% 
    pivot_wider(id_cols = c(name), names_from = c(year), names_prefix = "fdis_", values_from = FDis) %>% 
    rename(point = name) %>% 
    mutate(point = case_when(point == "ynp_2" ~ "ynp_02",
                             point == "ynp_1" ~ "ynp_01",
                             point == "ynp_4" ~ "ynp_04",
                             point == "ynp_3" ~ "ynp_03",
                             point == "ynp_6" ~ "ynp_06",
                             point == "ynp_5" ~ "ynp_05",
                             point == "ynp_7" ~ "ynp_07",
                             TRUE ~ point)) %>% 
    filter(point %in% d$point) %>% 
    arrange(point)
  
  ## average patch size
  pat <- filter(covs$pat, buffer == buff_width, decay == decay_rate, stat == "mean", win == 27) %>% 
    pivot_wider(id_cols = c(point), names_from = c(ref_yr), names_prefix = "pat_", values_from = log_ha) %>% 
    mutate(point = case_when(point == "ynp_2" ~ "ynp_02",
                             point == "ynp_1" ~ "ynp_01",
                             point == "ynp_4" ~ "ynp_04",
                             point == "ynp_3" ~ "ynp_03",
                             point == "ynp_6" ~ "ynp_06",
                             point == "ynp_5" ~ "ynp_05",
                             point == "ynp_7" ~ "ynp_07",
                             TRUE ~ point)) %>% 
    filter(point %in% d$point) %>% 
    arrange(point)
  
  ## pull out elevation
  elev <- covs$elev %>% 
    mutate(point = case_when(point == "ynp_2" ~ "ynp_02",
                             point == "ynp_1" ~ "ynp_01",
                             point == "ynp_4" ~ "ynp_04",
                             point == "ynp_3" ~ "ynp_03",
                             point == "ynp_6" ~ "ynp_06",
                             point == "ynp_5" ~ "ynp_05",
                             point == "ynp_7" ~ "ynp_07",
                             TRUE ~ point))  %>% 
    filter(point %in% d$point) 
  
  ## put them all together
  d2 <- left_join(d, sev) %>%
    left_join(fri) %>% 
    left_join(pd) %>% 
    left_join(pat) %>%
    left_join(elev) %>% 
    mutate(sev_mn = ifelse(year == 2013, sev_2012, sev_2020),
           fri = ifelse(year == 2013, fri_2012, fri_2020),
           pat = ifelse(year == 2013, pat_2012, pat_2020),
           fdis = case_when(year == 2013 ~ fdis_2013,
                            year == 2021 ~ fdis_2021,
                            year == 2022 ~ fdis_2022),
           area = str_sub(point, 1, 3),
           year = as.factor(year),
           fdis_s = (fdis - mean(fdis)) / sd(fdis),
           pat_s = (pat - mean(pat)) / sd(pat),
           elev_s = (elevation - mean(elevation)) / sd(elevation),
           sev_s = (sev_mn - mean(sev_mn)) / sd(sev_mn),
           fri_s = (fri - mean(fri)) / sd(fri),
           pt = as.integer(as.factor(point)),
           area = str_sub(plot, 1,3),
           park = ifelse(area == "aaw", "aaw", "ynp"),
           pk = as.integer(as.factor(park))) %>% 
    dplyr::select(point, year:soil, elevation:pk) %>% 
    ## We have one plot with missing soil data drop 
    filter(!is.na(soil))
  
  ## run a model
  m <- brm(n ~ fdis_s + I(fdis_s^2) +
             sev_s + I(sev_s^2) +
             fri_s + I(fri_s^2) +
             pat_s + I(pat_s^2) +
             elev_s + I(elev_s^2) +
             (1 | soil + pt),
                data = d2, family = 'poisson', chains = 3,
                control = list(adapt_delta = 0.99))
  
  ## Calculate waic
  m.waic <- waic(m)$estimates %>% 
    as.data.frame()
  
  m.loo <- loo(m)$estimates %>% 
    data.frame()
  
  if(!is.null(n.fold))
  {
    m.kfold <- kfold(m, K = n.fold)
    kfoldic <- m.kfold$estimates %>% 
      data.frame()
    
    ## bundle data with model
    out <- list(d = d2,
                m = m,
                waic = m.waic["waic", "Estimate"],
                waic_se = m.waic["waic", "SE"],
                looic = m.loo["looic", "Estimate"],
                looic_se = m.loo["looic", "SE"],
                kfoldic = kfoldic["kfoldic", "Estimate"],
                kfoldic_se = kfoldic["kfoldic", "SE"],
                taxa = "plants",
                decay_rate = decay_rate,
                buff_width = buff_width,
                tws = tws)
    
  } else
  {
  
  ## bundle data with model
  out <- list(d = d2,
              m = m,
              taxa = "plants",
              decay_rate = decay_rate,
              buff_width = buff_width,
              tws = tws)
  }


    
    if(save_model) {

      write_rds(out, paste0("models/plant_rich_d", decay_rate, "_w", buff_width, 
                            '_s', tws[1], 'f', tws[2], 'p', tws[3],  
                            "_full.rds"), 
                compress = "bz") 
    }
  
  ## Save waic score
  if(save_waic){
    out.waic <- data.frame(taxa = "plants", decay_rate = decay_rate, buff_width = buff_width, 
                           tws = paste(tws, collapse = " "),
                           waic = m.waic["waic", "Estimate"],
                           waic_se = m.waic["waic", "SE"],
                           looic = m.loo["looic", "Estimate"],
                           looic_se = m.loo["looic", "SE"],
                           kfoldic = kfoldic["kfoldic", "Estimate"],
                           kfoldic_se = kfoldic["kfoldic", "SE"])
    
    ## append to existing dataframe
    ap <- ifelse(file.exists('Data/Results/plant_full_waic.csv'), T, F)
    write_csv(out.waic, file = "Data/Results/plant_full_waic.csv", append = ap)
  }

  
  capture.output(summary(m), file = paste0("models/summaries/plant_rich_d", decay_rate, "_w", buff_width, 
                                                '_s', tws[1], 'f', tws[2], 'p', tws[3], 
                                                "_full_summary.txt"))
  
  if(return) {
    return(out)
  }
  
}
