# Copyright (C) 2025 Zoi Gutschke
#
# These scripts are part of the research materials for the MSc thesis:
#   Go with the flow: Improving Pollution Source Identification and Localization 
#   with Graph-Based Spatial Autocorrelation in Surface Waters
#
# These materials are licensed under the GNU General Public License (GPL) version 3 or later.
# You are free to redistribute and/or modify them under the terms of the GPL.
#
# This material is provided without any warranty; without even the implied warranty
# of merchantability or fitness for a particular purpose. See the GNU General Public License
# for more details: https://www.gnu.org/licenses/gpl-3.0.html
#
# You should have received a copy of the GNU General Public License along with this material.
# If not, see https://www.gnu.org/licenses/.

# GRS82000 - GRS Thesis
# Script 02: Create time windows
#
# This script creates n time windows.

############################################################################## #
#### 1 Preparation  ############################################################
############################################################################## #

# Load data and add uncertainty and method detection limit (mdl)
mw_df <- dbGetQuery(con, "
                    SELECT *
                    FROM s1_mv_oi m 
                    LEFT JOIN s0_uncertainties_raw h ON m.parametercode = h.param
                    WHERE NOT m.waterbeheerdernaam = 'Wetterskip Fryslan'")

# From data exploration
mo_n0 <- c("NL34_1100", "NL34_1101", "NL34_5009", "NL34_5012", "NL34_6109", "NL34_6402", "NL34_5010", # measurement objects with 0 samples in at least one heavy metal
           "NL33_6640")

# Filter out locations with too few measurements + 
#   calculate concentrations and uncertainties
mw_df_sel <- mw_df %>% 
  filter(!(meetobjectcode %in% mo_n0)) %>%
  mutate(numeriekewaarde = case_when(
          !is.na(limietsymbool) ~ numeriekewaarde / 2,
          is.na(limietsymbool) ~ numeriekewaarde),
        uncertainty = case_when(
          !is.na(limietsymbool) ~ (5/6) * mdl,
          is.na(limietsymbool)  ~ sqrt( (0.5 * mdl)^2 + (error * numeriekewaarde)^2 )))

# Standardize
mw_df_z <- mw_df_sel %>%
  select(meetobjectcode, parametercode, numeriekewaarde, begindatum) %>% 
  group_by(meetobjectcode, parametercode, begindatum) %>%
  mutate(set_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(meetobjectcode, begindatum, set_id),      
              names_from = parametercode, 
              values_from = numeriekewaarde  
              ) %>% 
  mutate(across(
    .cols = -c(meetobjectcode, begindatum, set_id),
    .fns = ~ scale(.)[,1] 
  )) %>% 
  drop_na() 

# Find measurement objects which have at least 12 samples per year
mo_n12 <- mw_df_sel %>% 
  mutate(ym = floor_date(begindatum, "month")) %>% 
  group_by(meetobjectcode, parametercode) %>% 
  mutate(
    n_months = n_distinct(ym)) %>% 
  ungroup() %>% 
  filter(n_months >= 12) %>% 
  distinct(meetobjectcode) %>% 
  pull(meetobjectcode)

############################################################################## #
#### 2 Split creation ##########################################################
############################################################################## #

#### 2.1 Set-up variables ######################################################
days <- seq(as.Date("2023-01-01"), as.Date("2023-12-31"), by = "1 day")
n_days <- length(days)
metric <- "mean" # or "sum"

#### 2.2 Create function to compute metric for a given contiguous interval #####
compute_interval_cost <- function(start_idx, end_idx) {
  start_date <- days[start_idx]
  end_date   <- days[end_idx]
  
  if (end_idx < start_idx) stop("end < start in cost function")
  
  # Select rows between these dates
  sel_z  <- mw_df_z$begindatum >= start_date & mw_df_z$begindatum <= end_date
  
  if (sum(sel_z) == 0) return(Inf) # no samples
  
  df_z_sel <- mw_df_z[sel_z,]
  df_z_sel_combs <- df_z_sel %>% 
    filter(meetobjectcode %in% mo_n12) %>% 
    distinct(meetobjectcode)
  
  if (nrow(df_z_sel_combs) != length(mo_n12)) return(Inf) # missing object codes
  
  vars <- sapply(df_z_sel[4:ncol(df_z_sel)], var, na.rm = TRUE)
  
  if (metric == "mean") {
    return(mean(vars, na.rm = TRUE))
  } else {
    return(sum(vars, na.rm = TRUE))
  }
}

#### 2.3 Precompute cost matrix and save it ####################################
cost_matrix <- matrix(Inf, nrow = n_days, ncol = n_days)
for (i in 1:n_days) {
  for (j in i:n_days) {
    cost_matrix[i, j] <- compute_interval_cost(i, j)
  }
}

#### 2.4 dynamic programming (DP) table for up to max_splits ###################
max_splits <- 12
dp <- matrix(Inf, nrow = max_splits, ncol = n_days)
prev_split <- matrix(NA_integer_, nrow = max_splits, ncol = n_days)

# Base case
for (d in 1:n_days) dp[1, d] <- cost_matrix[1, d]

# Fill DP
for (k in 2:max_splits) {
  for (d in k:n_days) {
    best_val <- Inf
    best_m <- NA
    for (m in (k-1):(d-1)) {
      test_val <- dp[k-1, m] + cost_matrix[m+1, d]
      if (test_val < best_val) {
        best_val <- test_val
        best_m <- m
      }
    }
    dp[k, d] <- best_val
    prev_split[k, d] <- best_m
  }
}

#### 2.7 Add columns for each split configuration ##############################
get_split_dates <- function(k) {
  splits_idx <- integer(k + 1) # creates a vector for all splits
  splits_idx[k + 1] <- n_days # assigns 365 as last split
  end_day <- n_days # defines last day
  kk <- k # why?
  while (kk > 0) {
    start_day <- if (kk == 1) 1 else prev_split[kk, end_day] + 1; start_day
    splits_idx[kk] <- start_day
    end_day <- start_day - 1
    kk <- kk - 1
  }
  sort(days[splits_idx])
}

for (k in c(4, 6, splits_n)) {
  split_dates_k <- get_split_dates(k)
  
  # Tag each row of mw_df_sel with group number for this k
  mw_df_sel[[paste0("split_", k)]] <- findInterval(
    mw_df_sel$begindatum, split_dates_k, rightmost.closed = TRUE
  )
  
  # Same for mw_df_z
  mw_df_z[[paste0("split_", k)]] <- findInterval(
    mw_df_z$begindatum, split_dates_k, rightmost.closed = TRUE
  )
}

#### 2.8 Save ##################################################################
dbWriteTable(con, "s2_mv_windows", mw_df_sel, overwrite = T, row.names = FALSE)
dbWriteTable(con, "s2_mv_z_windows", mw_df_z, overwrite = T, row.names = FALSE)

#### 2.9 Clean up ##############################################################
remove(cost_matrix, dp, mw_df, mw_df_sel, mw_df_z, prev_split, best_m, best_val,
       d, days, i, j, k, m, max_splits, metric, mo_n0, mo_n12, n_days,
       split_dates_k, test_val, compute_interval_cost, get_split_dates)
