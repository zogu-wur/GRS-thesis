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
# Script 04: Calculating spatial autocorrelation and preparing PMF
#
# This script calculates the global and local Moran's I and Geary's C, and 
#   creates the PMF input

############################################################################## #
#### 1 Preparation #############################################################
############################################################################## #

# Load data
mw_df_all <- dbGetQuery(con,"SELECT * FROM s2_mv_windows")
mw_df_z <- dbGetQuery(con,"SELECT * FROM s2_mv_z_windows")
meetobjecten_oi <- dbGetQuery(con,"SELECT * FROM s1_mo_oi")

districts_sf <- st_read(con, "s0_hea_districts_raw") %>% select(district)
districts_nzv <- st_read(con, "s0_nzv_districts_raw")

# Pre-processing
mw_df_z_unique <- mw_df_z %>%
  group_by(meetobjectcode, begindatum, split_12) %>%
  summarise(
    across(all_of(metal_oi), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop")

districts_sf <- districts_sf %>% 
  rbind(districts_nzv)

############################################################################## #
#### 2 Prepare regular PMF #####################################################
############################################################################## #

pmf_uncertainty <- function(df, goal, i = 1, splits_n = 12){
  
  # aggregate concentrations and uncertainties
  aggregated_df <- df %>% 
    group_by(meetobjectcode, parametercode) %>% 
    summarise(
      mean_c = mean(numeriekewaarde),
      mean_u = max(uncertainty)
    )
  
  if (goal == "concentration"){
    from <- "mean_c"
  } else if (goal == "uncertainty"){
    from <- "mean_u"
  }
  
  # change format
  df_sel <- aggregated_df %>% 
    select(meetobjectcode, parametercode, all_of(from)) %>% 
    pivot_wider(names_from = parametercode,
                values_from = all_of(from)) %>% 
    arrange(meetobjectcode)
  
  # create file paths to save and save
  fp_csv <- paste0(pmf_path, "pmf_split_", i, "_", goal, ".csv")
  fp_txt <- paste0(pmf_path, "pmf_split_", i, "_", goal, ".txt")
  
  write.csv(df_sel, file = fp_csv, row.names = FALSE)
  write.table(df_sel, file = fp_txt, sep = "\t", row.names = FALSE)
}

for (i in 1:splits_n){
  
  mw_df <- mw_df_all %>% 
    filter(split_12 == i)
  
  for (j in c("concentration", "uncertainty")){
    pmf_uncertainty(mw_df, j, i)
  }
}

############################################################################## #
#### 3 Global SAC ##############################################################
############################################################################## #

params <- metal_oi

global_autocorrelation <- function(net, params, i, permutations = 199) {
  # Add edge weights
  net <- net %>%
    activate("edges") %>%
    mutate(weight = edge_length())
  
  # Compute network distance and weights
  d_cm <- st_network_cost(net, weights = "weight")
  d_cm <- drop_units(d_cm)
  
  w_cm <- 1 / d_cm
  diag(w_cm) <- 0
  w_cm[!is.finite(w_cm)] <- 0
  
  res_ls <- list()
  
  for (j in params){
    node_sf <- net %>% activate("nodes") %>% st_as_sf()
    values <- node_sf[[j]]
    n <- length(values)
    mean_x <- mean(values)
    var_x <- var(values)
    W <- sum(w_cm)
    
    # Global Moran's I
    numerator_I <- sum(w_cm * outer(values - mean_x, values - mean_x))
    I_global <- (n / W) * (numerator_I / sum((values - mean_x)^2))
    I_global <- as.numeric(I_global)
    
    # Global Geary's C
    numerator_C <- sum(w_cm * outer(values, values, function(xi, xj) (xi - xj)^2))
    C_global <- ((n - 1) / (2 * W)) * (numerator_C / sum((values - mean_x)^2))
    C_global <- as.numeric(C_global)
    
    # Permutation test
    I_perms <- numeric(permutations)
    C_perms <- numeric(permutations)
    for (p in 1:permutations) {
      permuted <- sample(values)
      num_I_perm <- sum(w_cm * outer(permuted - mean_x, permuted - mean_x))
      I_perms[p] <- (n / W) * (num_I_perm / sum((permuted - mean_x)^2))
      
      num_C_perm <- sum(w_cm * outer(permuted, permuted, function(xi, xj) (xi - xj)^2))
      C_perms[p] <- ((n - 1) / (2 * W)) * (num_C_perm / sum((permuted - mean_x)^2))
    }
    
    # p-values
    p_I <- (sum(I_perms >= as.numeric(I_global)) + 1) / (permutations + 1)
    p_C <- (sum(C_perms <= as.numeric(C_global)) + 1) / (permutations + 1)
    
    # combine results
    res_df <- as.data.frame(list(
      parametercode = j,
      morans_I = I_global,
      morans_I_p = p_I,
      gearys_C = C_global,
      gearys_C_p = p_C
    ))
    
    res_ls <- append(res_ls, list(res_df))
    
  }
  
  df_res <- do.call(rbind, res_ls)
  df_res$split <- i
  
  return(df_res)
}

# Run
df_ls <- list()

for(i in 1:splits_n){
  
  fp_nodes <- paste0("temp/net_split_", i, "_nodes.gpkg")
  fp_edges <- paste0("temp/net_split_", i, "_edges.gpkg")
  
  nodes <- st_read(fp_nodes)
  edges <- st_read(fp_edges)
  
  window_i_c <- mw_df_z_unique %>% 
    filter(split_12 == i)
  
  net <- sfnetwork(nodes, edges, directed = T) %>%
    activate("nodes") %>%
    left_join(window_i_c)

  df_res <- global_autocorrelation(net, metal_oi, i)
  df_ls <- append(df_ls, list(df_res))
}

df_res <- do.call(rbind, df_ls)

# 1. Prepare data – create a mask for significant values based on p-value columns
df_long <- df_res %>%
  mutate(sig_morans = morans_I_p <= 0.05,
         sig_gearys = gearys_C_p <= 0.05) %>%
  pivot_longer(cols = c("morans_I", "gearys_C"),
               names_to = "stat_type",
               values_to = "stat_value") %>%
  mutate(significant = case_when(
    stat_type == "morans_I" ~ sig_morans,
    stat_type == "gearys_C" ~ sig_gearys
  ),
  stat_value_masked = ifelse(significant, stat_value, NA))

dbWriteTable(con, "s4_globalsac", df_long, overwrite = TRUE)

############################################################################## #
#### 5 Local SAC ###############################################################
############################################################################## #

calculate_local_autocorrelation_with_significance <- function(net, parameters, split_i, permutations = 199, 
                                                              metric = "moransI") {
  
  # Add edge weights
  net <- net %>%
    activate("edges") %>%
    mutate(weight = edge_length()) %>% 
    activate("nodes") %>% 
    arrange(district)
  
  districts <- net %>% 
    mutate(id = row_number()) %>% 
    pull(district)
  
  # Calculate distance matrix and spatial weights
  d_cm <- st_network_cost(net, weights = "weight")
  d_cm <- drop_units(d_cm)
  
  # Create a comparison matrix showing whether each (row, col) share the same district
  district_mask <- outer(districts, districts, FUN = "==")
  
  # Mask distances: set to NA if not in the same district
  d_cm[!district_mask] <- 0
  
  # Transfer into weights
  w_cm <- 1 / d_cm
  diag(w_cm) <- 0
  w_cm[!is.finite(w_cm)] <- 0
  
  # Get node data
  node_sf <- net %>% activate("nodes") %>% st_as_sf()
  
  # Loop over each parameter
  for (param in parameters) {
    cat("Processing:", param, "\n")
    
    x <- node_sf[[param]]
    n <- length(x)
    x_mean <- mean(x)
    x_var <- var(x)
    
    # Observed values
    moransI_obs <- numeric(n)
    gearysC_obs <- numeric(n)
    
    for (i in seq_len(n)) {
      v <- x[i] # xi
      w <- as.vector(w_cm[i, ]) # weights of i
      
      moransI_obs[i] <- ((v - x_mean) / x_var) * sum(w * (x - x_mean), na.rm = TRUE)
      gearysC_obs[i] <- (1 / x_var) * sum(w * (v - x)^2, na.rm = TRUE)
      }
    
    # Permutations
    moransI_perm <- matrix(NA, nrow = n, ncol = permutations)
    gearysC_perm <- matrix(NA, nrow = n, ncol = permutations)
    
    for (p in 1:permutations) {
      x_perm <- sample(x)
      x_perm_mean <- mean(x_perm)
      x_perm_var <- var(x_perm)
      
      for (j in seq_len(n)) {
        v <- x_perm[j]
        w <- as.vector(w_cm[j, ])
        
        moransI_perm[j, p] <- ((v - x_perm_mean) / x_perm_var) * sum(w * (x_perm - x_perm_mean), na.rm = TRUE)
        gearysC_perm[j, p] <- (1 / x_perm_var) * sum(w * (v - x_perm)^2, na.rm = TRUE)
      }
    }
    
    # P-values (two-sided)
    moransI_p <- rowMeans(sweep(moransI_perm, 1, moransI_obs, function(perms, obs) abs(perms) >= abs(obs)), na.rm = TRUE)
    gearysC_p <- rowMeans(sweep(gearysC_perm, 1, gearysC_obs, function(perms, obs) abs(perms - 1) >= abs(obs - 1)), na.rm = TRUE)
    
    # Add to sf
    node_sf[[paste0(param, "_local_moransI")]] <- moransI_obs
    node_sf[[paste0(param, "_local_moransI_p")]] <- moransI_p
    node_sf[[paste0(param, "_local_gearysC")]] <- gearysC_obs
    node_sf[[paste0(param, "_local_gearysC_p")]] <- gearysC_p
    
    value_col <- paste0(param, "_local_", ifelse(metric == "moransI", "moransI", "gearysC"))
    p_col <- paste0(value_col, "_p")
    
    val_min <- min(node_sf[[value_col]], na.rm = TRUE)
    val_max <- max(node_sf[[value_col]], na.rm = TRUE)
    
    scale_factor <- 1000
    sf_df <- node_sf %>%
      mutate(
        val_scaled = !!sym(value_col) * scale_factor,
        sig = (!!sym(p_col) <= 0.05) & !is.na(!!sym(p_col))
      )
    
    sig_count <- sum(sf_df$sig, na.rm = TRUE)
    total_count <- nrow(sf_df)
    
    # Check if all points are non-significant or all fo fill_color is NA
    if(all(!sf_df$sig)) {
      message("No points are significant (p <= 0.05) for ", param)
      sf_df <- sf_df %>%
        mutate(fill_color = NA_real_)  # all greyed out
      subtitle_text <- paste0("No significant points at p ≤ 0.05; Value range: ",
                              signif(val_min, 3), " to ", signif(val_max, 3), " (before scaling)")
    } else {
      sf_df <- sf_df %>%
        mutate(fill_color = ifelse(sig, val_scaled, NA_real_))
      subtitle_text <- paste0(sig_count, " of ", total_count, " points significant at p ≤ 0.05; Value range: ",
                              signif(val_min, 3), " to ", signif(val_max, 3), " (before scaling)")

      p <- ggplot(sf_df) +
        geom_sf(aes(fill = fill_color), size = 3, shape = 21) +
        scale_fill_distiller(
          palette = "RdBu", direction = 1,
          na.value = "grey80",
          name = paste0("Local ", ifelse(metric == "moransI", "Moran's I", "Geary's C"), " × ", scale_factor),
          labels = scales::label_number(accuracy = 0.01)
        ) +
        labs(
          title = paste("Local", ifelse(metric == "moransI", "Moran's I", "Geary's C"), "for", param),
          subtitle = subtitle_text
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "right",
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)
        )

      p

      file_name <- paste0("plots/sac/local_", metric, "_split_", split_i, "_", param, ".png")
      ggsave(file_name, p, width = 8, height = 6, dpi = 300, bg = "white")

      message("Saved: ", file_name)
    }
    
  }
  
  # Build relevant column names dynamically
  value_cols <- paste0(metal_oi, "_local_", metric)
  p_cols <- paste0(metal_oi, "_local_", metric, "_p")
  
  res_df <- sf_df %>% 
    select(meetobjectcode, contains(metric)) %>% 
    mutate(split = split_i)
  
  #### CALCULATE CHANGED UNCERTAINTY
  value_df <- sf_df %>% 
    as.data.frame() %>% 
    select(meetobjectcode, contains(metric)) %>% 
    arrange(meetobjectcode)  %>%
    rowwise() %>%
    mutate(
      across(
        all_of(value_cols),
        ~ {
          param <- str_remove(cur_column(), paste0("_", metric, "$"))
          pval <- get(paste0(param, "_", metric, "_p"))
          newval <- if (metric == "gearysC") {
            ifelse(pval <= 0.05, 2 - .x, 1)
          } else {
            ifelse(pval <= 0.05, 1 + .x, 1)
          }
          newval <- ifelse(newval < 0.01, 0.01, newval)
          newval <- ifelse(newval > 1.99, 1.99, newval)
          newval
        },
        .names = "{.col}_tmp"
      )
    ) %>%
    select(-all_of(c(value_cols, p_cols))) %>%
    rename_with(
      ~ str_remove(.x, paste0("_local_", metric, "_tmp$")),
      ends_with("_tmp")
    )
  
  # create file paths to load
  fp_csv <- paste0(pmf_path, "pmf_split_", split_i, "_uncertainty.csv")
  u_df <- read.csv(fp_csv)
  
  v_n <- unique(value_df$meetobjectcode)
  u_n <- unique(u_df$meetobjectcode)
  
  print(u_n[!(u_n %in% v_n)])
  
  # Multiply matching columns and preserve meetobjectcode:
  result_df <- value_df %>%
    select(meetobjectcode) %>%
    bind_cols(value_df[metal_oi] * u_df[metal_oi])
  
  # to save
  fp_csv <- paste0(pmf_path, "pmf_split_", split_i, "_uncertainty_", metric, ".csv")
  fp_txt <- paste0(pmf_path, "pmf_split_", split_i, "_uncertainty_", metric, ".txt")
  
  write.csv(result_df, file = fp_csv, row.names = FALSE)
  write.table(result_df, file = fp_txt, sep = "\t", row.names = FALSE)
 
  return(res_df)
}

# Run
res_ls <- list()

# for (i in 1:splits_n){
for (i in splits){
 
  fp_nodes <- paste0("temp/net_split_", i, "_nodes.gpkg")
  fp_edges <- paste0("temp/net_split_", i, "_edges.gpkg")
  
  nodes <- st_read(fp_nodes)
  edges <- st_read(fp_edges)
  
  # add district information
  nodes_d <- st_join(nodes, districts_sf) %>% 
    mutate(district = case_when(
      meetobjectcode == "NL33_1638" ~ "Oost",
      meetobjectcode %in% c("NL33_1103", "NL33_4106") ~ "Noord",
      meetobjectcode %in% c("NL34_1100", "NL34_1101") ~ "Drents Plateau met beekdalen",
      meetobjectcode == "NL34_5012" ~ "Zuidelijk Westerkwartier en benedenlopen beken",
      .default = district
    ))
  
  window_i_c <- mw_df_z_unique %>% 
    filter(split_12 == i)
  
  net <- sfnetwork(nodes_d, edges, directed = T) %>%
    activate("nodes") %>%
    left_join(window_i_c)
  b_ls <- list()
  for (m in c("moransI", "gearysC")){
    res <- calculate_local_autocorrelation_with_significance(net, metal_oi, i, metric = m)
    b_ls <- append(b_ls, list(res))
  }
  
  t <- b_ls[[1]] %>% 
    left_join(st_drop_geometry(b_ls[[2]]), by = c("meetobjectcode", "split"))

  res_ls <- append(res_ls, list(t))
}

# Run and save 
res_df <- do.call(rbind, res_ls)
st_write(res_df, con, layer= "s4_localsac", delete_layer = TRUE)

# Summarize for figures
morans_cols <- names(res_df)[grep("local_moransI$", names(res_df))]
morans_p_cols <- names(res_df)[grep("local_moransI_p$", names(res_df))]
gearys_cols <- names(res_df)[grep("local_gearysC$", names(res_df))]
gearys_p_cols <- names(res_df)[grep("local_gearysC_p$", names(res_df))]

# Convert to long format for Moran's I
morans_long <- res_df %>%
  select(meetobjectcode, all_of(morans_cols), all_of(morans_p_cols)) %>%
  st_drop_geometry() %>%  
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = all_of(morans_cols), names_to = "parameter", values_to = "moransI") %>%
  mutate(parameter = sub("_local_moransI$", "", parameter)) %>%
  left_join(
    res_df %>%
      st_drop_geometry() %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(cols = all_of(morans_p_cols), names_to = "parameter_p", values_to = "moransI_p") %>%
      mutate(parameter_p = sub("_local_moransI_p$", "", parameter_p)),
    by = c("row_id", "parameter" = "parameter_p")
  ) %>%
  select(row_id, meetobjectcode = meetobjectcode.x, parameter, moransI, moransI_p)

# Similar pivot for Geary's C
gearys_long <- res_df %>%
  select(meetobjectcode, geom, all_of(gearys_cols), all_of(gearys_p_cols)) %>%
  st_drop_geometry() %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = all_of(gearys_cols), names_to = "parameter", values_to = "gearysC") %>%
  mutate(parameter = sub("_local_gearysC$", "", parameter)) %>%
  left_join(
    res_df %>%
      st_drop_geometry() %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(cols = all_of(gearys_p_cols), names_to = "parameter_p", values_to = "gearysC_p") %>%
      mutate(parameter_p = sub("_local_gearysC_p$", "", parameter_p)),
    by = c("row_id", "parameter" = "parameter_p")
  ) %>%
  select(row_id, meetobjectcode = meetobjectcode.x, parameter, gearysC, gearysC_p)

# Join Moran's I and Geary's C data, restoring geometry by joining with unique meetobjectcode and geom
geom_lookup <- res_df %>% 
  select(meetobjectcode) %>% 
  distinct(meetobjectcode, .keep_all = TRUE)

summary_df <- morans_long %>%
  inner_join(gearys_long, by = c("row_id", "meetobjectcode", "parameter")) %>%
  left_join(geom_lookup, by = "meetobjectcode")

# Count for Moran's I significant positive spatial autocorrelation ("strong" = moransI > 0.5 & significant p < 0.05)
morans_counts <- summary_df %>%
  filter(!is.na(moransI), !is.na(moransI_p)) %>%
  mutate(moransI_sigStrong = (moransI > 0.5 & moransI_p < 0.05)) %>%
  group_by(meetobjectcode, parameter) %>%
  summarise(moransI = sum(moransI_sigStrong), .groups = "drop")

# Count for Geary's C significant strong positive spatial autocorrelation ("strong" = gearysC < 0.5 & significant p < 0.05)
gearys_counts <- summary_df %>%
  filter(!is.na(gearysC), !is.na(gearysC_p)) %>%
  mutate(gearysC_sigStrong = (gearysC < 0.5 & gearysC_p < 0.05)) %>% 
  group_by(meetobjectcode, parameter) %>%
  summarise(gearysC = sum(gearysC_sigStrong), .groups = "drop")

# Combine counts and geometry
counts_df <- morans_counts %>%
  full_join(gearys_counts, by = c("meetobjectcode", "parameter")) %>%
  left_join(geom_lookup, by = "meetobjectcode") %>% 
  st_as_sf(sf_column_name = "geom", crs = st_crs(28992)) %>% 
  rename_with(tolower)

# Save
st_write(counts_df, con, layer= "s4_localsac_counts", delete_layer = TRUE)

#### Clean
rm(b_ls, counts_df, df_long, df_ls, df_res, districts_nzv, districts_sf, edges,
   gearys_counts, gearys_long, geom_lookup, meetobjecten_oi, morans_counts,
   morans_long, mw_df, mw_df_all, mw_df_z, mw_df_z_unique, net, nodes, nodes_d,
   res, res_df, res_ls, summary_df, t, window_i_c, fp_edges, fp_nodes, gearys_cols,
   gearys_p_cols, i, j, m, morans_cols, morans_p_cols, params, 
   calculate_local_autocorrelation_with_significance, global_autocorrelation,
   pmf_uncertainty)
   