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
# Script 06: Visualizations
#
# This script creates the figures of the final report. Heading numbers are linked
#   to the result chapters of the report.

############################################################################## #
#### 1 Preparation #############################################################
############################################################################## #

# Load
globalSAC_df <- dbGetQuery(con, "SELECT * FROM s4_globalsac")
localSAC_counts_df <- st_read(con, "s4_localsac_counts")
localSAC_df <- st_read(con, "s4_localsac")
mw_df <- dbGetQuery(con, "SELECT * FROM s2_mv_windows")
districts_sf <- st_read(con, "s0_hea_districts_raw") %>% select(district)
districts_nzv <- st_read(con, "s0_nzv_districts_raw")
df_fl <- st_read(con, "s5_fragments")

fm_df_models <- dbGetQuery(con, "SELECT * FROM s0_factor_mapping_intra_split")
fm_df_splits <- dbGetQuery(con, "SELECT * FROM s0_factor_mapping_inter_split")

# Pre-processing
districts_sf <- districts_sf %>% 
  rbind(districts_nzv)
districts_sf <- districts_sf %>% 
  st_transform("EPSG:3857")

# Prepare basemap
bbox <- st_bbox(districts_sf)
bbox_sf <- st_sf(st_as_sfc(bbox)) %>%
  st_buffer(200) %>% 
  st_transform("EPSG:3857")
set_defaults(map_service = "esri", map_type = "world_light_gray_base")
basemap_magick(bbox_sf)

############################################################################## #
#### 2 Methods #################################################################
############################################################################## #

#### 2.3.1 Concentration means over time per heavy metal #######################

# Calculate mean concentration for each parametercode and split_12
mean_df <- mw_df %>%
  group_by(parametercode, split_12) %>%
  summarise(
    mean_value = mean(numeriekewaarde, na.rm = TRUE),
    sd_value = sd(numeriekewaarde, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
p <- ggplot(mean_df, aes(x = split_12, y = mean_value, group = 1)) +
  geom_ribbon(aes(ymin = ifelse((mean_value - sd_value) < 0, 0, (mean_value-sd_value)), ymax = mean_value + sd_value), 
              fill = "#E69F00", alpha = 0.3) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_point(size = 2, color = "#D55E00") +
  labs(
    x = "Time [splits]",
    y = "Mean Concentration [µg/l] ± SD"
  ) +
  facet_wrap(~parametercode, ncol = 4, nrow = 2, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
  theme_minimal(); p

# Save
fp <- paste0(fp_save, "c_over_time.png")
ggsave(fp, p, width = 150, height = 100, units = "mm", dpi = 200)

############################################################################## #
#### 3 Spatial Autocorrelation #################################################
############################################################################## #

#### 3.1 Heatmap of global SAC #################################################

plot_global_sac <- function(df, c){
  
  # create a color scale: blue (negative), orange (positive), grey (not significant)
  if (grepl("c", c, fixed = T)){
    stat <- "gearys_C"
      
    name <- "Geary's C"
    
    heatmap_palette <- c("#E69F00", "white")
    breaks <- c(0, 0.5, 1)
    labels <- c("0", "0.5", "1")
    limits <- c(0,1)
  } else {
    stat <- "morans_I"
    
    name <- "Moran's I"
    
    heatmap_palette <- c("white", "#E69F00")
    breaks <- c(0, 0.5, 1)
    labels <- c("0", "0.5", "1")
    limits <- c(0,1)
  }
  
  df_f <- df %>% 
    filter(stat_type == stat)
  
  # plot
  p <- ggplot(df_f, aes(x=factor(split), y=fct_rev(parametercode), fill=stat_value_masked)) +
    geom_tile(color="grey60", size=0.2, na.rm=FALSE) +
    scale_fill_gradientn(colors=heatmap_palette, 
                         limits=limits, 
                         breaks=breaks, 
                         labels=labels,
                         na.value="grey80", 
                         name=name) +
    theme_minimal(base_size=13) +
    labs(x="Split", y="Parameter", title=paste("Global", name)) +
    theme(axis.text.x = element_text(hjust=1)); p
  
  # save
  fp <- paste0(fp_save, "global_sac_", c, ".png")
  ggsave(fp, p, width = 140, height = 70, units = "mm", dpi = 200)
  
  message("Saved: ", fp)
}

# Plot and save
plot_global_sac(globalSAC_df, "moransi")
plot_global_sac(globalSAC_df, "gearysc")

#### 3.1 Summary maps of local SAC #############################################

# sources: 
#   https://stackoverflow.com/questions/31876236/maps-with-facet-wrap-in-ggplot2
#   https://stackoverflow.com/questions/3681647/ggplot-how-to-increase-spacing-between-faceted-plots
#   https://jakob.schwalb-willmann.de/basemaps/
#   https://www.datanovia.com/en/blog/how-to-change-ggplot-facet-labels/
#   https://stackoverflow.com/questions/17825201/is-there-a-way-to-increase-the-height-of-the-strip-text-bar-in-a-facet

  
plot_local_sac_sum <- function(df, c){
  
  # change CRS of input and ensure c is numeric
  df_tr <- st_transform(df, "EPSG:3857") %>% 
    mutate(!!c := as.numeric(.data[[c]]))
  
  # define things for the plot
  if (grepl("c", c, fixed = T)){
    name <- "Local Geary's C"
  } else {
    name <- "Local Moran's I"
  }
  title = paste("Number of times a point showed significant strong positive spatial autocorrelation for", name)
  
  # plot
  p <- ggplot() + 
    basemap_gglayer(bbox_sf)+
    scale_fill_identity() +
    coord_sf() +
    theme_map() +
    ggtitle(title)+
    geom_sf(data=districts_sf, fill = NA, col = "skyblue", linewidth =0.5)+
    geom_sf(data=df_tr, aes(col = .data[[c]]), size = 2)+
    scale_color_gradient(
      name = "Count",
      low = "#F0E442", high = "#D55E00",
      limits = c(0, 12), breaks = 0:12
    )+
    theme(legend.position="right",
          panel.spacing = unit(1, "lines"),
          strip.background = element_rect(
            fill="white", linetype = "blank"),
          strip.text.x = element_text(
            size = 14, lineheight=3.0
          ),
          legend.key.height = unit(2, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          plot.title = element_text(size=14, face = "bold", hjust = 0.5),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          plot.caption = element_text(size = 10))+
    facet_wrap(~parameter, nrow=2, ncol=4)+
    labs(caption = "Esri (2025)"); p
  
  # save
  fp <- paste0(fp_save, "local_sac_sum_", c, ".png")
  ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)
  
  message("Saved: ", fp)
}

# Plot and save
plot_local_sac_sum(localSAC_counts_df, "gearysc")
plot_local_sac_sum(localSAC_counts_df, "moransi")

############################################################################## #
#### 3 Source identification ###################################################
############################################################################## #

#### 3.2 PMF table #############################################################

prep_pmf_table <- function(fp, s = 7){
  
  df <- tibble(
    file = fp,
    u = str_extract(fp, "^(gearys|normals)") |> str_remove("s$"),
    split = as.numeric(str_extract(fp, "(?<=y|l)s\\d+") |> str_remove("s")),
    k = as.numeric(str_extract(fp, "(?<=k)\\d+"))
  ) %>% 
    select(-file)
  
  #### PMF fps
  if (grepl("bs100", fp)){
    start_fp <- paste0(pmf_path, "results100/")
    bs <- "bs100"
  } else {
    start_fp <- paste0(pmf_path, "results/")
    bs <- ""
  }
  base_fp <- paste0(start_fp, fp)
  diag_fp <- paste0(start_fp, df$u, "s", df$split, "k", df$k, bs, "_diagnostics.xlsx")
  error_fp <- paste0(start_fp, df$u, "s", df$split, "k", df$k, bs, "_BaseErrorEstimationSummary.xlsx")
  
  # read PMF
  base_df_rc <- read_xlsx(base_fp, sheet = "Run Comparison")
  base_df_r <- read_xlsx(base_fp, sheet = "Residuals", col_names = F)
  diag_df_i <- read_xlsx(diag_fp, sheet = "Input", skip = 2)
  diag_df_br <- read_xlsx(diag_fp, sheet = "Base Runs")
  diag_df_bo <- read_xlsx(diag_fp, sheet = "Bootstrap Runs")
  diag_df_ees <- read_xlsx(error_fp)
  
  # High scaled residuals
  start <- which(diag_df_br[[1]] == "Scaled residuals beyond 3 standard deviations (IDs by species):")
  end   <- which(diag_df_br[[1]] == "Scaled residuals beyond 3 standard deviations (species by ID):")
  
  r <- nrow(diag_df_br %>% 
    slice((start+2):(end-2)))
  
  # Minimum S/N-ratio
  sn_min <- diag_df_i %>% 
    summarise(min_sn = min(`S/N`, na.rm = TRUE)) %>% 
    pull(min_sn)
  
  # Number of samples
  start <- which(base_df_r[[1]] == "Residuals from Base Run #1 (Convergent Run)")
  end <- which(base_df_r[[1]] == "Scaled Residuals from Base Run #1 (Convergent Run)")
  
  n <- nrow(base_df_r %>% 
                 slice((start+2):(end-2)))
    
  # BS %
  start <- which(diag_df_bo[[1]] == "Mapping of bootstrap factors to base factors:")
  end <- which(diag_df_bo[[1]] == "Q(Robust) Percentile Report:")
  
  # Map k to Excel column
  col_letter <- LETTERS[1 + df$k]
  
  diag_df_bo <- read_xlsx(diag_fp, sheet = "Bootstrap Runs", range = paste0("B", (start+2), ":", col_letter, (end-1)))
  
  bs_r <- diag_df_bo %>%
    summarise(
      bs_r = sum(across(everything()), na.rm = TRUE)
    ) %>% 
    pull(bs_r)
  
  bs_s <- diag_df_bo %>%
    { 
      df <- .
      # factor columns = all except "Unmapped"
      diag_cols <- setdiff(names(df), "Unmapped")
      
      # total bootstrap runs (sum of everything)
      bs_r <- sum(unlist(df), na.rm = TRUE)
      
      # sum of diagonal entries: convert factor cols to matrix and take diag()
      mat <- as.matrix(select(df, all_of(diag_cols)))
      diag_sum <- sum(diag(mat), na.rm = TRUE)
      
      # success rate (diagonal / total)
      success_rate <- diag_sum / bs_r
      
      tibble(bs_r = bs_r,
             diag_sum = diag_sum,
             success_rate = success_rate) %>% 
        pull(success_rate)
    }

  # R2
  start <- which(diag_df_br[[2]] == "Species")[1]
  
  df_r <- read_xlsx(diag_fp, sheet = "Base Runs", 
                           range = paste0("B", (start+1), ":H", (start+s))) %>%
    summarise(
      below_05 = sum(`r^2` < 0.5, na.rm = TRUE),
      below_07 = sum(`r^2` < 0.7, na.rm = TRUE)
    )
  
  # DISP
  start <- which(diag_df_ees[[1]] == "DISP Diagnostics:")
  
  disp_q <- diag_df_ees %>% 
    slice((start+3)) %>% 
    select(3) %>% 
    pull() %>% 
    as.numeric()
  
  disp_c <- diag_df_ees %>% 
    slice((start+4)) %>% 
    select(3:ncol(diag_df_ees)) %>% 
    mutate_if(is.character,as.numeric) %>% 
    rowSums(na.rm = T)
  
  # BS-DISP
  start <- which(diag_df_ees[[1]] == "BS-DISP Diagnostics:")
  
  if (length(start)==0){
    bs_disp_r <- NA
    bs_disp_c <- NA
  } else {
    bs_disp_r <- diag_df_ees %>% 
      slice((start+1)) %>% 
      select(3) %>% 
      pull() %>% 
      as.numeric()
    
    bs_disp_c <- diag_df_ees %>% 
      slice((start+8)) %>% 
      select(3:ncol(diag_df_ees)) %>% 
      mutate_if(is.character,as.numeric) %>% 
      rowSums(na.rm = T)
  }
  
  # Q
  best_run <- base_df_rc %>%  
    slice(1) %>% 
    pull(2) %>% 
    str_split(pattern = " ") %>% 
    pluck(1) %>% 
    pluck(7) %>% 
    as.numeric()
  
  start <-  which(diag_df_br[[1]] == "Base run summary table:")
  df_q <- diag_df_br %>% 
    slice((start+best_run+1))
  
  # Add to df
  df2 <- df %>% 
    mutate(s = s,
           sn_min = sn_min,
           n = n,
           r = r,
           r2_less05 = df_r$below_05,
           r2_less07 = df_r$below_07,
           bs_s = bs_s,
           bs_r = bs_r,
           disp_q_less1 = disp_q < 1,
           disp_c = disp_c,
           bs_disp_r = bs_disp_r,
           bs_disp_c = bs_disp_c,
           df = n*s - k*(n+s),
           Qexp = df,
           Qrobust = as.numeric(df_q$'...3'),
           Qtrue = as.numeric(df_q$'...4'),
           Qdelta = as.numeric(df_q$'...7'),
           Qdelta2 = Qtrue/Qexp
    )
  
  return(df2)
}

# repeat for bs = 50 and bs = 100
f_ls50 <- list.files(paste0(pmf_path, "results/"), pattern = "base"); f_ls50
f_ls100 <- list.files(paste0(pmf_path, "results100/"), pattern = "base"); f_ls100

df_bs50 <- lapply(f_ls50, prep_pmf_table)
df_bs50 <- bind_rows(df_bs50) %>% 
  mutate(Qdelta1 = 1-Qdelta)
dbWriteTable(con, "s6_bs50_table", df_bs50, overwrite = TRUE)

df_bs100 <- lapply(f_ls100, prep_pmf_table)
df_bs100 <- bind_rows(df_bs100) %>% 
  mutate(Qdelta1 = 1-Qdelta)
dbWriteTable(con, "s6_bs100_table", df_bs100, overwrite = TRUE)

#### 3.2 Q vs k ################################################################

df_bs50_l <- df_bs50 %>%
  mutate(
    Model = recode(u,
                   "normal" = "Regular",
                   "geary" = "SAC-corrected"),
    Model = factor(Model, levels = c("Regular", "SAC-corrected"))
  )

# Custom colors for the models, matching new labels/order
model_colors <- c("Regular" = "#0072B2", "SAC-corrected" = "#E69F00")

p <- ggplot(df_bs50_l, aes(x = k, y = Qdelta, color = Model)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  xlab("Number of factors (k)") +
  ylab("Qtrue/ Qexpected") +
  scale_color_manual(
    name = "Model",
    values = model_colors
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 28),     
    axis.text  = element_text(size = 24),     
    legend.title = element_text(size = 28), 
    legend.text = element_text(size = 24)      
  ); p

fp <- paste0(fp_save, "k_vs_Q.png")
ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)

#### 3.3 Pretty PMF table ######################################################

pretty_table <- function(df, type = "small"){
  
  df2 <- df %>% 
    select(-c(s, bs_r, disp_q_less1, bs_disp_r, df, Qrobust, Qdelta2, Qdelta1)) %>% 
    mutate(u = str_to_title(u),
           sn_min = round(sn_min, 2),
           bs_s = round(100*bs_s, 2),
           Qtrue = round(Qtrue, 2),
           Qdelta = round(Qdelta, 2)) 
    
  if (type == "small"){
    df3 <- df2 %>% 
      select("Dataset" = u, "Split" = split, "Factors (k)" = k, "Minimum S/N" = sn_min,
             n, "High residuals" = r, "r2 < 0.5" = "r2_less05", "r2 < 0.7" = "r2_less07",
             "Bootstrap success rate" = bs_s, "DISP Factor swaps in first dQmax" = disp_c,
             Qexp, Qtrue, "Qtrue/ Qexp" = Qdelta)
      
      fp = paste0(fp_save, "pmf_bs50.xls")
      write.xlsx(df3, fp)
  } else {
    df3 <- df2 %>% 
      select("Dataset" = u, "Split" = split, "Factors (k)" = k, "Minimum S/N" = sn_min,
             n, "High residuals" = r, "r2 < 0.5" = "r2_less05", "r2 < 0.7" = "r2_less07",
             "Bootstrap success rate" = bs_s, "DISP Factor swaps in first dQmax" = disp_c,
             "BS-DISP Factor swaps in first dQmax" = bs_disp_c, Qexp, Qtrue, "Qtrue/ Qexp" = Qdelta)
    
    fp = paste0(fp_save, "pmf_bs100.xls")
    write.xlsx(df3, fp)
  }

}

pretty_table(df_bs50)

#### 3.4 Stacked bar plots #####################################################

folder_fp <- paste0(pmf_path, "results100/")

reorder_by_mapping <- function(df, mapping) {
  df_long <- df %>%
    pivot_longer(cols = starts_with("F"), names_to = "old_factor", values_to = "value") %>%
    # Join mapping to get new_factor for the model
    left_join(mapping, by = c("model", "old_factor")) %>%
    # Replace old_factor by new_factor; for rows without a mapping, keep old_factor
    mutate(new_factor = if_else(is.na(new_factor), old_factor, new_factor)) %>%
    select(-old_factor) %>%
    # Pivot wider by the new_factor
    pivot_wider(names_from = new_factor, values_from = value) %>%
    # Keep all original columns besides F columns (Parameter, model)
    select(Parameter, F1, F2, F3, F4, F5, model)
  
  return(df_long)
}

plot_sankey <- function(s, mapping_models, mapping_splits){  
  
  df_ls <- list()
  
  for (m in c("normals", "gearys")){
    
    # create fp-lists
    base_ls <- list.files(folder_fp, pattern = "base")
    
    # create fps
    base_fp <-  paste0(folder_fp, base_ls[grepl(paste0(m, s), base_ls)])
  
    # load data
    base_rc_df <- read_xlsx(base_fp, sheet = "Run Comparison")
    base_p_df <- read_xlsx(base_fp, sheet = "Profiles")
    
    # determine best run
    best_run <- base_rc_df %>% 
      slice(1) %>% 
      pull(2) %>% 
      str_split(pattern = " ") %>% 
      pluck(1) %>% 
      pluck(7) %>% 
      as.numeric()
    
    # get information of best run
    row_ids <- base_p_df %>% 
      mutate(row_id = row_number()) %>%
      filter(.[[1]] == paste0("Factor Profiles (% of species sum) from Base Run #", best_run," (Convergent Run)")) %>% 
      pull(row_id)
    
    sel_df <- base_p_df %>% 
      slice((row_ids+3):(row_ids+9)) %>%
      select(2:7) %>% 
      rename_with(~ "Parameter", 1) %>%
      mutate(Parameter = str_remove_all(Parameter, '\\"')) %>%
      rename_with(~ paste0("F", 1:5), 2:6) %>% 
      mutate(model = m)
    
    df_ls <- append(df_ls, list(sel_df))
  }
  
  df <- bind_rows(df_ls)
  # re-order factors
  df_reordered_models <- df %>%
    filter(model == "gearys") %>%
    reorder_by_mapping(mapping_models) %>%
    bind_rows(df %>% filter(model != "gearys"))
  df_reordered_splits <- reorder_by_mapping(df_reordered_models, mapping_splits)
  
  # find difference
  df_dif <- df_reordered_splits %>%
    # 1. Pivot longer to gather F1 to F5 factors into one column
    pivot_longer(cols = starts_with("F"), names_to = "Factor", values_to = "Value") %>%
    # 2. Pivot wider to get separate columns for each model's value
    pivot_wider(names_from = model, values_from = Value) %>%
    # 3. Calculate difference (e.g., gearys - normals)
    mutate(Difference = as.numeric(gearys) - as.numeric(normals)) %>%
    # Optional: arrange to see large differences on top
    arrange(desc(abs(Difference)))
  
  # Example Wong color palette (only 5 used)
  wong_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
  
  model_order <- c("normals", "gearys")
  factor_order <- c("F1", "F2", "F3", "F4", "F5")
  
  df_all_long <- df_reordered_splits %>%
    pivot_longer(cols = starts_with("F"), names_to = "Factor", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = c("As", "B", "Ba", "Co", "Se", "U", "Zn")),
           Value = as.numeric(Value)) %>%
    mutate(
      # Parameter = factor(Parameter, levels = param_order),
      model = factor(model, levels = model_order),
      Factor = factor(Factor, levels = factor_order)
    )
  
  # Clustered + stacked bar chart
  p <- ggplot(df_all_long, aes(x = model, y = Value, fill = Factor)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = wong_colors) +
    labs(
      x = "Normal model (left) vs SAC-corrected model (right) per Parameter",
      y = "Contribution fo factors [%]",
      fill = "Factor",
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x=element_blank(),
      axis.text.y = element_text(size = 24), 
      axis.title.y = element_text(size = 28), 
      axis.title.x = element_text(size = 28), 
      strip.text = element_text(size = 28),   
      legend.title = element_text(size = 28),
      legend.text = element_text(size = 24),                
    ) + 
    facet_wrap(~Parameter, ncol = 7, strip.position="bottom"); p
  
  fp <- paste0(fp_save, "param_bar_split", s, ".png")
  ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)
  message("Saved: ", fp)

  
  df_bar <- df_dif %>% 
    filter(Difference >= 0) %>%  
    group_by(Parameter) %>% 
    summarize(p_sum = sum(Difference)) %>%
    mutate(
      Parameter = factor(Parameter, levels = rev(sort(unique(Parameter))))
    )
  
  return(df_bar)
}

bar_ls <- list()

for (s in splits){
  
  # Pivot the mapping longer to join it easily
  mapping_models <- fm_df_models %>%
    filter(split == s) %>% 
    pivot_longer(cols = starts_with("F"),
                 names_to = "old_factor",
                 values_to = "new_factor")
  mapping_splits <- fm_df_splits %>%
    filter(split == s) %>% 
    pivot_longer(cols = starts_with("F"),
                 names_to = "new_factor",
                 values_to = "old_factor")
  
  t <- plot_sankey(s, mapping_models, mapping_splits)
  bar_ls <- append(bar_ls, list(t))
}

df_t <- bind_rows(bar_ls) %>% 
  mutate(
  Parameter = factor(Parameter, levels = sort(unique(Parameter)))
  )

wong_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(df_t, aes(x=p_sum, y = Parameter, fill = Parameter))+
  geom_boxplot() +
  theme_minimal()+
  scale_fill_manual(values = wong_colors) +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("Difference between models [%]")+
  theme(
    axis.title = element_text(size = 28),      # axis titles
    axis.text  = element_text(size = 24),      # axis tick labels
    legend.title = element_text(size = 28),    # legend title
    legend.text = element_text(size = 24)      # legend items
  ); p
fp <- paste0(fp_save, "param_dif_all.png")
ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)

############################################################################## #
#### 5.3 Source localization ###################################################
############################################################################## #

#### 5.3 Area comparison per split #############################################

plot_area_comparison <- function(df, m){
  
  df <- st_transform(df, "EPSG:3857") %>% 
    filter(model == m) %>% 
    select(split)
  
  p <- ggplot()+
    basemap_gglayer(bbox_sf)+
    scale_fill_identity() +
    coord_sf() +
    theme_map() +
    geom_sf(data=df, fill = "#D55E00", alpha = 0.25)+
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)); p
  # save
  fp <- paste0(fp_save, "area_", m, ".png")
  ggsave(fp, p, width = 200, height = 210, units = "mm", dpi = 300)
  
  message("Saved: ", fp)
}

# Run and save
for (m in c("normals", "gearys")){
  plot_area_comparison(df_fl, m)
}

#### 5.3 Environmental make-up pers plit #######################################

# Pre-processing
df_mu <- df_fl %>% 
  select(split, model, factor, ends_with("_perc"), ends_with("_percent_area")) %>% 
  st_drop_geometry() %>% 
  group_by(split, model) %>% # previously also factor 
  summarize(across(ends_with("_perc"),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "{.col}"),
            across(ends_with("_percent_area"),
                   ~ mean(.x, na.rm = TRUE),
                   .names = "{.col}")) %>%
  rename_with(~ sub("_perc$", "", .x), ends_with("_perc"),
              ~ sub("_percent_area$", "", .x), ends_with("_percent_area"))

# Soil

df_soil <- df_mu %>% 
  select(split, model, starts_with("soil_"))

df_soil_n <- df_soil %>% 
  filter(model == "normals") %>% 
  arrange(split)
df_soil_g <- df_soil %>% 
  filter(model == "gearys") %>% 
  arrange(split)

t <- abs(df_soil_n[3:ncol(df_soil_n)] - df_soil_g[3:ncol(df_soil_g)])

df_long <- t %>% 
  rename(
    marine_clay = soil_zeekleigronden,
    other = soil_other,
    calcareous_sandy_soil = soil_kalkhoudende_zandgronden,
    thick_earth_soil = soil_dikke_eerdgronden,
    calcareous_loose_sandy_soil = soil_kalkloze_zandgronden,
    moorish_soil = soil_moerige_gronden,
    podzol_soil = soil_podzolgronden,
    peat_soil = soil_veengronden,
    boulder_clay_soil = soil_keileemgronden,
    defined_associations_soil = soil_gedefinieerde_associaties
  ) %>% 
  pivot_longer(cols = everything(),
               names_to = "soil_type",
               values_to = "value")

# Calculate range per soil_type
range_order <- df_long %>%
  group_by(soil_type) %>%
  summarise(range_value = max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) %>%
  arrange(desc(range_value)) %>%
  pull(soil_type)

# Reorder soil_type factor based on decreasing range
df_long <- df_long %>%
  mutate(soil_type = factor(soil_type, levels = range_order))

# Create a named vector for labels with nicer formatting
label_map <- c(
  marine_clay = "Marine Clay",
  other = "Other",
  calcareous_sandy_soil = "Calcareous Sandy Soil",
  thick_earth_soil = "Thick Earth Soil",
  calcareous_loose_sandy_soil = "Calcareous Loose Sandy Soil",
  moorish_soil = "Moorish Soil",
  podzol_soil = "Podzol Soil",
  peat_soil = "Peat Soil",
  boulder_clay_soil = "Boulder Clay Soil",
  defined_associations_soil = "Defined Associations Soil"
)
# Create boxplot
p <- ggplot(df_long, aes(x = soil_type, y = value, fill = "#D55E00")) +
  coord_flip()+
  geom_boxplot() +
  scale_x_discrete(labels = label_map) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size = 28),
    axis.text.y = element_text(size = 24),       # y-axis tick labels
    axis.title.y = element_text(size = 28),  # y-axis title
    axis.title.x = element_text(size = 28),
    legend.position = "none"
  )+
  labs(x = "Soil Type", y = "Difference [%]"); p

# save
fp <- paste0(fp_save, "area_content_soil.png")
ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)
message("Saved: ", fp)

# Land use

df_lu <- df_mu %>% 
  select(split, model, starts_with("lu_"))

df_lu_n <- df_lu %>% 
  filter(model == "normals") %>% 
  arrange(split)
df_lu_g <- df_lu %>% 
  filter(model == "gearys") %>% 
  arrange(split)

t <- abs(df_lu_n[3:ncol(df_lu_n)] - df_lu_g[3:ncol(df_lu_g)])
t2 <- t[, names(sort(colSums(t[3:45]), decreasing = T)[1:10])]

df_long <- t2 %>% 
  rename(
    grains = lu_granen_percent_area,
    grass_in_primary_built_up_area = lu_gras_in_primair_bebouwd_gebied_percent_area,
    fresh_water = lu_zoet_water_percent_area,
    # maize = lu_maïs_percent_area,
    other_crops = lu_overige_landbouwgewassen_percent_area,
    # raised_bog = lu_hoogveen_percent_area,
    leafy_forest = lu_loofbos_percent_area,
    beets = lu_bieten_percent_area,
    # grass_in_secondary_built_up_area = lu_gras_in_secundair_bebouwd_gebied_percent_area,
    forest_in_raised_bog_area = lu_bos_in_hoogveengebied_percent_area,
    grass_in_primary_built_up_area = lu_gras_in_primair_bebouwd_gebied_percent_area,
    potatoes = lu_aardappelen_percent_area,
    other_marsh_vegetation = lu_overige_moeras_vegetatie_percent_area,
    reed_vegetation = lu_rietvegetatie_percent_area
    # buildings_in_primary_built_up_area = lu_bebouwing_in_primair_bebouwd_gebied_percent_area
  ) %>% 
  pivot_longer(cols = everything(),
               names_to = "landuse_type",
               values_to = "value")

# Calculate range per soil_type
range_order <- df_long %>%
  group_by(landuse_type) %>%
  summarise(range_value = max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) %>%
  arrange(desc(range_value)) %>%
  pull(landuse_type)

# Reorder soil_type factor based on decreasing range
df_long <- df_long %>%
  mutate(soil_type = factor(landuse_type, levels = range_order))

# Create a named vector for labels with nicer formatting
label_map <- c(
  grains = "Grains",
  grass_in_primary_built_up_area = "Grass in Primary Built-up Area",
  fresh_water = "Fresh Water",
  maize = "Maize",
  other_crops = "Other Crops",
  raised_bog = "Raised Bog",
  leafy_forest = "Leafy Forest",
  beets = "Beets",
  grass_in_secondary_built_up_area = "Grass in Secondary Built-up Area",
  buildings_in_primary_built_up_area = "Buildings in Primary Built-up Area",
  forest_in_raised_bog_area = "Forest in Raised Bog Area",
  potatoes = "Potatoes",
  reed_vegetation = "Reed Vegetation",
  other_marsh_vegetation = "Other Marsh Vegetation"
)

# Create boxplot
p <- ggplot(df_long, aes(x = landuse_type, y = value, fill = "#D55E00")) +
  coord_flip()+
  geom_boxplot() +
  scale_x_discrete(labels = label_map) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size = 28),
    axis.text.y = element_text(size = 24),       # y-axis tick labels
    axis.title.y = element_text(size = 28),  # y-axis title
    axis.title.x = element_text(size = 28),
    legend.position = "none"
  )+
  labs(x = "Landuse Type", y = "Difference [%]"); p

# save
fp <- paste0(fp_save, "area_content_landuse.png")
ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)
message("Saved: ", fp)

# Functional use

df_fu <- df_mu %>% 
  select(split, model, starts_with("fu_"))

df_fu_n <- df_fu %>% 
  filter(model == "normals") %>% 
  arrange(split)
df_fu_g <- df_fu %>% 
  filter(model == "gearys") %>% 
  arrange(split)

t <- abs(df_fu_n[3:ncol(df_fu_n)] - df_fu_g[3:ncol(df_fu_g)])
t2 <- t[, names(sort(colSums(t[3:ncol(t)]), decreasing = T)[1:10])]

df_long <- t2 %>% 
  rename(
    national_park = fu_nationaal_park,
    nature_reserve = fu_natuurgebied,
    holiday_park = fu_vakantiepark,
    harbor = fu_haven,
    # bungalow_park = fu_bungalowpark,
    # cemetery = fu_begraafplaats,
    greenhouse_area = fu_kassengebied,
    # wind_farm = fu_windturbinepark,
    sports_complex = fu_sportterrein._sportcomplex,
    building_complex = fu_gebouwencomplex,
    park = fu_park,
    production_facility = fu_productie.installatie,
    ice_rink = fu_ijsbaan
  ) %>% 
  pivot_longer(cols = everything(),
               names_to = "fu_type",
               values_to = "value")

# Calculate range per soil_type
range_order <- df_long %>%
  group_by(fu_type) %>%
  summarise(range_value = max(value, na.rm = TRUE) - min(value, na.rm = TRUE)) %>%
  arrange(desc(range_value)) %>%
  pull(fu_type)

# Reorder soil_type factor based on decreasing range
df_long <- df_long %>%
  mutate(soil_type = factor(fu_type, levels = range_order))

# Create a named vector for labels with nicer formatting
label_map <- c(
  national_park = "National Park",
  nature_reserve = "Nature Reserve",
  holiday_park = "Holiday Park",
  harbor = "Harbor",
  bungalow_park = "Bungalow Park",
  cemetery = "Cemetery",
  greenhouse_area = "Greenhouse Area",
  wind_farm = "Wind Farm",
  sports_complex = "Sports Complex",
  building_complex = "Building Complex",
  park = "Park",
  production_facility = "Production Facility",
  ice_rink = "Ice Rink"
)
# Create boxplot
p <- ggplot(df_long, aes(x = fu_type, y = value, fill = "#D55E00")) +
  coord_flip()+
  geom_boxplot() +
  scale_x_discrete(labels = label_map) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(size = 28),
    axis.text.y = element_text(size = 24),       # y-axis tick labels
    axis.title.y = element_text(size = 28),  # y-axis title
    axis.title.x = element_text(size = 28),
    legend.position = "none"
  )+
  labs(x = "Functional Uses", y = "Difference [%]"); p
# save
fp <- paste0(fp_save, "area_content_function.png")
ggsave(fp, p, width = 297, height = 210, units = "mm", dpi = 300)
message("Saved: ", fp)
