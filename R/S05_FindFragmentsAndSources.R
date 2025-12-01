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
# Script 05: Identify fragments and localize sources

# This script identifies water network fragments of high-concentration-drop edges 
#   and looks for potential sources around them. 

############################################################################## #
#### 0 Preparation #############################################################
############################################################################## #

# Load
flowdirection_hea <- st_read(con, query = "SELECT * FROM s1_flow")
mo_oi <- st_read(con, query = "SELECT * FROM s1_mo_oi")
rwzis_sf <- st_read(con, query = "
                 SELECT * 
                 FROM s0_rwzis_raw e
                 JOIN s1_research_area r ON ST_Within(e.geometry, r.geom)
                 WHERE e.statusfunctioneren = 'http://data.gwsw.nl/1.6/totaal/InGebruik'
                 ") 

soil_sf <- st_read(con, query = "SELECT * FROM s0_soil_raw") %>% 
  group_by(mainsoilclassification) %>%
  summarize() 

lu_r <- pgGetRast(con, "s0_landuse10_raw")

fu_sf <- st_read(con, query = "SELECT * FROM s0_top10nl_functional_raw") %>% 
  group_by(typefunctioneelgebied) %>%
  summarize() 

lgn_lut <- dbGetQuery(con, "SELECT * FROM s0_landuse10_lut")

# Pre-processing
net_hea <- as_sfnetwork(flowdirection_hea, directed = T) %>% 
  activate("edges") 


# Make sure codes are integers and in same domain as raster values
lgn_lut$code <- as.integer(lgn_lut$code)

# Build levels: data.frame with columns value (code) and label
levels_df <- data.frame(value = lgn_lut$code, label = as.character(lgn_lut$lgn2023))

# Assign levels to the raster (this creates factor levels)
levels(r_f) <- list(levels_df)

cell_area <- res(lu_r)[1] * res(lu_r)[2]

# Globals
d_up <- 1500
d_buffer <- 200
k <- 5 # number of profiles in PMF - the same for all my cases

############################################################################## #
#### 1 Create function #########################################################
############################################################################## #

find_fragments_and_sources <- function(i, metric, k, d_up, d_buffer){  
  
  # Read pmf data
  pmf_fp <- paste0(pmf_path, "results100/", metric,i, "k", k, "bs100_base.xlsx")
  
  # Find best run
  pmf_df <- read_xlsx(pmf_fp, sheet = "Run Comparison")
  best_run <- pmf_df %>% 
    slice(1) %>% 
    pull(2) %>% 
    str_split(pattern = " ") %>% 
    pluck(1) %>% 
    pluck(7) %>% 
    as.numeric()
  
  nr_f <- k
  
  # Create a list of expressions for each factor column dynamically
  factor_exprs <- lapply(1:nr_f, function(i) {
    col_name <- paste0("F", i)
    expr <- rlang::expr(
      rowSums(across(ends_with(!!paste0("_", col_name)), as.numeric), na.rm = TRUE)
    )
    rlang::set_names(list(expr), col_name)
  }) |> unlist(recursive = FALSE)
  
  # read pmf data
  pmf_fp <- paste0(pmf_path, "results100/", metric, i, "k", k, "bs100_sourcecontributions.xls")
  pmf_df <- read_xls(pmf_fp)
  
  nodes_fp <- paste0("temp/net_split_", i , "_nodes.gpkg")
  nodes <- st_read(nodes_fp) %>% 
    select(meetobjectcode)
  n_l <- length(unique(nodes$meetobjectcode))
  
  # get all the rows belonging to the best run
  row_ids <- pmf_df %>%
    mutate(row_id = row_number()) %>%
    filter(grepl("Lowest Q Run", .[[1]])) %>%
    pull(row_id)
  # select those rows
  sel_df <- pmf_df %>% 
    slice((row_ids[best_run]+2):(row_ids[best_run] + ((n_l+3)*7+1)))
  
  # find columns with the right information
  row4 <- as.character(sel_df[2, ])
  # find columns where row 4 contains "N" or "Mean"
  cols_to_keep <- which(grepl("N|Mean", row4))
  # always keep the first column
  cols_to_keep <- unique(c(1, cols_to_keep))
  # subset dataframe and drop na
  filtered_sel_df <- sel_df[, cols_to_keep] 
  
  # find rows of metals
  row_ids <- filtered_sel_df %>%
    mutate(row_id = row_number()) %>%
    filter(grepl(paste(metal_oi, collapse = "|"), .[[1]])) %>%
    pull(row_id)

  for (j in 1:length(row_ids)){ # per metal
    
    # find metal
    metal <- metal_oi[j]
    
    start <- filtered_sel_df %>% 
      mutate(row_id = row_number(),
             location = str_remove_all(`****Source Contribution Statistics by Site****`, '"')) %>%
      filter(location == metal) %>% 
      pull(row_id)
    
    # select those rows
    sel_df <- filtered_sel_df %>% 
      slice(start:(start+n_l+1))
    
    # define column names
    col_names <- c("meetobjectcode", paste0(metal, "_F", c(1:nr_f)))
    col_ids <- c(1, seq(from = 3, to = (3*nr_f), by =3))
    
    # reshape/ select/ bla the df
    df_metal <- sel_df %>% 
      slice(3:n()) %>% 
      select(all_of(col_ids)) 
    names(df_metal) <- col_names
    df_metal <- df_metal %>% 
      mutate(meetobjectcode = gsub('"', '', meetobjectcode)) %>% 
      mutate(across(-meetobjectcode, ~ {
        x <- as.numeric(.)
        ifelse(x < 0, 0, x)
      }))
    
    nodes <- nodes %>% 
      left_join(df_metal)
  }
  
  # Then apply mutate with dynamic expressions
  nodes_factors <- nodes %>%
    mutate(!!!factor_exprs) %>%
    select(meetobjectcode, all_of(paste0("F", 1:nr_f)))
  
  # reshape to long format for easy faceting
  df_long <- nodes_factors %>%
    st_drop_geometry() %>%              
    pivot_longer(cols = starts_with("F"),
                 names_to = "Factor",
                 values_to = "Value")
  
  # calculate mean and sd per factor
  stats <- df_long %>%
    group_by(Factor) %>%
    summarise(
      mean_val = mean(Value, na.rm = TRUE),
      sd_val   = sd(Value, na.rm = TRUE)
    ) %>%
    mutate(line1 = mean_val,
           line2 = mean_val + sd_val,
           line3 = mean_val + 2 * sd_val)
  
  # Map factors
  nodes_long <- nodes_factors %>%
    pivot_longer(
      cols = starts_with("F"),
      names_to = "Factor",
      values_to = "Value"
    )
  
  edges <- net_hea %>%
    activate("edges") %>%
    st_as_sf() %>%
    st_geometry()
  
 
  # List of factors
  factors <- paste0("F", c(1:nr_f))
  f <- "F1"
  f_ls <- list()
  
  for (f in factors){

    nodes_long_sel <- nodes_long %>%
      filter(Factor == f)
    
    #### Find sources ##########################################################
    
    border <- as.numeric(stats[stats$Factor == f,]["line2"])
    
    problem_points <- nodes_factors %>% 
      filter(!!sym(f) >= border)
    
    sf_ls2 <- list()

    for (p in 1:nrow(problem_points)){
      
      poi <- problem_points[p,] 
      
      # Find the nearest network edge to p7
      net_hea <- net_hea %>%
        activate("nodes")
      nearest_node <- st_nearest_feature(poi, net_hea)
      
      # Measure distance between point and nearest edge
      net_hea_n <- net_hea %>% 
        slice(nearest_node)
      dist <- st_distance(poi, net_hea_n)#; dist
      
      if (as.numeric(dist) == 0){
        net_hea_poi <- net_hea
      } else {
        net_hea_poi <- st_network_blend(net_hea, poi) 
      }
    
      l_e <- net_hea_poi %>% 
        activate("edges") %>% 
        mutate(weight = edge_length()) %>% 
        arrange(desc(weight)) %>% 
        slice(1) %>% 
        pull(weight) %>% 
        as.numeric()
      
      net_hea_poi <- net_hea_poi %>% 
        activate("nodes") %>% 
        st_filter(poi, .predicate = st_is_within_distance, dist = l_e+d_up) 
      
      lines_net <- net_hea_poi %>% 
        activate("edges") %>% 
        st_as_sf()
      
      lines_sf <- lines_net %>% 
        filter(birectional == 1)
      reversed_duplicates <- st_reverse(lines_sf)
      
      duplicates_sf <- rbind(lines_net, reversed_duplicates)
      net_hea_poi <- as_sfnetwork(duplicates_sf)
      
      net_hea_poi <- net_hea_poi %>% 
        activate("edges") %>% 
        mutate(weight = edge_length()) %>% 
        activate("nodes") %>% 
        mutate(id = row_number())
      
      # Get the ID of the nearest node to poi (poi is already snapped, but this is safe & simple)
      poi_id <- st_nearest_feature(poi, st_as_sf(net_hea_poi, "nodes"))
      
      t <- net_hea_poi %>% 
        activate("edges") %>% 
        st_as_sf() %>% 
        filter(to == poi_id)
      
      if (nrow(t) == 0){
        sf_fv <- poi
      } else {
        netd <- net_hea_poi %>%
          activate("nodes") %>% 
          mutate(d_to_poi = node_distance_to(poi_id, weights = weight))
        
        directed_clipped <- netd %>%
          activate("edges") %>%
          mutate(
            d_v   = .N()$d_to_poi[to],                             
            w_e   =  as.numeric(edge_length()),
            Lkeep = pmax(0, pmin(w_e, d_up - d_v)), 
            frac_keep = Lkeep / w_e,
            start     = pmax(0, 1 - frac_keep))%>% 
          filter(Lkeep > 0)
        
        edge_df <- directed_clipped %>% 
          activate("edges") %>% 
          st_as_sf()
        
        sf_ls <- list()
        
        for (m in 1:nrow(edge_df)){
          r <- edge_df[m,]
          if (r$frac_keep == 1){
            sf_ls <- append(sf_ls, list(edge_df[m,]))
          } else {
            start <- r$start
            sf_l <- lwgeom::st_linesubstring(r, start, 1)
            sf_ls <- append(sf_ls, list(sf_l))
          }
          
        }
        
        sf_fv <- bind_rows(sf_ls)
      }
      
      sf_fv_buf <- st_buffer(sf_fv, d_buffer)
      st_fv_union <- st_union(sf_fv_buf)
      
      t <- st_as_sf(st_fv_union)
      
      # add information 
      t$split <- i
      t$model <- metric
      t$size <- st_area(t)
      
      #### Intersect ####
      
      # RWZIs
      count_rwzis <- lengths(st_intersects(st_fv_union, rwzis_sf))
      t$rwzis <- count_rwzis
      
      # Soil
      
      # Intersect soil polygons with union polygon
      soil_in_union <- st_intersection(soil_sf, st_fv_union)
      
      # Calculate area of each intersection
      soil_in_union$area <- st_area(soil_in_union)
      
      # Calculate total area (should be same for all) of the union polygon
      fv_area <- st_area(st_fv_union)
      
      # Summarize area by soil class
      summary <- soil_in_union %>%
        mutate(perc = 100 * area/ fv_area) %>% 
        st_drop_geometry()
      
      if (as.numeric(sum(summary$perc)) != 100){
        others <- data.frame("other", (fv_area-sum(summary$area)), (100-as.numeric(sum(summary$perc))))
        names(others) <- names(summary)
        
        summary <- rbind(summary, others)
        }
      
      summary <- summary %>%
        pivot_wider(
          names_from = mainsoilclassification,
          values_from = c(area, perc),
          names_glue = "soil_{tolower(mainsoilclassification)}_{.value}"
        ) %>% 
        rename_with(~ gsub(" ", "_", .x))
      
      t <- cbind(t, summary)
      
      # Functional Uses
      
      # Intersect soil polygons with union polygon
      fu_in_union <- st_intersection(fu_sf, st_fv_union)
      
      # Calculate area of each intersection
      fu_in_union$area <- st_area(fu_in_union)
      
      if (nrow(fu_in_union) != 0){
        summary <- fu_in_union %>% 
          mutate(perc = 100 * area/ fv_area) %>% 
          st_drop_geometry() %>%
          pivot_wider(
            names_from = typefunctioneelgebied ,
            values_from = c(area, perc),
            names_glue = "fu_{tolower(typefunctioneelgebied )}_{.value}"
          )%>% 
          rename_with(~ gsub(" ", "_", .x))
        
        t <- cbind(t, summary)
      }
      
      # Land use
      
      st_fv_v <- vect(st_fv_union)
      
      # Crop and mask raster to polygon extent
      lu_crop <- crop(lu_r, st_fv_v)
      lu_mask <- mask(lu_crop, st_fv_v)
      
      # Count frequency of each raster class inside polygon
      class_freq <- freq(lu_mask, digits=0)
      
      # Calculate class areas
      class_freq <- as.data.frame(class_freq)
      class_freq$area_m2 <- class_freq$count * cell_area
      
      # Calculate total polygon area
      poly_area <- expanse(st_fv_v, unit = "m") 
      
      # Calculate percent area per class
      class_freq$percent_area <- 100 * class_freq$area_m2 / poly_area
      
      class_freq2 <- class_freq %>% 
        select(-c(layer,count)) %>% 
        pivot_wider(
          names_from = value ,
          values_from = c(area_m2, percent_area),
          names_glue = "lu_{tolower(value)}_{.value}"
        ) %>% 
        rename_with(~ gsub(" ", "_", .x)) %>% 
        rename_with(~ str_remove_all(., "\\(|\\)"))
      
      t <- cbind(t, class_freq2)
      
      sf_ls2 <- append(sf_ls2, list(t))
    }
    
    combined_df <- bind_rows(sf_ls2) %>%
      mutate(across(-c(model, x), as.numeric),
             factor = f,
             across(everything(), ~replace_na(.x, 0)))
      
    f_ls <- append(f_ls, list(combined_df))
  }
  
  test2 <- bind_rows(f_ls) %>%
    mutate(across(-c(model, x, factor), as.numeric),
           across(everything(), ~replace_na(.x, 0)))
  
  st_write(test2, paste0("temp/factor_locations_s", i, "_", metric, ".gpkg"), delete_layer = T)
  
  return(test2)
}

# Run
k_list <- list()

for (split in splits){
  for (metric in c("normals", "gearys")){
    
    df <- find_fragments_and_sources(split, metric, k, d_up, d_buffer)
    k_list <- append(k_list, list(df))
  }
}

t_df <- bind_rows(k_list) %>%
  mutate(across(-c(model, x, factor), as.numeric),
         across(everything(), ~replace_na(.x, 0)))

st_write(t_df, con, layer = "s5_fragments", delete_layer = T)

#### Clean
rm(df, flowdirection_hea, fu_sf, k_list, lu_r, mo_oi, net_hea, rwzis_sf, soil_sf,
   t_df, cell_area, d_buffer, d_up, k, metric, split, find_fragments_and_sources,
   lgn_lut)

