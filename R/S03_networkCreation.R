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
# Script 03: Create networks
#
# This script creates the networks.

############################################################################## #
#### 1 Preparation  ############################################################
############################################################################## #

# Load data
flowdirection <- st_read(con, query = "SELECT * FROM s1_flow")
meetobjecten_oi <- st_read(con, query = "SELECT * FROM s1_mo_oi")
mw_df_z <- dbGetQuery(con, "SELECT * FROM s2_mv_z_windows")

# Create base graph
net <- as_sfnetwork(flowdirection, directed = T)

############################################################################## #
#### 2 Create graphs ###########################################################
############################################################################## #

createNetwork <- function(net, pois, split_i){
  # snap pois to network
  net_pois <- st_network_blend(net, pois)
  
  # clean network (remove multi edges and loops)
  net_clean <- net_pois %>% 
    activate("edges") %>% 
    mutate(weight = edge_length()) %>% 
    arrange(weight) %>% 
    filter(!edge_is_loop(),
           !edge_is_multiple()) 
  
  net_smooth <- net_clean
  net_sf <- net_smooth %>% 
    activate("edges") %>% 
    st_as_sf()
  
  lines_sf <- net_smooth %>% 
    activate("edges") %>% 
    st_as_sf() %>% 
    filter(birectional == 1)
  # duplicates <- net_lines[net_lines$birectional == 1,]
  reversed_duplicates <- st_reverse(lines_sf)
  
  duplicates_sf <- rbind(net_sf, reversed_duplicates)
  net_lines <- as_sfnetwork(duplicates_sf) %>% 
    activate("nodes") %>% 
    mutate(id = row_number())
  
  # Find indices of nearest nodes.
  nearest_nodes <- st_nearest_feature(pois, net_lines)
  pois$id <- nearest_nodes
  
  pois_df <- as.data.frame(pois) %>% select(-geometry)
  net_lines_pois <- net_lines %>% 
    activate("nodes") %>% 
    left_join(pois_df)
  
  net_lines_pois_id <- net_lines_pois %>% 
    activate("nodes") %>% 
    filter(!is.na(meetobjectcode)) %>% 
    pull(id)
  
  # remove pseudo nodes to speed up
  net_smooth2 <- convert(net_lines_pois, to_spatial_smooth, protect = net_lines_pois_id) %>% 
    activate("edges") %>% 
    mutate(weight = edge_length()) %>% 
    activate("nodes") %>% 
    mutate(id = row_number())
  
  net_id <- net_smooth2
  
  # define pois id
  pois_id <- net_id %>%
    activate("nodes") %>% 
    filter(!is.na(meetobjectcode)) %>%
    pull(id)
  
  # calculate cost matrix
  cm <- st_network_cost(net_id, from = pois_id, to = pois_id, weights = "weight")
  
  # create empty list
  con_ls <- list()
  iso_ls <- list()
  
  for (i in 1:nrow(cm)){ # per poi
    
    r <- as.vector(cm[i,])
    n <- sum(!is.infinite(r))-1; n
    
    # define neighbour ids in cm
    neighbours <- which(r > 0 & r < Inf)
    neighbours_id <- pois_id[neighbours]
    
    if (n>1){ ### SCENARIO 1: MORE THAN 1 NEIGHBOUR
      # define neighbour ids in network
      for (j in 1:length(neighbours_id)){ # for every neighbour
        # print(paste("j",j))
        node_sel <- neighbours_id[-j]#; node_sel
        
        # create temporary network that removes all neighbors except for j
        net_temp <- net_id %>%
          activate("nodes") %>%
          filter(!(id %in% node_sel))#; net_temp
        
        # CORRECTION FOR REMOVAL
        c_poi <- 0
        c_neighbour <- 0
        if (any(pois_id[i] > node_sel)){
          c_poi <- sum(pois_id[i] > node_sel) #-1 for the neighbour in observation
        }
        if (any(neighbours_id[j] > node_sel)){
          c_neighbour <- sum(neighbours_id[j] > node_sel)
        }
        
        # check whether there is a shortest path possible
        path <- st_network_paths(net_temp, from = (pois_id[i]-c_poi), to = (neighbours_id[j]-c_neighbour))
        
        # check if path exists
        if (length(path[1, "node_paths"][[1]][[1]]) > 0){
          k_ls <- list()
          for (k in path$edge_paths[[1]]){
            t3 <- net_temp %>%
              activate("edges") %>%
              slice(k) %>%
              activate("nodes") %>%
              filter(!node_is_isolated()) %>%
              activate("edges") %>%
              st_as_sf()
            k_ls <- append(k_ls, list(t3))
          }
          
          k_sf <- do.call(rbind, k_ls)
          
          k_coords <- unique(st_coordinates(k_sf)[,c("X", "Y")])
          k_sf2 <- st_sfc(st_linestring(k_coords), crs = "EPSG:28992")
          
          con_ls <- append(con_ls, list(k_sf2))
        }
      }
    } else if (n == 1) { ### SCENARIO 2: ONLY 1 NEIGHBOUR
      
      path <- st_network_paths(net_id, from = pois_id[i], to = neighbours_id)
      
      x_ls <- list()
      for (x in path$edge_paths[[1]]){
        t3 <- net_id %>%
          activate("edges") %>%
          slice(x) %>%
          activate("nodes") %>%
          filter(!node_is_isolated()) %>%
          activate("edges") %>%
          st_as_sf()
        x_ls <- append(x_ls, list(t3))
      }
      
      x_sf <- do.call(rbind, x_ls)
      # 
      x_coords <- unique(st_coordinates(x_sf)[,c("X", "Y")])
      x_sf2 <- st_sfc(st_linestring(x_coords), crs = "EPSG:28992")
      
      con_ls <- append(con_ls, list(x_sf2))
      
    } else if(n == 0){
      
      cm_test <- st_network_cost(net_id, from = pois_id[-i], to = pois_id[i], weights = "weight")
      cm_num <- as.numeric(cm_test)
      
      if (length(cm_num) == sum(is.infinite(cm_num))){
        iso <- net_id %>%
          activate("nodes") %>%
          slice(pois_id[i]) %>%
          pull(meetobjectcode)
        
        iso_ls <- append(iso_ls, list(iso))
      }
      
    }
  }
  
  # combine all linestrings, turn into network and add id
  all_sf <- do.call(c, con_ls)
  
  net_all <- as_sfnetwork(all_sf, directed = T)
  net_all_id <- net_all %>% 
    activate("nodes") %>% 
    mutate(id = row_number())
  
  # remove isolated nodes from pois 
  pois_con <- pois %>% 
    filter(!(meetobjectcode %in% unlist(iso_ls)))
  pois_iso <- pois %>% 
    filter(meetobjectcode %in% unlist(iso_ls))
  st_geometry(pois_iso) <- "x"
  
  # Find indices of nearest nodes.
  nearest_nodes <- st_nearest_feature(pois_con, net_all_id)
  # re-order pois
  pois_ordered <- pois_con %>% 
    mutate(id = nearest_nodes)%>% 
    as.data.frame() %>%
    select(-geometry)
  # join information
  net_con <- net_all_id %>% 
    activate("nodes") %>% 
    left_join(pois_ordered)
  
  # extract nodes and edges
  nodes <- st_as_sf(activate(net_con, "nodes"))
  edges <- st_as_sf(activate(net_con, "edges"))
  
  # add new node(s)
  nodes_new <- bind_rows(nodes, pois_iso)
  
  # rebuild the sfnetwork
  net_con_new <- sfnetwork(nodes_new, edges, directed = is_directed(net_con))
  
  fp <- paste0("temp/net_split_", split_i, "_nodes.gpkg")
  net_con_new %>%
    activate("nodes") %>%
    st_as_sf() %>%
    st_write(fp, delete_layer = T)
  message("Saved: ", fp)
  
  fp <- paste0("temp/net_split_", split_i, "_edges.gpkg")
  net_con_new %>%
    activate("edges") %>%
    st_as_sf() %>%
    st_write(fp, delete_layer = T)
  message("Saved: ", fp)
}

# Create function for parallel computing
run_this <- function(i){
  col_name <- paste0("split_", 12)
  
  # filter for all observations in split
  moc <- mw_df_z %>%
    filter(!!sym(col_name) == i)
  
  # get unique meetobjectcodes in split
  moi <- moc %>%
    distinct(meetobjectcode) %>% 
    pull(meetobjectcode)
  
  # get only meetobjectcodes within the water board and within the time split
  pois <- meetobjecten_oi %>%
    filter(meetobjectcode %in% moi)
  
  createNetwork(net, pois, i)
}

# Prepare parallel computing
ncores <- detectCores() - 4
cl <- makeCluster(ncores)

# Export necessary variables and packages to the cluster
clusterExport(cl, varlist = c("net", "mw_df_z", "meetobjecten_oi", "createNetwork"), envir = environment())
clusterEvalQ(cl, {
  library(dplyr)
  library(sfnetworks)
  library(stringr)
  library(tidygraph)
  library(sf)
  library(igraph)
  })

# Run with progress bar
results <- pbapply::pblapply(1:splits_n, run_this, cl = cl)

# Stop the cluster
stopCluster(cl)

#### Clean
rm(cl, flowdirection, meetobjecten_oi, mw_df_z, net, results, ncores, createNetwork,
   run_this)
