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
# Script 01: Pre-processing
#
# This script does some pre-processing and creates the research area.

############################################################################## #
#### 1 Create research area ####################################################
############################################################################## #

# Load data
provincies_sf <- st_read(con, "s0_provinces_raw")
waterschappen_sf <- st_read(con, "s0_waterboards_raw")
top10nlwater_sf <- st_read(con, "s0_top10nl_water_raw")

# Select provinces
groningen <- provincies_sf %>% filter(naam == "Groningen")
drenthe <- provincies_sf %>% filter(naam == "Drenthe")

# Create intersection between selected provinces and water boards 
nzv_sf <- waterschappen_sf %>%
  filter(naam == "Noorderzijlvest") %>%
  st_intersection(groningen %>% rbind(drenthe))

hea_sf <- waterschappen_sf %>%
  filter(naam == "Waterschap Hunze en Aa's") %>%
  st_intersection(groningen %>% rbind(drenthe))

ra_land <- rbind(nzv_sf, hea_sf) %>%
  st_union() %>%
  st_remove_holes()

# Remove sea water
water_sf <- top10nlwater_sf %>%
  filter(!typewater %in% c("waterloop","meer, plas")) %>%
  st_union() %>%
  st_remove_holes()

ra_sf <- st_difference(ra_land, water_sf)

# save
st_write(ra_sf, dsn = con, layer = "s1_research_area", delete_layer = T)

############################################################################## #
#### 2 Clean measurement values ################################################
############################################################################## #

# Load selected data
mw_df <- dbGetQuery(con,"
                    SELECT meetjaar, waterbeheerdercode, waterbeheerdernaam, meetobjectcode, begindatum, parametercode, 
                      parameteromschrijving, eenheidcode, eenheidomschrijving, limietsymbool, numeriekewaarde
                    FROM s0_mv_raw
                    WHERE grootheidcode = 'CONCTTE'
                      AND parametercode IN ('As', 'B', 'Ba', 'Co', 'Se', 'U', 'Zn')
                      AND NOT meetobjectcode IN ('NL33_6640', 'NL34_1133');
                    ")

# save
dbWriteTable(con, "s1_mv_oi", mw_df, overwrite = T, row.names = FALSE)

############################################################################## #
#### 3 Clean measurement objects ###############################################
############################################################################## #

# Load selected data
mo_sf <- st_read(con, query = "
                 SELECT mo.waterbeheerdercode, mo.waterbeheerdernaam, mo.meetobjectcode, mo.geometry
                 FROM s0_mo_raw mo
                 JOIN s1_research_area r ON ST_Within(mo.geometry, r.geom)
                  OR mo.meetobjectcode = 'NL33_1638'
                 WHERE mo.meetobjectcode IN (SELECT DISTINCT meetobjectcode FROM s1_mv_oi)
                 ")

# save
st_write(mo_sf, dsn = con, layer = "s1_mo_oi", delete_layer = T)

############################################################################## #
#### 4 Clean water networks ####################################################
############################################################################## #

#### 4.1 Water board HEA #######################################################

# Load data
flowdirection_hea <- st_read(con, query = "
  SELECT hea_naam_w, objectiden, ovksoort, geometry
  FROM s0_hea_flow_raw
")

# Reverse edges
hea_still <- c("Reitdiep", "Eemskanaal", "Oude Eemskanaal", "Zeehavenkanaal", "Balkenhaven", "Damsterhaven", "Schildmeer",
               "De Groeve", "Verbindingskanaal Lalleweer", "Termunterzijldiep", "Afwateringskanaal van Duurswold",
               'Vispassage "de Bult"', "Opdiep", "Slochterdiep (ook Rengersdiep, gegraven 1694)")

reverse_edges <- flowdirection_hea %>% 
  filter(!hea_naam_w %in% hea_still,
         !(hea_naam_w == "Oosterhornhaven" & objectiden == "W0646_OAF-E-00070"),
         !(hea_naam_w == "Winschoterdiep" & objectiden %in% c("W0646_OAF-E-00150", "W0646_OAF-E-00771", "W0646_OAF-E-00267",
                                                              "W0646_OAF-E-00159", "W0646_OAF-E-00535", "W0646_OAF-E-00534",
                                                              "W0646_OAF-E-00268", "W0646_OAF-E-00772", "W0646_OAF-E-00525",
                                                              "W0646_OAF-E-00526", "W0646_OAF-E-00527", "W0646_OAF-E-00912",
                                                              "W0646_OAF-E-00542", "W0646_OAF-E-00174", "W0646_OAF-E-00575",
                                                              "W0646_OAF-E-00574", "W0646_OAF-E-00576")),
         !(hea_naam_w == "Oosterhaven" & objectiden == "W0646_OAF-E-00463"),
         !(hea_naam_w == "Oude Winschoterdiep" & objectiden %in% c("W0646_OAF-E-00151", "W0646_OAF-E-00768", "W0646_OAF-E-00770")),
         !(hea_naam_w == "Westerwoldsche Aa" & objectiden %in% c("W0646_OAF-E-00221", "W0646_OAF-E-00217", "W0646_OAF-E-00218",
                                                                 "W0646_OAF-E-00272", "W0646_OAF-E-00216", "W0646_OAF-E-00538",
                                                                 "W0646_OAF-E-00682", "W0646_OAF-E-00209", "W0646_OAF-E-00271",
                                                                 "W0646_OAF-E-00090", "W0646_OAF-E-00432", "W0646_OAF-E-00520",
                                                                 "W0646_OAF-E-00198", "W0646_OAF-E-00544", "W0646_OAF-E-00578",
                                                                 "W0646_OAF-E-00577", "W0646_OAF-E-00197", "W0646_OAF-E-00511")),
         !(hea_naam_w == "Pekel Aa" & objectiden %in% c("W0646_OAF-E-00175", "W0646_OAF-E-00176", "W0646_OAF-E-00579")),
         !objectiden %in% c("W0646_OAF-E-00629", "W0646_OAF-E-00546", "W0646_OAF-E-00631", "W0646_OAF-E-00182",
                            "W0646_OAF-E-00180", "W0646_OAF-E-00405", "W0646_OAF-E-00181"),
         ovksoort %in% c("Boezemvak", "af- en aanvoervak")) %>% 
  pull(objectiden)

flowdirection_hea <- flowdirection_hea %>% 
  mutate(birectional = ifelse(objectiden %in% reverse_edges, 1, 0)) %>% 
  select(-c(hea_naam_w, objectiden, ovksoort))

#### 4.2 Water baord NZV #######################################################

# Load
flowdirection_nzv <- st_read(con, query = "
  SELECT code, ws_zomeraf, ws_winte_1, geometry
  FROM s0_nzv_flow_raw
")

# Fix wrong edge
fix_edge <- flowdirection_nzv %>% filter(code == "OAF016164")
fix_edge$geometry <- st_reverse(fix_edge$geometry)
flowdirection_nzv <- flowdirection_nzv %>%
  filter(code != "OAF016164") %>%
  bind_rows(fix_edge)

# Reverse edges
reverse_edges <- flowdirection_nzv %>%
  filter(ws_zomeraf == 2 | ws_winte_1 == 2) %>%
  pull(code)

flowdirection_nzv <- flowdirection_nzv %>%
  mutate(birectional = ifelse(code %in% reverse_edges, 1, 0)) %>%
  select(-c(ws_zomeraf, ws_winte_1, code))

#### 4.3 Combine water boards ##################################################

# Load connections
con_p <- st_read(con, query ="
                 SELECT *
                 FROM s0_connections_raw")
con_lut <- dbGetQuery(con,  "
                            SELECT *
                            FROM s0_connections_lut_raw")

# Make networks
hea_net <- as_sfnetwork(flowdirection_hea, directed = T)
nzv_net <- as_sfnetwork(flowdirection_nzv, directed = T)

# Create function to join points
add_points <- function(net, ws){
  
  # filter points for net
  hea_p <- con_p %>% 
    filter(ws == !!ws)
  
  # find nearest node
  nearest_nodes <- st_nearest_feature(hea_p, net)
  hea_p$id <- nearest_nodes

  # add information
  net_p <- net %>% 
    activate("nodes") %>% 
    mutate(id = row_number()) %>% 
    left_join(st_drop_geometry(hea_p)) %>% 
    slice(hea_p$id) %>% 
    st_as_sf() %>% 
    select(-id, ws) %>% 
    arrange(pair)
  
  return(net_p)
}

hea_p <- add_points(hea_net, "hea")
nzv_p <- add_points(nzv_net, "nzv")

# Combine
points_sf <- bind_rows(hea_p, nzv_p)

lines_sf <- points_sf %>%
  arrange(pair, desc(start)) %>% 
  group_by(pair) %>%
  summarise(do_union = FALSE, 
            geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING") %>% 
  left_join(con_lut)

reverse_edges <- lines_sf %>%
  filter(ws_zomeraf == 2 | ws_winte_1 == 2) %>%
  pull(pair)

lines_sf2 <- lines_sf %>%
  mutate(birectional = ifelse(pair %in% reverse_edges, 1, 0)) %>%
  select(-c(ws_zomeraf, ws_winte_1, pair))

lines_f <- rbind(lines_sf2, 
                 flowdirection_hea, 
                 flowdirection_nzv)

# Save
st_write(lines_f, dsn = con, layer = "s1_flow")

#### Clean 
rm(con_lut, con_p, drenthe, fix_edge, flowdirection_hea, flowdirection_nzv,
   groningen, hea, hea_net, hea_p, hea_sf, lines_f, lines_sf, lines_sf2, mo_sf,
   mw_df, nzv_net, nzv_p, nzv_sf, points_sf, provincies_sf, ra_land, ra_sf,
   top10nlwater_sf, water_sf, waterschappen_sf, hea_still, reverse_edges, 
   add_points)
