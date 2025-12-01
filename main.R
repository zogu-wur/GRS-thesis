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
# Main Script
#
# The main script is used to prepare the environment and run all other scripts.
#   Between the scripts 05 and 06, the user needs to manually execute the 
#     Positive Matrix Factorization (PMF) using EPA PMF 5.0 (more details can be
#     found in the report).
#   Uncomment heading 2 after having run the PMF.
#   Fill in DB-details to connect to DB at 1.2.

############################################################################## #
#### 1 Preparation #############################################################
############################################################################## #

#### 1.1 Load libraries ########################################################

library(DBI)
library(RPostgres)
library(tidyverse)
library(parallel)
library(sfnetworks)
library(sf)
library(readxl)
library(nngeo)
library(terra)
library(rpostgis)

library(ggthemes)
library(basemaps)
library(ggalluvial)
library(patchwork)
library(xlsx)
library(shadowtext)
library(tidygraph)
library(igraph)
library(stringr)
library(geojsonsf)
library(units)
library(lwgeom)

#### 1.2 Connect to DB #########################################################
con <- dbConnect(
  RPostgres::Postgres(),
  dbname = ,
  host = ,
  port = ,
  user = ,
  password = 
)

#### 1.3 Define constants ######################################################

splits <- c(2,4,6,8,10,12)

# From data exploration
splits_n <- 12
metal_oi <- c("As", "B", "Ba", "Co", "Se", "U", "Zn")

# Create folder if needed
ifelse(!dir.exists("temp"),
       dir.create("temp"),
       "Directory Exists")

ifelse(!dir.exists("PMF"),
       dir.create("PMF"),
       "Directory Exists")

ifelse(!dir.exists("Figures"),
       dir.create("Figures"),
       "Directory Exists")

# Define file paths
pmf_path <- "PMF/"
fp_save <- "Figures/"

############################################################################## #
#### 1 Run scripts pre PMF #####################################################
############################################################################## #

source("R/S01_.R") # Pre-processing of data
source("R/S02_.R") # Create time windows
source("R/S03_.R") # Network creation
source("R/S04_.R") # Spatial autocorrelation

############################################################################## #
#### 2 Run scripts post PMF ####################################################
############################################################################## #

# source("R/S05_.R") # Locate fragments and their content
# source("R/S06_.R") # Visualisations

# #### Clean
# disconnect from DB
# dbDisconnect(con)