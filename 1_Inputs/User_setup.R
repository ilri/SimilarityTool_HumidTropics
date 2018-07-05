#######################################################################
## CRP sites similarity mapping: parameters
#######################################################################
## Contact details:   marius.gilbert@gmail.com 
## Latest updated:    21/9/2016 updated by Catherine Pfeifer c.pfeifer@cgiar.org
## Version:           v1.0
## Description:       This script takes a series of covariates in the form
##                    of a list of asc files, then extract the pixel values
##                    contained in a shape file and map the similarity to the 
##                    pixel
#######################################################################
library("raster")
library("maptools")
library("rgdal")
library("dismo")
library("biomod2")

memory.limit(size = 400000)
# Sets the general path
#setwd("C:/Users/trobinson/Dropbox/IIASA" )
setwd("D:/Dropbox/similarity analysis tool box - HT" )


# 1.0 General name of the run
###############################################
P_myOutName = "HT base run"

# 1.1 Sets the bounding box of the study region
###############################################
P_XMin = 20
P_XMax = 52
P_YMin = -10
P_YMax = 17

# 1.2 Sets the list of predictors to be used 
# for similarity mapping
###############################################
P_myRasterList = "1_Inputs/Variables.csv"

# Folder where workspace images should be stored
###############################################
P_myExFolder = "2_Workspace/"

# Should the extraction be redone if the file is 
# already present ?
###############################################
P_RedoExtraction = 1

# Folder where outputs should be stored
###############################################
P_myOutData = "3_Outputs/"
  
# 1.3 What is the name of the Shape file containing the area
# from which similarity should be mapped ?
###############################################
P_myDepShape = "1_Inputs/StudyArea/test_eth.shp"
P_NPoints = 1000

# 1.4 Path of the general shape file used 
# for reporting maps
###############################################
P_myGenShape = "1_Inputs/MapBackground/country.shp"

# 1.5 runs the extraction in the sub-area
###############################################
source("0_Code/1_CreateWorkSpace.r")
source("0_Code/2_Similarity_Mahalanobis_v3.r")


