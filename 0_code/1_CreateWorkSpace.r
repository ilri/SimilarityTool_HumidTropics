##########################################################
# Extraction script
##########################################################
# Description: Crop all predictors with a fixed bounding
# box and stores the outputs in a specified folder
##########################################################
# Add respective weighting to the output
##########################################################
  sdate = date()
# Sets the extent of the cropping
  G_Extent = extent(P_XMin, P_XMax, P_YMin, P_YMax)
# Subset to fields to be extracted only
  G_myPredTable = read.csv(P_myRasterList)
  G_myPredTable$PATH = as.character(G_myPredTable$PATH)
  G_myPredTable$VARNAME = as.character(G_myPredTable$VARNAME)
  G_myPredTable$EXPATH = paste(P_myExFolder,"/",G_myPredTable$VARNAME,".tif", sep = "")
  myEPredTable = subset(G_myPredTable, EXTRACTION == 1)
  myMaskTable = subset(G_myPredTable, Mask == 1)
# Creates the mask
  mask <- raster(myMaskTable$PATH[1])>0
  maskv = values(mask)
  mask = setValues(mask,ifelse(maskv == 1,1,NA))
# If the output folder does not exists, creates it
  dir.create(P_myExFolder, showWarnings = FALSE)
# Starts the extractions
  print("#############################################################")    
  print(paste("## Extractions started on ", sdate))
  print("#############################################################")    
# Crop all images in the stack
  if (P_RedoExtraction == 1){
  print("## Info: Re-extracting all imagery")    
    for (i in 1:nrow(myEPredTable)){      
        myRaster = raster(myEPredTable$PATH[i])       
        myRaster@crs = CRS("+proj=longlat +datum=WGS84")
        names(myRaster) = myEPredTable$VARNAME[i]
        myRaster = mask(myRaster, mask)
            
        if (extent(myRaster)@xmin!= G_Extent@xmin | extent(myRaster)@xmax!= G_Extent@xmax | extent(myRaster)@ymin!= G_Extent@ymin | extent(myRaster)@ymax!= G_Extent@ymax) {
        myCr = crop(myRaster, G_Extent)
        writeRaster(myCr,filename = myEPredTable$EXPATH[i], format = "GTiff", overwrite = T )
        #writeGDAL(as(myCr, "SpatialGridDataFrame"), fname = myEPredTable$EXPATH[i], drivername = "RST") 
        } else {
        writeRaster(myRaster,filename = myEPredTable$EXPATH[i], format = "GTiff", overwrite = T )
        #writeGDAL(as(myRaster, "SpatialGridDataFrame"), fname = myEPredTable$EXPATH[i], drivername = "RST") 
        }
        
        myStr = paste("## ",myEPredTable$VARNAME[i], " extracted to: ",myEPredTable$EXPATH[i])     
        print(myStr)
   }
  } else {
  print("## Info: extracting only missing files")    
    for (i in 1:nrow(myEPredTable)){
        Test = file.access(myEPredTable$EXPATH[i], mode = 4)[[1]]
        if (Test == -1) {
        myRaster = raster(myEPredTable$PATH[i])   
        myRaster@crs = CRS("+proj=longlat +datum=WGS84")
        names(myRaster) = myEPredTable$VARNAME[i]
        myRaster = mask(myRaster, mask)
  
        if (extent(myRaster)@xmin!= G_Extent@xmin | extent(myRaster)@xmax!= G_Extent@xmax | extent(myRaster)@ymin!= G_Extent@ymin | extent(myRaster)@ymax!= G_Extent@ymax) {
          myCr = crop(myRaster, G_Extent)
        writeRaster(myCr,filename = myEPredTable$EXPATH[i], format = "GTiff", overwrite = T )
#          writeGDAL(as(myCr, "SpatialGridDataFrame"), fname = myEPredTable$EXPATH[i], drivername = "RST") 
        } else {
        writeRaster(myRaster,filename = myEPredTable$EXPATH[i], format = "GTiff", overwrite = T )
#          writeGDAL(as(myRaster, "SpatialGridDataFrame"), fname = myEPredTable$EXPATH[i], drivername = "RST") 
        }
        
        myStr = paste("## ",myEPredTable$VARNAME[i], " extracted to: ",myEPredTable$EXPATH[i])     
        print(myStr)
      }
   }
  }
 
  edate = date()
  print("#############################################################")
  print(paste("Extractions Finished on ", edate))
  print("#############################################################")
  print("## Generating extraction pdf report                       ")
  
  myGlobalShape = readShapePoly(P_myGenShape)
  myLocalShape = readShapePoly(P_myDepShape)
  
  # If the output folder does not exists, creates it
  dir.create(P_myOutData, showWarnings = FALSE)
  myPDFPath = paste(P_myOutData,"/0_StudyArea.pdf", sep = "")
  pdf(myPDFPath)
  plot(myGlobalShape, main = "Extraction area", xlim = c(P_XMin,P_XMax), ylim = c(P_YMin,P_YMax), axes = T)
  plot(myLocalShape, add = T, col = "grey")
  rect(P_XMin,P_YMin,P_XMax,P_YMax, density = 5, col = "red", angle = 45,border = "red", lwd = 1.5)
  for (i in 1:nrow(myEPredTable)){
  myR = raster(myEPredTable$EXPATH[i])
  plot(myR, main = myEPredTable$VARNAME[i])
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  }
  dev.off()

  plot(myGlobalShape, main = "Extraction area", xlim = c(P_XMin,P_XMax), ylim = c(P_YMin,P_YMax), axes = T)
  plot(myLocalShape, add = T, col = "grey")

  # Plot Bounding box
  rect(P_XMin,P_YMin,P_XMax,P_YMax,
     density = 10, col = "red", angle = 45,border = "red", lwd = 1.5)

print("## Report is finished                       ")
print("#############################################################")


