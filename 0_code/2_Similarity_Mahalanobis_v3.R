#######################################################################
## CRP sites similarity mapping
#######################################################################
## Contact details:   marius.gilbert@gmail.com 
## Latest updated:    02/07/2014
## Version:           v1.1
## Description:       This script takes a series of covariates in the form
##                    of a list of asc files, then extract the pixel values
##                    contained in a shape file and map the similarity to the 
##                    pixel
#######################################################################

ReclassQuant = function(myTempR){
  mySeq = seq(0,1,0.001)
  myVals = values(myTempR)
  myQuantile = quantile(myVals,mySeq, na.rm = T, names = F, type = 8) 
  myQuantile = myQuantile + seq_along(myQuantile) * .Machine$double.eps*1000
  myVL = length(mySeq)
  myReMatrix = cbind(myQuantile[-myVL],myQuantile[-1],mySeq[-myVL])
  myReMatrix = as.matrix(myReMatrix)
  myReMatrix[1,1] = min(myVals, na.rm = T)
  myReMatrix[nrow(myReMatrix),2] = max(myVals, na.rm = T)
  myReclassGrid = reclassify(myTempR,rcl = myReMatrix, right = T, include.lowest = T)
  plot(myReclassGrid)
  return(myReclassGrid)
}
 
mess<-function(X,V,full=TRUE){
messi<-function(p,v){
    niinf<-length(which(v<=p))
    f<-100 * niinf / length(v)
    if(f==0) simi<-100*(p-min(v))/(max(v)-min(v))
    if(0<f & f<=50) simi<-2*f
    if(50<=f & f<100) simi<-2*(100-f)
    if(f==100) simi<-100*(max(v)-p)/(max(v)-min(v))
    return(simi)
    }
E<-extract(x=X,y=1:ncell(X))
r_mess<-X
for (i in 1:(dim(E)[2])) {
    e<-data.frame(E[,i]) ; v<-V[,i]
    r_mess[[i]][]<-apply(X=e, MARGIN=1, FUN=messi, v=v)
    }
rmess<-r_mess[[1]]
E<-extract(x=r_mess,y=1:ncell(r_mess[[1]]))
rmess[]<-apply(X=E, MARGIN=1, FUN=min)
if(full==TRUE) {
    out <- addLayer(r_mess,rmess)
    layerNames(out)<-c(layerNames(X),"mess")
}
if(full==FALSE) out <- rmess
return(out)
}




  print("#############################################################")    
  print(paste("## Similarity mapping ",date()))
  print("#############################################################")    

# LOads the list of files to be considered in the analysis
G_myPredTable = read.csv(P_myRasterList)
G_myPredTable$PATH = as.character(G_myPredTable$PATH)
G_myPredTable$VARNAME = as.character(G_myPredTable$VARNAME)
G_myPredTable$EXPATH = paste(P_myExFolder,"/",G_myPredTable$VARNAME,".tif", sep = "")
myUSETable = subset(G_myPredTable, USE == 1)
myMaskTable = subset(G_myPredTable, Mask == 1)

# Extract values from the Shape file
covs = brick(as.list(myUSETable$EXPATH))

# sets the continental mask
mask <- raster(myMaskTable$EXPATH[1])>0
maskv = values(mask)
mask = setValues(mask,ifelse(maskv == 1,1,NA))

# Masks the brick
covs = mask(covs, mask)

# Apply the weighting
  for (i in 1:nrow(myUSETable)){
  covs[[i]] = scale(covs[[i]])*myUSETable$WEIGHT[i] 
  } 

# Create a dataframe with all predictor values
covs.vals = getValues(covs)


##############################################
## PCA and components visualisation
##############################################
print("## Info: Running PCA")    
 # Runs a PCA on the file
pca <- princomp(na.omit(getValues(covs)), scaled = F)
pca.vals <- predict(pca, getValues(covs))

# rescale these between 0 and 1
pca.vals <- apply(pca.vals, 2, function(x) {
  x <- x - min(x, na.rm = TRUE)
  x <- x / max(x, na.rm = TRUE)
  return (x)
  })

# stick them back in the rasterbrick
pca_brick <- setValues(covs, pca.vals)

# plot the first three axes in RGB
#plotRGB(pca_brick * 255, maxpixels = 1e+9)
#summary(pca)

##############################################
## Extract predictor variables in the sites
##############################################
  # loads the area Shape  
  myLocalShape = readShapePoly(P_myDepShape)
  myPtDF = spsample(myLocalShape, n = P_NPoints, type = "random")

  plot(covs[[1]])
  plot(myPtDF, add = T)

  myDFR = extract(covs, myPtDF)
  myDFP = extract(pca_brick, myPtDF)
  myDFR = na.omit(myDFR)
  myDFP = na.omit(myDFP)
  myDFRMn = colMeans(myDFR)
  myDFPMn = colMeans(myDFP)

# Estimates the Euclidian distance to pixel values
###################################################
  myDR = values(mask)*0
  for (i in 1:nrow(myUSETable)){
  myDR = myDR + (myDFRMn[i]-covs.vals[,i])^2 
  } 
  myDR = setValues(mask,myDR^0.5)

# Estimates the Euclidian distance to PCA-transformed values
###################################################
 print("## Info: Estimating Euclidian distance")

  myDP = values(mask)*0
  for (i in 1:nrow(myUSETable)){
  myDP = myDP + (myDFPMn[i]-pca.vals[,i])^2
  }
  myDP = setValues(mask,myDP^0.5)

# Mahalanobis distance with Raw observations
##############################
 print("## Info: Estimating Mahalanobis distance")
   myPtMx <- as.matrix(myDFR) # Create matrix with raw values for point locations 
   myAllMx <- as.matrix(covs.vals)  # Create matrix with raw values for all pixels
   myPtCov <- cov(na.omit(myPtMx)) # Covariance matrix for all pixels with data
   myMahaV <- mahalanobis(myAllMx, colMeans(myPtMx), myPtCov)
   myMahaP = setValues(mask, myMahaV)

# Bioclim distance
##############################
print("## Info: Bioclim distance")
  bc = bioclim(covs,myPtDF)
  pbc = predict(covs, bc, ext=extent(mask))

# Site similarity
##############################
  mess.out = mess(covs,myDFR, full = F)
  mess.out = mess.out * mask
  mess.pos = (mess.out >=0)*mess.out

# Apply quantile transform to distance metrics
###################################################
#1
  myEucPCAQ = 1-ReclassQuant(myDP)
#2
  myEucRawQ = 1-ReclassQuant(myDR)
#3
  myMahaRawQ = 1-ReclassQuant(myMahaP)
#
  myMESSGen = ReclassQuant(mess.out)
  myMESSPos = ReclassQuant(mess.pos)


# Plot the outputs
###################################################
  print("## Info: Generating pdf report                       ")
  myGlobalShape = readShapePoly(P_myGenShape)
  dir.create(P_myOutData, showWarnings = FALSE)
  myPDFPath = paste(P_myOutData,"/1_Similarity.pdf", sep = "")
  pdf(myPDFPath)
  #par(mfrow=c(2,1))
  plot(pca)
  plot(mask, "PCA in RGB", legend = F)
  plotRGB(pca_brick * 255, maxpixels = 1e+9, main = "PCA in RGB", axes = T, add = T)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  box()
  #par(mfrow=c(1,1))
  
  plot(myEucPCAQ)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  title("Euclidian similarity (quantiles; PCA-transformed covariates)")

  plot(myMahaRawQ)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  title("Mahalanobis similarity (quantiles; raw covariates)")
  
  plot(myMESSGen)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  title("Multivariate Env. Similarity Surfaces (quantiles; raw covariates)")

  plot(myMESSPos)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  title("Multivariate Env. Similarity Surfaces (only positives;)")

  plot(mess.pos)
  plot(myGlobalShape, add = T)
  plot(myLocalShape, add = T, border = "red")
  title("Multivariate Env. Similarity Surfaces (only positives; rescaled)")

  dev.off()
  print("#############################################################")    


myEuclPath = paste(P_myOutData,"/EuclSim.tif", sep = "")
myMahaPath = paste(P_myOutData,"/MahaSim.tif", sep = "")
myMessPath = paste(P_myOutData,"/MessSim.tif", sep = "")
writeRaster(myDP,myEuclPath, format = "GTiff", overwrite = T)
writeRaster(myMahaP,myMahaPath, format = "GTiff", overwrite = T)
writeRaster(mess.out,myMessPath, format = "GTiff", overwrite = T)

data<-stack(myEucPCAQ,myMahaRawQ, myMESSPos)
names(data) <- c("Euclidian", "Mahalanobis", "MESS")
my.at<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
my.colorkey<-list(at=my.at,
                  labels= list(
                    labels= c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'),
                    at=my.at))
mytheme <- rasterTheme(region = brewer.pal(9, "YlGnBu"))

library(RColorBrewer)
library(rasterVis)
tiff( '3_Outputs/highresmap.tif' ,  width=1072*4.5, height=372*4.5, res=300, units='px')
levelplot(data, par.settings=mytheme, margin=FALSE,at=my.at,
          colorkey=my.colorkey)+ layer(sp.polygons(myGlobalShape,col='grey'))+layer(sp.polygons(myLocalShape,col='red'))
dev.off()

