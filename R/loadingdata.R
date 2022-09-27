#' Load and prepare data for analysis.
#'
#' @param fileClinicalDataPath is the path of the file containing the clinical data of the subjects for all the models.
#' @param folderImagingDataPath is the path of the folders storing the image data e.g. "/home/yourname/folder-path/".
#' @param phenoType is the imaging phenotype e.g. "1" for S2S, "2" for Curvature.
#' @param nPoints is the number of points (vertices) in the mesh, this will change for each organ e.g. for kidney left: 4380.
#' @param organ is the organ segmentation e.g. liver, spleen, kidney_left etc.
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar read.table
#' @export
#' @examples
#' \dontrun{
#' model <- loadingdata(fileClinicalDataPath, folderImagingDataPath, phenoType, nPoints, organ)
#' }

loadingdata <- function(fileClinicalDataPath, folderImagingDataPath, phenoType, nPoints, organ){

  fileNames <- c(paste0("/", organ, "_registration_output/txt/", organ, "_mask_signeddistances.txt"),
                 paste0("/", organ, "_registration_output/txt/", organ, "_mask_curvature.txt")
  )
  
  
  foldersNames <- list.dirs(folderImagingDataPath,
                            full.names = FALSE,
                            recursive = F)
  clinicalData <- as.data.frame(read.csv(fileClinicalDataPath))
  nPatients <- nrow(clinicalData)
  
  iDel <- 1
  noK <- c()
  
  for(iP in 1:nPatients){
    filePath <- paste(folderImagingDataPath,
                      foldersNames[grep(clinicalData[iP, 1], foldersNames)[1]],
                      fileNames[phenoType],  sep = "")
    setTxtProgressBar(txtProgressBar(style=3),iP/nPatients)
    if(!file.exists(filePath[1])){
      noK[iDel] <- iP
      iDel <- iDel + 1
    }
  }
  
  
  if(length(noK)>0) clinicalData <- clinicalData[-noK,]
  nPatients <- nrow(clinicalData)
  
  Y <- matrix(0, ncol = nPoints, nrow = nPatients)
  Xc <- matrix(0, ncol = nPoints, nrow = nPatients)
  Yc <- matrix(0, ncol = nPoints, nrow = nPatients)
  Zc <- matrix(0, ncol = nPoints, nrow = nPatients)
  print(dim(Y))
  
  # READ ATLAS DATA
  # Y response matrix
  for(iP in 1:nPatients){
    filePath <- paste(folderImagingDataPath,
                      foldersNames[grep(clinicalData[iP, 1], foldersNames)[1]],
                      fileNames[phenoType],  sep = "")
    setTxtProgressBar(txtProgressBar(style=3),iP/nPatients)
    # Imaging Data Path for the iP patient
    if (file.exists(filePath)){  #if the folder exist
      rDataFrame <- lapply(file.path(filePath), function(x) {
        tryCatch(read.table(x), error=function(e) NULL)
      })
      if(!is.null(rDataFrame[[1]]$V4)){
        Xc[iP,] <- as.vector(rDataFrame[[1]]$V1) # take the 1st row of the mesh - X coordinate
        Yc[iP,] <- as.vector(rDataFrame[[1]]$V2) # take the 2nd row of the mesh - Y coordinate
        Zc[iP,] <- as.vector(rDataFrame[[1]]$V3) # take the 3rd row of the mesh - Z coordinate
        Y[iP,] <- as.vector(rDataFrame[[1]]$V4)  # take the 4th row of the mesh - 3D mesh-derived phenotype
      }
    }
  }

  model <- list("X" = clinicalData, "Y" = Y, "X_coordinate" = Xc, "Y_coordinate" = Yc, "Z_coordinate" = Zc)
  
  return(model)
}

