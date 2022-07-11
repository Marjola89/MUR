###################################################################################
################# Mass Univariate Analysis for UKBB data ##########################
###################################################################################
###################################################################################
#' Mass Univariate Analysis.
#'
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V.
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix  Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
#' @param mesh_Coordinates is the matrix from the 3D model mesh of the template.
#' @param phenotype is the imaging phenotype e.g. 1 for `S2S`, 2 for `Curvature`.
#' @param organ is the organ segmentation e.g. liver, spleen, kidney_left etc.
#' @param scale_range is the range of the independent variables to be scaled e.g. c(1:10).
#' @param nPermutations is the number of permutation testing e.g. 1000.
#' @param extract_range is the range of the independent variables to extract the beta coefficients and p-values after MUR analysis e.g. c(2, 5:8).
#' @param output_dir is the path of the output directory were the output with the coefficients and p-values for each covariate are saved as a ```.txt``` file.
#' @export
#' @examples mur_analysis <- mass_univariate_analysis(inputClinical, Y, A, NNmatrix, mesh_Coordinates, phenotype, organ, scale_range, nPermutations, extract_range, output_dir)

mass_univariate_analysis <- function(inputClinical, Y, A, NNmatrix, mesh_Coordinates, phenotype, organ, scale_range, nPermutations, extract_range, output_dir){

    library(data.table)
    library(multtest)
    library(mutools3D)

    # CLINICAL DATA MATRIX
    # NCOL = N COVARIATES UNDER STUDY
    # NROW = N PATIENTS
    head(inputClinical)

    X <- data.matrix(inputClinical[, -c(1)])

    X <- cbind(1, X)
    print(dim(X))
    print(head(X))

    # IMAGING DATA MATRIX
    # NCOL = N POINTS ON THE ATLAS
    # NROW = N PATIENTS
    # DATA PRE-PROCESSING
    print(dim(Y))

    # NUMBER OF CORES TO USE
    # All core detected here

    # set the column names of the organ mesh coordinates
    colnames(mesh_Coordinates) <- c("x", "y", "z")

    extractNames <- colnames(X)

    # Run the 3D Mass univariate regression analysis
    print(extract_range)

    for (iEx in extract_range){

        start.time <- Sys.time()
        extract <- iEx
        print(extract)
        result <- murq(X, Y, extract)

        # MULTIPLE TESTING CORRECTION
        corrected <- mt.rawp2adjp(result[, 3], proc = c("BH"), na.rm = FALSE)
        pvalueADJ5tsbh <- array(dim = length(result[, 3]))
        BHpvalues <- corrected$adjp[order(corrected$index),][, 2]


        meshCoordinates <- cbind(mesh_Coordinates,99999)
        # PRINT OUTPUT
        meshCoordinates[, 4] <- result[,  1]
        write.table(meshCoordinates, paste(output_dir, extractNames[extract], "_beta_", phenotype, "_", organ, ".txt",sep = ""), col.names = FALSE, row.names = FALSE)
        meshCoordinates[, 4] <- result[, 3]
        write.table(meshCoordinates, paste(output_dir, extractNames[extract], "_pvalues_", phenotype, "_", organ, ".txt",sep = ""), col.names = FALSE, row.names = FALSE)
        meshCoordinates[, 4] <- BHpvalues
        write.table(meshCoordinates, paste(output_dir, extractNames[extract], "_BHpvalues_", phenotype, "_", organ, ".txt",sep = ""), col.names = FALSE, row.names = FALSE)

        signif <- permFL_fast(X, Y, extract, as.matrix(A), as.matrix(NNmatrix), nPermutations, E = 0.5, H = 2)
        sign <- signif$pvalues # get p-values

        # MULTIPLE TESTING CORRECTION
        pfdr5TSBH <- mt.rawp2adjp(sign[, 1], proc=c("BH"), na.rm = FALSE)
        pvalueADJ5tsbh <- array(dim = length(sign[, 1]))
        BHpvaluesTFCE <- pfdr5TSBH$adjp[order(pfdr5TSBH$index),][, 2]

        # PRINT OUTPUT
        meshCoordinates[, 4] <- sign[, 1]
        write.table(meshCoordinates, paste(output_dir, extractNames[extract], "_pvaluesTFCE_", phenotype, "_", organ, ".txt", sep=""), col.names = FALSE, row.names = FALSE)
        meshCoordinates[, 4] <- BHpvaluesTFCE
        write.table(meshCoordinates, paste(output_dir, extractNames[extract], "_BHpvaluesTFCE_", phenotype, "_", organ, ".txt", sep=""), col.names = FALSE, row.names = FALSE)

        end.time <- Sys.time()
        time.taken<-end.time - start.time
        print(time.taken)

    }
}
# END

