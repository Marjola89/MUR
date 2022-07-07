# Mass Univariate Regression Analysis (MUR)


The **MUR** package implements R functions for mass univariate analysis of three-dimensional phenotypes. 

## Installation

Currently this package is only available on GitHub and requires R (>= 3.3.2). For installation, the easiest way is using `devtools`. If you don’t have `devtools` installed, please start by typing:

```r
install.packages("devtools")
library(devtools)
```

Then, to install the package just use `install_github`function.

```r
install_github("Marjola89/MUR")
```

## Code
### Produce Matrices

First prepare matrices X and Y in .rds files e.g `S2S_kidney_left.rds` for the model that will be used in Mass Univariate Analysis code.

```bash
model <- loadingdata(fileClinicalDataPath, folderImagingDataPath, indVar, nPoints, organ)
```
Input parameters:
* `fileClinicalDataPath`: Set the path of the file containing the clinical data of the subjects for all the models e.g. `"/home/yourname/clinicaldata-path/data/ClinicalData.csv"`.
* `folderImagingDataPath`: Set the path of the folders storing the image data e.g. `"/home/yourname/folder-path/"`.
* `indVar`: Define the independent variables of the model e.g. `c("individual_id","Age","Sex","Race","BMI","WHR")`.
* `nPoints`: Set the number of points (vertices) in the mesh, this will change for each organ e.g. for kidney left: `4380`.
* `organ`: Set organ segmentation e.g. `liver`, `spleen`, `kidney_left` etc.

Output parameters:
* `model`: a list including the design matrix `X`, the imaging matrix `Y`, the matrix with x-coordinates `Xc`, the matrix with y-coordinates `Yc`and the matrix with z-coordinates `Zc`.

Notes:

* It is essential to saved the output as:
```bash
phenoType <- c(1, 2)
# 1 = S2S, 2 = Curvature

phenoTypeNames <- c("S2S","Curvature")

for(iM in 1:length(phenoType)){
  model <- loadingdata(fileClinicalDataPath, folderImagingDataPath, indVar, phenoType[iM], nPoints)

  clinicalData <- model[[1]]
  Y <- model[[2]]
  Xc <- model[[3]]
  Yc <- model[[4]]
  Zc <- model[[5]]

  saveRDS(X, paste(data_dir, phenoTypeNames[phenoType[iM]], "_", organ, "_clinicalData.rds", sep = ""))
  saveRDS(Y, paste(data_dir, phenoTypeNames[phenoType[iM]], "_", organ, ".rds", sep = ""))
  saveRDS(Xc, paste(data_dir, phenoTypeNames[phenoType[iM]], "_", organ, "_Xcoordinate.rds", sep = ""))
  saveRDS(Yc, paste(data_dir, phenoTypeNames[phenoType[iM]], "_", organ, "_Ycoordinate.rds", sep = ""))
  saveRDS(Zc, paste(data_dir, phenoTypeNames[phenoType[iM]], "_", organ, "_Zcoordinate.rds", sep = ""))

}
```
* `data_dir`: is the path of the data directory  were the matrices X and Y for the model will be saved.

### Code for computing the NN list and the voronoi area
This is needed for the Mass Univariate Analysis. Usage:
```bash 
NN_Area <- create_nnlist_area(data_dir, organ, dpl) 
```
Input parameters:
* `data_dir`: Set the path of the data directory were the input matrices are saved.
* `organ`: Set organ segmentation e.g. `liver`, `spleen`, `kidney_left` etc.
* `dpl` : A list describing the structure of the dummy points to be added to the data being triangulated. e.g. ```list(ndx = 1.5, ndy = 1.5)```. 
  If the argument is ```NULL``` then no dummy point is added to the data. 
  ndx: The x-dimension of a rectangular grid; if either ndx or ndy is null, no grid is constructed. 
  ndy: The y-dimension of the aforementioned rectangular grid.

Notes:

* It is essential to save the output as a `.rds` file in the `data_dir` as `<organ>NNmatrix.rds` and  `<organ>area.rds` e.g. `kindey_leftNNmatrix.rds` and `kindey_leftarea.rds`.

* Alternative way to create NN matrices and area is using ```computeNNmatrix``` and ```computeVertAreas``` from the R package [`mutools3D`](https://github.com/UK-Digital-Heart-Project/mutools3D).

### Mass univariate regression analysis Code

Computes the associations between the 3D mesh-derived phenotype (e.g. signed distances) and other clinical variables using a linear regression framework.

Run the Mass Univariate Analysis with TFCE using the system's cores.

```bash
mur_analysis <- mass_univariate_analysis(data_dir, phenotype, organ, scale_range, nPermutations, data_range, extract_range, output_dir)
```
This code is dependent on three functions (`~/function/murq.R`, `~/function/permFL_fast.R` and `~/function/TFCE.R`) which are called within the code.

Input parameters:
* `data_dir`: Set the path of the data directory were the matrices X and Y for the model are saved.
* `phenotype`: Define the imaging phenotype e.g. `1` for S2S, `2` for Curvature.
* `organ`: Set organ segmentation e.g. `liver`, `spleen`, `kidney_left` etc.
* `scale_range`: Range of the independent variables to be scaled e.g. `c(1:10)`.
* `nPermutations`: Number of permutation testing e.g. `1000`.
* `data_range`: Range of the independent variables to be included in the model e.g. `c(1:8)`.
* `extract_range`: Range of the independent variables to extract the beta coefficients and p-values after MUR analysis e.g. `c(2, 5:8)`.
* `output_dir`: Set the path of the output directory were the output with the coefficients and p-values for each covariate are saved as a ```.txt``` file.

Notes:

* It is essential to have saved the mesh coordinates of the organ template as a `.csv` in the `data_dir` as `<organ>_meshCoordinates.csv` e.g. `kindey_left_meshCoordinates.csv`.


Sub-functions for computing the mass univariate regression analysis

* `murq.R` for the mass univariate regression analysis across the mesh vertices.
  ```bash
  result <- murq(X, Y, extract)```
* `permFL_fast.R` for the permutation the data and applying TFCE.
   ```bash 
   TFCEresults <- permFL_fast(X, Y, extract, A, NNmatrix, nPermutations, E = 0.5, H = 2)
   ```
   Input parameters:
   * `X` is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
   * `Y` is the imaging matrix. Number of rows = N. Number of columns = V.
   * `extract` is an array expressing which covariates in X you want to extract.
   * `A` a V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
   * `NNmatrix` Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
   * `nPermutations` number of permutations in the permutation test, default is 1000.
   * `E` is the TFCE parameter, by default fixed to 0.5.
   * `H` is the TFCE parameter, by default fixed to 2.

Notes:

* `TFCE.R` function is used from the R package [`mutools3D`](https://github.com/UK-Digital-Heart-Project/mutools3D).
* `murq.R` and `permFL_fast.R` functions are similar to the functions `mur.R` and `permFL.R` of the `mutools3D` package however, they are enhanced to allow process of big datasets .
