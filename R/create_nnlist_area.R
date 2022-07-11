#' Compute NN list and voronoi area.
#'
#' @param data_dir is the path of the data directory were the matrices X and Y for the model are saved.
#' @param organ is the organ segmentation e.g. liver, spleen, kidney_left etc.
#' @param dpl is a list describing the structure of the dummy points to be added to the data being triangulated.
#' e.g. list(ndx = 5, ndy = 5). If the argument is 'NULL' then no dummy point is added to the data.
#' ndx: The x-dimension of a rectangular grid; if either ndx or ndy is null, no grid is constructed.
#' ndy: The y-dimension of the aforementioned rectangular grid.
#' @export
#' @examples NN_Area <- create_nnlist_area(data_dir, organ, dpl)

create_nnlist_area <- function(data_dir, organ, dpl){

    library(data.table)
    library(Rvcg)
    library(deldir)
    # load the organ mesh coordinates
    mesh_Coordinates <- read.csv(paste(data_dir,"data/", organ, "_meshCoordinates.csv", sep = ""))
    data_mesh <- as.matrix(mesh_Coordinates)
    colnames(data_mesh) <- NULL
    NN_data_mesh <- apply(data_mesh, 2, as.numeric)

    # surface reconstruction providing a triangular data if class mesh3D
    NN_data_mesh_new <- vcgBallPivoting(NN_data_mesh)

    # compute the NN list
    edges <- vcgGetEdge(NN_data_mesh_new)
    NNmatrix <- edges[, 1:2]
    colnames(NNmatrix) <- c("v", "NN")

    # compute the voronoi area of the organ
    z <- deldir(data_mesh, plot = FALSE, dpl = dpl, digits = 2)
    Area <- z$summary$del.area[1:nrow(mesh_Coordinates)]
    Area[is.na(Area)] <- 0

    NN_Area <- list("NNmatrix" = NNmatrix, "Area" = Area)
    return(NN_Area)
}
# END
