#' @title Helper functions for the processing of .loom and .h5ad files into
#' relational database tables
#'
#' @rdname sce
#'
#' @param loom_filepath character() absolute path to downloaded .loom file
#'
#' @importFrom LoomExperiment import
#'
#' @return SingleCellLoomExperiment object
loom_to_sce <- function(loom_filepath) {
    loom_sce <- LoomExperiment::import(loom_filepath,
                                       type ="SingleCellLoomExperiment")
    loom_sce
}

#' @rdname sce
#'
#' @param h5ad_filepath character() absolute path to downloaded .h5ad file
#'
#' @importFrom zellkonverter readH5AD
#'
#' @return SingleCellExperiment object
h5ad_to_sce <- function(h5ad_filepath) {
    h5ad_sce <- zellkonverter::readH5AD(h5ad_filepath)
    h5ad_sce
}

#' @rdname sce
#'
#' @description convert a SummarizedExperiment assay to a
#' sparse matrix representation
#'
#' @param sce_obj SummarizedExperiment (or subclass) object
to_sparse_matrix <- function(matrix_obj){
    UseMethod(".to_sparse_matrix")
}

#' @rdname sce
#'
#' @description  function to change SingleCellLoomExperiment object's assay
#'  matrix into sparse representation (`dgCMatrix`)
#'
#' @param dgCMatrix sparse matrix
#'
#' @importFrom methods as
#'
#' @return dgCMatrix sparse matrix representation
.to_sparse_matrix.dgCMatrix <- function(c_matrix) {
    sparse_matrix <- as(c_matrix, "dgCMatrix")
    sparse_matrix
}

#' @rdname sce
#'
#' @description  function to change SingleCellLoomExperiment object's assay
#'  matrix into sparse representation (`dgRMatrix`)
#'
#' @param dgRMatrix sparse matrix
#'
#' @importFrom methods as
#'
#' @return dgCMatrix sparse matrix representation
.to_sparse_matrix.dgRMatrix <- function(r_matrix) {
    ## from: https://github.com/theislab/zellkonverter/issues/55
    sparse_matrix <- as(as(r_matrix, "CsparseMatrix"), "dgCMatrix")
    sparse_matrix
}


#' @rdname sce
#'
#' @description  function to change SingleCellLoomExperiment object's assay
#'  matrix into sparse representation (`DelayedMatrix`)
#'
#' @param DelayedMatrix sparse matrix
#'
#' @importFrom methods as
#'
#' @return dgCMatrix sparse matrix representation
.to_sparse_matrix.DelayedMatrix <- function(d_matrix) {
    sparse_matrix <- as(d_matrix, "dgCMatrix")
    sparse_matrix
}

#' @rdname sce
#'
#' @description  function to change SingleCellLoomExperiment object's assay
#'  matrix into sparse representation (`matrix`)
#'
#' @param matrix sparse matrix
#'
#' @importFrom methods as
#'
#' @return dgCMatrix sparse matrix representation
.to_sparse_matrix.matrix <- function(m_matrix) {
    sparse_matrix <- as(m_matrix, "dgCMatrix")
    sparse_matrix
}


#' @rdname sce
#'
#' @description convert a sparse matrix representation of the assay into a
#' reformatted tibble
#'
#' @param sparse_matrix dgCMatrix or dgRMatrix sparse matrix representation
sparse_mtx_to_assay_tbl <- function(sparse_matrix){
    ## extract file assay
    UseMethod(".sparse_matrix_to_assay_tbl")
}

#' @rdname sce
#'
#' @description convert a dgCMatrix sparse matrix representation of the assay
#' into a reformatted tibble
#'
#' @param sparse_matrix dgCMatrix sparse matrix representation
#'
#' @importFrom tibble tibble
#' @importFrom Matrix summary
#'
#' @return tibble with row_index, col_index, and value of each non-zero entry
#' in the assay matrix
## dgCMatrix implementation
.sparse_matrix_to_assay_tbl.dgCMatrix <- function(sparse_matrix) {
    stopifnot(
        ## sparse_matrix must inherit from class dgCMatrix
        `'sparse_matrix =' must inherit from class 'dgCMatrix'` =
            inherits(sparse_matrix, "dgCMatrix")
    )

    ## handling of sparse matrices depends on the initial file type
    ## help: https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    ##dim(sparse_matrix)
    ##length(sparse_matrix@i)
    ##length(sparse_matrix@p)
    ##length(sparse_matrix@x)

    ## increment each element of row_index by one to match R's 1-based indexing
    row_index <- sparse_matrix@i + 1
    values <- sparse_matrix@x

    ## diff between each consecutive integer in the sequence
    ## sparse_matrix@p gives number of non-zero values in each column
    ## as the sequence is made by reading the matrix
    ## bottom to top, then left to right

    n_per_column <-  diff(sparse_matrix@p)
    col_index <-  rep(seq_along(n_per_column), n_per_column)
    ##str(loom_col_index)

    # tibble is the row, column, and value for the entire matrix
    assay_tbl <- tibble(row_index, col_index, values)
    assay_tbl
}

## #' @rdname sce
## #'
## #' @description convert a dgRMatrix sparse matrix representation of the assay
## #' into a reformatted tibble
## #'
## #' @param sparse_matrix dgRMatrix sparse matrix representation
## #'
## #' @importFrom tibble tibble
## #' @importFrom Matrix summary
## #'
## #' @return tibble with row_index, col_index, and value of each non-zero entry
## #' in the assay matrix
## ## dgRMatrix implementation
## .sparse_matrix_to_assay_tbl.dgRMatrix <- function(sparse_matrix) {
##     stopifnot(
##         ## sparse_matrix must inherit from class dgRMatrix
##         `'sparse_matrix =' must inherit from class 'dgRMatrix'` =
##             inherits(sparse_matrix, "dgRMatrix")
##     )

##     ## handling of sparse matrices depends on the initial file type
##     ## same as with .loom matrices, just compressed row wise vs column wise

##     ## increment each element of col_index by one to match R's 1-based indexing
##     col_index <- sparse_matrix@j + 1
##     values <- sparse_matrix@x

##     ## diff between each consecutive integer in the sequence
##     ## sparse_matrix@p gives number of non-zero values in each row

##     n_per_column <-  diff(sparse_matrix@p)
##     row_index <-  rep(seq_along(n_per_column), n_per_column)
##     ##str(h5ad_col_index)

##     # tibble is the row, column, and value for the entire matrix
##     assay_tbl <- tibble(row_index, col_index, values)
##     assay_tbl
## }

#' @rdname sce
#'
#' @description function to change SingleCellExperiment or
#' SingleCellLoomExperiment object's `rowData` into a gene annotation table
#'
#' @param sce SingleCellExperiment or SingleCellLoomExperiment object
#'
#' @importFrom dplyr %>% as_tibble
#' @importFrom tibble add_column
#' @importFrom SummarizedExperiment rowData
#'
#' @return tibble of gene annotations
sce_rowdata_to_gene_tbl <- function(sce) {
    stopifnot(
        ## sce must be a SingleCellExperiment object
        `'sce =' must inherit from class 'SingleCellExperiment'` =
        inherits(sce, "SingleCellExperiment")
    )

    ## get number of genes to use as row index
    num_genes <- dim(sce)[1]

    # if rowData exists, start table with that
    if(dim(rowData(sce))[2] > 0) {
        gene_tbl <- dplyr::as_tibble(rowData(sce))
        ## add row index
        gene_tbl <- gene_tbl %>%
            tibble::add_column(row_index = 1:num_genes, .before = 1)
    ## else, start table with just row index
    } else {
        gene_tbl <- dplyr::as_tibble(list(row_index = 1:num_genes))
    }
    # if rownames exists, add
    if(length(rownames(sce)) > 0) {
        gene_tbl %>%
            tibble::add_column(rnames = rownames(sce))
    }
    ##print("Gene table looks like:")
    ##print(head(gene_tbl))
    gene_tbl
}

#' @rdname sce
#'
#' @description function to change SingleCellExperiment or
#' SingleCellLoomExperiment object's `colData` into a cell annotation table
#'
#' @param sce SingleCellExperiment or SingleCellLoomExperiment object
#'
#' @importFrom dplyr %>% as_tibble
#' @importFrom tibble add_column
#' @importFrom SummarizedExperiment colData
#'
#' @return tibble of cell annotations
sce_coldata_to_cell_tbl <- function(sce) {
    stopifnot(
        ## sce must be a SingleCellExperiment object
        `'sce =' must inherit from class 'SingleCellExperiment'` =
            inherits(sce, "SingleCellExperiment")
    )

    ## get number of cells to use as row index
    num_cells <- dim(sce)[2]

    # if colData exists, start table with that
    if(dim(colData(sce))[2] > 0) {
        cell_tbl <- dplyr::as_tibble(colData(sce))
        ## add col index
        cell_tbl <- cell_tbl %>%
            tibble::add_column(col_index = 1:num_cells, .before = 1)
        ## else, start table with just col index
    } else {
        cell_tbl <- dplyr::as_tibble(list(col_index = 1:num_cells))
    }
    # if colnames exists, add
    if(length(colnames(sce)) > 0) {
        cell_tbl %>%
            tibble::add_column(cnames = colnames(sce))
    }

    ##print("Cell tibble looks like:")
    ##print(head(cell_tbl))
    cell_tbl
}
