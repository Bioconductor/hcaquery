.AVAILABLE_TABLES <- c(
    "genes_tbl",
    "cells_tbl",
    "assays_tbl",
    "experiment_overviews"
)

#' @title Functions for querying and interacting with the database
#'
#' @rdname query
#'
#' @name available_tables
#' @description funtion returning a list of database table available for
#' querying
#'
#' @examples \dontrun{available_tables()}
#'
#' @export
available_tables <- function() {
    .AVAILABLE_TABLES
}
#'
#' @rdname query
#'
#' @name table_description
#' @description function returning a description of a specified database table
#' @param table_name character() name of database table to be described
#'
#' @importFrom dplyr %>% copy_to mutate across add_row tbl collect filter slice_sample
#' @importFrom tibble tibble
#' @importFrom tidyselect vars_select_helpers
#' @importFrom hca .is_scalar_character
#' @importFrom utils head
#'
#' @examples \dontrun{table_description("genes_tbl")}
#'
#' @export
table_description <- function(table_name) {
    stopifnot(
        ## table_name must be one of the permissible table in the database
        `table_name must be one of those returned by available_tables()` =
            .is_scalar_character(table_name) &&
            (table_name %in% available_tables())
    )

    ## connect to database
    db_connection <- .database_connection()

    tbl_choice <- tbl(db_connection, table_name) %>% head() %>% collect()

    # column names
    col_names <- colnames(tbl_choice)
    # data_types
    dtypes <- sapply(tbl_choice, class)
    # taking a row of the data to bee used as an example
    ex_row <- unname(t(slice_sample(tbl_choice, n = 1))[,1])

    table_description_tbl <- tibble::tibble(columns = col_names,
                                               data_types = dtypes,
                                               example = ex_row)

    table_description_tbl
}

#' @rdname query
#'
#' @name hca_tables
#' @description function returning instances of all queryable tables
#'
#' @importFrom dplyr %>% tbl
#'
#' @examples \dontrun{hca_tables()}
#'
#' @export
hca_tables <- function(){
    ## connect to database
    db_connection <- .database_connection()

    experiments_tbl <- db_connection %>%
        tbl("experiment_overviews")
    assays_tbl <- db_connection %>%
        tbl("assays_tbl")
    genes_tbl <- db_connection %>%
        tbl("genes_tbl")
    cells_tbl <- db_connection %>%
        tbl("cells_tbl")

    hca_tbls <- list(experiments_tbl = experiments_tbl,
                     assays_tbl = assays_tbl,
                     genes_tbl = genes_tbl,
                     cells_tbl = cells_tbl)
}

#' @rdname query
#'
#' @name hca_sql_query
#' @description function to execute a SQL statement against tables in the
#' database
#'
#' @param sql_statement character() SQL statement to be interpolated
#' and executed
#'
#' @export
hca_sql_query <- function(sql_statement){
    ## connect to database
    db_connection <- .database_connection()
    ## security considerations
    ## see: https://shiny.rstudio.com/articles/sql-injections.html

}

#' @rdname query
#'
#' @name hca_file_gene_query
#' @description function to query a file for a specific subset of genes
#'
#' @importFrom dplyr %>% tbl filter collect right_join left_join mutate
#' @importFrom Matrix sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @param genes character() genes of interest
#' @param file_ident character() ID of file to query
#'
#' @export
hca_file_gene_query <- function(genes = character(),
                                file_ident = character()){
    ## connect to database
    db_connection <- .database_connection()

    row_data_all <-
        tbl(db_connection, "genes_tbl") %>%
        filter(file_id == file_ident) %>%
        filter(gene %in% genes)

    row_data <-
        row_data_all %>%
        select(-c("row_index"))

    row_idx <-
        row_data_all %>%
        select(row_index)

    col_data_all <-
        tbl(db_connection, "cells_tbl") %>%
        filter(file_id == file_ident)

    col_data <-
        col_data_all %>%
        select(-c("col_index"))

    col_idx <-
        col_data_all %>%
        select(col_index)

    assay_idx <-
        tbl(db_connection, "assays_tbl") %>%
        right_join(row_idx) %>%
        left_join(col_idx) %>%
        collect()

    ## need to reset the row_index
    assay_idx <- assay_idx %>%
        dplyr::mutate(row_index_reset = row_index - min(row_index) + 1)

    exp_metadata <- tbl(db_connection, "experiment_overviews") %>%
                    filter(file_id == file_ident) %>%
                    collect()

    metadata <- list(donor_organism.genus_species = exp_metadata$donor_organism_genus_species,
                     expression_data_type = exp_metadata$expression_data_type,
                     library_preparation_protocol.library_construction_approach = exp_metadata$library_construction_approach,
                     pipeline_version = exp_metadata$pipeline_version,
                     specimen_from_organism.organ = exp_metadata$specimen_from_organism_organ)

    row_idx_dims <- row_idx %>% collect() %>% dim()
    col_idx_dims <- col_idx %>% collect() %>% dim()

    row_idx_length <- row_idx_dims[1]
    col_idx_length <- col_idx_dims[1]

    ## reconstruct a SingleCellExperiment object
    assay_matrix <- sparseMatrix(i = assay_idx$row_index_reset,
                                 j = assay_idx$col_index,
                                 x = assay_idx$values,
                                 dims = c(row_idx_length, col_idx_length))

    sce <- SingleCellExperiment(assays = assay_matrix,
                                colData = col_data,
                                rowData = row_data,
                                metadata = metadata)
}

