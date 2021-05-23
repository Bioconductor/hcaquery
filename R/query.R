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
#' @examples available_tables("genes_tbl")
#'
#' @export
available_tables <- function() {
    .AVAILABLE_TABLES
}
#'
#' @rdname query
#'
#' @param table_name character() name of database table to be described
#'
#' @importFrom DBI dbExistsTable dbDisconnect
#' @importFrom dplyr %>% copy_to mutate across add_row tbl collect filter
#' @importFrom tibble tibble
#' @importFrom tidyselect vars_select_helpers
#' @importFrom hca .is_scalar_character
#'
#' @examples table_description("genes_tbl")
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
