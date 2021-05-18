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
#' @importFrom DBI dbConnect dbExistsTable dbDisconnect
#' @importFrom RPostgres Postgres
#' @importFrom dplyr %>% copy_to mutate across add_row tbl collect filter %in%
#' @importFrom tibble tibble
#' @importFrom tidyselect vars_select_helpers
#' @importFrom hca .is_scalar_character
#' @importFrom getPass getPass
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
    ## gathering user credentials
    hcauser <- Sys.getenv("HCA_USER")
    if(is.null(hcauser) || hcauser == ""){
        hcauser <- readline(prompt="Database username: ")
    }
    ## print(paste("hcauser is: ", hcauser))


    hcapassword <- Sys.getenv("HCA_PASSWORD")
    if(is.null(hcapassword) || hcapassword == ""){
        hcapassword <- getPass(msg = "Database password: ",
                               noblank = TRUE,
                               forcemask = FALSE)
    }
    ## print(paste("hcapassword is: ", hcapassword))

    ## connect to database
    print("Establishing database connection...")
    db_connection <- DBI::dbConnect(RPostgres::Postgres(),
                                    host = "localhost",
                                    dbname = "bioc_hca",
                                    user = hcauser,
                                    port = 5432,
                                    password = hcapassword
    )

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
