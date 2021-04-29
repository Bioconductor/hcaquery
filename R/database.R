#' @title Populating database with data output of single cell experiments
#'
#' @rdname database
#'
#' @param file_tbl tbl_hca tibble of files as returned by hca::files()
#'
#' @importFrom dplyr select left_join
#' @importFrom tibble add_column
#' @importFrom hca filters files_download projects
#'
#' @export
files_to_db <- function(file_tbl = NULL) {
    stopifnot(
        ## input must be of class files_tbl_hca
        inherits(file_tbl, "files_tbl_hca")
    )

    ## download files and return location of files
    file_locations <- files_download(file_tbl)
    file_tbl_aug <- file_tbl |>
        tibble::add_column(file_locations = file_locations)

    ## generate a tibble of file locations and associated projectIds
    project_titles <- file_tbl_aug$projectTitle
    file_filter <- filters(projectTitle = list(is = project_titles))
    test_proj <- projects(file_filter)
    file_and_projIds <- file_tbl_aug |>
        left_join(test_proj, by = "projectTitle") |>
        select("fileId", "name", "projectTitle", "file_locations", "projectId")

    ## apply to each pair of file path and project ID
    mapply(
        .single_file_to_db,
        file_and_projIds$file_locations,
        file_and_projIds$fileId,
        file_and_projIds$projectTitle,
        file_and_projIds$projectId
    )
}

#' @rdname database
#'
#' @param file_path character() location of experiment output file
#' @param fileId character() unique identifier of file
#' @param projectTitle character() title of file's associated project
#' @param projectId character() unique identifier of the file's associated project
#'
#' @importFrom DBI dbConnect dbExistsTable
#' @importFrom RPostgres Postgres
#' @importFrom rstudioapi askForPassword
#' @importFrom dplyr copy_to mutate across add_row tbl collect
#' @importFrom tibble tibble
#' @importFrom tools file_ext
#' @importFrom tidyselect vars_select_helpers
#' @importFrom hca .is_scalar_character
.single_file_to_db <- function(file_path, fileId, projectTitle, projectId) {
    stopifnot(
        ## file_path must be a non-null character vector
        `'file_path =' must be a non-null character vector` =
            .is_scalar_character(file_path),
        ## fileId must be a non-null character vector
        `'fileId =' must be a non-null character vector` =
            .is_scalar_character(fileId),
        ## projectTitle must be a non-null character vector
        `'projectTitle =' must be a non-null character vector` =
            .is_scalar_character(projectTitle),
        ## projectId must be a non-null character vector
        `'projectId =' must be a non-null character vector` =
            .is_scalar_character(projectId)
    )
    ## connect to database
    con <- DBI::dbConnect(RPostgres::Postgres(),
                          host = "localhost",
                          dbname = "bioc_hca",
                          user = "hca_user",
                          port = 5432,
                          password = rstudioapi::askForPassword("Database password")
    )

    ## first, check to see if file already exists in the database as not to
    ## duplicate data
    overview_table_exists <- con |> DBI::dbExistsTable("experiment_overviews")
    file_exists_in_db <- FALSE

    if(overview_table_exists){
        existing_experiments_tbl <- con |>
                                    tbl("experiment_overviews") |>
                                    collect()
        files_available <- existing_experiments_tbl$fileId
        if(!is.null(files_available) && fileId %in% files_available){
            file_exists_in_db <- TRUE
        }
    } else {
        existing_experiments_tbl <- tibble(fileId = character(),
                                           projectId = character(),
                                           projectTitle = character())
    }

    ## if file does not already exist in the database, proceed with adding it
    if(!file_exists_in_db){
        file_ext <- tools::file_ext(file_path)
        sce <- switch(file_ext,
                      "loom" = loom_to_sce(file_path),
                      "h5ad" = h5ad_to_sce(file_path))
        sparse_matrix <- switch(file_ext,
                                "loom" = loom_sparse_matrix(file_path),
                                "h5ad" = h5ad_sparse_matrix(file_path))

        assay_tbl <- switch(file_ext,
                            "loom" = loom_sparse_mtx_to_assay_tbl(sparse_matrix),
                            "h5ad" = h5ad_sparse_mtx_to_assay_tbl(sparse_matrix))

        gene_tbl <- sce_rowdata_to_gene_tbl(sce)

        cell_tbl <- sce_coldata_to_cell_tbl(sce)

        ## if any column in any table is of type "raw" i.e. byte data
        ## conversion is needed
        assay_tbl_recast <- assay_tbl |>
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical))

        gene_tbl_recast <- gene_tbl |>
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical))

        cell_tbl_recast <- cell_tbl |>
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical))

        ## figure out how we want to name tables
        dplyr::copy_to(con, assay_tbl_recast,
                       paste(c(fileId, "assay"), collapse = "_"),
                       temporary = FALSE)
        dplyr::copy_to(con, gene_tbl_recast,
                       paste(c(fileId, "gene"), collapse = "_"),
                       temporary = FALSE)
        dplyr::copy_to(con, cell_tbl_recast,
                       paste(c(fileId, "cell"), collapse = "_"),
                       temporary = FALSE)

        ## add details to overview table
        existing_experiments_tbl <- existing_experiments_tbl |>
                                    add_row(fileId = fileId,
                                            projectId = projectId,
                                            projectTitle = projectTitle)
        dplyr::copy_to(con, existing_experiments_tbl,
                       name = "experiment_overviews",
                       temporary = FALSE,
                       overwrite = TRUE)
    }
}
