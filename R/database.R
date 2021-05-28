.database_user <- local({
    ## `local()` creates an environment 'local' to the function, so
    ## that `database_user` is available for reading and writing for
    ## the duration of the session. We use this so that we only need
    ## to enter the username only once per session.

    database_user <- "" # user not yet provided

    function(user = "") {
        stopifnot(.is_scalar_character(user))

        ## determine user from the argument, or the system variable,
        ## or a previous value, or prompt, in that order
        if (!nzchar(user))
            user <- Sys.getenv("HCA_USER")
        if (!nzchar(user))
            user <- database_user
        if (!nzchar(user)) {
            user <- readline(prompt="Database username: ")
            if (!nzchar(user))
                stop("no username entered")
        }

        ## remember the user as database_user for next time, and
        ## return the user (without displaying)
        database_user <<- user
    }
})

#' @importFrom getPass getPass
.database_password <- local({
    ## `local()` creates an environment 'local' to the function, so
    ## that `database_user` is available for reading and writing for
    ## the duration of the session. We use this so that we only need
    ## to enter the username only once per session.

    database_password <- "" # password not yet provided

    function(password = "") {
        stopifnot(.is_scalar_character(password))

        ## determine password from the argument, or the system variable,
        ## or a previous value, or prompt, in that order
        if (!nzchar(password))
            password <- Sys.getenv("HCA_PASSWORD")
        if (!nzchar(password))
            password <- database_password
        if (!nzchar(password)) {
            password <- getPass(msg = "Database password: ",
                                noblank = TRUE,
                                forcemask = FALSE)
            if (is.null(password))
                stop("no password entered")
        }

        ## remember the password as database_password for next time,
        ## and return the password (without displaying)
        database_password <<- password
    }
})

#' @importFrom DBI dbConnect
#' @importFrom RPostgres Postgres
.database_connection <- function(user = "", password = "") {
    stopifnot(
        .is_scalar_character(user),
        .is_scalar_character(password)
    )

    ## gathering user credentials
    hcauser <- .database_user(user)
    hcapassword <-  .database_password(password)

    message("Establishing database connection...")
    dbConnect(
        RPostgres::Postgres(),
        host = "localhost",
        dbname = "bioc_hca",
        user = hcauser,
        port = 5432,
        password = hcapassword
    )
}

#' @title Populating database with data output of single cell experiments
#'
#' @rdname database
#'
#' @param file_tbl tbl_hca tibble of files as returned by hca::files()
#'
#' @importFrom dplyr %>% select left_join
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
    file_tbl_aug <- file_tbl %>%
        tibble::add_column(file_locations = file_locations)

    ## generate a tibble of file locations and associated projectIds
    project_titles <- file_tbl_aug$projectTitle
    file_filter <- filters(projectTitle = list(is = project_titles))
    test_proj <- projects(file_filter)
    file_and_projIds <- file_tbl_aug %>%
        left_join(test_proj, by = "projectTitle") %>%
        select("fileId", "version", "name", "projectTitle",
               "file_locations", "projectId")

    ## apply to each pair of file path and project ID
    mapply(
        .single_file_to_db_compact,
        file_and_projIds$file_locations,
        file_and_projIds$fileId,
        file_and_projIds$version,
        file_and_projIds$projectTitle,
        file_and_projIds$projectId
    )
}

#' @rdname database
#'
#' @param file_path character() location of experiment output file
#' @param fileId character() unique identifier of file
#' @param version character() file version
#' @param projectTitle character() title of file's associated project
#' @param projectId character() unique identifier of the file's
#' associated project
#'
#' @importFrom DBI dbExistsTable dbDisconnect
#' @importFrom dplyr %>% copy_to mutate across add_row tbl collect filter bind_rows
#' @importFrom tibble tibble
#' @importFrom tools file_ext
#' @importFrom tidyselect vars_select_helpers
#' @importFrom hca .is_scalar_character
.single_file_to_db_compact <- function(file_path,
                                        fileId,
                                        version,
                                        projectTitle,
                                        projectId) {
    stopifnot(
        ## file_path must be a non-null character vector
        `'file_path =' must be a non-null character vector` =
            .is_scalar_character(file_path),
        ## fileId must be a non-null character vector
        `'fileId =' must be a non-null character vector` =
            .is_scalar_character(fileId),
        ## version must be a non-null character vector
        `'version =' must be a non-null character vector` =
            .is_scalar_character(version),
        ## projectTitle must be a non-null character vector
        `'projectTitle =' must be a non-null character vector` =
            .is_scalar_character(projectTitle),
        ## projectId must be a non-null character vector
        `'projectId =' must be a non-null character vector` =
            .is_scalar_character(projectId)
    )

    ## connect to database
    db_connection <- .database_connection()

    ## first, check to see if file already exists in the database as not to
    ## duplicate data
    overview_table_exists <- db_connection %>%
        DBI::dbExistsTable("experiment_overviews")
    file_exists_in_db <- FALSE

    if(overview_table_exists){
        existing_experiments_tbl <- db_connection %>%
            tbl("experiment_overviews")
        files_available <- existing_experiments_tbl$file_id
        ## print("Files currently in db are: ")
        ## print(files_available)
        if(!is.null(files_available) && fileId %in% files_available){
            ## get instance of file already in database and check version
            existing_version <- existing_experiments_tbl %>%
                filter(file_id == fileId) %>%
                select(version) %>%
                as.character()
            if(existing_version == version){
                message(fileId, " with version ",
                        version, " already exists in the data")
                file_exists_in_db <- TRUE
            } else { # will need to add file to existing tables
                existing_assays_tbl <- db_connection %>%
                    tbl("assays_tbl")
                existing_genes_tbl <- db_connection %>%
                    tbl("genes_tbl")
                existing_cells_tbl <- db_connection %>%
                    tbl("cells_tbl")
            }
        } else {# will need to add file to existing tables
            existing_assays_tbl <- db_connection %>%
                tbl("assays_tbl")
            existing_genes_tbl <- db_connection %>%
                tbl("genes_tbl")
            existing_cells_tbl <- db_connection %>%
                tbl("cells_tbl")
        }
    } else {
        message("This is the first experiment to be added to the database")
        existing_experiments_tbl <- tibble(file_id = character(),
                                           version = character(),
                                           project_id = character(),
                                           project_title = character(),
                                           donor_organism_genus_species = character(),
                                           expression_data_type = character(),
                                           library_construction_approach = character(),
                                           pipeline_version = character(),
                                           specimen_from_organism_organ = character())

        existing_assays_tbl <- tibble(row_index = numeric(),
                                      col_index = numeric(),
                                      values = numeric())

        existing_genes_tbl <- tibble(row_index = integer(),
                                     gene = character(),
                                     antisense_reads = numeric(),
                                     duplicate_reads = numeric(),
                                     ensembl_ids = character(),
                                     fragments_per_molecule = numeric(),
                                     fragments_with_single_read_evidence = numeric(),
                                     gene_names = character(),
                                     genomic_read_quality_mean = numeric(),
                                     genomic_read_quality_variance = numeric(),
                                     genomic_reads_fraction_bases_quality_above_30_mean = numeric(),
                                     genomic_reads_fraction_bases_quality_above_30_variance = numeric(),
                                     molecule_barcode_fraction_bases_above_30_mean = numeric(),
                                     molecule_barcode_fraction_bases_above_30_variance = numeric(),
                                     molecules_with_single_read_evidence = numeric(),
                                     n_fragments = numeric(),
                                     n_molecules = numeric(),
                                     n_reads = numeric(),
                                     noise_reads = numeric(),
                                     number_cells_detected_multiple = numeric(),
                                     number_cells_expressing = numeric(),
                                     perfect_molecule_barcodes = numeric(),
                                     reads_mapped_exonic = numeric(),
                                     reads_mapped_intronic = numeric(),
                                     reads_mapped_multiple = numeric(),
                                     reads_mapped_uniquely = numeric(),
                                     reads_mapped_utr = numeric(),
                                     reads_per_fragment = numeric(),
                                     reads_per_molecule = numeric(),
                                     spliced_reads = numeric(),
                                     file_id = character())

        existing_cells_tbl <- tibble(col_index = integer(),
                                     cell_id = character(),
                                     antisense_reads = integer(),
                                     cell_barcode_fraction_bases_above_30_mean = numeric(),
                                     cell_barcode_fraction_bases_above_30_variance = numeric(),
                                     cell_names = character(),
                                     duplicate_reads = integer(),
                                     emptydrops_FDR = numeric(),
                                     emptydrops_IsCell = logical(),
                                     emptydrops_Limited = logical(),
                                     emptydrops_LogProb = numeric(),
                                     emptydrops_PValue = numeric(),
                                     emptydrops_Total = integer(),
                                     fragments_per_molecule = numeric(),
                                     fragments_with_single_read_evidence = integer(),
                                     genes_detected_multiple_observations = integer(),
                                     genomic_read_quality_mean = numeric(),
                                     genomic_read_quality_variance = numeric(),
                                     genomic_reads_fraction_bases_quality_above_30_mean = numeric(),
                                     genomic_reads_fraction_bases_quality_above_30_variance = numeric(),
                                     input_id = character(),
                                     molecule_barcode_fraction_bases_above_30_mean = numeric(),
                                     molecule_barcode_fraction_bases_above_30_variance = numeric(),
                                     molecules_with_single_read_evidence = integer(),
                                     n_fragments = integer(),
                                     n_genes = integer(),
                                     n_mitochondrial_genes = integer(),
                                     n_mitochondrial_molecules = integer(),
                                     n_molecules = integer(),
                                     n_reads = integer(),
                                     noise_reads = integer(),
                                     pct_mitochondrial_molecules = numeric(),
                                     perfect_cell_barcodes = integer(),
                                     perfect_molecule_barcodes = integer(),
                                     reads_mapped_exonic = integer(),
                                     reads_mapped_intergenic = integer(),
                                     reads_mapped_intronic = integer(),
                                     reads_mapped_multiple = integer(),
                                     reads_mapped_too_many_loci = integer(),
                                     reads_mapped_uniquely = integer(),
                                     reads_mapped_utr = integer(),
                                     reads_per_fragment = numeric(),
                                     reads_unmapped = integer(),
                                     spliced_reads = integer(),
                                     file_id = character())
    }

    ## if file does not already exist in the database, proceed with appending
    ## it to the gene annotation and cell annotation tables
    ## https://dplyr.tidyverse.org/reference/bind.html
    if(!file_exists_in_db) {
        file_ext <- tools::file_ext(file_path)
        sce <- switch(file_ext,
                      "loom" = loom_to_sce(file_path),
                      "h5ad" = h5ad_to_sce(file_path))
        sparse_matrix <- switch(file_ext,
                                "loom" = loom_sparse_matrix(file_path),
                                "h5ad" = h5ad_sparse_matrix(file_path))

        new_assay_tbl <- sparse_mtx_to_assay_tbl(sparse_matrix)

        new_gene_tbl <- sce_rowdata_to_gene_tbl(sce)

        new_cell_tbl <- sce_coldata_to_cell_tbl(sce)

        ## if any column in any table is of type "raw" i.e. byte data
        ## conversion is needed
        new_assay_tbl_recast <- new_assay_tbl %>%
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical)) %>%
            add_column(file_id = fileId)

        new_gene_tbl_recast <- new_gene_tbl %>%
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical)) %>%
            add_column(file_id = fileId) %>%
            dplyr::rename(gene = Gene)

        new_cell_tbl_recast <- new_cell_tbl %>%
            #mutate(across(where(is.raw), ~ rawToChar(.x, multiple = T)))
            mutate(across(tidyselect::vars_select_helpers$where(is.raw),
                          as.logical)) %>%
            add_column(file_id = fileId) %>%
            dplyr::rename(cell_id = CellID)

        ## appending to existing tables
        fin_assays_tbl <- existing_assays_tbl %>%
            collect() %>%
            bind_rows(new_assay_tbl_recast)
        dplyr::copy_to(db_connection, fin_assays_tbl,
                       name = "assays_tbl",
                       temporary = FALSE,
                       overwrite = TRUE)

        fin_genes_tbl <- existing_genes_tbl %>%
            collect() %>%
            bind_rows(new_gene_tbl_recast)
        dplyr::copy_to(db_connection, fin_genes_tbl,
                       name = "genes_tbl",
                       temporary = FALSE,
                       overwrite = TRUE)

        fin_cells_tbl <- existing_cells_tbl %>%
            collect() %>%
            bind_rows(new_cell_tbl_recast)
        dplyr::copy_to(db_connection, fin_cells_tbl,
                       name = "cells_tbl",
                       temporary = FALSE,
                       overwrite = TRUE)

        ## add details to overview table
        fin_experiments_tbl <- existing_experiments_tbl %>%
            collect() %>%
            add_row(file_id = fileId,
                    version = version,
                    project_id = projectId,
                    project_title = projectTitle,
                    donor_organism_genus_species = metadata(sce)[['donor_organism.genus_species']],
                    expression_data_type = metadata(sce)[['expression_data_type']],
                    library_construction_approach = metadata(sce)[['library_preparation_protocol.library_construction_approach']],
                    pipeline_version = metadata(sce)[['pipeline_version']],
                    specimen_from_organism_organ = metadata(sce)[['specimen_from_organism.organ']])
        dplyr::copy_to(db_connection, fin_experiments_tbl,
                       name = "experiment_overviews",
                       temporary = FALSE,
                       overwrite = TRUE)

        DBI::dbDisconnect(db_connection)
    }
}
