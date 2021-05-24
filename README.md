## Files to Database
- functions will take either a `.loom` file or `.h5ad` file and create a triplet
of database tables, one per experiment
    1. `<project>_gene_annotations` for row-wise gene annotations of an assay
    matrix
    2. `<project>_cell_annotations` for column-wise cell annotations of an
    assay matrix
    3. `<project>_assay_values` for the row index, column index, and cell
    value for all ***non-zero*** values of the assay matrix
- this makes use of a dockerized PostgreSQL database (see `docker-compose.yaml`
for details)
- some useful aliases:
```
alias hcadbup='docker-compose -f docker-compose.yaml up -d && docker-compose -f docker-compose.yaml logs -f'

alias hcadbdownv='docker-compose -f docker-compose.yaml down -v'

alias hcadb='docker exec -ti bioc-hca-db psql -U hca_user bioc_hca'
```

- use `\dn` to list database schemas

### Password Management
- User must provide a username and password for connecting to the database
- The program will first check to see if the appropriate environment variables
are set in `.Renviron` (or your `.Renviron` of choice)
for `HCA_USER` and `HCA_PASSWORD`
- If they are not, the user will be prompted to input their username
and password

### Examples
```
devtools::load_all()
library(hca)
library(dplyr)
```
- loading the two smallest `.loom` files
```
loom_filter <- hca::filters(fileFormat = list(is = c("loom")))
loom_tbl <- hca::files(filters = loom_filter,
                            size = 2, sort = "fileSize", order = "asc")
files_to_db(loom_tbl)
```

- loading `.loom` files produced by a single project, and processed by the HCA
```{r}
filters <- filters(
    fileFormat = list(is = "loom"),
    fileSource = list(is = "DCP/2 Analysis")
)
files <- files(filters) # how many? 46 files

files |>  # number and total size per project
    group_by(projectTitle) |>
    summarize(n = n(), GB = sum(size) / (1024^3)) |>
    arrange(GB)

filters <- filters(
    ## a particular project, with 3 files, 1/2 GB
    projectId = list(is = "88ec040b-8705-4f77-8f41-f81e57632f7d"),
    fileFormat = list(is = "loom"),
    fileSource = list(is = "DCP/2 Analysis")
)

loom_tbl <- files(filters)
files_to_db(loom_tbl)
```

- loading the
```
h5ad_filter <- hca::filters(fileFormat = list(is = c("h5ad")))
h5ad_tbl <- hca::files(filters = h5ad_filter,
                            size = 2, sort = "fileSize", order = "asc")

files_to_db(h5ad_tbl)
```
