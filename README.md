# `hcaquery`
- This package, in conjunction with the  [**hca package**](https://bioconductor.org/packages/release/bioc/html/hca.html), allows
for exploration of the Human Cell Atlas' data, made available by their API.

- Its core functionality is the generation of relational database
tables containing data from `.loom` and `.h5ad` files

- This is done by using the `hca` package to obtain `.loom` and `.h5ad` files
from the HCA API, each of which can be parse into a
[`SingleCellExperiment`](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) object, providing us with access to
(1) cell annotations, (2) gene annotations, (3) assay matrices,
and (4) file metadata.

- Each of these data categories in turn corresponds to a table in the relational
database schema:
    1. **genes_tbl**
    2. **cells_tbl**
    3. **assays_tbl**
    4. **experiment_overviews**

### Examples (*in development*)
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
