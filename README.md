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
