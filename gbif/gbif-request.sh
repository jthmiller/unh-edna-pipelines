




## To get a download key
curl --include --user userName:PASSWORD --header "Content-Type: application/json" --data @query.json https://api.gbif.org/v1/occurrence/download/request


## A download key is returned. Querying that download key shows the download information, including the download link and DOI once the download is ready. Run this repeatedly, until you see SUCCEEDED, replacing the key with the key for your download
## GBIF.org (29 January 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.paw44d
curl -Ss https://api.gbif.org/v1/occurrence/download/0001005-130906152512535

## You can then download the resulting file:
curl --location --remote-name https://api.gbif.org/occurrence/download/request/0001005-130906152512535.zip
