# bin/coloc_susie.R needs both coloc (for coloc.susie()/coloc.abf()) and
# susieR (for susie_rss()) in the same image, plus optparse/data.table like
# every other custom container in this repo. Neither the existing
# ghcr.io/farissham/coloc:5.2.3 nor ghcr.io/farissham/susie image bundles both,
# so this is its own dedicated container rather than trying to extend either —
# same rationale as docker/coloc/Dockerfile and docker/susie/Dockerfile: build
# our own with everything the script actually needs, versions pinned to match
# what's already used elsewhere in this repo (coloc 5.2.3 from CRAN, susieR
# 0.14.2 via remotes — same as docker/susie/Dockerfile).
FROM rocker/r-ver:4.4.0

RUN install2.r --error \
        data.table \
        optparse \
        coloc \
        remotes && \
    Rscript -e 'remotes::install_version("susieR", version="0.14.2", repos="https://cran.r-project.org")'
