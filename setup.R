# Dependency checks only. Installation is an explicit separate call.
dependencies <- read.csv("dependencies.csv",stringsAsFactors=FALSE)
check_dependencies <- function() {
  data.frame(package=dependencies$package,
    installed=vapply(dependencies$package,requireNamespace,logical(1),quietly=TRUE),
    row.names=NULL)
}
install_cran_dependencies <- function(repos="https://cloud.r-project.org") {
  wanted <- dependencies$package[dependencies$source=="CRAN"]
  missing <- wanted[!vapply(wanted,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) install.packages(missing,repos=repos)
}

install_bioconductor_dependencies <- function() {
  wanted <- dependencies$package[dependencies$source=="Bioconductor"]
  missing <- wanted[!vapply(wanted,requireNamespace,logical(1),quietly=TRUE)]
  if(!length(missing)) return(invisible(NULL))
  if(!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install(missing,ask=FALSE,update=FALSE)
}

# Install the pinned upstream GenomicSEM release for the checked native interfaces.
install_genomicsem <- function() {
  if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
  remotes::install_github(
    "GenomicSEM/GenomicSEM@0a63ac0ea01b61d28bd17e4a204e0fa561ce5040",
    upgrade="never")
}
