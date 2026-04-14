r <- getOption("repos")
r["CRAN"] <- "https://cran.r-project.org"
options(repos = r)

# Install xgboost 1.7.8.1 from CRAN archive
install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.8.1.tar.gz",
                 repos = NULL, type = "source")

# Install bambu 3.5.1 from Bioconductor archive
# (3.5.1 is not on conda; install via remotes from the Bioc git tag)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_version("bambu",
                         version = "3.5.1",
                         repos = BiocManager::repositories())
