r <- getOption("repos")
r["CRAN"] <- "https://cran.r-project.org"
options(repos = r)
install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.8.1.tar.gz",
                 repos = NULL, type = "source")
