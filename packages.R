# packages.R - packages not available via conda
install.packages('BiocManager', repos='https://cloud.r-project.org')
BiocManager::install(c('MSstatsLiP', 'EmpiricalBrownsMethod'), ask=FALSE)
devtools::install_github('yaowuliu/ACAT')
