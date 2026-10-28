.libPaths(unique(c('mp10_pathology_regions_v2/Rlib','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
options(timeout=600)
for(pkg in c('SpaCET','Banksy')) install.packages(paste0('mp10_pathology_regions_v2/sources/',pkg),repos=NULL,type='source',lib='mp10_pathology_regions_v2/Rlib')
for(pkg in c('SpaCET','Banksy')) {library(pkg,character.only=TRUE);cat(pkg,as.character(packageVersion(pkg)),'\n')}
