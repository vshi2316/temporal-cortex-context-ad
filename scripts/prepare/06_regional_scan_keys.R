# Recover the regional-analysis scan identifiers from retained genetic frames.
if(!exists("study_path",mode="function")) source("config.R",encoding="UTF-8")
library(data.table)
result_root <- study_path("project","results")
frozen <- fread(file.path(result_root,"13H55_transition_ordered_structural_chain_release","01_frozen_frame.csv"))
base <- as.data.table(readRDS(file.path(result_root,"13F_R_only_structural_analysis","04_ADNI_R_only_longitudinal_structural_frame.rds")))
frozen[,`:=`(RID=as.character(RID),mri_date=as.IDate(mri_date))]
base[,`:=`(RID=as.character(RID),mri_date=as.IDate(as.character(mri_date)),
           primary_complete=is.finite(suppressWarnings(as.numeric(ad_vulnerable_index))))]
setorderv(base,c("RID","mri_date","primary_complete","IMAGEUID"),c(1L,1L,-1L,-1L),na.last=TRUE)
base <- base[, .SD[1L],by=.(RID,mri_date)]
keys <- merge(frozen[,.(RID,mri_date)],base[,.(RID,mri_date,main_imageuid=as.character(IMAGEUID))],
              by=c("RID","mri_date"),all.x=TRUE,sort=FALSE)
stopifnot(nrow(keys)==4659L,uniqueN(keys$RID)==695L,!anyNA(keys$main_imageuid),
          !anyDuplicated(keys,by=c("RID","mri_date")))
keys[,in_main_model:=TRUE]
dest <- study_path("project","input/analysis_definitions")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
fwrite(keys,file.path(dest,"scan_keys.csv"))
