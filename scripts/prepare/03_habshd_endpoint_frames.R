# Endpoint-specific HABS-HD preparation only
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

project_dir <- study_path("project")
habs_root <- study_path("habs_release")
freeze_dir <- file.path(project_dir, "results", "12B_ADNI_HABSHD_structural_MRI_feature_freeze")
out_dir <- file.path(project_dir, "results", "12E_HABSHD_external_structural_MRI_validation")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

need <- function(p){if(!file.exists(p))stop("Missing input: ",p,call.=FALSE);p}
num <- function(x){z=suppressWarnings(as.numeric(x));z[z<=-9000]=NA_real_;z}
zfun <- function(x){x=num(x);s=sd(x,na.rm=TRUE);if(!is.finite(s)||s==0)return(rep(NA_real_,length(x)));(x-mean(x,na.rm=TRUE))/s}
has_var <- function(x) length(unique(x[!is.na(x)])) >= 2L
normalize_id <- function(x){
  z <- trimws(as.character(x))
  n <- suppressWarnings(as.numeric(z))
  z[is.finite(n)] <- format(n[is.finite(n)], scientific=FALSE, trim=TRUE, digits=15)
  z <- sub("\\.0+$", "", z)
  z[z %in% c("", "NA", "N/A", ".")] <- NA_character_
  z
}
clean_names <- function(x){x=trimws(as.character(x));x=gsub("[^A-Za-z0-9]+","_",x);x=gsub("^_+|_+$","",x);make.unique(x,sep="_")}
log_msg <- function(...)message(sprintf("[%s] [INFO] %s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),paste0(...)))

month_num <- function(x){
  y=suppressWarnings(as.integer(x)); bad=!is.finite(y)
  if(any(bad)){m=match(tolower(substr(trimws(as.character(x[bad])),1,3)),tolower(month.abb));y[bad]=m}
  y
}
make_date <- function(year,month){y=suppressWarnings(as.integer(year));m=month_num(month);as.IDate(sprintf("%04d-%02d-15",y,m))}

read_inclusions <- function(path, pattern="inclusion|included"){
  ss=excel_sheets(path);keep=ss[grepl(pattern,ss,ignore.case=TRUE)&!grepl("excl|pending",ss,ignore.case=TRUE)]
  if(!length(keep))stop("No inclusion sheets: ",basename(path),call.=FALSE)
  rbindlist(lapply(keep,function(s){x=suppressWarnings(as.data.table(read_excel(path,sheet=s,col_types="text",na=c("","NA","N/A","."),.name_repair="minimal")));setnames(x,clean_names(names(x)));x[,source_sheet:=s];x}),fill=TRUE)
}
find_col <- function(d, patterns, required=TRUE){
  for(p in patterns){hit=grep(p,names(d),value=TRUE,ignore.case=TRUE);if(length(hit))return(hit[1])}
  if(required)stop("Column not found: ",paste(patterns,collapse=" | "),call.=FALSE)
  NA_character_
}
has_var <- function(x)length(unique(x[!is.na(x)]))>=2L

log_msg("Step 12E started: HABS-HD locked external structural MRI validation.")
ct=as.data.table(readRDS(need(file.path(freeze_dir,"01_HABSHD_frozen_cortical_features.rds"))))
ct[,participant_id:=normalize_id(participant_id)]

# Recover MRI month/year from the source release without changing frozen indices.
ct_source=read_inclusions(need(file.path(habs_root,"HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")))
idc=find_col(ct_source,c("^Med_ID$","^MED_ID$"));tpc=find_col(ct_source,c("^Time_point_MRI$","^Time_point$"));mc=find_col(ct_source,c("^MRI_month$"));yc=find_col(ct_source,c("^MRI_year$"))
ct_dates=unique(ct_source[,.(participant_id=normalize_id(get(idc)),mri_timepoint=as.character(get(tpc)),mri_date=make_date(get(yc),get(mc)))])
ct=merge(ct,ct_dates,by=c("participant_id","mri_timepoint"),all.x=TRUE,sort=FALSE)

# HABS-HD hippocampus and whole-brain volume files.
hip=read_inclusions(need(file.path(habs_root,"HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx")))
hid=find_col(hip,c("^Med_ID$"));hmo=find_col(hip,c("^MRI_month$"));hyr=find_col(hip,c("^MRI_year$"));hl=find_col(hip,c("^L_hippocampus$"));hr=find_col(hip,c("^R_hippocampus$"));hicv=find_col(hip,c("^ICV_HD$"))
hip[,`:=`(participant_id=normalize_id(get(hid)),hip_date=make_date(get(hyr),get(hmo)),hippocampus_icv_ratio=(num(get(hl))+num(get(hr)))/num(get(hicv)))]
hip=hip[is.finite(hippocampus_icv_ratio),.(participant_id,hip_date,hippocampus_icv_ratio)]

brain=read_inclusions(need(file.path(habs_root,"HABS_HD_RP_8_Whole_Brain_Volume_Data_BraskieLab_20260501.xlsx")))
bid=find_col(brain,c("^Med_ID$"));bmo=find_col(brain,c("^MRI_month$"));byr=find_col(brain,c("^MRI_year$"));bv=find_col(brain,c("^BrainSegNotVent$"));bicv=find_col(brain,c("^ICV_HD$"))
brain[,`:=`(participant_id=normalize_id(get(bid)),brain_date=make_date(get(byr),get(bmo)),brain_icv_ratio=num(get(bv))/num(get(bicv)))]
brain=brain[is.finite(brain_icv_ratio),.(participant_id,brain_date,brain_icv_ratio)]

# Amyloid and Tau PET releases.
amy=read_inclusions(need(file.path(habs_root,"HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx")))
aid=find_col(amy,c("^Med_ID$"));amo=find_col(amy,c("FBB_PET_month"));ayr=find_col(amy,c("FBB_PET_year"));acl=find_col(amy,c("Centiloid_Units_GlobalAmyloid","Global_FBB_SUVR"))
amy[,`:=`(participant_id=normalize_id(get(aid)),amyloid_date=make_date(get(ayr),get(amo)),amyloid_value=num(get(acl)))]
amy=amy[is.finite(amyloid_value),.(participant_id,amyloid_date,amyloid_value)]

tau=read_inclusions(need(file.path(habs_root,"HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx")))
tid=find_col(tau,c("^Med_ID$"));tmo=find_col(tau,c("PI2620_PET_month"));tyr=find_col(tau,c("PI2620_PET_year"));tv=find_col(tau,c("Median_Medial_Temporal_PI2620_SUVR","Medial_Temporal_PI2620_SUVR"))
tau[,`:=`(participant_id=normalize_id(get(tid)),tau_date=make_date(get(tyr),get(tmo)),tau_value=num(get(tv)))]
tau=tau[is.finite(tau_value),.(participant_id,tau_date,tau_value)]

# Cognitive memory composite and global CDR.
zs=fread(need(file.path(habs_root,"zscores_12Aug2026.csv")),showProgress=FALSE);zs[,participant_id:=normalize_id(MED_ID)];zs[,cog_date:=make_date(VISIT_YEAR,VISIT_MONTH)]
memvars=intersect(c("LM1_AB_ZSCORE","LM2_AB_ZSCORE","SEVLT_T1235_ZSCORE","SEVLT_DR_ZSCORE"),names(zs))
for(v in memvars)zs[,(v):=num(get(v))]
zs[,memory_composite:=rowMeans(.SD,na.rm=TRUE),.SDcols=memvars];zs[!is.finite(memory_composite),memory_composite:=NA_real_]
zs=zs[is.finite(memory_composite),.(participant_id,cog_date,memory_composite)]

cdr=fread(need(file.path(habs_root,"consensuscdr_12Aug2026.csv")),showProgress=FALSE);cdr[,participant_id:=normalize_id(MED_ID)];cdr[,cdr_date:=make_date(VISIT_YEAR,VISIT_MONTH)];cdr[,cdr_sum:=num(CDR_SUM)];cdr=cdr[is.finite(cdr_sum),.(participant_id,cdr_date,cdr_sum)]

gen=as.data.table(read_excel(need(file.path(habs_root,"RP_HD_7_Genomics.xlsx")),col_types="text"));setnames(gen,clean_names(names(gen)));gid=find_col(gen,c("^Med_ID$"));ga=find_col(gen,c("APOE4_Positivity"));gen[,participant_id:=normalize_id(get(gid))];gen[,APOE4_positive:=fifelse(tolower(as.character(get(ga)))%in%c("1","yes","positive","pos","true"),1,fifelse(tolower(as.character(get(ga)))%in%c("0","no","negative","neg","false"),0,num(get(ga))))];gen=unique(gen[is.finite(APOE4_positive),.(participant_id,APOE4_positive)])

nearest <- function(index,candidate,date_col,prefix,max_days){
  a=copy(index);b=copy(candidate);a[,row__:=.I];setnames(b,date_col,"candidate_date__")
  z=merge(a,b,by="participant_id",allow.cartesian=TRUE,all.x=TRUE);z[,gap__:=abs(as.integer(mri_date-candidate_date__))];setorder(z,row__,gap__);z=z[,head(.SD,1),by=row__]
  cand=setdiff(names(b),"participant_id");z[!is.finite(gap__)|gap__>max_days,(cand):=NA]
  setnames(z,cand,paste0(prefix,cand));z[,c("row__","gap__"):=NULL];z
}

linked=nearest(ct,hip,"hip_date","hip_",548L)
linked=nearest(linked,brain,"brain_date","brain_",548L)
linked=nearest(linked,amy,"amyloid_date","amyloid_",548L)
linked=nearest(linked,tau,"tau_date","tau_",548L)
linked=nearest(linked,zs,"cog_date","cog_",365L)
linked=nearest(linked,cdr,"cdr_date","cdr_",365L)
linked=merge(linked,gen,by="participant_id",all.x=TRUE,sort=FALSE)
setorder(linked,participant_id,mri_date)
index=linked[,head(.SD,1),by=participant_id]
index_amyloid=linked[is.finite(amyloid_amyloid_value),head(.SD,1),by=participant_id]
index_tau=linked[is.finite(tau_tau_value),head(.SD,1),by=participant_id]
index_hippocampus=linked[is.finite(hip_hippocampus_icv_ratio),head(.SD,1),by=participant_id]
index_brain=linked[is.finite(brain_brain_icv_ratio),head(.SD,1),by=participant_id]
index_memory=linked[is.finite(cog_memory_composite),head(.SD,1),by=participant_id]
index_cdr=linked[is.finite(cdr_cdr_sum),head(.SD,1),by=participant_id]
index_apoe=linked[is.finite(APOE4_positive),head(.SD,1),by=participant_id]

outcomes=c("global_cortical_thickness","ad_vulnerable_index","frontoposterior_gradient","spatial_dispersion","mean_absolute_asymmetry","hip_hippocampus_icv_ratio","brain_brain_icv_ratio")
actual_outcomes=c("global_cortical_thickness","ad_vulnerable_index","frontoposterior_gradient","spatial_dispersion","mean_absolute_asymmetry","hip_hippocampus_icv_ratio","brain_brain_icv_ratio")
frames=list(index=index,index_amyloid=index_amyloid,index_tau=index_tau,index_hippocampus=index_hippocampus,index_brain=index_brain,index_memory=index_memory,index_cdr=index_cdr,index_apoe=index_apoe)
frames=lapply(frames,function(q){
  for(v in actual_outcomes)if(v%in%names(q))q[,paste0(v,"_z"):=zfun(get(v))]
  q[,age_z:=zfun(age)];q[,scanner_factor:=factor(scanner)];q
})
index=frames$index;index_amyloid=frames$index_amyloid;index_tau=frames$index_tau
index_hippocampus=frames$index_hippocampus;index_brain=frames$index_brain
index_memory=frames$index_memory;index_cdr=frames$index_cdr;index_apoe=frames$index_apoe

saveRDS(frames,file.path(out_dir,"06_HABSHD_endpoint_specific_validation_frames.rds"),compress="xz")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"99_sessionInfo.txt"))
