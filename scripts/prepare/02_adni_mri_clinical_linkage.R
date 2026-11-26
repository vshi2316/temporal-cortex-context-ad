# ADNI MRI, polygenic scores and visit linkage
if (!exists("study_path", mode="function")) source("config.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

project_dir <- study_path("project")
adni_root <- study_path("cohorts", "ADNI/data")
out_dir <- file.path(project_dir,"results","12C_ADNI_common_structural_indices_multimodal_linkage")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
need <- function(p){if(!file.exists(p))stop("Missing input: ",p,call.=FALSE);p}
log_msg <- function(...)message(sprintf("[%s] [INFO] %s",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),paste0(...)))
num <- function(x)suppressWarnings(as.numeric(x))
date0 <- function(x)suppressWarnings(as.IDate(x))

fs_file <- need(file.path(adni_root,"MRI","MRI_FreeSurfer.csv"))
dict_file <- need(file.path(adni_root,"MRI","MRI_FreeSurfer_Data_Dictionary.csv"))
clinical_file <- need(file.path(adni_root,"Clinical","ADNIMERGE_03Jan2026.csv"))
pgs_file <- need(file.path(project_dir,"results","10A_ADNI_gMAV_PGS_assembly","02_ADNI_gMAV_gAD_PGS_longitudinal_clinical.tsv"))
hip_file <- need(file.path(adni_root,"MRI","MRI_Hippocampus_ASHS.csv"))
wmh_file <- need(file.path(adni_root,"MRI","UCD_WMH_03Jan2026.csv"))
tau_file <- need(file.path(adni_root,"PET","PET_Tau_Berkeley.csv"))
amy_file <- need(file.path(adni_root,"PET","PET_Amyloid_Berkeley.csv"))

log_msg("Step 12C started: common structural indices and multimodal linkage; no association model is fitted.")
fs=fread(fs_file,showProgress=FALSE); dd=fread(dict_file,showProgress=FALSE)
setnames(dd,toupper(names(dd)))
if(!all(c("FLDNAME","TEXT")%in%names(dd)))stop("Unexpected FreeSurfer dictionary format")

# Resolve the 68 Desikan-Killiany thickness-average variables from dictionary text.
map=dd[grepl("Thickness Average of",TEXT,ignore.case=TRUE) & FLDNAME%in%names(fs)]
map[,hemisphere:=fifelse(grepl("Right",TEXT,ignore.case=TRUE),"R",fifelse(grepl("Left",TEXT,ignore.case=TRUE),"L",NA_character_))]
map[,region:=tolower(gsub("[^A-Za-z]","",sub(".*(Right|Left)","",TEXT,ignore.case=TRUE)))]
map=map[hemisphere%in%c("L","R") & nzchar(region)]
if(nrow(map)<60)stop("Resolved fewer than 60 FreeSurfer thickness variables: ",nrow(map),call.=FALSE)
fwrite(map[,.(field=FLDNAME,description=TEXT,hemisphere,region)],file.path(out_dir,"01_ADNI_ST_anatomical_mapping.csv"),bom=TRUE)

for(v in map$FLDNAME)set(fs,j=v,value=num(fs[[v]]))
region_sets=list(
 AD_vulnerable=c("entorhinal","fusiform","inferiortemporal","middletemporal","temporalpole","parahippocampal"),
 frontal=c("caudalmiddlefrontal","rostralmiddlefrontal","superiorfrontal","lateralorbitofrontal","medialorbitofrontal","frontalpole"),
 posterior=c("precuneus","posteriorcingulate","isthmuscingulate","inferiorparietal","superiorparietal","lateraloccipital")
)
fields_for=function(regions)map[region%in%regions,FLDNAME]
rowmean=function(d,v){x=as.matrix(d[,v,with=FALSE]);z=rowMeans(x,na.rm=TRUE);z[rowSums(is.finite(x))==0]=NA_real_;z}
all_t=map$FLDNAME
fs[,global_cortical_thickness:=rowmean(.SD,all_t)]
fs[,ad_vulnerable_index:=rowmean(.SD,fields_for(region_sets$AD_vulnerable))]
fs[,frontal_index:=rowmean(.SD,fields_for(region_sets$frontal))]
fs[,posterior_index:=rowmean(.SD,fields_for(region_sets$posterior))]
fs[,frontoposterior_gradient:=frontal_index-posterior_index]
fs[,spatial_dispersion:=apply(as.matrix(.SD[,all_t,with=FALSE]),1,sd,na.rm=TRUE)]

pairs=merge(map[hemisphere=="L",.(region,L=FLDNAME)],map[hemisphere=="R",.(region,R=FLDNAME)],by="region")
L=as.matrix(fs[,pairs$L,with=FALSE]);R=as.matrix(fs[,pairs$R,with=FALSE])
fs[,mean_absolute_asymmetry:=rowMeans(abs(L-R)/pmax(1e-6,(L+R)/2),na.rm=TRUE)]

# QC is categorical text. Exclude explicit failures only; missing QC is retained and flagged.
qc_text=tolower(trimws(as.character(fs$OVERALLQC)))
fs[,qc_explicit_fail:=grepl("fail|exclude|poor|unusable",qc_text)]
fs[,qc_missing:=is.na(OVERALLQC)|trimws(as.character(OVERALLQC))==""]
fs[,mri_date:=date0(EXAMDATE)]
mri=fs[qc_explicit_fail == FALSE,.(RID=num(RID),PTID,VISCODE,VISCODE2,IMAGEUID,mri_date,FIELD_STRENGTH,FSVER,OVERALLQC,qc_missing,global_cortical_thickness,ad_vulnerable_index,frontal_index,posterior_index,frontoposterior_gradient,spatial_dispersion,mean_absolute_asymmetry)]
saveRDS(mri,file.path(out_dir,"02_ADNI_common_cortical_indices.rds"),compress="xz")

# Longitudinal PGS/clinical file supplies ancestry, PCs, APOE and cognition.
pgs=fread(pgs_file,showProgress=FALSE);pgs[,RID:=num(RID)];pgs[,clinical_date:=date0(EXAMDATE)]
clin=fread(clinical_file,showProgress=FALSE);clin[,RID:=num(RID)];clin[,clinical_date:=date0(EXAMDATE)]

nearest_by_rid <- function(index, candidate, candidate_date, prefix, max_days=365L) {
  a=copy(index);b=copy(candidate)
  a[,row_id__:=.I]; setnames(b,candidate_date,"candidate_date__")
  z=merge(a,b,by="RID",allow.cartesian=TRUE,all.x=TRUE)
  z[,gap__:=abs(as.integer(mri_date-candidate_date__))]
  setorder(z,row_id__,gap__);z=z[,head(.SD,1),by=row_id__]
  z[!is.finite(gap__)|gap__>max_days,(setdiff(names(b),"RID")):=NA]
  candidate_columns <- setdiff(names(b),"RID")
  setnames(z,candidate_columns,paste0(prefix,candidate_columns))
  standardized_date <- paste0(prefix,"candidate_date__")
  if (standardized_date %in% names(z)) setnames(z,standardized_date,paste0(prefix,"date"))
  z[,c("row_id__","gap__"):=NULL];z
}

pgs_keep=intersect(c("RID","clinical_date","age_bl","sex","education","DX","primary_european_ancestry","APOE4","gMAV_z_eur","gMAV_noAPOE_z_eur","gAD_noAPOE_z_eur",paste0("PC",1:10),"MMSE","ADAS13","CDRSB","FAQ","MOCA","mPACCdigit","mPACCtrailsB"),names(pgs))
linked=nearest_by_rid(mri,pgs[,..pgs_keep],"clinical_date","clinical_",365L)

# Link amyloid and tau PET within 18 months; keep raw endpoint fields for 12D freezing.
amy=fread(amy_file,showProgress=FALSE);amy[,RID:=num(RID)];amy[,pet_date:=date0(SCANDATE)]
tau=fread(tau_file,showProgress=FALSE);tau[,RID:=num(RID)];tau[,pet_date:=date0(SCANDATE)]
amy_keep=intersect(c("RID","pet_date","TRACER","AMYLOID_STATUS","CENTILOIDS","SUMMARY_SUVR"),names(amy))
tau_keep=intersect(c("RID","pet_date","TRACER","META_TEMPORAL_SUVR","CTX_ENTORHINAL_SUVR"),names(tau))
linked=nearest_by_rid(linked,amy[,..amy_keep],"pet_date","amyloid_",548L)
linked=nearest_by_rid(linked,tau[,..tau_keep],"pet_date","tau_",548L)

# Link hippocampal volume and WMH within one year.
hip=fread(hip_file,showProgress=FALSE);hip[,RID:=num(RID)];hip[,endpoint_date:=date0(EXAMDATE)];hip[,hippocampus_icv_ratio:=(num(LEFT_HIPP_VOL)+num(RIGHT_HIPP_VOL))/num(ICV)]
wmh=fread(wmh_file,showProgress=FALSE);wmh[,RID:=num(RID)];wmh[,endpoint_date:=date0(EXAMDATE)];wmh[,wmh_brain_ratio:=num(TOTAL_WMH)/num(TOTAL_BRAIN)]
linked=nearest_by_rid(linked,hip[,.(RID,endpoint_date,hippocampus_icv_ratio)],"endpoint_date","hip_",365L)
linked=nearest_by_rid(linked,wmh[,.(RID,endpoint_date,wmh_brain_ratio)],"endpoint_date","wmh_",365L)

saveRDS(linked,file.path(out_dir,"03_ADNI_MRI_PGS_ATVN_linked.rds"),compress="xz")
fwrite(head(linked,1000),file.path(out_dir,"03B_ADNI_linked_preview_first1000.csv"),bom=TRUE)

availability=data.table(
 endpoint=c("MRI","PGS/clinical","amyloid PET","tau PET","hippocampus","WMH"),
 n_rows=c(nrow(linked),sum(!is.na(linked$clinical_date)),sum(!is.na(linked$amyloid_date)),sum(!is.na(linked$tau_date)),sum(is.finite(linked$hip_hippocampus_icv_ratio)),sum(is.finite(linked$wmh_wmh_brain_ratio))),
 n_people=c(uniqueN(linked$RID),uniqueN(linked[!is.na(clinical_date),RID]),uniqueN(linked[!is.na(amyloid_date),RID]),uniqueN(linked[!is.na(tau_date),RID]),uniqueN(linked[is.finite(hip_hippocampus_icv_ratio),RID]),uniqueN(linked[is.finite(wmh_wmh_brain_ratio),RID]))
)
fwrite(availability,file.path(out_dir,"04_multimodal_linkage_flow.csv"),bom=TRUE)

gate=data.table(gate="ADNI_COMMON_STRUCTURAL_INDICES_MULTIMODAL_LINKAGE",resolved_thickness_variables=nrow(map),paired_bilateral_regions=nrow(pairs),MRI_rows=nrow(linked),MRI_people=uniqueN(linked$RID),PGS_people=availability[endpoint=="PGS/clinical",n_people],amyloid_people=availability[endpoint=="amyloid PET",n_people],tau_people=availability[endpoint=="tau PET",n_people],WMH_people=availability[endpoint=="WMH",n_people],raw_MRI_required=FALSE,status=if(nrow(map)>=60&&availability[endpoint=="PGS/clinical",n_people]>=500)"READY_FOR_12D" else "REQUIRES_REPAIR")
fwrite(gate,file.path(out_dir,"05_12C_gate.csv"),bom=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"99_sessionInfo.txt"),useBytes=TRUE)
log_msg("Step 12C completed: ",gate$status)
