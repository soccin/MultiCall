x11 = function (...) grDevices::x11(...,type='cairo')

format_maf<-function(maf) {
    maf %>%
        select(1,5,6,9,11,13,16,HGVSp_Short,t_vaf,t_alt_count,n_vaf,n_depth,CALLERS,FILTERS) %>%
        mutate(Tumor_Sample_Barcode=fct_reorder(Tumor_Sample_Barcode,as.numeric(gsub(".*RR","",Tumor_Sample_Barcode)))) %>%
        arrange(Tumor_Sample_Barcode,desc(t_vaf),Start_Position)
}


require(tidyverse)

callerOrder=c("mutect2", "freebayes", "strelka", "vardict")

mm=fs::dir_ls("int",rec=T,regex=".rda") %>% map(readRDS) %>% bind_rows %>% type_convert

QUAL=mm$vcf_qual
QUAL=as.numeric(ifelse(QUAL==".",-1,QUAL))
QUAL[QUAL<0]=NA
mm$QUAL=QUAL
ps=min(QUAL[QUAL>0],na.rm=T)

mm=mm %>%
    mutate(t_vaf=t_alt_count/t_depth,n_vaf=n_alt_count/n_depth) %>%
    mutate(FILTER=case_when(
                    CALLER=="freebayes" & QUAL>15000 & t_vaf>5*n_vaf ~ "PASS",
                    T ~ FILTER
                )
    ) %>%
    mutate(CALLER=factor(CALLER,levels=callerOrder)) %>%
    arrange(ETAG,CALLER,Tumor_Sample_Barcode)

val1Cols=colnames(mm)[setdiff(which(apply(mm,2,\(x)len(unique(x)))==1),1:20)]

m1=mm %>%
    select(-any_of(val1Cols)) %>%
    mutate(STATUS=ifelse(is.na(STATUS),"",STATUS)) %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(
        N.CALL=n(),N.PASS=sum(FILTER=="PASS"),
        CALLERS=paste(CALLER,collapse=";"),
        FILTERS=paste(FILTER,collapse="|"),
        STATUS=paste(STATUS,collapse="|")
        ) %>%
    ungroup

T_VAF_CUT=0.05

maf0=m1 %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

maf1=m1 %>%
    filter(N.PASS>=1) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

mafHC=m1 %>%
    filter(N.PASS>=2) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

maf10=m1 %>%
    filter(N.PASS>=1 & N.CALL>=2) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

maf2=m1 %>%
    filter(N.CALL>=2) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

maf12=m1 %>%
    filter(N.CALL>=2 | N.PASS>=1) %>%
    filter(t_vaf >= T_VAF_CUT & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))

