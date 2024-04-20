require(tidyverse)

mm=fs::dir_ls("int",rec=T,regex=".rda") %>% map(readRDS) %>% bind_rows %>% type_convert

dd=mm %>%
    select(Hugo_Symbol,ETAG,Tumor_Sample_Barcode,CALLER,FILTER,vcf_qual,matches("^[tn]_")) %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(N.CALL=n(),N.PASS=sum(FILTER=="PASS"),CALLS=paste(CALLER,collapse=";")) %>%
    ungroup


qual=as.numeric(dd$vcf_qual[dd$CALLER=="freebayes"])

#dd %>% mutate(PASS=FILTER=="PASS") %>% count(CALLER,PASS) %>% spread(CALLER,n)
tbl=dd %>% mutate(PASS=FILTER=="PASS") %>% count(CALLER,PASS) %>% group_by(CALLER) %>% mutate(PCT=n/sum(n)) %>% select(-n) %>% spread(CALLER,PCT)

qualCut=quantile(qual,c(.975))

callerOrder=c("mutect2", "freebayes", "strelka")

QUAL=mm$vcf_qual
QUAL=as.numeric(ifelse(QUAL==".",-1,QUAL))
QUAL[QUAL<0]=NA
mm$QUAL=QUAL
ps=min(QUAL[QUAL>0],na.rm=T)

m1=mm %>%
    mutate(t_vaf=t_alt_count/t_depth,n_vaf=n_alt_count/n_depth) %>%
    mutate(CALLER=factor(CALLER,levels=callerOrder)) %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(N.CALL=n(),N.PASS=sum(FILTER=="PASS"),CALLS=paste(CALLER,collapse=";")) %>%
    ungroup

x11 = function (...) grDevices::x11(...,type='cairo')

m2=m1 %>%
    mutate(FILTER=case_when(
                    CALLER=="freebayes" & QUAL>1 & t_vaf>5*n_vaf & t_alt_count>=8 ~ "PASS",
                    T ~ FILTER
                )
    ) %>%
    group_by(ETAG,Tumor_Sample_Barcode) %>%
    mutate(N.CALL=n(),N.PASS=sum(FILTER=="PASS"),CALLS=paste(CALLER,collapse=";")) %>%
    ungroup

maf0=m2 %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short)) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T)

mafHC=m2 %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short)) %>%
    filter(N.PASS>=2) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T)

maf10=m2 %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short)) %>%
    filter(N.PASS>=1 & N.CALL>=2) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T)

maf1=m2 %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short)) %>%
    filter(N.CALL>=2) %>%
    distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T)
