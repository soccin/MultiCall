x11 = function (...) grDevices::x11(...,type='cairo')

format_maf<-function(maf) {
    maf %>%
        select(1,5,6,9,11,13,16,HGVSp_Short,t_vaf,t_alt_count,n_vaf,n_depth,CALLERS,FILTERS) %>%
        mutate(Tumor_Sample_Barcode=fct_reorder(Tumor_Sample_Barcode,as.numeric(gsub(".*RR","",Tumor_Sample_Barcode)))) %>%
        arrange(Tumor_Sample_Barcode,desc(t_vaf),Start_Position)
}


require(tidyverse)

callerOrder=c("mutect2", "freebayes", "strelka", "vardict")

mm=fs::dir_ls("int",rec=T,regex=".rda") %>%
    map(readRDS) %>%
    bind_rows %>%
    type_convert(locale=locale(grouping_mark=""))

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


#
# All Events
#

maf0=m1 %>% distinct(ETAG,Tumor_Sample_Barcode,.keep_all=T) %>%
    arrange(Chromosome,Start_Position,Tumor_Sample_Barcode)

#
# Remove silent events and population events (dbSNP)
#

maf0f=maf0 %>%
    filter(!is.na(HGVSp_Short) & !grepl("=$",HGVSp_Short))


maf1=maf0f %>%
    filter(t_vaf >= 0.05 & t_alt_count >= 8 & t_depth >= 20 & t_vaf > 5*n_vaf) %>%
    filter(N.PASS>=1)

maf3=maf1 %>% filter(N.PASS>=2)

tbl3=maf3 %>%
    select(Sample=Tumor_Sample_Barcode,Gene=Hugo_Symbol,
        Type=Variant_Classification,dbSNP=dbSNP_RS,Alteration=HGVSp_Short,
        MAF=t_vaf,t_depth,t_alt_count,n_depth,n_alt_count,CALLERS,FILTERS)

class(tbl3$MAF)="percentage"

library(openxlsx)
wb=createWorkbook()
addWorksheet(wb,sheetName="HCEvents")
addWorksheet(wb,sheetName="HighSensMAF")
writeDataTable(wb,sheet=1,tbl3,tableStyle="none",withFilter=F)
setColWidths(wb,sheet=1,cols=1:ncol(tbl3),widths="auto")
writeDataTable(wb,sheet=2,maf1,tableStyle="none",withFilter=F)

projNo=grep("^Proj",strsplit(getwd(),"/")[[1]],value=T)
if(len(projNo)==0) {
    projNo=""
}

rFile=cc(projNo,"mutationReport","MusVarV1.xlsx")
rDir="post/reports"
fs::dir_create(rDir)

saveWorkbook(wb,file.path(rDir,rFile),overwrite=T)

mFile=cc(projNo,"MergedUnFiltered","MAF.txt")
write_tsv(maf0,file.path(rDir,mFile))

