suppressPackageStartupMessages({
    require(tidyverse)
})

argv=commandArgs(trailing=T)
manifest=read_csv(argv,progress=F,show_col_types = FALSE)
normals=manifest %>% filter(status==0) %>% select(patient,normal=sample)
tumors=manifest %>% filter(status==1) %>% select(patient,tumor=sample)
pairs=full_join(normals,tumors,by="patient") %>% select(patient,normal,tumor)
cat(format_tsv(pairs,col_names=F))
