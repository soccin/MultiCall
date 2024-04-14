suppressPackageStartupMessages({
    require(tidyverse)
})

argv=commandArgs(trailing=T)

mafs=map(argv,read_tsv,col_types=cols(.default="c"),comment="#") %>%
  bind_rows %>%
  type_convert %>%
  mutate(Chromosome=factor(Chromosome,levels=c(1:19,"X","Y","MT"))) %>%
  arrange(Chromosome,Start_Position,CALLER) %>%
  mutate(ETAG=paste0(
            Chromosome,":",Start_Position,":",End_Position,":",
            Reference_Allele,":",Tumor_Seq_Allele2
            )
        )

saveRDS(mafs,file.path(dirname(argv[1]),"merge.maf.rda"),compress=T)
write_tsv(mafs,file.path(dirname(argv[1]),"merge.maf.tsv.gz"))
