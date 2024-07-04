suppressPackageStartupMessages({
    require(tidyverse)
})

argv=commandArgs(trailing=T)

manifestCols=cols(
  patient = col_character(),
  sample = col_character(),
  status = col_double(),
  lane = col_character(),
  fastq_1 = col_character(),
  fastq_2 = col_character()
)


manifest=read_csv(argv[1],progress=F,col_types = manifestCols) %>%
    select(patient,sample,status) %>%
    distinct

if(len(argv)==1) {

    normals=manifest %>% filter(status==0) %>% select(patient,normal=sample)
    tumors=manifest %>% filter(status==1) %>% select(patient,tumor=sample)
    pairs=full_join(normals,tumors,by="patient") %>% select(patient,normal,tumor)
    cat(format_tsv(pairs,col_names=F))

} else if(len(argv)==2) {

    bams=fs::dir_ls(argv[2],recur=T,regex=".recal.cram$")
    bams=tibble(sample=gsub(".recal.(cram)","",basename(bams)),bam=bams)

    manifest=left_join(manifest,bams,by="sample")

    normals=manifest %>% filter(status==0) %>% select(patient,normal=sample,normalBam=bam)
    tumors=manifest %>% filter(status==1) %>% select(patient,tumor=sample,tumorBam=bam)

    pairs=left_join(normals,tumors,by="patient") %>% select(normalBam,tumorBam)

    cat(format_tsv(pairs,col_names=F))

} else {

    cat("\n usage:\n")
    cat("    getSarekPairs.R sarek_input.csv\n")
    cat("    getSarekPairs.R sarak_input.csv path/to/bams\n\n")
    quit()

}
