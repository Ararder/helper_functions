library(googledrive)
library(googlesheets4)
library(tidyverse)
library(lubridate)
library(data.table)
library(here)
library(fs)
library(glue)
library(broom)


ldsc_munge <- function(infile, N, B = "B") {
  out <- path(path_dir(infile), "ldsc")
  
  if(missing(N)) {
    glue(
      "munge_sumstats.py ",
      "--sumstats {infile} ",
      "--out {out} ",
      "--snp RSID ",
      "--a1 EffectAllele ",
      "--a2 OtherAllele ",
      "--signed-sumstats {B},0 ",
      "--merge-alleles /nas/depts/007/sullilab/shared/bin/ldsc/w_hm3.snplist ",
      "--chunksize 500000"
    )
  } else {
    glue(
      "munge_sumstats.py ",
      "--sumstats {infile} ",
      "--out {out} ",
      "--snp RSID ",
      "--a1 EffectAllele ",
      "--a2 OtherAllele ",
      "--N {N} ",
      "--signed-sumstats B,0 ",
      "--merge-alleles /nas/depts/007/sullilab/shared/bin/ldsc/w_hm3.snplist ",
      "--chunksize 500000"
    )
    
  }
  
}
ldsc_intercept <- function(infile) {
  inf <- path(path_dir(infile), "ldsc.sumstats.gz")
  out <- path(path_dir(infile), "ldsc_h2")
  
  glue(
    "ldsc.py ",
    "--h2 {inf} ",
    "--ref-ld-chr /nas/depts/007/sullilab/shared/bin/ldsc/eur_w_ld_chr/ ",
    "--w-ld-chr /nas/depts/007/sullilab/shared/bin/ldsc/eur_w_ld_chr/ ",
    "--out {out}"
    
  )
}
colname_getter <- function(path) {
  colnames(fread(path, nrows=20))
}
phenotype_getter <- function(path) {
  path_split(path)[[1]][[9]]
}
col_checker <- function(vec) {
  n = "N" %in% vec
  b = "B" %in% vec
  z = "Z" %in% vec
  tibble(N = n, B = b, Z =z)
}
parse_ldsc_log <- function(path){
  pheno = phenotype_getter(path)
  numbers <- read_tsv(path) %>% 
    slice(43:45) %>% 
    pull(1) %>% 
    str_extract(., "\\d..\\d\\d")
  
  tibble(phenotype = pheno, mean_chi2=numbers[1], lambda_gc = numbers[2], max_chi2 = numbers[3])
}
parse_ldsc_h2 <- function(path){
  pheno <- phenotype_getter(path)
  numbers <- read_tsv(path) %>% 
    slice(24,25,26,27,28) %>% 
    pull(1) %>% 
    # str_extract(., " \\d.\\d\\d\\d") %>% 
    str_extract(., " \\d.\\d{1,4}") %>%
    str_remove_all(" ") %>% 
    as.numeric()
  
  tibble(phenotype = pheno,
         obs_h2 = numbers[1], lambda_gc = numbers[2],
         mean_chi2 = numbers[3], intercept = numbers[4],
         ratio = numbers[5]
  )
}
