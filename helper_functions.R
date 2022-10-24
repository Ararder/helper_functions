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
      "--N-cas-col CaseN ",
      "--N-con-col ControlN "
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
ldsc_getter <- function(path) {
  p = phenotype_getter(path)
  numbers <-
    read_tsv(path) %>%
    slice(24,25:28) %>%
    pull(1) %>%
    str_extract(" \\d\\.\\d{1,3}") %>%
    as.numeric()
  
  h2_se <-
    read_tsv(path) %>%
    slice(24) %>%
    str_extract("\\(\\d.\\d{1,6}") %>%
    str_remove("\\(") %>%
    as.numeric()
  tibble(
    phenotype = p,
    obs_h2 = numbers[1],
    obs_h2_se = h2_se,
    lambda_gc = numbers[2],
    mean_chi2 = numbers[3],
    intercept = numbers[4],
    ratio = numbers[5],
    
  )
}

make_cleansumstats_job <- function(dir) {
  metafile <- path(dir, "meta.txt")
  filename <- path_file(dir)
  glue(
    "/nas/depts/007/sullilab/shared/gwas_sumstats/cleansumstats/cleansumstats.sh ",
    "-i {metafile} ",
    "-d /nas/depts/007/sullilab/shared/gwas_sumstats/cleansumstats/out_dbsnp ",
    "-k /nas/depts/007/sullilab/shared/gwas_sumstats/cleansumstats/out_1kgp ",
    "-o {dir}/{filename}",
  )  
}
col_args <- function(colsname) {
  
  code <-
    c("col_CHR: CHR","col_POS: BP","col_SNP: SNP","col_BETA: BETA",
      "col_EffectAllele: A1","col_OtherAllele: A2","col_P: P", "col_N: N",
      "col_EAF: FREQ", "col_SE: SE", "col_Z: Z", "col_INFO: INFO", "col_CaseN: NCAS",
      "col_ControlN: NCON")
  
  possible <- c("CHR","BP","SNP","BETA","A1","A2",
                "P","N","FREQ","SE","Z","INFO", "NCAS",
                "NCON")
  
  vec <- vector("character", length = 12)
  for(i in seq_along(possible)) {
    if(any(possible[i] %in% colsname)) {
      vec[i] <- code[i]
    }
  }
  vec[vec != ""]
}

col_args2 <- function(colsname) {
  
  code <-
    c("col_CHR: CHR","col_POS: BP","col_SNP: ID","col_BETA: EFFECT",
      "col_EffectAllele: REF","col_OtherAllele: ALT","col_P: P", "col_N: N",
      "col_EAF: FREQ", "col_SE: SE", "col_Z: Z", "col_INFO: INFO", "col_CaseN: NCAS",
      "col_ControlN: NCON")
  
  possible <- c("CHR","BP","SNP","BETA","A1","A2",
                "P","N","FREQ","SE","Z","INFO", "NCAS",
                "NCON")
  
  vec <- vector("character", length = 12)
  for(i in seq_along(possible)) {
    if(any(possible[i] %in% colsname)) {
      vec[i] <- code[i]
    }
  }
  vec[vec != ""]
}

construct_metadata_file <- function(path, model="logistic") {
  
  
  test_folder <- "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"
  name <- path_ext_remove(path_ext_remove(path_ext_remove(fs::path_file(path))))
  
  # create dir and copy over raw sumstat
  dir <- dir_create(path(test_folder, name))
  file_copy(path, dir)
  
  slurm_header <- c(
    "#!/bin/bash",
    "#SBATCH --ntasks=1",
    "#SBATCH --cpus-per-task 6",
    "#SBATCH --time=3:00:00",
    "#SBATCH --mem=40g",
    glue("#SBATCH --output={dir}/slurm-%j.out")
  )
  # create the header
  header <- c(
    "cleansumstats_metafile_kind: minimal",
    glue("path_sumStats: {fs::path_file(path)}"),
    glue("stats_Model: {model}")
  )
  sumstat_colnames <- colnames(fread(path, nrows=0))
  cols <- col_args(sumstat_colnames)
  
  meta_out <- paste0(dir, "/", "meta.txt")
  writeLines(c(header, cols), meta_out)
  
  sbatch_out <- paste0(dir, "/run.sh")
  sbatch_job <- make_cleansumstats_job(dir)
  writeLines(c(slurm_header,sbatch_job), sbatch_out)
  system(glue("sbatch {sbatch_out}"))
  
}


standard_checks <- function(path) {
  name <- path_file(path_dir(path))
  df <- fread(path)
  checks <- list()
  if("Z" %in% colnames(df)) {
    checks[["mean_z"]] = mean(df$Z)
    checks[["median_z"]] =  median(df$Z)
    checks[["zero_z"]] =  any(df$Z == 0)
    
  }
  
  if("B" %in% colnames(df)) {
    checks[["mean_b"]] = mean(df$B)
    checks[["median_b"]] =  median(df$B)
    checks[["zero_b"]] =  any(df$B == 0)
  }
  
  if("EAF" %in% colnames(df)) {
    checks[["EAF_range"]] = min(df$EAF) > 0 & max(df$EAF) < 1
  }
  
  if("N" %in% colnames(df)) {
    checks[["unique_n_values"]] = length(unique(df$N))
  }
  
  checks[["p_range"]] <- min(df$P) >= 0 & max(df$P) <= 1
  checks[["n_min_p"]] <- nrow(filter(df, P ==  min(P)))
  
  
  map2(checks, names(checks), ~tibble({{ .y }} := .x)) %>% 
    reduce(bind_cols) %>% 
    add_column(phenotype = name)
  
}
