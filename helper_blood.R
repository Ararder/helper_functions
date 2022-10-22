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
    c("col_CHR: CHR","col_POS: BP","col_SNP: ID","col_BETA: EFFECT",
      "col_EffectAllele: ALT","col_OtherAllele: REF","col_P: PVAL", "col_N: N",
      "col_EAF: ALT_FREQ", "col_SE: SE", "col_Z: Z", "col_INFO: INFO", "col_CaseN: NCAS",
      "col_ControlN: NCON")
  
  possible <- c("CHR_num","BP","ID","EFFECT","ALT","REF",
                "PVAL","N","ALT_FREQ","SE","Z","INFO", "NCAS",
                "NCON")
  
  vec <- vector("character", length = 12)
  for(i in seq_along(possible)) {
    if(any(possible[i] %in% colsname)) {
      vec[i] <- code[i]
    }
  }
  vec[vec != ""]
}

construct_metadata_file <- function(path, model="linear mixed-model") {
  
  
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