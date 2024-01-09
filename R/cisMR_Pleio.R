#' Cis Multivariate Mendelian Randomization (cisMVMR) Analysis Using PRS
#'
#' This function performs a comprehensive Cis Multivariate Mendelian Randomization analysis. It starts by processing GWAS summary data for an outcome and multiple exposures. PRS-CSx is then utilized to estimate the Polygenic Risk Scores (PRS) based on these summary statistics and genotype data from the UK Biobank. Finally, it performs MVMR to investigate the causal effects of the exposures on the outcome.
#'
#' @param outcomefile Path to the outcome file containing GWAS summary statistics for the outcome of interest. It should include columns for SNP identifiers, chromosome index, base pair positions, effect alleles, other alleles, Z-scores, and P-values.
#' @param CHR The chromosome of interest for the analysis.
#' @param BPcenter The central base pair position of the genomic region in the chromosome of interest.
#' @param BPtol A threshold for the furthest SNP from the central base pair, used to define the genomic region of interest.
#' @param prscsxpath Path to PRSCSx directory
#' @param plinkpath Path to Plink software
#' @param conda_env The name of the environment of PRSCSx in conda
#' @param ref_dir The path of refenece panels used in PRSCSx
#' @param bfile The PLINK  bed files used in prediction.
#' @param indMR The IDs of individual in UK Biobank used to perform MVMR analysis.
#' @param pleiotropy A vector of SNP used to account for direct effect.
#' @importFrom data.table fread
#' @importFrom dplyr `%>%` select
#' @importFrom glue glue
#' @importFrom reticulate use_condaenv
#' @return An object of class lm representing the linear model fit
#' @examples
#' MR_PRS(outcomefile = "path/to/outcomefile.csv",
#'        CHR = 1,
#'        BPcenter = 15000000,
#'        BPtol = 1500000,
#'        exposures_list = list("path/to/exposure1.csv", "path/to/exposure2.csv"),
#'        prscsxpath = "/path/to/PRScsx/",
#'        plinkpath = "/path/to/Plink/",
#'        conda_env = "prscsx_env",
#'        indMR = c(1, 2, 3, 4))
#' @export


cisMR_Pleio=function(outcomefile, CHR, BPcenter, BPtol, prscsxpath, plinkpath, conda_env, ref_dir, bfile, indMR, pleiotropy){
  temp_dir <- tempfile()
  dir.create(temp_dir)
  prs_file_dir <- file.path(temp_dir, "prs_file")
  temporary_file_dir <- file.path(temp_dir, "temporary_file")
  dir.create(prs_file_dir)
  dir.create(temporary_file_dir)

  print("Step 1: Processing the data")
  outcome=fread(outcomefile)%>%as.data.frame(.)
  SNP=c(outcome$SNP[which(outcome$CHR==CHR&abs(outcome$BP-BPcenter)<BPtol)])%>%unique(.)
  outcome=outcome[which(outcome$SNP%in%SNP),]%>%dplyr::select(SNP,A1,A2,BETA=Zscore,P,N)
  Noutcome=round(median(outcome$N))
  write.table(outcome,glue("{temporary_file_dir}/outcome.txt"),row.names=F,quote=F,sep="\t")

  print("Step 2: Estimation of PRS using PRSCS")
  use_condaenv(conda_env, required = TRUE)
  setwd(prscsxpath)
  system(glue("python PRScsx.py --ref_dir={ref_dir} --bim_prefix={bfile} --sst_file={temporary_file_dir}/outcome.txt --n_gwas={Noutcome} --pop=EUR --out_dir={prs_file_dir} --out_name=outcome --chrom={CHR}"))

  print("Step 3: Calculation of PRS using PLINK")
  setwd(plinkpath)
  system(glue("./plink --bfile {bfile} --score {prs_file_dir}/outcome_EUR_pst_eff_a1_b0.5_phiauto_chr{CHR}.txt 2 4 6 header sum --out {prs_file_dir}/outcome"))
  write.table(pleiotropy,file = "{prs_file_dir}/snplist.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  system(glue("./plink --bfile {bfile} --extract {prs_file_dir}/snplist.txt --recode --out {prs_file_dir}/mydata_selected_snps"))
  Ple=fread(glue("{prs_file_dir}/mydata_selected_snps.ped"))%>%as.data.frame(.)
  names(Ple)[7:ncol(Ple)]=pleiotropy
  Ple=Ple[which(Ple$ID%in%indMR),]

  print("Step 4: Performing MVMR using predicted scores")
  pred_outcome=fread(glue("{prs_file_dir}/outcome.profile"))
  PRSCS=data.frame(ID=pred_outcome$FID,outcome=pred_outcome$SCORESUM)
  PRSCS=PRSCS[which(PRSCS$ID%in%indMR),]
  PRSCS[,"outcome"]=PRSCS[,"outcome"]/sd(PRSCS[,"outcome"])
  PRSCS=cbind(PRSCS,Ple)
  predictors_string <- paste(pleiotropy, collapse = " + ")
  full_formula_string <- paste0("outcome", " ~ ", predictors_string)
  full_formula <- as.formula(full_formula_string)
  fitegger <- lm(full_formula, data = PRSCS)

  return(A=list(fitegger=fitegger))
}
