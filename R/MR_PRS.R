#' Multivariate Mendelian Randomization (MVMR) Analysis Using PRS
#'
#' This function performs a comprehensive Multivariate Mendelian Randomization analysis. It starts by processing GWAS summary data for an outcome and multiple exposures. PRS-CSx is then utilized to estimate the Polygenic Risk Scores (PRS) based on these summary statistics and genotype data from the UK Biobank. Finally, it performs MVMR to investigate the causal effects of the exposures on the outcome.
#'
#' @param outcomefile Path to the outcome file containing GWAS summary statistics for the outcome of interest. It should include columns for SNP identifiers, chromosome index, base pair positions, effect alleles, other alleles, Z-scores, and P-values.
#' @param CHR The chromosome of interest for the analysis.
#' @param BPcenter The central base pair position of the genomic region in the chromosome of interest.
#' @param BPtol A threshold for the furthest SNP from the central base pair, used to define the genomic region of interest.
#' @param eQTL_list A list of paths to exposure files. Exposure file corresponding to each path should have the same column as outcome file.
#' @param prscsxpath Path to PRScsx directory
#' @param plinkpath Path to Plink software
#' @param conda_env The name of the environment of PRSCSX in conda
#' @param indMR The IDs of individual in UK Biobank used to perform MVMR analysis.
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


MR_PRS=function(outcomefile, CHR, BPcenter, BPtol, eQTL_list, prscsxpath, plinkpath, conda_env, indMR){

print("Step 1: Processing the data")
NAM=names(eQTL_list)
GENE <- sapply(strsplit(NAM, "_", fixed = TRUE), `[`, 1)
eQTL_data=list()
NeQTL=c(1:length(NAM))
for(i in 1:length(NAM)){
A=fread(eQTL_list[[i]])
A=A[which(A$GeneSymbol==GENE[i]),]
A=A%>%dplyr::select(SNP,A1,A2,BETA=Zscore,P,N)
eQTL_data[[i]]=A
NeQTL[i]=median(A$N)
}

outcome=fread(outcomefile)%>%as.data.frame(.)
SNP=c(outcome$SNP[which(outcome$CHR==CHR&abs(outcome$BP-BPcenter)<BPtol)])%>%unique(.)
outcome=outcome[which(outcome$SNP%in%SNP),]%>%dplyr::select(SNP,A1,A2,BETA=Zscore,P,N)
Noutcome=median(outcome$N)

write.table(outcome,glue("~/MR_PRS/temporary_file/outcome.txt"),row.names=F,quote=F,sep="\t")
for(i in 1:length(NAM)){
write.table(eQTL_data[[i]],glue("~/MR_PRS/temporary_file/{NAM[i]}.txt"),row.names=F,quote=F,sep="\t")
}

print("Step 2: Estimation of PRS using PRSCS")
use_condaenv(conda_env, required = TRUE)
setwd(prscsxpath)
system(glue("python PRScsx.py --ref_dir=/mnt/rstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/ukbb --bim_prefix=/mnt/vstor/SOM_EPBI_XXZ10/njl96/data/ukbb/hapmap3/CHR{CHR} --sst_file=/home/yxy1234/MR_PRS/temporary_file/outcome.txt --n_gwas={Noutcome} --pop=EUR --out_dir=/home/yxy1234/MR_PRS/prs_file --out_name=outcome --chrom={CHR}"))
for(i in 1:length(NAM)){
system(glue("python PRScsx.py --ref_dir=/mnt/rstor/SOM_EPBI_XXZ10/yxy1234/LDpanel/ukbb --bim_prefix=/mnt/vstor/SOM_EPBI_XXZ10/njl96/data/ukbb/hapmap3/CHR{CHR} --sst_file=/home/yxy1234/MR_PRS/temporary_file/{NAM[i]}.txt --n_gwas={NeQTL[i]} --pop=EUR --out_dir=/home/yxy1234/MR_PRS/prs_file --out_name={NAM[i]} --chrom={CHR}"))
}

print("Step 3: Calculation of PRS using PLINK")
setwd(plinkpath)
system(glue("./plink --bfile /mnt/rstor/SOM_EPBI_XXZ10/njl96/data/ukbb/hapmap3/CHR{CHR} --score /home/yxy1234/MR_PRS/prs_file/outcome_EUR_pst_eff_a1_b0.5_phiauto_chr{CHR}.txt 2 4 6 header sum --out /home/yxy1234/MR_PRS/prs_file/outcome"))
for(i in 1:length(NAM)){
system(glue("./plink --bfile /mnt/rstor/SOM_EPBI_XXZ10/njl96/data/ukbb/hapmap3/CHR{CHR} --score /home/yxy1234/MR_PRS/prs_file/{NAM[i]}_EUR_pst_eff_a1_b0.5_phiauto_chr{CHR}.txt 2 4 6 header sum --out /home/yxy1234/MR_PRS/prs_file/{NAM[i]}"))
}

print("Step 4: Performing MVMR using predicted scores")
pred_outcome=fread("~/MR_PRS/prs_file/outcome.profile")
PRSCS=data.frame(ID=pred_outcome$FID,outcome=pred_outcome$SCORESUM)
for(i in 1:length(NAM)){
PRSCS[[NAM[i]]]=fread(glue("~/MR_PRS/prs_file/{NAM[i]}.profile"))$SCORESUM
}
PRSCS=PRSCS[which(PRSCS$ID%in%indMR),]
predictors_string <- paste(NAM, collapse = " + ")
full_formula_string <- paste0("outcome", " ~ ", predictors_string,"-1")
full_formula <- as.formula(full_formula_string)
fit <- lm(full_formula, data = PRSCS)

return(fit)
}
