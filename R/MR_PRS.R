#' Multivariate Mendelian Randomization (MVMR) Analysis Using PRS
#'
#' This function performs a comprehensive Multivariate Mendelian Randomization analysis. It starts by processing GWAS summary data for an outcome and multiple exposures. PRS-CSx is then utilized to estimate the Polygenic Risk Scores (PRS) based on these summary statistics and genotype data from the UK Biobank. Finally, it performs MVMR to investigate the causal effects of the exposures on the outcome.
#'
#' @param outcomefile Path to the outcome file containing GWAS summary statistics for the outcome of interest. It should include columns for SNP identifiers, chromosome index, base pair positions, effect alleles, other alleles, Z-scores, and P-values.
#' @param exposure_list A list of paths to exposure files. Exposure file corresponding to each path should have the same column as outcome file.
#' @param prscsxpath Path to PRSCSx directory
#' @param plinkpath Path to Plink software
#' @param conda_env The name of the environment of PRSCSx in conda
#' @param ref_dir The path of refenece panels used in PRSCSx
#' @param bfile A list of PLINK bed files used in prediction.
#' @param indMR The IDs of individual in UK Biobank used to perform MVMR analysis.
#' @param intercept If an intercept is included in the MVMR analysis.
#' @importFrom data.table fread setnames
#' @importFrom dplyr `%>%` select
#' @importFrom glue glue
#' @importFrom reticulate use_condaenv
#' @return An object of class lm representing the linear model fit
#' @export


MR_PRS=function(outcomefile, exposure_list, prscsxpath, plinkpath, conda_env, ref_dir, bfile, indMR, intercept=F){
  temp_dir <- tempfile()
  dir.create(temp_dir)
  prs_file_dir <- file.path(temp_dir, "prs_file")
  temporary_file_dir <- file.path(temp_dir, "temporary_file")
  dir.create(prs_file_dir)
  dir.create(temporary_file_dir)

  use_condaenv(conda_env, required = TRUE)
  NAM=names(exposure_list)
  exposure_data=list()
  for(i in 1:length(NAM)){
    A=fread(exposure_list[[i]])%>%setDT(.)
    A=A%>%dplyr::select(SNP,CHR,A1,A2,BETA=Zscore,P,N)
    A$P=as.numeric(A$P)
    A$BETA=as.numeric(A$BETA)
    exposure_data[[i]]=A
  }
  outcome=fread(outcomefile)%>%setDT(.)
  outcome=outcome%>%dplyr::select(SNP,CHR,A1,A2,BETA=Zscore,P,N)
  outcome$P=as.numeric(outcome$P)
  outcome$BETA=as.numeric(outcome$BETA)

  for(CHR_target  in 1:22){
    print(glue("CHR{CHR_target }: Processing the data"))
    outcome1=outcome[CHR==CHR_target]
    Noutcome=median(outcome1$N)
    outcome1$CHR=outcome1$N=NULL
    write.table(outcome1,glue("{temporary_file_dir}/outcome.txt"),row.names=F,quote=F,sep="\t")
    Nexposure=c(1:length(NAM))
    for(i in 1:length(NAM)){
      A=exposure_data[[i]]
      A=A[CHR==CHR_target]
      Nexposure[i]=median(A$N)
      A$CHR=A$N=NULL
      write.table(A,glue("{temporary_file_dir}/{NAM[i]}.txt"),row.names=F,quote=F,sep="\t")
    }

    print("CHR{CHR_target}: Estimation of PRS using PRSCS")
    setwd(prscsxpath)
    system(glue("python PRScsx.py --ref_dir={ref_dir} --bim_prefix={bfile[CHR_target]} --sst_file={temporary_file_dir}/outcome.txt --n_gwas={Noutcome} --pop=EUR --out_dir={prs_file_dir} --out_name=outcome --chrom={CHR_target}"))
    for(i in 1:length(NAM)){
      system(glue("python PRScsx.py --ref_dir={ref_dir} --bim_prefix={bfile[CHR_target]} --sst_file={temporary_file_dir}/{NAM[i]}.txt --n_gwas={Nexposure[i]} --pop=EUR --out_dir={prs_file_dir} --out_name={NAM[i]} --chrom={CHR_target}"))
    }

    print("Step 3: Calculation of PRS using PLINK")
    setwd(plinkpath)
    system(glue("./plink --bfile {bfile[CHR_target]} --score {prs_file_dir}/outcome_EUR_pst_eff_a1_b0.5_phiauto_chr{CHR_target}.txt 2 4 6 header sum --out {prs_file_dir}/outcome_{CHR_target}"))
    for(i in 1:length(NAM)){
      system(glue("./plink --bfile {bfile[CHR_target]} --score {prs_file_dir}/{NAM[i]}_EUR_pst_eff_a1_b0.5_phiauto_chr{CHR_target}.txt 2 4 6 header sum --out {prs_file_dir}/{NAM[i]}_{CHR_target}"))
    }

  }

  print("Step 4: Performing MVMR using predicted scores")
  PRSCS=combine_chrwise_profiles(prs_file_dir=prs_file_dir,NAM=c("outcome",NAM),num_chr=22)
  PRSCS=PRSCS[which(PRSCS$ID%in%indMR),]
  if(intercept==T){
    PRSCS[,"outcome"]=PRSCS[,"outcome"]/sd(PRSCS[,"outcome"])
    for(i in 1:length(NAM)){
      PRSCS[,NAM[i]]=PRSCS[,NAM[i]]/sd(PRSCS[,NAM[i]])
    }
    predictors_string <- paste(NAM, collapse = " + ")
    full_formula_string <- paste0("outcome", " ~ ", predictors_string)
    full_formula <- as.formula(full_formula_string)
    fitjoint <- lm(full_formula, data = PRSCS)
    sumdata=as.data.frame(summary(fitjoint)$coefficient[-1,])
    A=matrix(0,length(NAM),2)
    for(i in 1:length(NAM)){
      full_formula_string <- paste0("outcome", " ~ ", NAM[i])
      full_formula <- as.formula(full_formula_string)
      fit=lm(full_formula,data=PRSCS)
      A[i,]=c(summary(fit)$coefficient[2,1:2])
    }
    sumdata$cor=A[,1];sumdata$corse=A[,2]
    sumdata$pratt=sumdata[,1]*sumdata[,"cor"]
    sumdata$prattse=sqrt(sumdata[,1]^2*sumdata[,"corse"]^2+sumdata[,"cor"]^2*sumdata[,2]^2+(1-summary(fitjoint)$r.squared)/nrow(outcome)*sumdata[,"pratt"])
  }

  if(intercept==F){
    PRSCS[,"outcome"]=PRSCS[,"outcome"]/sqrt(sum(PRSCS[,"outcome"]^2))
    for(i in 1:length(NAM)){
      PRSCS[,NAM[i]]=PRSCS[,NAM[i]]/sqrt(sum(PRSCS[,NAM[i]])^2)
    }
    predictors_string <- paste(NAM, collapse = " + ")
    full_formula_string <- paste0("outcome", " ~ ", predictors_string,"-1")
    full_formula <- as.formula(full_formula_string)
    fitjoint <- lm(full_formula, data = PRSCS)
    sumdata=as.data.frame(summary(fitjoint)$coefficient)
    A=matrix(0,length(NAM),2)
    for(i in 1:length(NAM)){
      full_formula_string <- paste0("outcome", " ~ ", NAM[i],"-1")
      full_formula <- as.formula(full_formula_string)
      fit=lm(full_formula,data=PRSCS)
      A[i,]=c(summary(fit)$coefficient[1:2])
    }
    sumdata$cor=A[,1];sumdata$corse=A[,2]
    sumdata$pratt=sumdata[,1]*sumdata[,"cor"]
    sumdata$prattse=sqrt(sumdata[,1]^2*sumdata[,"corse"]^2+sumdata[,"cor"]^2*sumdata[,2]^2+(1-summary(fitjoint)$r.squared)/nrow(outcome)*sumdata[,"pratt"])
  }
  unlink(temp_dir, recursive = TRUE)
  fitjoint$summarydata=sumdata
  return(fitjoint)
}
