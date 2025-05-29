setwd('~/Desktop/VGP_10_diseases/results/')
rr = function(x,digit=10) return(round(x,digit))

eval_single_PRS = function(data_df, pheno = "trait", prs_name, covar_list, isbinary=F, liabilityR2=F, alpha=0.05, regression_output=F) { # nolint

  colnames(data_df)[which(colnames(data_df)==pheno)] = "trait"
  
  if (isbinary & !liabilityR2) {
    data_df$trait = as.numeric(data_df$trait)
    formula = as.formula(paste0(pheno, " ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
    model_full = glm(formula, data=data_df, family="binomial")
    r_full = suppressWarnings(logLik(model_full, REML=FALSE))[1]
    
    formula = as.formula(paste0(pheno, " ~ ", paste0(covar_list, collapse="+")))
    model_null = glm(formula, data=data_df, family="binomial")
    r_null = suppressWarnings(logLik(model_null, REML=FALSE))[1]
    
    m = r_full
    n = r_null
    N = nobs(model_full)
    cs = 1 - exp(-2/N * (m - n))
    nk = cs/(1 - exp(2/N * n))
    R2 = nk
    
  } else {
    data_df$trait = as.numeric(data_df$trait)
    
    formula = as.formula(paste0("trait ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
    model_full = lm(formula, data=data_df)
    r_full = summary(model_full)$r.squared
    
    formula = as.formula(paste0("trait ~ ", paste0(covar_list, collapse="+")))
    model_null = lm(formula, data=data_df)
    r_null = summary(model_null)$r.squared
    
    N = nobs(model_full)
    R2 = r_full - r_null		
  }
  if (isbinary & liabilityR2) {
    N = nrow(data_df)
    K = mean(data_df$trait, na.rm=T)
    R2 = R2 * K * (1-K) / (dnorm(qnorm(p=1-K, lower.tail=T))^2)
  }
  
  NCP = N * R2 / (1-R2)
  power = 1-pnorm(qnorm(1-alpha/2)-NCP^0.5) + pnorm(qnorm(alpha/2)-NCP^0.5)
  
  vv = (4*R2*(1-R2)^2 *(N-2)^2) / ((N^2-1)*(N+3))
  
  se = sqrt(vv)
  lower_r2 = R2 - 1.97*se
  upper_r2 = R2 + 1.97*se
  pval = pchisq((R2/se)^2, df=1, lower.tail=F)
  
  r2_out = paste0(rr(R2,3), " (", rr(lower_r2,3), "-", rr(upper_r2,3), ")")
  
  if (regression_output) {
    return(data.frame(
      pgs=prs_name, 
      R2=R2, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval=pval, power=power, 
      coef_regression=coef(summary(model_full))[2,1], 
      se_regression=coef(summary(model_full))[2,2], 
      pval_regression=coef(summary(model_full))[2,4]))
  } else {
    return(data.frame(pgs=prs_name, R2=R2, R2_out=r2_out, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval=pval, power=power))
  }
}

eval_multiple_PRS = function(data_df, pgs_list, covar_list, liabilityR2=F, alpha=0.05, isbinary=F, ncores=1, regression_output=F, pheno="trait") {
  
  colnames(data_df)[which(colnames(data_df)==pheno)] = "trait"	
  
  writeLines("Eval all PGS")
  if (isbinary) {
    writeLines("Case - control numbers:")
    print(table(data_df$trait))
  }
  
  idx = which(!pgs_list %in% colnames(data_df))
  if (length(idx)>0) {
    writeLines(paste0(length(idx), " scores not found or sum equal to 0 in data"))
    pgs_list = pgs_list[-idx]
  }
  
  pred_acc_test = NULL
  for (prs_i in 1:length(pgs_list)) {
    
    if (prs_i %% 100 == 0) {
      writeLines(paste0("Evaluated ",  prs_i, " scores"))
    }
    
    prs_name = pgs_list[prs_i]
    pred_acc_test_tmp = eval_single_PRS(data_df=data_df, pheno="trait", prs_name=prs_name, covar_list=covar_list, isbinary=isbinary, liabilityR2=liabilityR2, alpha=alpha, regression_output=regression_output)
    pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp)
  }
  
  pred_acc_test = pred_acc_test[order(pred_acc_test$R2, decreasing=T),]
  
  return(pred_acc_test)
}

# covariate
cov = read.delim('./target_data/covar_list.tsv', sep='')
diseases = c('breast_cancer', 'colorectal_cancer', 'gastric_cancer', 'pd', 'ckd',
              'osteoporosis', 'cad', 'hyperlipidemia', 'osteoarthritis')

# Evaluate PRS for each disease with covar_list: PC1-10, Age and Sex
for (d in diseases) {
  prs_df = read.delim(sprintf('./results/%s/PRS.score.%s.samples.age.txt', d, d))
  prs_df = merge(prs_df, cov)
  pgs_list = c(
    grep("score", names(prs_df), value = TRUE),
    grep("PGS", names(prs_df), value = TRUE),
    grep("prscs", names(prs_df), value = TRUE),
    grep("PRSice", names(prs_df), value = TRUE)
  )
  covar_list = c(paste0("PC", 1:10), "Age", "Sex")
  
  prs_acc = eval_multiple_PRS(data_df=prs_df, 
                              pgs_list=pgs_list, 
                              covar_list=covar_list, 
                              liabilityR2=T, alpha=0.05, 
                              isbinary=T, regression_output=T, pheno="PHENO")
  print(c("Disease:", d))
  print(head(prs_acc, 5))
  write.table(prs_acc, sprintf('./%s/PRS.acc.age.txt', d), sep='\t', row.names = FALSE)
}

# Evaluate PRS for each disease with covar_list: PC1-10 and Sex
for (d in diseases) {
  prs_df = read.delim(sprintf('./results/%s/PRS.score.%s.samples.no.age.txt', d, d))
  prs_df = merge(prs_df, cov)
  pgs_list = c(
    grep("score", names(prs_df), value = TRUE),
    grep("PGS", names(prs_df), value = TRUE),
    grep("prscs", names(prs_df), value = TRUE),
    grep("PRSice", names(prs_df), value = TRUE),
    grep("SBayesRC", names(prs_df), value = TRUE)
  )
  covar_list = c(paste0("PC", 1:10), "Sex")
  
  prs_acc = eval_multiple_PRS(data_df=prs_df, 
                              pgs_list=pgs_list, 
                              covar_list=covar_list, 
                              liabilityR2=T, alpha=0.05, 
                              isbinary=T, regression_output=T, pheno="PHENO")
  write.table(prs_acc, sprintf('./%s/PRS.acc.no.age.txt', d), sep='\t', row.names = FALSE)
}




