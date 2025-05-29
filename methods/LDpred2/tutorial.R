setwd('~/Desktop/VGP_10_diseases/methods/LDpred2/')
library(bigsnpr)

# read bed, bim, fam to generate .rds and .bk file
tryCatch({
  snp_readBed("../../target_data/eas.1kgp.hapmap3.annotate.bed")
  }, error = function(err) {
    print(err)
  })

obj.bigSNP <- snp_attach("../../target_data/eas.1kgp.hapmap3.annotate.rds")
hapmap3 <- readRDS('./map_hm3_plus.rds')

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <- nb_cores()
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".", ncores = NCORES)

# read GWAS
sumstats <- bigreadr::fread2("../../sum_stats_data/cad/gwas.BBJ.CAD/gwas.LDpred2.txt")
sumstats <- sumstats[sumstats$rsid %in% hapmap3$rsid, ]
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE, match.min.prop=0.001)

# Calculate LD score
tmp <- tempfile(tmpdir = "tmp_data")

print("Calculate LD score in chr:")
for (chr in 1:22) {
  print(chr)
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(G, 
                   ind.col = ind.chr2, 
                   size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], 
                   ncores = as.integer(NCORES/4))
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, 
                               sample_size = n_eff, blocks = NULL))
ldsc_h2_est <- ldsc[["h2"]]
print(paste0("ldsc_h2 = ", ldsc_h2_est))

# Perform LDpred2-grid
set.seed(1)
print("Run LDpred2-grid...")
h2_seq <- round(ldsc_h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 0.5, length.out = 15), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = as.integer(NCORES/4))

# Perform LDpred2-auto
set.seed(1)
print("Run LDpred2-auto...")
multi_auto <- snp_ldpred2_auto(
  corr, 
  df_beta, 
  h2_init = ldsc_h2_est,
  vec_p_init = seq_log(5e-5, 0.5, length.out = 30), 
  ncores = as.integer(NCORES/4),
  use_MLE = FALSE,
  allow_jump_sign = FALSE,
  shrink_corr = 0.95)

range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# evaluate PRS
samples = read.delim('../../results/breast_cancer/samples.txt')
prs_list = cbind(target.fam, pred_auto, pred_grid)
colnames(prs_list) = c("IID", "auto", paste0("grid", 1:90))
prs_list = merge(samples, prs_list)

eval_single_PRS = function(data_df, pheno = "trait", prs_name, covar_list, isbinary=F, 
                           liabilityR2=F, alpha=0.05, regression_output=F) {
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

eval_multiple_PRS = function(data_df, pgs_list, covar_list, liabilityR2=F, alpha=0.05, 
                             isbinary=F, ncores=1, regression_output=F, pheno="trait") {
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
    pred_acc_test_tmp = eval_single_PRS(data_df=data_df, pheno="trait", prs_name=prs_name,
                                        covar_list=covar_list, isbinary=isbinary, liabilityR2=liabilityR2, 
                                        alpha=alpha, regression_output=regression_output)
    pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp)
  }
  
  pred_acc_test = pred_acc_test[order(pred_acc_test$R2, decreasing=T),]
  
  return(pred_acc_test)
}

cov = read.delim('../../target_data/target.data.impute.eigenvec', sep=' ')
prsice2 = read.delim('../../results/breast_cancer/PRSice2.best', sep="")[c("IID", "PRS")]
prs_list = merge(prs_list, cov)
prs_list = merge(prs_list, prsice2)
num_cols = c("auto", paste0("grid", 1:90), paste0("PC", 1:10), "PRS")
prs_list[num_cols] = sapply(prs_list[num_cols], as.numeric)

eval_prs = eval_multiple_PRS(data_df = prs_list, pgs_list = c("auto", paste0("grid", 1:90), "PRS"),
                            covar_list = paste0("PC", 1:10), liabilityR2 = F, alpha = 0.05,
                            isbinary = T, pheno = "PHENO", regression_output = T)



