# Load required libraries
# install.packages(c("pROC", "dplyr"))
setwd('~/Desktop/VGP_10_diseases/results/')
library(pROC)
library(dplyr)

# Function to calculate delta AUC
calculate_delta_auc <- function(data, y_col, covariate_cols, prs_col) {
  # Ensure y is a factor
  data[[y_col]] <- as.factor(data[[y_col]])

  # Null model formula and PRS model formula
  null_formula <- as.formula(paste(y_col, "~", paste(covariate_cols,
                                                     collapse = " + ")))
  prs_formula <- as.formula(paste(y_col, "~ scale(", prs_col, ") + ", paste(c(covariate_cols),
                                                    collapse = " + ")))

  # Fit null model and PRS model
  null_model <- glm(null_formula, data = data, family = binomial())
  prs_model <- glm(prs_formula, data = data, family = binomial())

  # Calculate AUC for null model
  null_probs <- predict(null_model, type = "response")
  null_roc <- roc(data[[y_col]], null_probs)
  null_auc <- auc(null_roc)

  # Calculate AUC for PRS model
  prs_probs <- predict(prs_model, type = "response")
  prs_roc <- roc(data[[y_col]], prs_probs)
  prs_auc <- auc(prs_roc)

  # Calculate delta AUC
  delta_auc <- prs_auc - null_auc

  # Return results
  return(data.frame(
    pgs = prs_col,
    null_auc = null_auc,
    prs_auc = prs_auc,
    delta_auc = delta_auc
  ))
}

multiple_calculate_delta_auc <- function(data, y_col, covariate_cols, prs_cols) {
  res_full = NULL
  for (prs_col in prs_cols) {
    tryCatch({
      res_temp = calculate_delta_auc(data=data, y_col = y_col,
                                     covariate_cols = covariate_cols,
                                     prs_col = prs_col)
      res_full = rbind(res_full, res_temp)
    }, error = function(e) {
      print(e)
    })
  }

  return(res_full[order(res_full$delta_auc, decreasing=T),])
}

# covariate
cov = read.delim('./target_data/covar_list.tsv', sep='\t')
covar_list = c(paste0("PC", 1:10), "Sex", "Age")
diseases = c('breast_cancer', 'colorectal_cancer', 'gastric_cancer', 'pd', 'ckd',
             'osteoporosis', 'cad', 'hyperlipidemia', 'osteoarthritis')

# run for all
for (d in diseases) {
  data = read.delim(sprintf('./results/%s/PRS.score.%s.samples.age.txt', d, d))
  data = merge(data, cov)
  pgs_list = c(
    grep("score", names(data), value = TRUE),
    grep("PGS", names(data), value = TRUE),
    grep("prscs", names(data), value = TRUE),
    grep("PRSice", names(data), value = TRUE)
  )
  res = multiple_calculate_delta_auc(data=data, y_col = "PHENO",
                                     covariate_cols = covar_list,
                                     prs_cols = pgs_list)
  print(head(res, 10))
  write.table(res, sprintf('./results/%s/PRS.delta.AUC.age.txt', d), sep='\t', row.names = FALSE)
}

# run without age
covar_list = c(paste0("PC", 1:10), "Sex")
# run for all
for (d in diseases) {
  data = read.delim(sprintf('./results/%s/PRS.score.%s.samples.txt', d, d))
  data = merge(data, cov)
  pgs_list = c(
    grep("score", names(data), value = TRUE),
    grep("PGS", names(data), value = TRUE),
    grep("prscs", names(data), value = TRUE),
    grep("PRSice", names(data), value = TRUE)
  )
  res = multiple_calculate_delta_auc(data=data, y_col = "PHENO",
                                     covariate_cols = covar_list,
                                     prs_cols = pgs_list)
  print(head(res, 10))
  write.table(res, sprintf('./results/%s/PRS.delta.AUC.no.age.txt', d), sep='\t', row.names = FALSE)
}