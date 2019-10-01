gset120340 <- readRDS("gset120340_exprs_df_t.rds")
gset21138 <- readRDS("gset21138_exprs_df_t.rds")
gset27383 <- readRDS("gset27383_exprs_df_t.rds")
gset98793 <- readRDS("gset98793_exprs_df_t.rds")
gset92538 <- readRDS("gset92538_exprs_df_t.rds")

gset120340_avg = t(apply(gset120340, 1, function(x) tapply(x, colnames(gset120340), mean)))
gset21138_avg = t(apply(gset21138, 1, function(x) tapply(x, colnames(gset21138), mean)))
gset27383_avg = t(apply(gset27383, 1, function(x) tapply(x, colnames(gset27383), mean)))
gset98793_avg = t(apply(gset98793, 1, function(x) tapply(x, colnames(gset98793), mean)))
gset92538_avg = t(apply(gset92538, 1, function(x) tapply(x, colnames(gset92538), mean)))


intersecting_columns = Reduce(intersect, list(colnames(gset120340_avg),colnames(gset21138_avg),
                                              colnames(gset27383_avg), colnames(gset98793_avg),
                                              colnames(gset92538_avg)))

gset120340_intersect = gset120340_avg[ , colnames(gset120340_avg) %in% intersecting_columns]
gset21138_intersect = gset21138_avg[ , colnames(gset21138_avg) %in% intersecting_columns]
gset27383_intersect = gset27383_avg[ , colnames(gset27383_avg) %in% intersecting_columns]
gset98793_intersect = gset98793_avg[ , colnames(gset98793_avg) %in% intersecting_columns]
gset92538_intersect = gset92538_avg[ , colnames(gset92538_avg) %in% intersecting_columns]

combined = rbind(gset120340_intersect, gset21138_intersect, gset27383_intersect, gset98793_intersect, gset92538_intersect)
saveRDS(combined, file="combined_expression_data.rds")

pheno27383 <- readRDS("pheno27383.rds")
pheno27383$gsm = rownames(pheno27383)
pheno21138 <- readRDS("pheno21138.rds")
pheno21138$gsm = rownames(pheno21138)
pheno120340 <- readRDS("pheno120340.rds")
pheno120340$gsm = rownames(pheno120340)
pheno98793 <- readRDS("pheno98793.rds")
pheno98793$gsm = rownames(pheno98793)
pheno92538 <- readRDS("pheno92538.rds")
pheno92538$gsm = rownames(pheno92538)

combined <- readRDS("combined_expression_data.rds")

combined_labels = (bind_rows(pheno27383, pheno21138, pheno120340, pheno98793, pheno92538))

combat_exprs = ComBat(dat = t(combined), batch = combined_labels$study, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
gset_exprs_df = as.data.frame(combat_exprs)
combat_combined = t(combat_exprs)

saveRDS(combat_combined, file="combat_combined.rds")

control = combined_labels[combined_labels$disease == 'Control',]
disease = combined_labels[!(combined_labels$disease == 'Control'),]
control_ids = control$gsm
disease_ids = disease$gsm

saveRDS(combined, file="combined.rds")
saveRDS(combined_control, file="combined_control.rds")
saveRDS(combined_disease, file="combined_disease.rds")

propsfeatures = props(as.data.frame(combat_combined_control), as.data.frame(combat_combined_disease))

#These props features were exported as a csv 