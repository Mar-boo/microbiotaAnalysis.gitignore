############################################################
# Longitudinal Microbiome Analysis — Simulated Data
# Alpha & Beta Diversity + Mixed Models
# Author: Maria Montagud (simulated example)
############################################################
# NOTE: All data are simulated for methodological demonstration


# -----------------------------
# 1. Libraries
# -----------------------------

librerias <- c("tidyverse", "vegan", "MASS", "leaps", "lme4", "lmerTest", "nlme", "ordinal", "car", "ggplot2", "dplyr", "ggpubr")

for (lib in librerias) {
  if (!lib %in% installed.packages()[, "Package"]) {
    install.packages(lib, dependencies = TRUE)
  }
  library(lib, character.only = TRUE)
}


set.seed(123)

# =============================
# 2. Simulate longitudinal metadata
# =============================
n_subjects <- 50
n_time <- 3

metadata <- expand.grid(
  id = factor(1:n_subjects),
  time = factor(c("T1", "T2", "T3"), ordered = TRUE)
)

metadata <- metadata %>%
  mutate(
    edad = round(rnorm(n(), 45, 6), 1),
    imc = round(rnorm(n(), 25, 3), 1),
    meno = factor(sample(c("Menopausia", "Perimenopausia", "No menopausia"),
                         n(), replace = TRUE)),
    alcohol = factor(sample(c(0, 1), n(), replace = TRUE)),
    oral = factor(sample(c(0, 1), n(), replace = TRUE)),
    tiempo = as.numeric(time)
  )


# =============================
# 3. Simulate microbiome abundance table
# =============================
n_taxa <- 25

otu_table <- matrix(
  rgamma(n_subjects * n_time * n_taxa, shape = 1.5, scale = 1),
  nrow = n_subjects * n_time,
  ncol = n_taxa
)

colnames(otu_table) <- paste0("Genus_", 1:n_taxa)
rownames(otu_table) <- paste0(metadata$id, "_", metadata$time)

# Relative abundance
otu_rel <- otu_table / rowSums(otu_table)

# =============================
# 4. Alpha diversity
# =============================
metadata$shannon_genus <- diversity(otu_rel, index = "shannon")

# Exploratory
summary(metadata$shannon_genus)
shapiro.test(metadata$shannon_genus)

# =============================
# 5. Linear Mixed Model (Alpha diversity)
# =============================
lmm_alpha <- lmer(
  shannon_genus ~ meno + imc + alcohol + tiempo + (1 | id),
  data = metadata
)

summary(lmm_alpha)

# Diagnostics
shapiro.test(residuals(lmm_alpha))
qqPlot(residuals(lmm_alpha))
plot(lmm_alpha)


# =============================
# 6. Beta diversity (Bray–Curtis)
# =============================
bray_dist <- vegdist(otu_rel, method = "bray")


# =============================
# 7. Ordination (PCoA)
# =============================
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  Sample = rownames(pcoa_res$points),
  PCoA1 = pcoa_res$points[,1],
  PCoA2 = pcoa_res$points[,2]
)

pcoa_df <- cbind(pcoa_df, metadata)

ggplot(pcoa_df, aes(PCoA1, PCoA2, color = meno)) +
     geom_point(alpha = 0.7) +
     theme_minimal() +
     labs(title = "PCoA – Bray-Curtis")


# =============================
# 7. PERMANOVA
# =============================
adonis_res <- adonis2(
  bray_dist ~ meno + tiempo + imc,
  data = metadata,
  permutations = 999
)

adonis_res

sample_n <- expand.grid(
  id = factor(1:n_subjects),
  time = factor(c("T1", "T2", "T3"), ordered = TRUE)
)

F_value <- adonis_res$F[1]
p_value <- adonis_res$`Pr(>F)`[1]
n_value <- sample_n


ggplot(pcoa_df, aes(PCoA1, PCoA2)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  
  geom_point(
    aes(shape = time, color = time),
    size = 3,
    alpha = 0.8
  ) +
  
  stat_ellipse(
    aes(group = time, color = time),
    type = "norm",
    linetype = 2,
    level = 0.68
  ) +
  
  theme_bw() +
  labs(
    title = "Bray–Curtis distance method",
    subtitle = paste0(
      "PERMANOVA: F = ", round(F_value, 3),
      ", p = ", round(p_value, 3),
      ", n = ", length(n_value)
    ),
    x = "PCoA 1",
    y = "PCoA 2",
    color = "Time",
    shape = "Time"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )


# =============================
# 8. Homogeneity of dispersion
# =============================
disp <- betadisper(bray_dist, metadata$meno)
anova(disp)
permutest(disp)


# =============================
# 10. Simulate CST clusters (multinomial outcome)
# =============================
metadata$DomCST2 <- factor(
  sample(c("Lactobacillus", "L_iners", "Anaerobic"),
         nrow(metadata), replace = TRUE),
  ordered = TRUE
)

# =============================
# 11. Ordinal mixed model (CLMM)
# =============================
clmm_model <- clmm(
  DomCST2 ~ meno + imc + tiempo + (1 | id),
  data = metadata,
  threshold = "symmetric",
  Hess = TRUE
)

summary(clmm_model)
confint(clmm_model)

# =============================
# 12. Model comparison
# =============================
clmm_null <- clmm(
  DomCST2 ~ 1 + (1 | id),
  data = metadata,
  threshold = "symmetric"
)

anova(clmm_null, clmm_model)


# =============================
# 13. Save data
# =============================


save(metadata, lmm_alpha, adonis_res, beta_long, clmm_model, clmm_null, pcoa_df, file = "Micro.RData")


############################################################
# End of script
############################################################






