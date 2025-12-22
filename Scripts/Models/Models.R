############################################################
# Models for Longitudinal Microbiome Analysis
# Alpha & Beta Diversity lineal models + Multinomial mixed model of cumulative link
# Author: Maria Montagud
############################################################



# -----------------------------
# 1. Alpha Model 
# -----------------------------

shannon_genus <- diversity(GenusAbun[, -212],
                               index = " shannon ") # Shannon - Wienner
# We calculate the average of the index

H2 <- as.factor(substring(GenusAbun[ ,212], 1, 5))
shannon <- tapply(shannon_genus, H2,
                    FUN = mean, na.rm = TRUE )
base$meno0m <- relevel(base$meno0m, ref = " Menopausia ") # reference category
mod1 <- lm(shannon ~ IMC0m + meno0m + oral0m ,
            data = base , na.action = na.omit )


# Weighted model
## Grouping of the BMI variable

levels(base$IMC0m)[levels(base$IMC0m) %in % c("<25","25-29.9")] <- "24-29.9"
## Weighted variable
iden <- substring(dta1$ID, 1, 5) # the time is identified
idn <- table(iden)
pip <- cbind(data.frame(idn, base))
mod1_pon <- lm(shannon ~ IMC0m + meno0m + oral0m,
                  data = base, na.action = na.omit,
                  weights = (pip$Freq))



# -----------------------------
# 2. Beta Model 
# -----------------------------


abund_table_norm <- decostand(physeq2@otu_table, "normalize")
BCGenus <- vegdist (abund_table_norm, "bray")
base_be$ meno0m <- relevel (base_be$ meno0m , ref = " Menopausia ")
mod4_be2 <- lm(BetaIntraM ~ alco6m + meno0m + X1.rel,
                  data = base_be)

# -----------------------------
# 3. Taxonomy Compositional Model 
# -----------------------------


abund_table_norm <- decostand(physeq2@otu_table, "normalize")
bc_dist <- vegdist(abund_table_norm, method = "bray")
## Cluster
cluster_b <- hclust(bc_dist, method = 'ward.D2')
plot(cluster_b, xlab = " Distancia Bray - Curtis ", cex = 0.5)
rect.hclust(cluster_b,
                k=3, border = "red")

## ACoP

PCoA <- cmdscale(bc_dist, eig = TRUE, k = 3) # applies to 3 dimensions

## It is plotted
plot(PCoA$points[ ,1], PCoA$points[ ,2], ylim = c( -0.5 , 0.5),
       xlab = paste("PCoA 1 (", explainedvar1, " %)", sep = ""),
       ylab = paste("PCoA 2 (", explainedvar2, " %)", sep = ""),
       pch = 5, cex = 1.0 , type = "n", cex.lab = 1.0 , cex.axis = 1.2 ,
       axes = FALSE)
# The axis are included
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2 , las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2 , las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# The points and names of the genes are included

points(PCoA$points[ ,1], PCoA$points[ ,2],
         pch = 19, cex = 3, bg = "blue", col = "blue")
text(PCoA$points[ ,1], PCoA$points[ ,2],
       labels = row.names(PCoA$points), cex = 0.56)
plot(fit, p.max = 0.05 , col = "red")

## We will look at each subject in each cluster.
# We transform the grouping into a factor.

CST <- as.factor(cutree(cluster_b, k = 3))
levels(CST) <- c("1", "2", "3")

## Then a database is created that allows us to know which ID belongs to each cluster and which gen.

cla <- data.frame(CST = CST, id = names(CST))
uni <- cbind(data.frame(GenusAbun, cla))
cat1 <- uni[CST ==1 ,]
DomCST1 <- predominantCST(cat1[,-c (212:214)])
cat2 <- uni[CST ==2 ,]
DomCST2 <- predominantCST(cat2[,-c (212:214)])
cat3 <- uni[CST ==3 ,]
DomCST3 <- predominantCST(cat3[,-c (212:214)])

## We create a variable

ide2 <- data.frame(DomCST2, ID = cat2$id)
data_comp <- left_join(dta1, ide2, by = c("ID"="ID"))
data_comp$DomCST2[is.na(data_comp$DomCST2)] <- 1
data_comp$DomCST2[data_comp$DomCST2!= 1] <- 3 # 3 refers to anaerobic bacteria.

ins3 <- cat3$id
idi3 <- data_comp$ID %in% ins3
data_comp$DomCST2[idi3] <- 2 # 2 refers to l.iners, which are lactobacillus iners.
data_comp$DomCST <- data_comp$DomCST2
data_comp$tiempo <- as.numeric(substring(data_comp$time, 7, 8))
# Model
data_comp$meno <- relevel(data_comp$meno, ref = "Menopausia")
fmm1.time <- clmm(DomCST ~ meno + tiempo + (1|id),
                      threshold = "symmetric",
                      data = data_comp, Hess = TRUE)





