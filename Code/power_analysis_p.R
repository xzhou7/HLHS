library(Seurat)
library(ggplot2)
library(dplyr)
library(pwr)



DefaultAssay(lhv.Endothelium) <- "RNA"

data <- as.data.frame(lhv.Endothelium@assays$RNA@data)

data[1:10, 1:10]

table(lhv.Endothelium$integrated_snn_res.0.2, lhv.Endothelium$condition)

mean1 <- lapply(data, FUN = mean)
sd1 <- lapply(data, FUN=sd)

mean1 <- as.data.frame(mean1)
sd1 <- as.data.frame(sd1)

table1 <- as.data.frame(cbind(t(mean1), t(sd1)))
metadata <- lhv.Endothelium[[]]
metadata1 <- cbind(metadata,table1)

means <- aggregate(V1 ~  condition + integrated_snn_res.0.2 , metadata1, mean)

sds1 <- metadata1 %>% group_by(group) %>% summarise_at("V2","sd")

metadata1$group <- paste(metadata1$integrated_snn_res.0.2, metadata1$condition)

p <- ggplot(metadata1, aes(x=condition, y=V1)) + geom_boxplot()
p <- p + facet_grid(.~integrated_snn_res.0.2)
p

#determine the singificance value 
p.adjust(2.45E-06, method = "bonferroni",n = 20445)

#determine the cohenD
cohens_mean <- function(mx, sdx, my, sdy, n1, n2) {
    md <- abs(mx - my)
    n1 <- n1 - 1
    n2 <- n2 - 1
    csd <- (n1 * (sdx^2)) + (n2 * (sdy^2))
    csd <- csd/(n1 + n2)
    csd <- sqrt(csd)
    
    cd <- md/csd
    print(cd)
}
#0
d1 <- cohens_mean(0.1582891,0.0174,0.1217646,0.0101,614,1775)

#1
d2 <- cohens_mean(0.1711772,0.0170,0.1283304,0.00933,243,786)

#2
d3 <- cohens_mean(0.1683286,0.0181,0.1418091,0.0121,404,522)

#3
d4 <- cohens_mean(0.1611230,0.0175,0.1270450,0.00983,225,488)

#4
d5 <- cohens_mean(0.1484054,0.0165,0.1356253,0.0130,286,76)

#5
d6 <- cohens_mean(0.1603662,0.0317,0.1204256,0.0442,63,95)

pwr.t.test(n=614, d = 2.947402, sig.level =2.45E-06)
pwr.t.test(n=243, d = 3.692644, sig.level =2.45E-06)
pwr.t.test(n=404, d = 1.766241, sig.level =2.45E-06)
pwr.t.test(n=225, d = 2.671902, sig.level =2.45E-06)
pwr.t.test(n=76, d = 0.8070916, sig.level =2.45E-06)
pwr.t.test(n=63, d = 1.005906, sig.level =2.45E-06)

c(1,1,1,1,0.53,0.75)
