#plot Do
library(stringr)
library(ggplot2)
library(GOplot)

setwd("~/Box/XinZhouFiles/Projects/ZXE7_SingleCell_HLHS/NewAnalysis/GOmap/")
DOUP <- read.csv("./GOUP.csv",header = T)
DOUP$count <- str_count(DOUP$Symbols, ',')

DODOWN <- read.csv("./GODOWN.csv",header = T)
DODOWN$count <- str_count(DODOWN$Symbols, ',')

DOUP$Cluster <- factor(DOUP$Cluster, levels=c(2,0,3,1))
DODOWN$Cluster <- factor(DODOWN$Cluster, levels=c(2,0,3,1))

unique(DOUP$Term)
unique(DODOWN$Term)

DOUP$Term <- factor(DOUP$Term,levels = c("GO:0071559", "GO:0008593", "GO:0001936", "GO:0043534",
                                         "GO:0003158", "GO:0048514", "GO:0030335","GO:0034329"))


DODOWN$Term <- factor(DODOWN$Term, levels = c("GO:0051092", "GO:0038061", "GO:0006096","GO:0001666", 
                                              "GO:0010948", "GO:0002218", "GO:0000715" ))

p.up <- ggplot(DOUP, aes(x=Cluster, y=Term)) + geom_point(aes(size=count, color = -LogP))
p.up <- p.up + theme(legend.position="bottom",legend.box = "vertical")
p.up <- p.up + scale_color_gradient(low="darkgrey", high="red")
p.up
ggsave2("./GOUP.pdf", p.up, width = 3, height = 5, dpi=300)

p.down <- ggplot(DODOWN, aes(x=Cluster, y=Term)) + geom_point(aes(size=count, color = -LogP))
p.down <- p.down + theme(legend.position="bottom",legend.box = "vertical")
p.down <- p.down + scale_color_gradient(low="darkgrey", high="blue")
p.down
ggsave2("./GODown.pdf", p.down, width = 3, height = 5, dpi=300)

p.up + p.down

DO1 <- subset(DODOWN, Cluster==1)


GOChord(DO1,space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
