# phylogeny and sequences

library(ape)
library(phangorn)

set.seed(2)
tre <- rphylo(5, 1, 0.9)

set.seed(123)
rootS <- sample(c('a', 'c', 'g', 't'), 10, replace = TRUE)
Q <- matrix(c(0.0, 0.2, 1.0, 0.2,
              0.2, 0.0, 0.0, 1.0, 
              1.0, 0.2, 0.2, 0.2, 
              0.2, 1.0, 0.2, 0.0), 
            byrow = TRUE, nrow = 4)
diag(Q) <- - rowSums(Q)
set.seed(123)
sdat <- simSeq(tre, l = length(rootS), rootseq = rootS, type = 'DNA', rate = 0.1, Q = Q)
s <- toupper(as.character(sdat))

pdf('phylogenetics/fig_phyloSeq1.pdf', width = 5, height = 4)

layout(matrix(1:2, nrow = 1), widths = c(1.5, 1))
par(mar = rep(0, 4))

plot(tre, show.tip.label = FALSE, edge.width = 4)
yy <- get("last_plot.phylo", envir = .PlotPhyloEnv)$yy[1:nrow(s)]
cols <- hsv(c(0, 0.5, 0.12, 0.7), c(0.8, 1, 1, 0.8), c(0.8, 0.8, 1, 0.8))
names(cols) <- c('A', 'C', 'G', 'T')

plot(1:nrow(s), type = 'n', axes = FALSE)
for(i in 1:nrow(s)) {
    for(j in 1:ncol(s)) {
        text(1 + strwidth('M', cex = 1.2) * j, yy[i], labels = s[i, j], 
             adj = c(0.5, 0.5), col = cols[s[i, j]], cex = 1.4)
    }
}

dev.off()
