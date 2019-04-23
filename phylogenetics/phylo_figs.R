# phylogeny and sequences ----

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

# unrooted phylo
pdf('phylogenetics/fig_phyloUnroot.pdf', width = 4, height = 4)
par(mar = rep(0, 4))
plot(tre, type = 'unrooted', show.tip.label = FALSE, edge.width = 3)
dev.off()


# likelihood maximization ----

library(MASS)
library(viridis)
x <- rnorm(1000)
y <- rnorm(length(x))
z <- kde2d(x, y)$z

colz <- viridis(100)
zFacetCenter <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
zFacetRange<-cut(zFacetCenter, length(colz))

pdf('phylogenetics/fig_lik_surface.pdf', width = 4, height = 4)
par(mar = rep(0, 4))
persp(z, 
      theta = -40, phi = 30, expand = 3, scale = FALSE,
      shade = NA, col = colz[zFacetRange], border = 'grey',
      zlim = c(0.01, 0.06), box=FALSE)
dev.off()

# human, cat, squirrel, bird phylo ----
mb <- read.tree(text = '(bird, (cat, (squirrel, human)));')

pdf('phylogenetics/fig_phylo_mamBird.pdf', width = 4, height = 4)
par(mar = rep(0, 4))
plot(mb, edge.width = 3)
dev.off()

# non-chronogram
mb$edge.length <- c(6, 8.25, 4, 0.25, 0.75, 1.75)
pdf('phylogenetics/fig_phylo_mamBird_noChrono.pdf', width = 4, height = 4)
par(mar = rep(0, 4))
plot(mb, edge.width = 3)
dev.off()


# chronogram
mbCalib <- read.tree('phylogenetics/mb.nwk')

pdf('phylogenetics/fig_phylo_mamBird_chrono.pdf', width = 4, height = 5)
layout(matrix(1:2, nrow = 2), heights = c(4, 2))
par(mar = c(0, 1, 0, 0))
plot(mbCalib, show.tip.label = FALSE, edge.width = 3)
xlim <- par('usr')[2:1]
par(mar = c(5, 1, 0, 0), xpd = NA)
plot(1, type = 'n', xlim = xlim, xaxs = 'i', axes = FALSE, xlab = '')
socorro::paleoAxis(1)
mtext('Millions of years ago', side = 1, line = 4)
dev.off()

# long branch attraction ----

library(ape)
library(phangorn)

set.seed(2)
tre <- rphylo(5, 1, 0.9)
tre <- bind.tree(read.tree(text = sprintf('(a:%s, b:0.2);', 0.2 + max(node.depth.edgelength(tre)))), tre, 2)
plot(drop.tip(tre, 't4'))

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

pdf('phylogenetics/fig_longBranch.pdf', width = 5, height = 4)
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
