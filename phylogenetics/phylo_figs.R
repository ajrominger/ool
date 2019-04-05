# phylogeny and sequences

library(ape)
library(phangorn)

set.seed(2)
tre <- rphylo(5, 1, 0.9)

set.seed(123)
rootS <- sample(c('a', 'c', 'g', 't'), 10, replace = TRUE)
rootS <- rep(c('a', 'c', 'g', 't'), 3)[1:10]
Q <- matrix(c(0.0, 0.2, 0.2, 1.0,
              0.2, 0.0, 1.0, 0.2, 
              0.2, 1.0, 0.0, 0.2, 
              1.0, 0.2, 0.2, 0.0), 
            byrow = TRUE, nrow = 4)
diag(Q) <- - rowSums(Q)

s <- simSeq(tre, l = length(rootS), rootseq = rootS, type = 'DNA', rate = 0.1, Q = Q)
rbind(rootS, as.character(s))

