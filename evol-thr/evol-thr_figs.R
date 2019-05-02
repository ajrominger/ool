# flower change through time ----

rtat <- function(xy, x0, y0, theta) {
    r <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
                nrow = 2)
    xy[, 1] <- xy[, 1] - x0
    xy[, 2] <- xy[, 2] - y0
    
    o <- xy %*% r
    o[, 1] <- o[, 1] + x0
    o[, 2] <- o[, 2] + y0
    
    return(o)
}

flwr <- function(x0, y0, r0, r1, col = 'blue') {
    s <- seq(0, 2 * pi, length.out = 50)
    
    p1 <- cbind(x0 + r1 + r1 * cos(s), y0 + 0.65 * r0 * sin(s))
    p2 <- cbind(x0 + 0.65 * r1 + 0.65 * r1 * cos(s), y0 + 0.5 * 0.65 * r0 * sin(s))
    np <- 12 + 1
    for(i in seq(0, 2 * pi, length.out = np)[-np]) {
        polygon(rtat(p1, x0, y0, i), col = col, lwd = 1)
        polygon(rtat(p2, x0, y0, i), col = hsv(0.15, 0.5, 1))
    }
    
    polygon(x0 + r0 * cos(s), y0 + r0 * sin(s), col = hsv(0.12, 0.7, 0.9), lwd = 1)
    
}


x0y0 <- function(n, lim) {
    r <- diff(lim) / 2
    prop <- pi * r^2 / (2 * r)^2
    print(prop)
    n <- ceiling(sqrt(n / prop))
    
    xiyi <- expand.grid(seq(lim[1], lim[2], length.out = n), seq(lim[1], lim[2], length.out = n))
    d <- sqrt(xiyi[, 1]^2 + xiyi[, 2]^2)
    
    return(as.matrix(xiyi[d <= r, ]))
}

pop <- function(n, cols, probs) {
    xy <- jitter(x0y0(n, c(-10, 10)))
    
    par(mar = rep(0, 4))
    plot(1, xlim = c(-12, 12), ylim = c(-12, 12), asp = 1, type = 'n', axes = FALSE)
    
    for(i in 1:nrow(xy)) {
        flwr(xy[i, 1], xy[i, 2], 0.4, 0.8, cols[sample(1:2, 1, prob = probs)])
    }
}

pdf('evol-thr/fig_pop-change.pdf', width = 16, height = 4)
layout(matrix(1:4, nrow = 1))
for(i in c(0, 0.25, 0.75, 1)) {
    pop(30, hsv(c(0.7, 0.12), c(0.7, 0.8), c(0.8, 0.9)), c(1 - i, i))
    s <- seq(0, 2 * pi, length.out = 100)
    lines(12.5 * cos(s), 12.5 * sin(s), lwd = 2)
}
dev.off()


# coalescence ----

xy <- expand.grid(1:5, 1:8)
n <- c(8, 6, 4, 2, 1)

e <- list(matrix(c(8, 7, 
                   7, 7, 
                   6, 6, 
                   5, 5, 
                   4, 5,
                   3, 4,
                   2, 3,
                   1, 2), byrow = TRUE, ncol = 2), 
          matrix(c(7, 6, 
                   6, 5, 
                   5, 5, 
                   4, 4, 
                   3, 4, 
                   2, 3), byrow = TRUE, ncol = 2), 
          matrix(c(6, 5, 
                   5, 4, 
                   4, 4, 
                   3, 4), byrow = TRUE, ncol = 2), 
          matrix(c(5, 5, 
                   4, 5), byrow = TRUE, ncol = 2))


coal <- function(ngen, xy, edges, bgcol = 'gray35') {
    par(mar = rep(0, 4), cex = 3, lwd = 2)
    plot(xy, col = bgcol, axes = FALSE)
    points(1, 5, pch = 16, col = 'gray35')
    
    if(ngen < 5) {
        for(i in 4:ngen) {
            x <- 5 - i + 1
            j <- (8 - n[i]) / 2
            y <- (1 + ceiling(j)):(8 - floor(j))
            
            segments(x, edges[[i]][, 1], x - 1, edges[[i]][, 2], col = 'gray35', lwd = 3)
            points(rep(x, length(y)), y, pch = 16, col = 'gray35')
        }
    }
}

for(k in 5:1) {
    pdf(sprintf('evol-thr/fig_coal%s.pdf', 6 - k), width = 6, height = 4)
    coal(k, xy, e, bgcol = 'transparent')
    dev.off()
}

pdf('evol-thr/fig_coal_K.pdf', width = 6, height = 4)
coal(1, xy, e, bgcol = 'gray35')
dev.off()


# coal with selection

eSelect <- e
eSelect[[1]][2, ] <- c(7, 6)
eSelect[[1]][6, ] <- c(3, 4)
eSelect[[1]][7, ] <- c(2, 4)
eSelect[[2]][4, ] <- c(4, 5)

pdf('evol-thr/fig_coal-select.pdf', width = 6, height = 4)
coal(1, xy, eSelect)
dev.off()


# drift and selection time series ----

wfDrift <- function(f0, N, tgen, s = 0) {
    ff <- rep(0, tgen)
    ff[1] <- f0
    
    if(s == 0) {
        p <- function(f) f / N
    } else {
        p <- function(f) {
            o <- f * s / N 
            if(o > 1) o <- 1
            
            return(o)
        }
    }
    for(i in 2:tgen) {
        ff[i] <- rbinom(1, N, p(ff[i - 1]))
    }
    
    return(ff / N)
}


# drift

B <- 1000
set.seed(1)
foo <- replicate(B, wfDrift(1, 100, 100))

pdf('evol-thr/fig_drift-timeseries.pdf', width = 6, height = 4)
par(mar = c(1, 3, 0, 0) + 0.5, mgp = c(2, 0.5, 0), cex = 1.4)
plot(foo[, which(foo[100, ] == 1)[1]], type = 'l', lwd = 2, 
     col = hsv(0.56, 1, 0.5, alpha = 0.8), 
     xlab = '', ylab = 'Frequency', xaxt = 'n')
mtext('Time', side = 1, line = 0.25, cex = 1.4)
points(foo[, which(colSums(foo > 0.4) > 0 & foo[100, ] == 0)[2]], type = 'l', lwd = 2,
       col = hsv(0.67, 0.7, 0.8, alpha = 0.8))
dev.off()


# selection

B <- 200
set.seed(2)
doo <- replicate(B, wfDrift(1, 100, 100, 1.08))

pdf('evol-thr/fig_selection-timeseries.pdf', width = 6, height = 4)
par(mar = c(1, 3, 0, 0) + 0.5, mgp = c(2, 0.5, 0), cex = 1.4)
plot(doo[, which.max(colSums(doo == 1))], type = 'l', lwd = 2, 
     col = hsv(0.02, 1, 0.5, alpha = 0.8), 
     xlab = '', ylab = 'Frequency', xaxt = 'n')
mtext('Time', side = 1, line = 0.25, cex = 1.4)
points(doo[, which.min(abs(colSums(doo == 1) - 30))], type = 'l', lwd = 2,
       col = hsv(0.1, 0.7, 0.8, alpha = 0.8))
dev.off()


# fitness landscape ----

