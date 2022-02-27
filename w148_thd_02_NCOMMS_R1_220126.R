#########################################################
# Accompanying code for Merker et al. 2022
# Jean-Philippe Rasigade
# AGPL3.0 license

#########################################################
# PACKAGES

# install THD package from github if required,
# devtools::install_github("rasigadelab/thd")

library(thd)
library(ape)
library(readxl)
library(data.table)
library(stringr)
library(lmerTest)
library(beeswarm)

# data.table::order uses C-locale which is inconsistent with base::order
# Force locale to C
Sys.setlocale("LC_ALL","C")

#########################################################
# Genomic data

# Load FASTA then build SNP distance matrix
dna <- read.dna("W148_720strains.fasta", "fasta")
dna_len <- 5264
dna_names <- sort(labels(dna))

# Hamming distance matrix
H <- dist.dna(dna, "raw", pairwise.deletion = T, as.matrix = T) * dna_len

#########################################################
# Metadata

# Load metadata
load("metadata.Rdata") # data.table format, variable name `d`
H <- H[ d$Names, d$Names ]
stopifnot(all(d$Names == dna_names))

############################################################
# COMPUTE UNWEIGHTED AND WEIGHTED THD

# Rationale: as long as we can estimate variations of sampling effort, these can be taken into account
# during THD computation by weighting the isolates in the dataset. E.g. isolates with perfect sampling have
# weight 1, but isolates with 50% sampling effort are random representatives of a group twice as large: they receive
# weight 2.

# Reimplementation of main THD function using weighted mean

# Probability mass function of right-truncated geometric distribution with parameter alpha and
# truncation limit p
dtgeom <- function(x, alpha, p) {
	k <- (1 - alpha) / (1 - alpha^(p+1)) * alpha^x
	k[!is.finite(k)] <- 1/p
	return(k)
}

# Geometric kernel density estimation, weighted version, provides weighted THD estimates
# See https://github.com/rasigadelab/thd for unweighted implementation
# See https://www.nature.com/articles/srep45326 for mathematical details
gkdew <- function(h, alpha = 0.5, p = max(h), weights, from, to, skipself = T) {
  if(missing(from) & missing(to)) {
    # Pairwise densities
    n <- nrow(h)
    k <- NULL 
    if(skipself) 
      k <- sapply(1:n, function(i) {
        weighted.mean(dtgeom(h[,i][-i], alpha, p), weights[-i])
      })
    else
      k <- sapply(1:n, function(i) {
        weighted.mean(dtgeom(h[,i], alpha, p), weights[-i])
      })
    return(k)		
  } else {
    # Densities from one group to another
    # Chekc whether from = to and skip diagonal
    if(missing(from)) stop("Argument 'to' cannot be provided alone")
    if(missing(to)) to <- from
    from[is.na(from)] <- FALSE
    to[is.na(to)]     <- FALSE
    skipdiag <- all(from == to)
    hh <- h[to, from] # Scan columns rather than rows for efficiency
    if(sum(to) == 1)   hh <- matrix(hh, 1, sum(from))
    if(sum(from) == 1) hh <- matrix(hh, 1, sum(to))
    n <- ncol(hh)
    kdefunc <- if(skipdiag) {
      function(i, hh, alpha, p) {weighted.mean(dtgeom(hh[,i][-i], alpha, p), weights[-i])}
    } else {
      function(i, hh, alpha, p) {weighted.mean(dtgeom(hh[,i], alpha, p), weights[-i])}
    }
    k <- sapply(1:n, kdefunc, hh = hh, alpha = alpha, p = p)
    return(k)
  }
}

# THD parameters

# Effective genome size (H37RV size with 10% masked)
effsize <- 3969000
# Mutation rate per site per year
mu <- 1.07e-7
# 10y timescale, recent epidemic
timescale <- 10
bandwidth <- thd.bandwidth(timescale, effsize, mu, 1/2)

# Unweighted THD
thd <- thd::thd(H, timescale, effsize, mu)
d[, thd := thd]

# Weighted THD
thdw <- gkdew(H, bandwidth, effsize, d$thd_weight)
d[, thdw := thdw]

# Comparison of weighted and unweighted THD. Differences are globally moderate
plot(thd, thdw)

d[, thd0 := thd(H, t = 0.5, m = effsize, mu = mu)]

plot(d$thd0, d$length, log = "x")

# Correlation of terminal branch lengths and THD as a function of timescale

xseq <- seq(0, 8, length.out = 100)
xx <- 2^(xseq)
yy <- sapply(xseq, function(t) cor(d$length, log(thd(H, t = 2^t, m = effsize, mu = mu)), method = "p"))

svg(file = "thd tbl correlation.svg", 6, 6)
plot(xx, yy, log = "x", type = "l", lwd = 2, col = "darkblue",
	 xlab = "THD timescale (years, log scale)", ylab = "Pearson correlation of log-THD and TBL")
dev.off()
############################################################
# THD MODELLING

# d$xdr2 is the resistance status
table(d$xdr2)
table(d$xdr2, d$possible_compensation)
fisher.test(d$xdr2, d$possible_compensation)

# Consider the basic interaction model:
summary(lm(thdw ~ genotypic_resistances:possible_compensation:xdr2, d))

############################################################
# CORRECTION FOR POPULATION STRUCTURE
# See https://onlinelibrary.wiley.com/doi/full/10.1111/eva.12991 for details

# Control variables:
# Control for population structure using relevant genetic principal components
# Control for country-level variations as a random effect

# Find relevant genetic PCs for control
mds <- cmdscale(H, k = 20)
stopifnot(all(rownames(mds) == dna_names))

# Extract most relevant PCs using stepwise regression
summary(step(lm(thdw ~ ., as.data.frame(mds))))

# PCs are 1,2,3,7; include them in data table
d[, mds1 := mds[,1]]
d[, mds2 := mds[,2]]
d[, mds3 := mds[,3]]
d[, mds7 := mds[,7]]

##############################################
# FIGURE 5 THD / COMPENSATION BOXPLOT

# Use display form of THD, x100
d[, thdw100 := thdw * 100]

# Figure 5 panel a

# Helper function to control color transparency
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Figure 5a
{
  svg("thd_boxplot.svg", 6, 6)
  par(mar = c(6,6,2,2))
  comp_colors <- c("darkturquoise", "lightcoral")
  
  boxplot(thdw100 ~ xdr2 + possible_compensation, d, col = rep(add.alpha(comp_colors, .25), each = 3), las = 2, outline = F, xaxt = "n", xlab = "", ylab = "THD success index x100", ylim = c(0, 0.025)*100)
  beeswarm(thdw100 ~ xdr2 + possible_compensation, d, method = "hex", add = T, pch = 19, cex = 0.6, col = rep(add.alpha(comp_colors, .75), each = 3), corral = "gutter")
  
  legend("topleft", bty = "n",
         title = "Possible compensatory mutations",
         legend = c("0", ">0"), fill = comp_colors, xjust = .5)
  axis(1, at = 1:6, labels = rep(c(
    "MDR", "Pre-XDR", "XDR"
  ), 2))
  dev.off()
}

# Pairwise comparison of THD distributions across groups using Mann-whitney U-test
pairwise.wilcox.test(thdw, paste(d$possible_compensation, d$xdr2), p.adjust.method = "none")

# Complement analysis with terminal branch lengths
hist(d$length)


# Tip length comparison
pairwise.wilcox.test(d$length, d$xdr2)

boxplot(length ~ xdr2, d)
boxplot(length ~ possible_compensation, d)

{
	svg("tbl_boxplot.svg", 6, 6)
	par(mar = c(6,6,2,2))
	comp_colors <- c("darkturquoise", "lightcoral")
	# comp_colors <- c(rgb(0,.7,.7,.5), rgb(1,.5,.5,.5))
	boxplot(length ~ xdr2 + possible_compensation, d, col = rep(add.alpha(comp_colors, .25), each = 3), 
			las = 2, outline = F, xaxt = "n", xlab = "", ylab = "Terminal branch length (subst/site)", 
			ylim = c(0, 0.004))
	beeswarm(length ~ xdr2 + possible_compensation, d, method = "hex", add = T, pch = 19, cex = 0.6, col = rep(add.alpha(comp_colors, .75), each = 3), corral = "gutter")
	# legend("topleft", bty = "n", legend = c("No compensatory mut", "With compensatory mut"), fill = comp_colors)
	legend("topleft", bty = "n",
		   title = "Possible compensatory mutations",
		   legend = c("0", ">0"), fill = comp_colors, xjust = .5)
	axis(1, at = 1:6, labels = rep(c(
		"MDR", "Pre-XDR", "XDR"
	), 2))
	dev.off()
}

# Figure 5 panel b
# Compare number of resistance mutations in isolates with and without compensatory mutation(s)
{
  svg("mutcomp_boxplot.svg", 3.5, 6)
  par(mfrow = c(1,1))
  par(mar = c(6,6,2,2))
  boxplot(genotypic_resistances ~ possible_compensation, d, col = add.alpha(comp_colors, .25), las = 2, outline = T, xaxt = "n", xlab = "Possible compensatory mutations", ylab = "Resistance mutations count", ylim = c(0, 12))
  axis(1, at = 1:2, labels = c("0", ">0"))
  dev.off()
}

wilcox.test(genotypic_resistances ~ possible_compensation, d)

######################################################
# Figure 5 panel c

# Show variations of THD as a function of the accumulation of resistance-conferring mutations
# in isolates with and without compensatory mutations. Show confidence bands of model after correction for
# population structure. Note that LMM with correction for the country of isolation is not used here
# as confidence bands cannot be reproduced.

mod <- lm(thdw100 ~ mds1 + mds2 + mds3 + mds7 + genotypic_resistances:possible_compensation:xdr2, d[])
summary(mod)

# Prediction DF
xx <- seq(0, 11.5, by = 0.25)
l <- length(xx)
predDf <- data.table(
	genotypic_resistances = xx
)
predDf[, mds1 := 0][, mds2 := 0][, mds3 := 0][, mds7 := 0]

# Make three panels for not-XDR, pre-XDR, XDR with fixed axes, and two separate bands for compensatory and no-compensatory (blue vs red)
# Shift with/without mut on the left and right for readability

displayShift <- .16
d[ , resmutDisplay := genotypic_resistances]
d[ possible_compensation == "no", resmutDisplay := resmutDisplay - displayShift]
d[ possible_compensation == "yes", resmutDisplay := resmutDisplay + displayShift]

# Define colors

comp_colors <- c("darkturquoise", "lightcoral")
band_alpha <- 0.25
dot_alpha <- 0.75

d[ , resmutCol := comp_colors[1]]
d[ possible_compensation == "yes", resmutCol := comp_colors[2]]

# Small routine to add prediction bands on plot
addBands <- function() {
	# Blue band
	predDf[, possible_compensation := "no"]
	band <- predict(mod, predDf, interval = "confidence")
	
	polygon(c(xx, rev(xx)), c(band[,2], rev(band[,3])), col = add.alpha(comp_colors[1], band_alpha), border = NA)
	lines(xx, band[,1], col = comp_colors[1], lwd = 2)
	
	# Red band
	predDf[, possible_compensation := "yes"]
	band <- predict(mod, predDf, interval = "confidence")
	
	polygon(c(xx, rev(xx)), c(band[,2], rev(band[,3])), col = add.alpha(comp_colors[2], band_alpha), border = NA)
	lines(xx, band[,1], col = comp_colors[2], lwd = 2)  
}

drawPanel <- function(xdrCase) {
	predDf[, xdr2 := xdrCase]
	
	xdr2filt <- d$xdr2 == xdrCase
	plot(thdw100 ~ resmutDisplay, d[xdr2filt], pch = 19, col = add.alpha(d$resmutCol[xdr2filt], dot_alpha),
		 cex = 1.5, xlim = range(d$resmutDisplay[xdr2filt], na.rm = TRUE) + c(-.5,+.5), ylim = c(0, 2.5), las = 1,
		 xlab = "Resistance mutations count",
		 ylab = "THD success index x100",
		 main = xdrCase) 
	
	addBands()
	if(xdrCase == "MDR") legend("topright", bty = "n",
		   title = "Possible compensatory mutations",
		   legend = c("0", ">0"), fill = comp_colors)
}

{
	svg("mutations_modelplots.svg", 9, 3)
	par(mfrow = c(1,3))
	par(mar = c(4,4,4,4))
	drawPanel("MDR")
	drawPanel("pre-XDR")
	drawPanel("XDR") 
	dev.off()
}

######################################################
# Supplementary Figure 5 panel c

# Same approach as above but using terminal branch lengths as outcome variable
mod <- lm(length ~ genotypic_resistances:possible_compensation:xdr2, d[])
summary(mod)

# Prediction DF
xx <- seq(0, 11.5, by = 0.25)
l <- length(xx)
predDf <- data.table(
	genotypic_resistances = xx
)

# Make three panels for not-XDR, pre-XDR, XDR with fixed axes, and two separate bands for compensatory and no-compensatory (blue vs red)
# Shift with/without mut on the left and right for readability

displayShift <- .16
d[ , resmutDisplay := genotypic_resistances]
d[ possible_compensation == "no", resmutDisplay := resmutDisplay - displayShift]
d[ possible_compensation == "yes", resmutDisplay := resmutDisplay + displayShift]

# Define colors

comp_colors <- c("darkturquoise", "lightcoral")
band_alpha <- 0.25
dot_alpha <- 0.75

d[ , resmutCol := comp_colors[1]]
d[ possible_compensation == "yes", resmutCol := comp_colors[2]]

# Small routine to add prediction bands on plot
addBands <- function() {
	# Blue band
	predDf[, possible_compensation := "no"]
	band <- predict(mod, predDf, interval = "confidence")
	
	polygon(c(xx, rev(xx)), c(band[,2], rev(band[,3])), col = add.alpha(comp_colors[1], band_alpha), border = NA)
	lines(xx, band[,1], col = comp_colors[1], lwd = 2)
	
	# Red band
	predDf[, possible_compensation := "yes"]
	band <- predict(mod, predDf, interval = "confidence")
	
	polygon(c(xx, rev(xx)), c(band[,2], rev(band[,3])), col = add.alpha(comp_colors[2], band_alpha), border = NA)
	lines(xx, band[,1], col = comp_colors[2], lwd = 2)  
}

drawPanel <- function(xdrCase) {
	predDf[, xdr2 := xdrCase]
	
	xdr2filt <- d$xdr2 == xdrCase
	plot(length ~ resmutDisplay, d[xdr2filt], pch = 19, col = add.alpha(d$resmutCol[xdr2filt], dot_alpha),
		 cex = 1.5, xlim = range(d$resmutDisplay[xdr2filt], na.rm = TRUE) + c(-.5,+.5), ylim = c(0, 0.0045), las = 1,
		 xlab = "Resistance mutations count",
		 ylab = "Terminal branch length",
		 main = xdrCase) 
	
	addBands()
	if(xdrCase == "MDR") legend("topright", bty = "n",
								title = "Possible compensatory mutations",
								legend = c("0", ">0"), fill = comp_colors)
}

{
	svg("tbl_mutations_modelplots.svg", 9, 3)
	par(mfrow = c(1,3))
	par(mar = c(4,4,4,4))
	drawPanel("MDR")
	drawPanel("pre-XDR")
	drawPanel("XDR") 
	dev.off()
}


#######################################################
# REGRESSION MODELS

# Base models

# Compensatory mutations and resistance status do not improve model fit when
# considered separately
summary(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + possible_compensation + (1 | country_of_isolation), d))
anova(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + xdr2 + (1 | country_of_isolation), d))

# Interaction model

# Interaction term improves model fit
anova(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + possible_compensation*xdr2 + (1 | country_of_isolation), d[]))
summary(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + possible_compensation*xdr2 + (1 | country_of_isolation), d[]))

# Resistance accumulation model

# Resistance accumuation negatively coorelates with THD
summary(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + genotypic_resistances + (1 | country_of_isolation), d))

# 3-way interaction model
summary(lmer(thdw ~ mds1 + mds2 + mds3 + mds7 + genotypic_resistances:possible_compensation:xdr2 + (1 | country_of_isolation), d[]))


#######################################################
# REGRESSION MODELS - TERMINAL BRANCH LENGTHS

# Base models

# Find relevant genetic PCs for control
mds <- cmdscale(H, k = 20)
stopifnot(all(rownames(mds) == dna_names))

# Extract most relevant PCs using stepwise regression
summary(step(lm(d$length ~ ., as.data.frame(mds))))

# PCs are 2, 11, 18, 9; include them in data table
d[, mds2 := mds[,2]]
d[, mds11 := mds[,11]]
d[, mds18 := mds[,18]]
d[, mds9 := mds[,9]]


# Compensatory mutations and resistance status do not improve model fit when
# considered separately
summary(lmer(length ~ mds2 + mds11 + mds18 + mds9 + possible_compensation + (1 | country_of_isolation), d))
anova(lmer(length ~ mds2 + mds11 + mds18 + mds9 + xdr2 + (1 | country_of_isolation), d))

# Interaction model

# Interaction term improves model fit
anova(lmer(length ~ mds2 + mds11 + mds18 + mds9 + possible_compensation*xdr2 + (1 | country_of_isolation), d[]))
summary(lmer(length ~ mds2 + mds11 + mds18 + mds9 + possible_compensation*xdr2 + (1 | country_of_isolation), d[]))

# Resistance accumulation model

# Resistance accumuation does not correlate with terminal branch length
summary(lmer(length ~ mds2 + mds11 + mds18 + mds9 + genotypic_resistances + (1 | country_of_isolation), d))

# 3-way interaction model
summary(lmer(length ~ mds2 + mds11 + mds18 + mds9 + genotypic_resistances:possible_compensation:xdr2 + (1 | country_of_isolation), d[]))







