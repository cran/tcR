### R code from vignette source 'tcrvignette.Rnw'

###################################################
### code chunk number 1: tcrvignette.Rnw:81-85 (eval = FALSE)
###################################################
## data(twa)
## head(twa[[1]])
## data(twb)
## head(twb[[1]])


###################################################
### code chunk number 2: tcrvignette.Rnw:89-90 (eval = FALSE)
###################################################
## ?genealphabets


###################################################
### code chunk number 3: tcrvignette.Rnw:119-122
###################################################
library(tcR)
data(twb)
head(twb[[1]])


###################################################
### code chunk number 4: tcrvignette.Rnw:153-159
###################################################
# Load the package.
library(tcR)
# Load additional packages for making this vignette.
# Load the twins data, provided with the package.
data(twb)
mitcr.stats(twb)


###################################################
### code chunk number 5: tcrvignette.Rnw:167-169
###################################################
                            # How many clones fill up approximately
clonal.proportion(twb, 25)  # the 25% of the sum of values in 'Read.count'?


###################################################
### code chunk number 6: tcrvignette.Rnw:176-179
###################################################
                          # What accounts a proportion of the top-10 clones' reads
top.proportion(twb, 10)   # to the overall number of reads?
vis.top.proportions(twb)  # Plot this proportions.


###################################################
### code chunk number 7: tcrvignette.Rnw:184-187
###################################################
                                # What is a proportion of sequences which
                                # have 'Read.count' <= 100 to the
tailbound.proportion(twb, 100)  # overall number of reads?


###################################################
### code chunk number 8: tcrvignette.Rnw:193-198
###################################################
imm.in <- get.inframes(twb) # Return all in-frame sequences from the 'twb'.

                            # Count the number of out-of-frame sequences
count.outframes(twb, 5000)  # from the first 5000 sequences.
head(freq.Vb(imm.in)[,2] / freq.Vb(twb)[,2])       # Compare V-usage between in-frames and all seq.


###################################################
### code chunk number 9: tcrvignette.Rnw:201-207
###################################################
imm.in <- get.frames(twb, 'in') # Similar to 'get.inframes(twb)'.

count.frames(twb[[1]], 'all')   # Just return number of rows.

flag <- 'out'
count.frames(twb, flag, 5000)   # Similar to 'count.outframes(twb, 5000)'.


###################################################
### code chunk number 10: tcrvignette.Rnw:215-218
###################################################
# Equivalent to freq.Vb(twb[[1]]) by default.
imm1.vs <- freq.segments(twb[[1]])
head(imm1.vs)


###################################################
### code chunk number 11: tcrvignette.Rnw:221-223
###################################################
imm.vs.all <- freq.segments(twb)  # Equivalent to freq.Vb(twb) by default.
imm.vs.all[1:10, 1:4]


###################################################
### code chunk number 12: tcrvignette.Rnw:226-228
###################################################
imm1.vj <- freq.segments.2D(twb[[1]])
imm1.vj[1:5, 1:5]


###################################################
### code chunk number 13: tcrvignette.Rnw:232-234
###################################################
# Put ".dodge = F" to get distinct plot for every data frame in the given list.
vis.J.usage(twb, .cast.freq = T, .main = 'twb J-usage dodge', .dodge = T)


###################################################
### code chunk number 14: tcrvignette.Rnw:237-238
###################################################
vis.J.usage(twb, .cast.freq = T, .main = 'twb J-usage column', .dodge = F, .ncol = 2)


###################################################
### code chunk number 15: tcrvignette.Rnw:241-242
###################################################
vis.V.usage(imm1.vs, .cast.freq = F, .main = 'twb[[1]] V-usage', .coord.flip = F)


###################################################
### code chunk number 16: tcrvignette.Rnw:249-251
###################################################
cmv <- data.frame(CDR3.amino.acid.sequence = c('CASSSANYGYTF', 'CSVGRAQNEQFF', 'CASSLTGNTEAFF', 'CASSALGGAGTGELFF', 'CASSLIGVSSYNEQFF'),
                  V.segments = c('TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1'), stringsAsFactors = F)


###################################################
### code chunk number 17: tcrvignette.Rnw:254-255
###################################################
cmv


###################################################
### code chunk number 18: tcrvignette.Rnw:258-285
###################################################
twb <- set.rank(twb)
# Case 1.
cmv.imm.ex <- 
  find.clonotypes(.data = twb[1:2], .targets = cmv[,1], .method = 'exact',
                  .col.name = c('Read.count', 'Total.insertions'),
                  .verbose = F)
head(cmv.imm.ex)

# Case 2.
# Search for CDR3 sequences with hamming distance <= 1
# to the one of the cmv$CDR3.amino.acid.sequence with
# matching V-segments. Return ranks of found sequences.
cmv.imm.hamm.v <- 
  find.clonotypes(twb[1:3], cmv, 'hamm', 'Rank', 
                  .target.col = c('CDR3.amino.acid.sequence',
                                  'V.segments'),
                  .verbose = F)
head(cmv.imm.hamm.v)

# Case 3.
# Similar to the previous example, except
# using levenshtein distance and the "Read.count" column.
cmv.imm.lev.v <- 
  find.clonotypes(twb[1:3], cmv, 'lev', 
                  .target.col = c('CDR3.amino.acid.sequence', 'V.segments'),
                  .verbose = F)
head(cmv.imm.lev.v)


###################################################
### code chunk number 19: tcrvignette.Rnw:292-298
###################################################
# data(twb)
# Compute summary space of clones, that occupy
# [0, .05) and [.05, 1] proportion.
clonal.space.homeostasis(twb, c(Low = .05, High = 1))
# Use default arguments:
clonal.space.homeostasis(twb[[1]])


###################################################
### code chunk number 20: tcrvignette.Rnw:310-323
###################################################
# Equivalent to intersect(twb[[1]]$CDR3.nucleotide.sequence,
#                         twb[[2]]$CDR3.nucleotide.sequence)
# or intersectCount(twb[[1]]$CDR3.nucleotide.sequence,
#                    twb[[2]]$CDR3.nucleotide.sequence)
# "n" stands for a "CDR3.nucleotide.sequence" column, "e" for exact match.
intersect(twb[[1]], twb[[2]], 'n0e')
# "a" stands for "CDR3.amino.acid.sequence" column.
# "v" means that intersect should also use the "V.segments" column.
intersect(twb[[1]], twb[[2]], 'ave')
# Works also on lists, performs all possible pairwise intersections.
intersect(twb, 'ave')
# Plot a heatmap of number of shared clonotypes.
vis.heatmap(intersect(twb, 'ave'), .title = 'twb - (ave)-intersection', .labs = '')


###################################################
### code chunk number 21: tcrvignette.Rnw:330-336
###################################################
# Get logic vector of shared elements, where
# elements are tuples of CDR3 nucleotide sequence and corresponding V-segment
imm.1.2 <- intersectLogic(twb[[1]], twb[[2]],
                           .col = c('CDR3.amino.acid.sequence', 'V.segments'))  
# Get elements which are in both twb[[1]] and twb[[2]].
head(twb[[1]][imm.1.2, c('CDR3.amino.acid.sequence', 'V.segments')])


###################################################
### code chunk number 22: tcrvignette.Rnw:342-344
###################################################
twb.top <- top.cross(.data = twb, .n = seq(500, 10000, 500), .verbose = F, .norm = T)
top.cross.plot(twb.top)


###################################################
### code chunk number 23: tcrvignette.Rnw:351-359
###################################################
# Evaluate the diversity of clones by the ecological diversity index.
sapply(twb, function (x) diversity(x$Read.count))
# Compute the diversity as inverse probability of choosing two similar clones.
sapply(twb, function (x) inverse.simpson(x$Read.count))
# Evaluate the skewness of clonal distribution.
sapply(twb, function (x) gini(x$Read.count))
# Compute diversity of repertoire using Chao index.
t(sapply(twb, function (x) chao1(x$Read.count)))


###################################################
### code chunk number 24: tcrvignette.Rnw:377-380
###################################################
cols <- c('CDR3.amino.acid.sequence', 'Read.count')
# Apply the Morisitas overlap index to the each pair of repertoires.
apply.symm(twb, function (x,y) morisitas.index(x[, cols], y[, cols]), .verbose = F)


###################################################
### code chunk number 25: tcrvignette.Rnw:394-404
###################################################
                              # Transform "0:100" to distribution with Laplace correction 
entropy(0:100, .laplace = 1)  # (i.e., add "1" to every value before transformation).
entropy.seg(twb)  # Compute entropy of V-segment usage for each data frame. Same to
                  # apply(freq.Vb(twb)[,-1], 2, entropy)
# Next expression is equivalent to the expression
# js.div(freq.Vb(twb[[1]])[,2], freq.Vb(twb[[2]])[,2], .norm.entropy = T)
js.div.seg(twb[[1]], twb[[2]], .verbose = F)
# Also works when input arguments are list of data frames.
imm.js <- js.div.seg(twb, .verbose = F) 
vis.radarlike(imm.js, .ncol = 2)


###################################################
### code chunk number 26: tcrvignette.Rnw:410-412
###################################################
pca.segments(twb)                       # Plot PCA results of V-segment usage.
class(pca.segments(twb, .do.plot = F))  # Return object of class "prcomp"


###################################################
### code chunk number 27: tcrvignette.Rnw:418-425
###################################################
# Compute shared repertoire of amino acid CDR3 sequences and V-segments
# which has been found in two or more people.
imm.shared <- shared.repertoire(.data = twb, .type = 'avc', .min.ppl = 2, .verbose = F)
head(imm.shared)
shared.representation(imm.shared)  # Number of shared sequences.
cosine.sharing(imm.shared)         # Compute cosing similarity on shared sequences.
# It seems like repetoires are clustering in three groups: (1,2), (3,4) and (5,6).


###################################################
### code chunk number 28: tcrvignette.Rnw:440-443 (eval = FALSE)
###################################################
## p1 <- vis.count.len(twb[[1]])
## p2 <- vis.number.count(twb[[1]])
## grid.arrange(p1, p2, ncol = 2)


###################################################
### code chunk number 29: tcrvignette.Rnw:462-464
###################################################
imm.pca <- pca.segments(twb, scale. = T, .do.plot = F)
vis.pca(imm.pca, list(AB = c(1,2), CD = c(3,4)))


###################################################
### code chunk number 30: tcrvignette.Rnw:471-473
###################################################
d <- kmer.profile(c('CASLL', 'CASSQ', 'CASGL'))
vis.logo(d)


###################################################
### code chunk number 31: tcrvignette.Rnw:481-484
###################################################
# data(twb)
twb.space <- clonal.space.homeostasis(twb)
vis.clonal.space(twb.space)


###################################################
### code chunk number 32: tcrvignette.Rnw:491-495
###################################################
# data(twb)
twb.shared <- shared.repertoire(twb, .head = 1000, .verbose = F)
G <- mutation.network(twb.shared)
G


###################################################
### code chunk number 33: tcrvignette.Rnw:500-512
###################################################
# data(twb)
# twb.shared <- shared.repertoire(twb, .head = 1000)
# G <- mutation.network(twb.shared)
G <- set.group.vector(G, "twins", list(A = c(1,2), B = c(3,4)))  # <= refactor this
get.group.names(G, "twins", 1)
get.group.names(G, "twins", 300)
get.group.names(G, "twins", c(1,2,3), F)
get.group.names(G, "twins", 300, F)
# Because we have only two groups, we can assign more readable attribute.
V(G)$twin.names <- get.group.names(G, "twins")
V(G)$twin.names[1]
V(G)$twin.names[300]


###################################################
### code chunk number 34: tcrvignette.Rnw:517-521
###################################################
# data(twb)
# twb.shared <- shared.repertoire(twb, .head = 1000)
# G <- mutation.network(twb.shared)
head(mutated.neighbours(G, 1)[[1]])


###################################################
### code chunk number 35: tcrvignette.Rnw:532-534
###################################################
head(get.kmers(twb[[1]]$CDR3.amino.acid.sequence, 100, .meat = F, .verbose = F))
head(get.kmers(twb[[1]], .meat = T, .verbose = F))


###################################################
### code chunk number 36: tcrvignette.Rnw:544-548
###################################################
revcomp(c('AAATTT', 'ACGTTTGGA'))
cbind(bunch.translate(twb[[1]]$CDR3.nucleotide.sequence[1:10]),
      twb[[1]]$CDR3.amino.acid.sequence[1:10])
gc.content(twb[[1]]$CDR3.nucleotide.sequence[1:10])


###################################################
### code chunk number 37: tcrvignette.Rnw:554-560
###################################################
codon.variants('LQ')
translated.nucl.sequences(c('LQ', 'CASSLQ'))
reverse.translation('LQ')
translated.nucl.sequences('LQ', 'XXXXXG')
codon.variants('LQ', 'XXXXXG')
reverse.translation('LQ', 'XXXXXG')


