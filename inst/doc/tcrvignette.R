### R code from vignette source 'tcrvignette.Rnw'

###################################################
### code chunk number 1: tcrvignette.Rnw:80-85 (eval = FALSE)
###################################################
## data(human.alphabets)
## V_ALPHA_ALPHABET
## J_ALPHA_ALPHABET
## V_BETA_ALPHABET
## J_BETA_ALPHABET


###################################################
### code chunk number 2: tcrvignette.Rnw:91-94 (eval = FALSE)
###################################################
## data(mouse.alphabets)
## V_BETA_ALPHABET
## J_BETA_ALPHABET


###################################################
### code chunk number 3: tcrvignette.Rnw:101-105 (eval = FALSE)
###################################################
## data(twa)
## head(twa[[1]])
## data(twb)
## head(twb[[1]])


###################################################
### code chunk number 4: tcrvignette.Rnw:136-138 (eval = FALSE)
###################################################
## startmitcr('raw/TwA1_B.fastq.gz', 'mitcr/TwA1_B.txt', .file.path = '~/data/',
##             pset = 'flex', level = 1, 'debug', .mitcr.path = '~/programs/', .mem = '8g')


###################################################
### code chunk number 5: tcrvignette.Rnw:141-143 (eval = FALSE)
###################################################
## startmitcr(.file.path = '~/data/raw', pset = 'flex', level = 1, 'debug',
##             .mitcr.path = '~/programs/', .mem = '8g')


###################################################
### code chunk number 6: tcrvignette.Rnw:147-153 (eval = FALSE)
###################################################
## # Parse file in "~/data/twb1.txt".
## twb1 <- parse.file("~/data/twb1.txt")
## # Parse files "~/data/twb1.txt" and "~/data/immdat2.txt".
## twb12 <- parse.file.list(c("~/data/twb1.txt", "~/data/twb2.txt"))
## # Parse all files in "~/data/".
## twb <- parse.folder("~/data/")


###################################################
### code chunk number 7: tcrvignette.Rnw:160-163
###################################################
library(tcR)
data(twb)
head(twb[[1]])


###################################################
### code chunk number 8: tcrvignette.Rnw:187-195
###################################################
# Load the package.
library(tcR)
# Load additional packages for making this vignette.
# Load the twins data, provided with the package.
data(twb)
# Load human alphabets of V-genes and J-genes, provided with the package.
data(human.alphabets)
mitcr.stats(twb)


###################################################
### code chunk number 9: tcrvignette.Rnw:203-205
###################################################
                            # How many clones fill up approximately
clonal.proportion(twb, 25)  # the 25% of the sum of values in 'Read.count'?


###################################################
### code chunk number 10: tcrvignette.Rnw:212-215
###################################################
                          # What accounts a proportion of the top-10 clones' reads
top.proportion(twb, 10)   # to the overall number of reads?
vis.top.proportions(twb)  # Plot this proportions.


###################################################
### code chunk number 11: tcrvignette.Rnw:220-223
###################################################
                                # What is a proportion of sequences which
                                # have 'Read.count' <= 100 to the
tailbound.proportion(twb, 100)  # overall number of reads?


###################################################
### code chunk number 12: tcrvignette.Rnw:229-234
###################################################
imm.in <- get.inframes(twb) # Return all in-frame sequences from the 'twb'.

                            # Count the number of out-of-frame sequences
count.outframes(twb, 5000)  # from the first 5000 sequences.
head(freq.Vb(imm.in)[,2] / freq.Vb(twb)[,2])       # Compare V-usage between in-frames and all seq.


###################################################
### code chunk number 13: tcrvignette.Rnw:237-243
###################################################
imm.in <- get.frames(twb, 'in') # Similar to 'get.inframes(twb)'.

count.frames(twb[[1]], 'all')   # Just return number of rows.

flag <- 'out'
count.frames(twb, flag, 5000)   # Similar to 'count.outframes(twb, 5000)'.


###################################################
### code chunk number 14: tcrvignette.Rnw:251-254
###################################################
# Equivalent to freq.Vb(twb[[1]]) by default.
imm1.vs <- freq.segments(twb[[1]])
head(imm1.vs)


###################################################
### code chunk number 15: tcrvignette.Rnw:257-259
###################################################
imm.vs.all <- freq.segments(twb)  # Equivalent to freq.Vb(twb) by default.
imm.vs.all[1:10, 1:4]


###################################################
### code chunk number 16: tcrvignette.Rnw:262-264
###################################################
imm1.vj <- freq.segments.2D(twb[[1]])
imm1.vj[1:5, 1:5]


###################################################
### code chunk number 17: tcrvignette.Rnw:268-270
###################################################
# Put ".dodge = F" to get distinct plot for every data frame in the given list.
vis.J.usage(twb, .cast.freq = T, .main = 'twb J-usage dodge', .dodge = T)


###################################################
### code chunk number 18: tcrvignette.Rnw:273-274
###################################################
vis.J.usage(twb, .cast.freq = T, .main = 'twb J-usage column', .dodge = F, .ncol = 2)


###################################################
### code chunk number 19: tcrvignette.Rnw:277-278
###################################################
vis.V.usage(imm1.vs, .cast.freq = F, .main = 'twb[[1]] V-usage', .coord.flip = F)


###################################################
### code chunk number 20: tcrvignette.Rnw:285-287
###################################################
cmv <- data.frame(CDR3.amino.acid.sequence = c('CASSSANYGYTF', 'CSVGRAQNEQFF', 'CASSLTGNTEAFF', 'CASSALGGAGTGELFF', 'CASSLIGVSSYNEQFF'),
                  V.segments = c('TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1', 'TRBV4-1'), stringsAsFactors = F)


###################################################
### code chunk number 21: tcrvignette.Rnw:290-291
###################################################
cmv


###################################################
### code chunk number 22: tcrvignette.Rnw:294-320
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
                  .target.col = c('CDR3.amino.acid.sequence', 'V.segments'),
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
### code chunk number 23: tcrvignette.Rnw:332-345
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
### code chunk number 24: tcrvignette.Rnw:352-358
###################################################
# Get logic vector of shared elements, where
# elements are tuples of CDR3 nucleotide sequence and corresponding V-segment
imm.1.2 <- intersectLogic(twb[[1]], twb[[2]],
                           .col = c('CDR3.amino.acid.sequence', 'V.segments'))  
# Get elements which are in both twb[[1]] and twb[[2]].
head(twb[[1]][imm.1.2, c('CDR3.amino.acid.sequence', 'V.segments')])


###################################################
### code chunk number 25: tcrvignette.Rnw:364-366
###################################################
twb.top <- top.cross(.data = twb, .n = seq(500, 10000, 500), .verbose = F, .norm = T)
top.cross.plot(twb.top)


###################################################
### code chunk number 26: tcrvignette.Rnw:373-379
###################################################
# Evaluate the diversity of clones by the ecological diversity index.
sapply(twb, function (x) diversity(x$Read.count))
# Compute the diversity as inverse probability of choosing two similar clones.
sapply(twb, function (x) inverse.simpson(x$Read.count))
# Evaluate the skewness of clonal distribution.
sapply(twb, function (x) gini(x$Read.count))


###################################################
### code chunk number 27: tcrvignette.Rnw:397-400
###################################################
cols <- c('CDR3.amino.acid.sequence', 'Read.count')
# Apply the Morisitas overlap index to the each pair of repertoires.
apply.symm(twb, function (x,y) morisitas.index(x[, cols], y[, cols]), .verbose = F)


###################################################
### code chunk number 28: tcrvignette.Rnw:414-424
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
### code chunk number 29: tcrvignette.Rnw:430-432
###################################################
pca.segments(twb)                       # Plot PCA results of V-segment usage.
class(pca.segments(twb, .do.plot = F))  # Return object of class "prcomp"


###################################################
### code chunk number 30: tcrvignette.Rnw:438-445
###################################################
# Compute shared repertoire of amino acid CDR3 sequences and V-segments
# which has been found in two or more people.
imm.shared <- shared.repertoire(.data = twb, .type = 'avc', .min.ppl = 2, .verbose = F)
head(imm.shared)
shared.representation(imm.shared)  # Number of shared sequences.
cosine.sharing(imm.shared)         # Compute cosing similarity on shared sequences.
# It seems like repetoires are clustering in three groups: (1,2), (3,4) and (5,6).


###################################################
### code chunk number 31: tcrvignette.Rnw:455-458
###################################################
p1 <- vis.count.len(twb[[1]])
p2 <- vis.number.count(twb[[1]])
grid.arrange(p1, p2, ncol = 2)


###################################################
### code chunk number 32: tcrvignette.Rnw:478-480
###################################################
imm.pca <- pca.segments(twb, scale. = T, .do.plot = F)
vis.pca(imm.pca, list(AB = c(1,2), CD = c(3,4)))


###################################################
### code chunk number 33: tcrvignette.Rnw:491-493
###################################################
head(get.kmers(twb[[1]]$CDR3.amino.acid.sequence, 100, .meat = F, .verbose = F))
head(get.kmers(twb[[1]], .meat = T, .verbose = F))


###################################################
### code chunk number 34: tcrvignette.Rnw:503-506
###################################################
revcomp(c('AAATTT', 'ACGTTTGGA'))
cbind(bunch.translate(twb[[1]]$CDR3.nucleotide.sequence[1:10]), twb[[1]]$CDR3.amino.acid.sequence[1:10])
gc.content(twb[[1]]$CDR3.nucleotide.sequence[1:10])


###################################################
### code chunk number 35: tcrvignette.Rnw:512-518
###################################################
codon.variants('LQ')
translated.nucl.sequences(c('LQ', 'CASSLQ'))
reverse.translation('LQ')
translated.nucl.sequences('LQ', 'XXXXXG')
codon.variants('LQ', 'XXXXXG')
reverse.translation('LQ', 'XXXXXG')


