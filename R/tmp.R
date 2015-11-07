# mat.list <- list()
# 
# mat1.hh <- mat1[1:8, 1:8]
# mat1.hi <- mat1[1:8, 9:27]
# mat1.ha <- mat1[1:8, 28:39]
# mat.list[['hh']] <- mat1.hh
# mat.list[['hi']] <- mat1.hi
# mat.list[['ha']] <- mat1.ha
# 
# mat1.ih <- mat1[9:27, 1:8]
# mat1.ii <- mat1[9:27, 9:27]
# mat1.ia <- mat1[9:27, 28:39]
# mat.list[['ih']] <- mat1.ih
# mat.list[['ii']] <- mat1.ii
# mat.list[['ia']] <- mat1.ia
# 
# mat1.ah <- mat1[28:39, 1:8]
# mat1.ai <- mat1[28:39, 9:27]
# mat1.aa <- mat1[28:39, 28:39]
# mat.list[['ah']] <- mat1.ah
# mat.list[['ai']] <- mat1.ai
# mat.list[['aa']] <- mat1.aa
# 
# # mat.dcov <- matrix(0, 9, 9)
# # row.names(mat.dcov) <- names(mat.list)
# # colnames(mat.dcov) <- names(mat.list)
# # for (i in 1:9) {
# #   for (j in 1:9) {
# #     if (nrow(mat.list[[i]]) == nrow(mat.list[[j]])) {
# #       mat.dcov[i, j] <- energy::dcor(mat.list[[i]], mat.list[[j]])
# #     }
# #   }
# # }
# # 
# # mat.dcov
# 
# 
# #===================
# 
# 
# h.size <- 8
# i.size <- 27 - 8 + 1
# a.size <- 39 - 28 + 1
# 
# 
# 
# # inner-group overlap Monte-Carlo Test
# # iterations:
# #   randomly re-assign groups
# #   compute mean or any other statistic for each group
# # find how often experimental mean is within resampled mean for each group,
# # i.e., p.value for each group
# overlapTest <- function (.mat, .group = c(rep('h', 8), rep('i', 27 - 8 + 1), rep('a', 39 - 28 + 1)), .fun = mean, .n = 5000) {
#   fun.results <- list()
#   groups.names <- unique(.group)
#   for (gn in .group) { }
#   fun.results[['h']] <- c()
#   fun.results[['i']] <- c()
#   fun.results[['a']] <- c()
#   # values <- .mat[upper.tri(.mat)]
#   values <- c(mat1.hh[upper.tri(mat1.hh)], mat1.ii[upper.tri(mat1.ii)], mat1.aa[upper.tri(mat1.aa)])
#   for (iter in 1:.n) {
#     values.tmp <- values
#     h1.indices <- sample(1:length(values.tmp), 8, F)
#     h1 <- values.tmp[h1.indices]
#     values.tmp <- values.tmp[-h1.indices]
#     
#     i1 <- sample(values.tmp, 20, F)
#     values.tmp <- values.tmp[!(values.tmp %in% i1)]
#     a1 <- values.tmp
#     
#     fun.results[['h']] <- c(fun.results[['h']], .fun(h1))
#     fun.results[['i']] <- c(fun.results[['i']], .fun(i1))
#     fun.results[['a']] <- c(fun.results[['a']], .fun(a1))
#   }
#   
#   res.list <- list()
#   res <- c(sum(.fun(mat1.hh[upper.tri(mat1.hh)]) < fun.results[['h']]) / .n,
#            sum(.fun(mat1.ii[upper.tri(mat1.ii)]) < fun.results[['i']]) / .n,
#            sum(.fun(mat1.aa[upper.tri(mat1.aa)]) < fun.results[['a']]) / .n)
#   names(res) <- names(fun.results)
#   res.list[['lower']] <- res
#   
#   res <- c(sum(.fun(mat1.hh[upper.tri(mat1.hh)]) > fun.results[['h']]) / .n,
#            sum(.fun(mat1.ii[upper.tri(mat1.ii)]) > fun.results[['i']]) / .n,
#            sum(.fun(mat1.aa[upper.tri(mat1.aa)]) > fun.results[['a']]) / .n)
#   names(res) <- names(fun.results)
#   res.list[['upper']] <- res
#   
#   res.list
# }
# # print(overlapTest(mat1))
# 
# # inter-group Monte-Carlo significancy test
# # ???