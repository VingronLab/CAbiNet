predir = './testdata/'
data <- readRDS(file.path(predir, "mini_lympho_example.rds"))


dists = as.matrix(dist(data))
saveRDS(dists, file.path(predir, 'handmade_Euclidean_dist.rds'))

test_that("Euclidean distance calculation by function 'calc_euclidean'. ", {
  eudist = readRDS(file.path(predir, 'handmade_Euclidean_dist.rds'))
  test = calc_euclidean(as.matrix(data))
  expect_equal(eudist, test)  
  })

test_that("The gene-cell association matrix is the transpose of cell-gene association matrix,",{
  ca <- suppressWarnings(APL::cacomp(data, princ_coords = 3, ntop = nrow(df)))
  expect_equal(calc_assR(ca, direction = "cells"), t(calc_assR(ca, direction = "genes")))
})

