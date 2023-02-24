test_that("check inputs", {
  expect_s3_class(TADkit::tad_1_10kb.bed, "data.frame")
  expect_s3_class(TADkit::tad_2_10kb.bed, "data.frame")
  expect_s4_class(TADkit::matrix_1_chr25_50kb, "dgCMatrix")
  expect_s4_class(TADkit::matrix_2_chr25_50kb, "dgCMatrix")
  expect_s4_class(TADkit::matrix_1_chr26_50kb, "dgCMatrix")
  expect_s4_class(TADkit::matrix_2_chr26_50kb, "dgCMatrix")
  expect_s3_class(TADkit::IS_1_10kb.bedgraph, "data.frame")
  expect_s3_class(TADkit::IS_2_10kb.bedgraph, "data.frame")
  expect_s3_class(TADkit::PC1_1_50kb.bedgraph, "data.frame")
  expect_s3_class(TADkit::PC1_2_50kb.bedgraph, "data.frame")
  })

