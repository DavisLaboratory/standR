test_that("plotRLE works with Eset", {
  library(ALL)
  data(ALL)

  #Default param
  p <- plotRLE(ALL)
  expect_silent(print(p))

  #color with column
  p <- plotRLE(ALL, color=age)
  expect_silent(print(p))

  #color with expression
  p <- plotRLE(ALL, color=age > 50)
  expect_silent(print(p))

  #multiple aesthetics
  p <- plotRLE(ALL, color=age > 50, shape=sex)
  expect_error(print(p), NA)

  #order with one column
  p <- plotRLE(ALL, ordannots=age)
  expect_silent(print(p))
  p <- plotRLE(ALL, ordannots=c(age))
  expect_silent(print(p))

  #order with multiple columns
  p <- plotRLE(ALL, ordannots=c(sex, age))
  expect_silent(print(p))

  #error when the column is not in the data
  p <- plotRLE(ALL, color=age2)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'age2' not found"
  )

  #will not error if the variable is present in the parent env
  age2 <- ALL$age
  p <- plotRLE(ALL, color=age2)
  expect_error(print(p), NA)

  #SCE
  sce = SingleCellExperiment::SingleCellExperiment(list('a' = matrix(rnorm(6e3), 3)))
  expect_silent(print(plotRLE(sce[, 1:10])))
  expect_silent(print(plotRLE(sce[, 1:100])))
  expect_silent(print(plotRLE(sce)))
})
