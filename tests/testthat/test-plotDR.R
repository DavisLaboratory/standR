test_that("plotPCA works with Eset", {
  library(ALL)
  data(ALL)

  #Default param
  p <- plotPCA(ALL)
  expect_silent(print(p))

  #color with column
  p <- plotPCA(ALL, color=age)
  expect_silent(print(p))

  #color with expression
  p <- plotPCA(ALL, color=age > 50)
  expect_silent(print(p))

  #multiple aesthetics
  p <- plotPCA(ALL, color=age > 50, shape=sex)
  expect_error(print(p), NA)

  #error when the column is not in the data
  p <- plotPCA(ALL, color=age2)
  expect_error(
    # We need to force eval here to throw the error
    print(p),
    "object 'age2' not found"
  )

  #will not error if the variable is present in the parent env
  age2 <- ALL$age
  p <- plotPCA(ALL, color=age2)
  expect_error(print(p), NA)
})
