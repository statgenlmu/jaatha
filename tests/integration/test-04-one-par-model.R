context("One Parameter Model")

test_that("test.oneParModel", {
    dm <- dm.createDemographicModel(5:6, 10)
    dm <- dm.addSpeciationEvent(dm, fixed = 1)
    dm <- dm.addMutation(dm, 1, 5)
    ss <- dm.simSumStats(dm, 2.5)
    expect_true(sum(ss$jsfs) > 0)
    jaatha <- Jaatha.initialize(dm, ss)
    jaatha <- Jaatha.initialSearch(jaatha, 20, 2)
    jaatha <- Jaatha.refinedSearch(jaatha, 1, 20, max.steps = 10)
})

