context('Simulation Inferface')

test_that("test.createSimProgram", {
    expect_true(exists(".jaatha"))
    expect_true(exists("sim_progs", envir = .jaatha))
    sim_prog_nr = length(.jaatha$sim_progs)
    expect_error(createSimProgram("name", "feature", sin, cos))
    expect_error(createSimProgram("name"))
    createSimProgram("test1", "feature", "sum.stat", sin)
    expect_equal(length(.jaatha$sim_progs), sim_prog_nr + 1)
    createSimProgram("test2", "feature", "sum.stat", sin, cos)
    expect_equal(length(.jaatha$sim_progs), sim_prog_nr + 2)
    createSimProgram("test2", "feature", "sum.stat", sin, cos, 
        tan)
    expect_equal(length(.jaatha$sim_progs), sim_prog_nr + 2)
})

