context("seqgen simulation interface")

test_that("test.F81Model", {
    if (!test_seqgen) return()
    set.seed(12)
    jsfs <- dm.simSumStats(dm.f81, c(1, 10))
    expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.GtrModel", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    set.seed(12)
    jsfs <- dm.simSumStats(dm.gtr, c(1, 10))
    expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.HkyModel", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    set.seed(12)
    jsfs <- dm.simSumStats(dm.hky, c(1, 10))
    expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.RateHeterogenity", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    set.seed(12)
    dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5, categories.number = 5)
    jsfs <- dm.simSumStats(dm.rh, c(1, 10, 1))
    expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.callSeqgen", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    opts <- c("seq-gen", " -mHKY", " -l", dm.getLociLength(dm.tt), 
        " -p", dm.getLociLength(dm.tt) + 1, " -q")
    dm.tt <- dm.addSummaryStatistic(dm.tt, "trees")
    ms.options <- jaatha:::generateMsOptions(dm.tt, c(1, 10))
    ms.file <- jaatha:::callMs(ms.options, dm.tt)
    seqgen.file <- callSeqgen(opts, ms.file)
    expect_true(file.exists(seqgen.file))
    expect_true(file.info(seqgen.file)$size != 0)
    unlink(ms.file)
    unlink(seqgen.file)
})

test_that("test.finalizeSeqgen", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    finalizeSeqgen = getSimProgram("seq-gen")$finalization_func
    dm.hky <- finalizeSeqgen(dm.hky)
    dm.f81 <- finalizeSeqgen(dm.f81)
    dm.gtr <- finalizeSeqgen(dm.gtr)
    expect_false(is.null(dm.hky@options[["seqgen.cmd"]]))
    expect_false(is.null(dm.hky@options[["tree.model"]]))
    expect_false(is.null(dm.f81@options[["seqgen.cmd"]]))
    expect_false(is.null(dm.f81@options[["tree.model"]]))
    expect_false(is.null(dm.gtr@options[["seqgen.cmd"]]))
    expect_false(is.null(dm.gtr@options[["tree.model"]]))
})

test_that("test.generateSeqgenOptions", {
    jaatha:::setJaathaVariable("seqgen.exe", "seq-gen")
    dm.hky@options$seqgen.cmd <- NULL
    opts <- jaatha:::generateSeqgenOptions(dm.hky, c(1, 10))
    opts <- strsplit(opts, " ")[[1]]
    expect_true(opts[1] == "seq-gen")
    expect_true("-l" %in% opts)
    expect_true("-p" %in% opts)
    expect_true("-z" %in% opts)
    expect_true("-q" %in% opts)
    expect_true("-mHKY" %in% opts)
    expect_true("-t" %in% opts)
    expect_true("-f" %in% opts)
    expect_true("-s" %in% opts)
})

test_that("test.generateTreeModel", {
    for (dm in c(dm.hky, dm.f81, dm.gtr)) {
        dm.ms <- dm.finalize(generateTreeModel(dm))
        sum.stats <- dm.simSumStats(dm.ms, c(1, 5))
        expect_false(is.null(sum.stats$file))
        expect_true(file.exists(sum.stats$file))
        unlink(sum.stats$file)
    }
})

test_that("test.seqgenMutationParameterNotLast", {
    dm.test <- dm.hky
    dm.test@parameters <- dm.test@parameters[3:1, ]
    cmd <- paste(generateSeqgenOptionsCmd(dm.test), collapse = "")
    expect_equal(grep("theta", cmd), 1)
    expect_equal(grep("tau", cmd), integer(0))
    expect_equal(grep("rho", cmd), integer(0))
})

test_that("test.seqgenSingleSimFunc", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: jar not found")
        return()
    }
    seqgenSingleSimFunc = getSimProgram("seq-gen")$sim_func
    expect_error(seqgenSingleSimFunc(dm.tt, c(1, 10)))
    set.seed(100)
    sum.stats <- seqgenSingleSimFunc(dm.hky, c(1, 10))
    expect_true(is.list(sum.stats))
    expect_true(is.array(sum.stats$jsfs))
    expect_true(sum(sum.stats$jsfs) > 0)
    set.seed(100)
    sum.stats2 <- seqgenSingleSimFunc(dm.hky, c(1, 10))
    expect_equal(sum.stats2$jsfs, sum.stats$jsfs)
})

test_that("test.seqgenWithMsms", {
    if ((!checkForSeqgen(FALSE)) | !checkForMsms(FALSE)) {
        warning("Can not test seqgen with msms: not found")
        return()
    }
    dm.selsq <- dm.addPositiveSelection(dm.f81, 100, 500, population = 1, 
        at.time = "0.1")
    dm.selsq <- dm.finalize(dm.selsq)
    expect_false(is.null(dm.selsq@options[["seqgen.cmd"]]))
    expect_false(is.null(dm.selsq@options[["tree.model"]]))
    set.seed(4444)
    sum.stats <- dm.simSumStats(dm.selsq, c(1, 5, 250))
    expect_false(is.null(sum.stats$jsfs))
    set.seed(4444)
    sum.stats2 <- dm.simSumStats(dm.selsq, c(1, 5, 250))
    expect_false(is.null(sum.stats2$jsfs))
    expect_equal(sum.stats2, sum.stats)
})

test_that("test.simulateFpcWithSeqgen", {
    if (!checkForSeqgen(FALSE)) {
        warning("Can not test seqgen: binary not found")
        return()
    }
    sum.stats <- dm.simSumStats(dm.sgfpc, c(1, 5))
    expect_false(is.null(sum.stats$fpc))
    expect_equal(sum(sum.stats$fpc), 5)
})

