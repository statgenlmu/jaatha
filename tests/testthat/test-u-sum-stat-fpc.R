context("FGC Summary Statistic")

test_that("test.calcFpcBreaks", {
    dm = calcFpcBreaks(dm.fpc, seg.sites)
    expect_false(is.null(dm@options[["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["fpc.breaks.far"]]))
    dm = calcFpcBreaks(dm.fpc, seg.sites, group = 1)
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.far"]]))
    dm = calcFpcBreaks(dm, seg.sites, group = 2)
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.far"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc.breaks.far"]]))
})

