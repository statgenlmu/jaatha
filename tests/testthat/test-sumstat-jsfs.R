context("SumStat JSFS")

test_that('Stat_JSFS works', {
  stats <- c(Stat_JSFS, Stat_JSFS_folded, Stat_JSFS_border, Stat_JSFS_smooth)
  for (stat in stats) {
    jsfs = stat$new(sum.stats.mig$seg.sites, dm.mig)
    expect_that(sum(jsfs$get_data()), is_more_than(0))
    expect_equal(jsfs$transform(sum.stats.mig), jsfs$get_data())
    
    # With groups
    jsfs = stat$new(sum.stats.grp$seg.sites.2, dm.grp, 2)
    expect_that(sum(jsfs$get_data()), is_more_than(0))
    expect_equal(jsfs$transform(sum.stats.grp), jsfs$get_data())
  }
})
