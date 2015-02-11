context("SumStat JSFS")

test_that('Stat_JSFS works', {
  stats <- list(Stat_JSFS, Stat_JSFS_folded, Stat_JSFS_border, Stat_JSFS_smooth)
  for (stat in stats) {
    jsfs = stat$new(sumstat_tt$seg.sites, dm_tt)
    expect_that(sum(jsfs$get_data()), is_more_than(0))
    expect_equal(jsfs$transform(sumstat_tt), jsfs$get_data())
    
    # With groups
    #jsfs = stat$new(sum.stats.grp$seg.sites.2, dm.grp, 2)
    #expect_that(sum(jsfs$get_data()), is_more_than(0))
    #expect_equal(jsfs$transform(sum.stats.grp), jsfs$get_data())
  }
})
