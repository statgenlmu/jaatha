context("SumStat JSFS")

test_that('Stat_JSFS works', {
  stats <- list(Stat_JSFS, Stat_JSFS_border, Stat_JSFS_smooth)
  for (stat in stats) {
    jsfs = stat$new(sumstat_tt$seg_sites, dm_tt,
                    coalsimr::get_summary_statistics(dm_tt)$jsfs)
    expect_that(sum(jsfs$get_data()), is_more_than(0))
    expect_equal(jsfs$transform(sumstat_tt), jsfs$get_data())
  }
})
