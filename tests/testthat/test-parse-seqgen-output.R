context("parsing seq-gen output")

test_that("parseSeqgenOutput works with a single file", {
  dm_tmp <- dm.createThetaTauModel(c(4,6), 2, 10)
  dm_tmp <- dm.addOutgroup(dm_tmp, "2*tau")
  dm_tmp <- dm.setMutationModel(dm_tmp, "HKY")

  seqgen_file <- tempfile('seqgen')
  cat(" 11 10
s11       AATTTTGCCT
s2        TTCCCAAGTT 
s4        TTCACAAGTG
s1        TTCCCAAGTG 
s3        TTCCTAAGTG
s5        TCGGAAGCAG
s7        TCGGAAGCAG
s6        CCGGAAGCCT
s8        GCGGAAGCCT
s9        CCGGCTGCAG
s10       CCTCAGGGCC
 11 10
s11       ATTGAACCGC
s5        GTATATTTAC
s9        GAATATGAAG
s6        CTATATTTAG
s8        CTAAATGAGG
s7        CTATATGAAC
s10       CTATATGAAC
s1        CCATACGATA
s2        CTTGACGGTA
s3        GCAGACGGTA
s4        GCTGATAATA", file = seqgen_file)
  
  seg_sites <- parseSeqgenOutput(list(seqgen_file), 
                                 dm.getSampleSize(dm_tmp), 
                                 dm.getLociNumber(dm_tmp))
  expect_is(seg_sites, 'list')
  expect_equal(length(seg_sites), 2)
    
  seg_sites_1 <- matrix(c(1,1,1,1,1,1,1,
                          1,1,1,1,1,1,0,
                          1,0,1,1,1,1,1,
                          1,1,1,1,1,1,1,
                          1,1,1,0,0,1,1,
                          1,1,1,0,0,0,0,
                          1,1,1,0,0,1,1,
                          1,1,1,0,0,0,0,
                          1,1,0,0,0,1,1,
                          0,1,1,0,1,0,1), 
                        10, 7, byrow=TRUE)
  attr(seg_sites_1, 'positions') <- c(2, 4:9)/9
  expect_equal(seg_sites[[1]], seg_sites_1)
  
  seg_sites_2 <- matrix(c(1,1,1,1,1,
                          0,0,0,1,1,
                          1,1,0,1,1,
                          1,0,0,1,1,
                          0,1,1,1,0,
                          0,1,1,1,1,
                          0,1,1,1,0,
                          0,1,1,0,1,
                          1,1,1,1,1,
                          0,1,1,1,0), 
                        10, 5, byrow=TRUE)
  attr(seg_sites_2, 'positions') <- c(1,2,3,8,9)/9
  expect_equal(seg_sites[[2]], seg_sites_2)
  
  # Muliple files
  seg_sites <- parseSeqgenOutput(list(seqgen_file, seqgen_file), 
                                 dm.getSampleSize(dm_tmp), 4)
  expect_is(seg_sites, 'list')
  expect_equal(length(seg_sites), 4)
  expect_equal(seg_sites[[1]], seg_sites_1)
  expect_equal(seg_sites[[2]], seg_sites_2)
  expect_equal(seg_sites[[3]], seg_sites_1)
  expect_equal(seg_sites[[4]], seg_sites_2)
})
