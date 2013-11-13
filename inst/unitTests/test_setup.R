# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# Setup for tests
dm.tt        <- dm.createThetaTauModel(11:12, 10)
sum.stats.tt <- dm.simSumStats(dm.tt, c(1, 5))
jaatha.tt    <- Jaatha.initialize(dm.tt, jsfs=sum.stats.tt) 

dm.mig        <- dm.addSymmetricMigration(dm.tt, 1, 5)
sum.stats.mig <- dm.simSumStats(dm.mig, c(1, 1, 5))
jaatha.mig    <- Jaatha.initialize(dm.mig, jsfs=sum.stats.mig) 
