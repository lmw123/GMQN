#author: Mengwei Li
#email: limengwei@big.ac.cn

options(stringsAsFactors = F)
library(minfi)

source("./script/betaEst2.R")
source("./script/blc2.R")
source("./script/BMIQ_1.4.R")
source("./script/gmqn.R")

#  test 450K 
rgset = read.metharray.exp("./450K_test_data")
mset = preprocessIllumina(rgset)
m = data.frame(getMeth(mset))
um = data.frame(getUnmeth(mset))

load("./script/probe_type_450k.RData")
beta = gmqn(as.numeric(m[,1]), as.numeric(um[,1]), row.names(m))


#  test 850K 
rgset = read.metharray.exp("./850K_test_data")
mset = preprocessIllumina(rgset)
m = data.frame(getMeth(mset))
um = data.frame(getUnmeth(mset))

load("./script/probe_type_850k.RData")
beta = gmqn(as.numeric(m[,1]), as.numeric(um[,1]), row.names(m))







