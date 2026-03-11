myCluster=parcluster('local');
myCluster.NumWorkers=12;
parpool(myCluster,myCluster.NumWorkers)
