% Script to initialize a parallel environment on the BU SCC cluster.

cluster = parcluster('local');
cluster.JobStorageLocation = getenv('TMPDIR');
nslots = str2double(getenv('NSLOTS'));
cluster.NumWorkers = nslots;
parpool(cluster, nslots);
