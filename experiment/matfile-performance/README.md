# matfile IO performance

prep\_series currently holds everything in memory, and waits until the end to write. This
design was a fix for super-slow incremental writes, and a means to enable parallelism.

We clearly need to change the approach in a way that:
1. Does not choke our performance
2. If possible, allows us to stay parallel
3. Does not add undue complexity
