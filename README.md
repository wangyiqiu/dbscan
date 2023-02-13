# Update

## Further development has been moved to a new repository [https://github.com/wangyiqiu/dbscan-python](https://github.com/wangyiqiu/dbscan-python)

# Overview

This repository hosts code for Theoretically Efficient and Practical Parallel DBSCAN ([link](https://arxiv.org/abs/1912.06255)), research work done at MIT, and presented at SIGMOD'20. This paper bridges the gap between theory and practice of parallel DBSCAN by presenting new parallel algorithms for Euclidean exact DBSCAN and approximate DBSCAN that match the work bounds of their sequential counterparts, and are highly parallel (polylogarithmic depth).

Our experiments on a 36-core machine with hyper-threading show that we outperform existing parallel DBSCAN implementations by up to several orders of magnitude, and achieve speedups by up to 33x over the best sequential algorithms.

# Code

* Please find the usage instructions in the ``code`` directory.

* For SIGMOD'20 reproducility, please find the instructions at [https://github.com/wangyiqiu/dbscan-reproducibility](https://github.com/wangyiqiu/dbscan-reproducibility).
