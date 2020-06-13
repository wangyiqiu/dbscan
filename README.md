# Overview
This repository contains code for our SIGMOD'20 paper: Theoretically Efficient and Practical Parallel DBSCAN (Full version: https://arxiv.org/abs/1912.06255). Currently it is the fastest sequential and parallel Euclidean DBSCAN code, and it processes the largest dataset in the DBSCAN literature just on a single multicore machine.

The main function can be found in ``DBSCANTime.C``, which calls the DBSCAN routine declared in ``DBSCAN.h`` and implemented in ``DBSCAN.C``. ``DBSCAN.C`` also contains the core subroutines for DBSCAN. By toggling the macro switches on the top of ``DBSCAN.C``, the user is able to switch among different implementations of DBSCAN described in our paper (by default, the fastest exact DBSCAN is enabled). Please see more details in the tutorial below on running our code.

If you use our work please also cite our paper:
    
    @inproceedings{wang2020theoretically,
      author = {Wang, Yiqiu and Gu, Yan and Shun, Julian},
      title = {Theoretically-Efficient and Practical Parallel DBSCAN},
      year = {2020},
      isbn = {9781450367356},
      publisher = {Association for Computing Machinery},
      address = {New York, NY, USA},
      url = {https://doi.org/10.1145/3318464.3380582},
      doi = {10.1145/3318464.3380582},
      booktitle = {Proceedings of the 2020 ACM SIGMOD International Conference on Management of Data},
      pages = {2555–2571},
      numpages = {17},
      keywords = {parallel algorithms, spatial clustering, DBScan},
      location = {Portland, OR, USA},
      series = {SIGMOD ’20}
    }
    
    @article{wang2019theoretically,
      title={Theoretically-Efficient and Practical Parallel DBSCAN},
      author={Wang, Yiqiu and Gu, Yan and Shun, Julian},
      journal={arXiv preprint arXiv:1912.06255},
      year={2019}
    }

# Tutorial

## Quick Start
The code can be compiled and run on a 64-bit Linux system with a g++ compiler.
* Clone the repository:

      git clone https://github.com/wangyiqiu/dbscan.git
        
* Navigate to the folder and compile with parallel execution:

      make clean; GCILK=1 make -j;

* Run using our example dataset ``tiny.txt`` (using ``Epsilon=9.9`` and ``Minpts=2``):

      ./DBSCAN -eps 9.9 -minpts 2 -o clusters.txt ./tiny.txt

The clustering result can then be read from clusters.txt. The generic command is:

    ./DBSCAN -eps <Epsilon> -minpts <Minpts> -o <Output-Path> <Dataset-Path>
    
To compile without parallelization and just run sequentially:

    make clean; make -j;

## Input Format

The input dataset should be a text file of line-delimited multi-dimensional points in floating point format. The first row should specify the dimension. Here's a small example dataset with five 3-dimensional points, which is also included in the repository (tiny.txt):

    dim 3
    10 20 30
    3.2 3 5
    3.1 3 5
    10 20 31
    30 0 1000
    
## Output Format

When the dataset above is run with ``Epsilon=9.9`` and ``Minpts=2``, we get output file with content:

    0 1
    1 0
    2 0
    3 1

The first column denotes the indices of the data points, for instance index ``2`` refers to data point ``3.1 3 5``. The second column denotes the cluster assignments of the data points, for instance, points ``10 20 30`` and ``10 20 31`` are in the same cluster, with id of ``1``. Note that point with index ``4`` does not have a cluster assignment for being a noise point.
