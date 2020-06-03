# Theoretically Efficient and Practical Parallel DBSCAN
This repository contains code for our SIGMOD'20 paper: Theoretically Efficient and Practical Parallel DBSCAN (Full version: https://arxiv.org/abs/1912.06255).

Currently it is the fastest sequential and parallel Euclidean DBSCAN code, and it processes the largest dataset in the DBSCAN literature just on a single multicore machine. If you use our work please cite our paper:

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
