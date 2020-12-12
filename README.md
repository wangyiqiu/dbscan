# Overview
This repository contains code for our SIGMOD'20 paper: Theoretically Efficient and Practical Parallel DBSCAN (Full version: https://arxiv.org/abs/1912.06255). Currently it is the fastest sequential and parallel Euclidean DBSCAN code, and it processes the largest dataset in the DBSCAN literature just on a single multicore machine.

The source code and instructions can be found in the ''c++'' subdirectory. For SIGMOD reproducibility, please refer to the ''reproducibility'' subdirectory.

*News: We are developing a Python version!*

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

