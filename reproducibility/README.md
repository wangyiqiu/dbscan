# Reproducibility (SIGMOD'20)
## Basic Info
### Source code
* Programming Language: C++
* Compiler Info: g++ 7.4
* Packages/Libraries Needed: CilkPlus (built-in to g++ 7.4)
### Datasets
Download the data sets from [Dropbox link](https://www.dropbox.com/sh/ehhv9thpuvb36jq/AADQowvv9FfQ8ZYdAPL9qJs1a?dl=0) into ``dbscan/reprodubility/datasets/`` directory.
### Hardware
Amazon Web Services EC2 c5.18xlarge with Ubuntu 18.04
* Processor: Intel Xeon Platinum 8124M (3.00GHz), 36 Cores in 2 Sockets
* Caches: L1 1.125 MiB, L2 18 MiB, L3 24.75 MiB
* Memory: 144 GiB DDR4-2666
* Secondary Storage: SSD

## Tutorial
### Compilation
First, clone the repository ``git clone https://github.com/wangyiqiu/dbscan.git``. Navigate to the source code folder (``dbscan/c++``). Compile the source with parallel execution enabled ``make clean; GCILK=1 PBBSIO=1 make -j``. Now the program can process data sets with the ''pbbs'' extension, for example ``./DBSCAN -eps 100 -minpts 10 ../reproducibility/datasets/2D_VisualVar_10M.pbbs``.

### Using our test script
Since our paper contains more than a thousand experiment data points, we provide test scripts that we used to run the experiments and generate the plots in our paper. Navigate to ``dbscan/reproducibility/scripts/``, create a symbolic link to the ``DBSCAN`` program that we compiled in the previous step ``ln -s ../../c++/DBSCAN``.

The main test script is ``genTest.py``, which will generate, run and record all the tests. First, create a subdirectory to hold raw terminal outputs in the current directory, ``mkdir outputs``. Then run the test script with a few parameters ``python3 genTest.py -m "reproducibility-dbscan-exact" -r 3 -p 72 -x``. The test script will automatically generate one or more text files prefixed with ''reproducibility-dbscan-exact'' in the ``outputs`` folder we just created.

We provide an additional script ``python3 ./parseOutput.py`` that parses the raw output into a format easier to manipulate. It will automaticaly parse all the raw output files ``./outputs/*.txt``, and generate a file called ``reproduced-dbscan.py``, containing running time decompositions for various tests that are easy to manipulate in Python. Specifically, the file contains a Python dictionary called ``plotData``, and it contains three keys, ``epsPlots``, ``minptsPlots`` and ``numprocsPlots``, whose values correspond to three sets of experiments respectively:
* experiments that test varying epsilon parameter;
* experiments that test varying minPts parameter;
* experiments that test scalability under chosen parameters.

Under each key, there will be a subdictionary keyed by the name of the data set, whose value is a subsubdictionary containing specific timing information. For example,
* to view how the program performs under varying thread-counts for data set ``2D_VisualSim_10M``, read the Python sub-dictionary by ``dataPoint = plotData['numprocsPlots']['2D_VisualSim_10M']``, where the thread-counts will be in the list ``dataPoint['x']``, and the respective timing will be in the list ``dataPoint['totaltime']``;
* to view how the program performs under varying $\epsilon$ for data set ``2D_VisualSim_10M``, read the Python sub-dictionary by ``dataPoint = plotData['epsPlots']['2D_VisualSim_10M']``, where the varying epsilon values will be in the list ``dataPoint['x']``, and the respective timing will be in the list ``dataPoint['totaltime']``.

With ``reproduced-dbscan.py``, it will be easy to verify the figures containing plots in our paper, either by plotting ``dataPoint['totaltime']`` against ``dataPoint['x']`` for each data set and each experiment; or by reading the numbers and comparing directly with our plots. For the convenience of verification, we provide all of our parsed data (similar to ``reproduced-dbscan.py`` above) and plotting scripts in ``./submitted/``. To re-generate the plots based on our data, navigate to directly``dbscan/reproducibility/scripts/submitted``, and run ``mkdir plots`` to make a new directly, followed by ``python3 paperPlots.py`` on a machine with ``matplotlib`` installed. The plots can be found in the ``plots`` folder just created.

### Testing variants

The tutorial above covers the replication of experiments for ``our-exact-bucketing`` algorithm, and testing other variants of our methods are similar -- the difference only lies in compiling the program. To switch a different method that we have included in our paper, open ``dbscan/c++/DBSCAN.C`` and modify the macros on the top the file by commenting/uncommenting them. We summarize which macro to uncomment for each respective method (and comment-out the rest):
* our-exact-qt: ``OUR_EXACT_QT``
* our-exact-qt-bucketing: ``OUR_EXACT_QT`` + ``USE_BUCKETING``
* our-exact: ``OUR_EXACT``
* our-exact-bucketing (default): ``OUR_EXACT`` + ``USE_BUCKETING``
* our-approx-qt: ``OUR_APPROX_QT``
* our-approx-qt-bucketing: ``OUR_APPROX_QT``+ ``USE_BUCKETING``
* our-approx: ``OUR_APPROX``
* our-approx-bucketing: ``OUR_APPROX``+ ``USE_BUCKETING``

Then, recompile the program to generate a new ``DBSCAN`` binary, and the rest of the steps are the same. We note that the script that parses the raw terminal output, ``parseOutput.py`` will parse all the files in the ``dbscan/reproducibility/scripts/outputs`` folder, whose contents need to be deleted before testing a different method.
