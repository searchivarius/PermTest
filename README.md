PermTest
========

Permutation algorithms for unadjusted pair-wise significance testing and multiple comparison adjustments
The code is released under the Apache License Version 2.0 http://www.apache.org/licenses/.

 EvalUtil:
 
    The test program itself: permtest. It accepts a matrix of performance scores (ERR, MAP, etc). 
    Each row of the matrix represent one retrieval method (called run in TREC terminology). 
    Column I represents performance scores for the I-th query.

 ConvScripts:
 
    Scripts to convert TREC output file (to the matrix format). Each script accepts a registry file,
    which lists names of the files with outputs of TREC utility. One file per retrieval method,
    or, in other words per TREC run.

 More detail (a working example):
 
    To see how it works 
    
    1) Compute the Eval util
    2) Go to the directory SampleData
    3) Run the shell script sample_run.sh
    4) Read the comments inside the script


 For technical/theoretical details see:
 
   Leonid Boytsov, Anna Belova, Peter Westfall, 2013, 
   Deciding on an Adjustment for Multiplicity in IR Experiments.
   In Proceedings of SIGIR 2013.

