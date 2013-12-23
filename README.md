PermTest
========

Permutation algorithms for unadjusted pair-wise significance testing and testing with adjustment for multiple comparisons. This is designed primarily for a standard IR evaluation, where one or more method is represented by vectors of real-value performance scores. Each vector element is a query-specific value of an effectiveness metric such as ERR, NDCG, or MAP. Additionally, we can compute p-values of the f-score and the accuracy for binary classification. In this case, all the values are 0s and 1s. The first row in the file represents ground truth labels.


The code is released under the Apache License Version 2.0 http://www.apache.org/licenses/.

 EvalUtil:
 
    The test program itself: permtest. It accepts a matrix of performance scores (ERR, MAP, etc). 
    Each row of the matrix represent one retrieval method (called run in TREC terminology). 
    Column I represents performance scores for the I-th query. 
    In the case of binary classification, all values are 0s and 1s. The first row represents ground truth labels.

 ConvScripts:
 
    Scripts to convert TREC output file (to the matrix format). Each script accepts a registry file,
    which lists names of the files, which contain an output of a TREC evalution utility, e.g., trec_eval. 
    Each such file should represent a single run.

 A working example:
 
    To see how it works 
    
    1) Compile the Eval util
    2) Go to the directory SampleData
    3) Run the shell script sample_run.sh
    4) Read the comments inside the script


 For technical/theoretical details see:
 
   Leonid Boytsov, Anna Belova, Peter Westfall, 2013, 
   Deciding on an Adjustment for Multiplicity in IR Experiments.
   In Proceedings of SIGIR 2013.
   
 If you use our software, please, consider citing this paper.

