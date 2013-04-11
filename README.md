PermTest
========

Permutation tests for unadjusted pair-wise significance testing and multiple comparison adjustments

 EvalUtil
    The binary utility permtest. It accepts a matrix of performance scores. Each row
    represent one method. Column i represents performance scores for the i-th query.

 ConvScripts
    Scripts to convert TREC output file (to the matrix format). Each script accepts a registry file,
    which lists names of the files with outputs of TREC utility. One file per retrieval method:
    or, per run in TREC terminology.
    

 For more details see:
   Leonid Boytsov, Anna Belova, Peter Westfall, 2013, 
   Deciding on an Adjustment for Multiplicity in IR Experiments.
   In Proceedings of SIGIR 2013.
