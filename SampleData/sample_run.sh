#!/bin/bash

# 1. Convert output of TREC evaluation utilities to matrix format.
#    The TREC utilities supported:
#       i) trec_eval      http://trec.nist.gov/trec_eval/
#       ii) gdeval        http://trec.nist.gov/data/web10.html
#
#    Metrics are classic IR evaluation metrics:
#
#    err    -  Expected Reciprocal Rank
#    ndcg   - Non-Discounted Cumulative Gain
#    P_10   - Precision at ten
#    Rprec  - R-precision
#
#    You can use the matrix FORMAT directly.
#    i) The number of rows in the matrix is equal to the number of methods
#    ii) The number of columns in the matrix is equal to the number of queries
#    iii) Row X, column Y represents a performance metric measured for query Y and system X
#    
#


# Sample retrieval algorithms
# id 0: BM25
# id 1: BM25 + PageRank
# id 2: BM25 + Morphology
# id 3: BM25 + proximity (based on closed-pair statistics)
for metric in err ndcg ; do
  ../ConvScriptsTREC/conv_gdeval.pl $metric registry.gdeval InputConvertedToMatrixFormat/matrix.$metric
done

for metric in map P_10 Rprec ; do
  ../ConvScriptsTREC/conv_treceval.pl $metric registry.trec_eval InputConvertedToMatrixFormat/matrix.$metric
done

# For speedy processing set the number of permutation to be small
# In a real test, you need 50-100K permuations for the significance level of 5%
PermNum=10000

# 2. Compute unadjusted p-values
#    Parameters to control the adjustment type: -a none
for metric in map err ndcg P_10 Rprec ; do
  ../EvalUtil/permtest -n $PermNum -a none InputConvertedToMatrixFormat/matrix.$metric pValuesUnadjusted/pvalues.$metric  
done


# 3. Compute p-values adjusted for multiplicity using the MaxT algorithm:
#    We compare methods against the baseline, which is BM25 in our example.
#    Parameters to control the adjustment type: -r maxt -a baseline -b $BaselineId

BaselineId=0
for metric in map err ndcg P_10 Rprec ; do
  ../EvalUtil/permtest -r maxt -n $PermNum -a baseline -b $BaselineId InputConvertedToMatrixFormat/matrix.$metric pValuesMaxT/pvalues.$metric  
done

# 4. Compute p-values adjusted for multiplicity using closed testing.
#    These p-values should be approximately the same as values computed by MaxT.
#    We compare methods against the baseline, which is BM25 in our example.
#    Parameters to control adjustment type: -r cls -a baseline -b $BaselineId

BaselineId=0
for metric in map err ndcg P_10 Rprec ; do
  ../EvalUtil/permtest -r cls -n $PermNum -a baseline -b $BaselineId InputConvertedToMatrixFormat/matrix.$metric pValuesClosedTest/pvalues.$metric  
done

# 5. Compute p-values adjusted for multiplicity using the MaxT algorithm.
#    This time around, we do ALL-pairwise comparisons.
#    Parameters to control adjustment type: -r max -a pairwise 
for metric in map err ndcg P_10 Rprec ; do
  ../EvalUtil/permtest -r maxt -n $PermNum -a pairwise InputConvertedToMatrixFormat/matrix.$metric pValuesMaxTPairWise/pvalues.$metric  
done

# 6. Printing some values for comparison
#    The output matrix format (matrix should be symmetric):
#    row X and column Y gives the p-value corresponding to the hypothesis that system X is different from Y
#
#    In the case of the BASELINE comparison, systems are compared only against the baseline.
#    Therefore, only the row $BaselineId and the column $BaselineId contain non-zero values. 
#    The p-values where neither column nor row is $BaselineId are zero (and should be ignored).

echo "Compare unadjusted p-values, and p-values adjusted for multiplicity"
echo "Unadjusted (ERR):"
cat pValuesUnadjusted/pvalues.err
echo ""
echo "Pairwise adjustments (ERR):"
cat pValuesMaxTPairWise/pvalues.err
