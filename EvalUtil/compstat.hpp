/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Leonid Boytsov, http://boytsov.info
 *
 * For details see:
 *   Leonid Boytsov, Anna Belova, Peter Westfall, 2013, 
 *   Deciding on an Adjustment for Multiplicity in IR Experiments.
 *   In Proceedings of SIGIR 2013.
 *
 */
#ifndef COMPSTAT_HPP
#define COMPSTAT_HPP

#include <vector>
#include "floattype.hpp"

struct CSampleMean {
    FloatType operator()(const std::vector<FloatType>& vVal1, const std::vector<FloatType>& vVal2) const;
};

/*
 * T-statistic
 */
struct CTStat {
    FloatType operator()(const std::vector<FloatType>& vVal1, const std::vector<FloatType>& vVal2) const;
};

/*
 * F-score for binary classification.
 */
struct CFScoreStat {
    CFScoreStat(const std::vector<FloatType>& labels) : labels(labels) {}
    FloatType operator()(const std::vector<FloatType>& vVal1, const std::vector<FloatType>& vVal2) const;

    const std::vector<FloatType> labels;
};

/*
 * Accuracy for binary classification.
 */
struct CAccuracyStat {
    CAccuracyStat(const std::vector<FloatType>& labels) : labels(labels) {}
    FloatType operator()(const std::vector<FloatType>& vVal1, const std::vector<FloatType>& vVal2) const;

    const std::vector<FloatType> labels;
};

FloatType CalcFscore(const std::vector<FloatType>& vVal, const std::vector<FloatType>& labels);
FloatType CalcAccuracy(const std::vector<FloatType>& vVal, const std::vector<FloatType>& labels);
FloatType CalcSTD(const std::vector<FloatType>& vVal);
FloatType CalcMean(const std::vector<FloatType>& vVal);

#endif
