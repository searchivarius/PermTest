/**
 * This is code is released under the
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

FloatType CalcSTD(const std::vector<FloatType>& vVal);
FloatType CalcMean(const std::vector<FloatType>& vVal);

#endif
