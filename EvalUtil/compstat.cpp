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
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <limits>

#include <cmath>
#include <cstdlib>

#include "compstat.hpp"

using namespace std;

FloatType 
CSampleMean::operator()(const vector<FloatType>& vVal1, const vector<FloatType>& vVal2) const
{
    unsigned    N = vVal1.size();

    if (N != vVal2.size()) {
        stringstream err;
        err << "Internal error, vector quantities are not equal: " << vVal1.size() << " vs. vVal2.size() ";

        throw runtime_error(err.str());
    }

    FloatType   fDelta = 0;

    for (unsigned i = 0; i < N; ++i) {
        fDelta += vVal1[i] - vVal2[i];
    }

    return fDelta / N;
}

FloatType 
CTStat::operator()(const vector<FloatType>& vVal1, const vector<FloatType>& vVal2) const
{
    unsigned    N = vVal1.size();

    if (N != vVal2.size()) {
        stringstream err;
        err << "Internal error, vector quantities are not equal: " << vVal1.size() << " vs. vVal2.size() ";

        throw runtime_error(err.str());
    }

    FloatType   fDelta = 0;

    for (unsigned i = 0; i < N; ++i) {
        fDelta += vVal1[i] - vVal2[i];
    }

    FloatType   fMean = fDelta / N;

    FloatType   fSqSum = 0;

    for (unsigned i = 0; i < N; ++i) {
        FloatType d = (vVal1[i] - vVal2[i] - fMean);

        fSqSum += d*d;
    }

    if (fSqSum < numeric_limits<FloatType>::min()) {
//return numeric_limits<FloatType>::max();
return 0;
        for (unsigned i = 0; i < N; ++i) {
            cerr << vVal1[i] << "\t";
        }
        cerr << endl;
        for (unsigned i = 0; i < N; ++i) {
            cerr << vVal2[i] << "\t";
        }
        cerr << endl;
        // TODO/TBD are going to process this?
        cerr << "The sum of squared differences is: " << fSqSum << endl;
        throw runtime_error("Cannot compute statistics in the presence of duplicate rows");
    }

    return (fMean) / sqrt(fSqSum/((N - 1)*N));
}

FloatType 
CalcMean(const vector<FloatType>& vVal)
{
    unsigned    N = vVal.size();

    FloatType   fSum = 0;

    for (unsigned i = 0; i < N; ++i) {
        fSum += vVal[i];
    }

    FloatType   fMean = fSum / N;
    return fMean;
}

FloatType 
CalcSTD(const vector<FloatType>& vVal)
{
    unsigned    N = vVal.size();

    FloatType   fSum = 0;

    for (unsigned i = 0; i < N; ++i) {
        fSum += vVal[i];
    }

    FloatType   fMean = fSum / N;

    FloatType   fSqSum = 0;

    for (unsigned i = 0; i < N; ++i) {
        FloatType d = (vVal[i] - fMean);

        fSqSum += d*d;
    }

    return (fMean) / sqrt(fSqSum/(N - 1));
}
