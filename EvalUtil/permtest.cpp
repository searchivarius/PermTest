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
#include <algorithm>
#include <limits>

#include "permtest.hpp"

void
CMultiplePerm::AddObservRow(const string& line) {
    vector<FloatType> CurrRow;

    if (!ReadVec(line, CurrRow)) {
        throw runtime_error("Parsing failed");
    }

    if (mObservQty) {
        if (mObservQty != CurrRow.size()) {
            stringstream err;
            err << "Error in line# " << (mMethQty + 1) << " the number of fields: " << CurrRow.size() <<
                   " does not match the number of fields in previous rows: " << mObservQty;
            throw runtime_error(err.str());
        }
    } else {
        mObservQty = CurrRow.size();
        mvvSourceData.resize(mObservQty);
    }
    ++mMethQty;
    FloatType sum = 0.0; 

    for (unsigned i = 0; i < mObservQty; ++i) {
        mvvSourceData[i].push_back(CurrRow[i]);
        sum += CurrRow[i];
    }
    mvSum.push_back(sum);
};


void 
CMultiplePerm::DoMeanSpeedOptim(size_t PermQty, ostream& OutFile,
                  bool bDoAdjustment, bool bOnlyBaselineAdj, const unsigned BaselineId)
{
    vector<vector<FloatType> >      vvDistMatr(mMethQty, vector<FloatType>(mMethQty));
    vector<vector<size_t> >         vvEqualOrGreaterQty(mMethQty, vector<size_t>(mMethQty));

    if (mMethQty < 2) {
        throw runtime_error("There should be at least 2 rows!");
    }

    if (BaselineId >= mMethQty) {
        cerr << "Baseline id: " << BaselineId << " is too high!" << endl;
        throw runtime_error("Baseline ID");
    }

    vector<CStatistics> vStat;

    if (mVerbosityLevel) {
        cout << "Means:" << endl;
    
        for (unsigned i = 0; i < mMethQty; ++i)  {
            cout << setw(5) << setprecision(3) << mvSum[i] / mObservQty << "\t";
        }
        cout << endl;
    }

    if (mVerbosityLevel) {
        vector<vector<FloatType> >      vvSourceDataByRows(mMethQty, vector<FloatType>(mObservQty));

        for (unsigned i = 0; i < mMethQty; ++i)  {
            for (unsigned j = 0; j < mObservQty; ++j) {
                vvSourceDataByRows[i][j] = mvvSourceData[j][i];
            }
        }

        cout << "STD dev:" << endl;

        for (unsigned i = 0; i < mMethQty; ++i)  {
            cout << setw(5) << setprecision(3) << CalcSTD(vvSourceDataByRows[i])  << "\t";
        }
        cout << endl;
    }

    if (mVerbosityLevel) {
        cout << "The statistics matrix:" << endl;
    }

    for (unsigned i = 0; i < mMethQty; ++i)  {
        for (unsigned j = 0; j < mMethQty; ++j)  {
            FloatType   Delta = vvDistMatr[i][j] = fabs(mvSum[i] - mvSum[j]);

            if (mVerbosityLevel) {
                cout << setw(5) << setprecision(3) << Delta / mObservQty << "\t" ; 
            }

            if (j < i && bDoAdjustment) {
                if (!bOnlyBaselineAdj || BaselineId == i || BaselineId == j) {
                    vStat.push_back(CStatistics(i, j, Delta / mObservQty));
                }
            }
        }

        if (mVerbosityLevel) {
            cout << endl;
        }
    }

    sort(vStat.begin(), vStat.end());

    if (mVerbosityLevel) {
        cout << "Statistics values:" << endl;
        for (unsigned i = 0; i < vStat.size(); ++i ) {
            cout << vStat[i] << " ";
        }
        cout << endl;
    }

    vector<FloatType>   vCurrSum(mMethQty);
    vector<unsigned>    mIndex(mMethQty);


    for (size_t nIter = 0; nIter < PermQty; ++nIter) {
        PrintIter(nIter);
        for (unsigned i = 0; i < mObservQty; ++i) {
            // Compute a random matrix permutation (before computing statistic values)
            for (unsigned k = 0; k < mMethQty; ++k)  {
                mIndex[k] = k;
            }
            const FloatType*    pData = &mvvSourceData[i][0];
            for (unsigned k = mMethQty; k; --k)  {
                unsigned kSwap = k > 1 ? GenRand() % k : 0;
                swap(mIndex[kSwap], mIndex[k - 1]);
                vCurrSum[mIndex[k - 1]] += pData[k - 1];
            }

        }
        if (bDoAdjustment) {
            FloatType       Delta = 0;

            for (unsigned p = 0; p < vStat.size(); ++p) {
                unsigned i = vStat[p].mMethI;
                unsigned j = vStat[p].mMethJ;
                Delta = max(Delta, fabs(vCurrSum[i] - vCurrSum[j]));

                if (Delta >= vvDistMatr[i][j]) {
                    vvEqualOrGreaterQty[i][j]++;
                    vvEqualOrGreaterQty[j][i]++;
                }
            }
        } else {
            // Update the empirical distribution
            // These are unadjusted p-values!!!
            for (unsigned i = 0; i < mMethQty; ++i)  {
                for (unsigned j = 0; j < i; ++j)  {
                    FloatType Delta = fabs(vCurrSum[i] - vCurrSum[j]);
                    if (Delta >= vvDistMatr[i][j]) {
                        vvEqualOrGreaterQty[i][j]++;
                        vvEqualOrGreaterQty[j][i]++;
                    }
                }
            }
        }

        fill(vCurrSum.begin(), vCurrSum.end(), 0);
    }
    cout << endl;

    vector<vector<FloatType> >      vvPVals(mMethQty, vector<FloatType>(mMethQty));
    for (unsigned i = 0; i < mMethQty; ++i)  {
        for (unsigned j = 0; j < mMethQty; ++j)  {
            vvPVals[i][j] = FloatType(vvEqualOrGreaterQty[i][j])/PermQty;
        }
    }

    
    if (bDoAdjustment) {
        // Enforce monotonicity
        FloatType   Prev = 0.0;
        for (unsigned p = vStat.size(); p ; --p) {
            unsigned i = vStat[p - 1].mMethI;
            unsigned j = vStat[p - 1].mMethJ;

            vvPVals[j][i] = vvPVals[i][j] = Prev = max(Prev, vvPVals[i][j]);
        }
    }

    OutputPVals(vvPVals, bDoAdjustment, OutFile);
}

void
CMultiplePerm::CheckMatrix(vector<vector<FloatType> >& vvPVals)
{
    if (vvPVals.size() != mMethQty) {
        stringstream err;

        err << "Internal error, the # of methods is: " << mMethQty << ", but the # of rows in the matrix is: " << vvPVals.size();
        throw runtime_error(err.str());
    }
    if (!vvPVals.size()) {
        throw runtime_error("Internal error, zero # of rows in the matrix");
    }
    if (vvPVals[0].size() != mMethQty) {
        stringstream err;

        err << "Internal error, the # of methods is: " << mMethQty << ", but the # of columns in the matrix is: " << vvPVals.size();
        throw runtime_error(err.str());
    }
}

void
CMultiplePerm::UpdatePVals(vector<vector<FloatType> >& vvPVals, MaskType MethMask, bool bOnlyBaselineAdj, unsigned BaselineId, FloatType pval)
{
    if (mVerbosityLevel) {
        cout << "Set p-value " << pval << " for : " << endl;
    }

    for (MaskType i1 = 0, Mask1 = 1; i1 < mMethQty; ++i1, Mask1 <<= 1) {
        if ((Mask1 & MethMask) == 0) continue;
        if (mVerbosityLevel) {
            cout << i1 << " ";
        }
        for (MaskType i2 = 0, Mask2 = 1; i2 < mMethQty; ++i2, Mask2 <<= 1) {
            if (bOnlyBaselineAdj && i1 != BaselineId && i2 != BaselineId) continue;
            if (i1 != i2 && (Mask2 & MethMask)) {
                FloatType p = max(vvPVals[i1][i2], pval);

                vvPVals[i1][i2] = vvPVals[i2][i1] = p;
            }
        }
    }

    if (mVerbosityLevel) {
        cout << endl;
    }
}

void
CMultiplePerm::OutputPVals(vector<vector<FloatType> >& vvPVals, bool bDoAdjustment, ostream& OutFile) 
{
    if (mVerbosityLevel) {
        cout << "The " << (bDoAdjustment ? "adjusted":"unadjusted") << " p-value matrix:" << endl;
    }

    for (unsigned i = 0; i < mMethQty; ++i)  {
        for (unsigned j = 0; j < mMethQty; ++j)  {
            float fPValue = vvPVals[i][j];

            if (mVerbosityLevel) {
                cout << setw(7) << setprecision(OutputPrecision) << fPValue << "\t" ; 
            }
            OutFile << setprecision(OutputPrecision) <<  fPValue << "\t" ; 
        }
        if (mVerbosityLevel) {
            cout << endl;
        }
        OutFile << endl;
    }
}
