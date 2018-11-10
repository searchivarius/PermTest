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
#ifndef PERMTEST_HPP
#define PERMTEST_HPP
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <queue>
#include <set>
#include <map>

#include <cmath>
#include <cstdlib>
#include <cctype>

#include "floattype.hpp"
#include "compstat.hpp"
#include "rand.hpp"

const unsigned OutputPrecision = 5;
const unsigned MaxClosedTestMethodQty = 12;

typedef unsigned long long MaskType;

using namespace std;

struct ToLower {
    char operator()(const char c) const {
        return tolower(c);
    }
};

inline void TransformToLower(std::string& Str) {
    ToLower Transformer;
    std::transform(Str.begin(), Str.end(), Str.begin(), Transformer);
}

template <typename T>
bool ReadVec(const string& line, vector<T>& v)
{
    v.clear();
    stringstream str(line);

    str.exceptions(ios::badbit);

    FloatType e;

    try {
        while (str >> e) {
            v.push_back(e);
        }
    } catch (const exception &e) {
        cerr << "Exception: " << e.what() << endl;
        cerr << "Failed to parse the line: '" << line << "'" << endl;
        return false;
    }
    return true;
}

class CMultiplePerm {
public:
    CMultiplePerm(unsigned VerbosityLevel = 0) : mVerbosityLevel(VerbosityLevel), mObservQty(0), mMethQty(0) {
        InitRand(0);
    }
    void AddObservRow(const string& line); 
    // A speed-optimized version for DoStat with StatObj = CSampleMean
    void DoMeanSpeedOptim(size_t PermQty, ostream& OutFile, 
                bool bDoAdjustment, bool bOnlyBaselineAdj, unsigned BaselineId);

    template <class TStatistic>

/*
 * Randomization (permutation) with adjustments for multiple comparisons.
 * The details of MaxT test can be found in:  
 *      1. S. Dudoit, Jp Schaffer, and Jc Boldrick. 
 *         Multiple hypothesis testing in microarray experiments. Statistical Science, 18(1):71–103, 2003.
 *
 *      2. WESTFALL, P. H. and YOUNG, S. S. (1993). Resampling-Based
 *         Multiple Testing: Examples and Methods for p-Value Adjustment.
 *         Wiley, New York. 
 */
   void DoStat(size_t PermQty, ostream& OutFile, 
                bool bDoAdjustment, bool bOnlyBaselineAdj, unsigned BaselineId,
                const TStatistic& StatObj);

/*
 *  The closed testing is a classic procedure due to:
 *      Marcus, R; Peritz, E; Gabriel, KR (1976). 
 *      "On closed testing procedures with special reference to ordered analysis of variance". 
 *      Biometrika 63: 655–660. doi:10.1093/biomet/63.3.655
 *
 *      NOTE: All closed testing procedures implemented in this file support only
 *            adjustment against the baseline.
 *      
 *  The closed testing procedure executes about 2^N permutation tests to verify all
 *  combination/intersection hypotheses that may include up to N methods.
 *  If you specify the maximum p-value this number can be reduced by 10-20%.
 *
 *  You can also approximate the closed-testing procedure by checking combinations of methods
 *  that contain at most M methods. The complexity of this approximate procedure is O(N^M).
 *  There is no guarantee that it provides a strong control of the Family Wise Error Rate (FWER).
 *  This approximation is less conservative, but it has a higher number of false positives.
 *
 */
    template <class TStatistic>
    void DoStatClosedTest(size_t PermQty, 
                unsigned MaxTestSubsetSize, FloatType MaxPValue,
                ostream& OutFile, 
                unsigned BaselineId,
                const TStatistic& StatObj);
    template <class TStatistic>
/*
 * This experimental procedure should control the FWER, but it may have less power than the classic closed testing method. 
 * It is a step-up closed-testing like procedure that stops when a certain number of hypotheses is tested.
 * Instead of specifying the exact # of tests, you can specify (approximately) the maximum size of the 
 * intersection hypotheses. If MaxApproxTestSubsetSize is M and the total # of methods is N, then
 * the amount of work is approximately O(M^N). You can change the amount of work by specifying the
 * ExtraWorkCoeff. If the value is larger than 1, more work is done.
 *
 * This procedure keeps track of the set of active hypotheses called "frontier". There is a subset of frontier
 * hypotheses that are tentatively rejected. When this procedure stops, it uses the Holm-Bonferroni method
 * to decide which tentatively rejected hypotheses should not be rejected. Non-rejection (i.e., "acceptance")
 * is simulated by assigning these hypotheses the p-value ImputedPValue.  We "accept" frontier hypotheses to
 * ensure that the p-value for the all-true-NULL hypothesis (the intersection/combination of all true NULL hypotheses )
 * is <= than MaxPValue
 */
    void DoStatClosedTestFrontier(size_t PermQty, 
                unsigned MaxApproxTestSubsetSize, FloatType ExtraWorkCoeff, 
                FloatType MaxPValue, FloatType ImputedPValue,
                ostream& OutFile, 
                unsigned BaselineId,
                const TStatistic& StatObj);
    unsigned GetMethQty() const { return mMethQty; }
protected:
/*
 * Randomization (permutation) test that verifies the complete NULL
 * hypothesis for a given subset of methods
 */
    template <class TStatistic>
    FloatType ComputeCompleteNullMethSubset(size_t PermQty, MaskType MethMask, const TStatistic& StatObj, unsigned BaselineId);
    void UpdatePVals(vector<vector<FloatType> >& vvPVals, MaskType MethMask, bool bOnlyBaselineAdj, unsigned BaselineId, FloatType pval);
    void OutputPVals(vector<vector<FloatType> >& vvPVals, bool bDoAdjustment, ostream& OutFile);
    void CheckMatrix(vector<vector<FloatType> >& vvPVals);
    unsigned GetMaxMethQty() const {
        return sizeof(MaskType) * 8;
    }
private:
    unsigned                    mVerbosityLevel;
    unsigned                    mObservQty;
    unsigned                    mMethQty;
    vector<vector<FloatType> >  mvvSourceData;
    vector<FloatType>           mvSum;

    void PrintIter(unsigned i) {
        if (mVerbosityLevel < 2) return;
        if (i % 1000 == 0) {
            cout << i << "\t";
            if (i % 10000 == 0) {
                cout << endl;
            }
            cout.flush();
        }
    }
};

struct CStatistics {
    unsigned    mMethI, mMethJ; // Computed for methods represented by indices mMethI and mMethJ
    FloatType   mT;
    CStatistics(unsigned MethI = 0, unsigned MethJ = 0, FloatType T = 0) :
                 mMethI(MethI), mMethJ(MethJ), mT(T) {}

    bool operator<(const CStatistics& That) const {
        return mT < That.mT; // Need to sort statistics value in ascending order
    }
};

inline ostream& operator<<(ostream& out, const CStatistics& obj) {
    return out << "["<< obj.mMethI << "," << obj.mMethJ << " " << setprecision(3) << obj.mT << "]";
}

inline unsigned 
GetSubsetQty(MaskType MethMask) {
    unsigned SubsetMethQty = 0;

    for (; MethMask; MethMask >>= 1) {
        if (1 & MethMask) {
            ++SubsetMethQty;
        }
    }

    return SubsetMethQty;
}

inline void
PrintSubset(MaskType MethMask) {
    unsigned SubsetMethQty = 0;

    cout << "[";
    for (unsigned i = 0; MethMask; ++i) {
        if (1 & MethMask) {
            cout << i << " ";
            ++SubsetMethQty;
        }
        MethMask >>= 1;
    }
    cout << "]";
}

template <class TStatistic>
void 
CMultiplePerm::DoStatClosedTest(size_t PermQty, 
            unsigned MaxTestSubsetSize, FloatType MaxPValue,
            ostream& OutFile, 
            unsigned BaselineId,
            const TStatistic& StatObj)
{
    vector<MaskType> vNonPromis;

    if (mMethQty < 2) {
        throw runtime_error("There should be at least 2 rows!");
    }

    if (BaselineId >= mMethQty) {
        cerr << "Baseline id: " << BaselineId << " is too high!" << endl;
        throw runtime_error("Baseline ID");
    }

    if (!MaxTestSubsetSize) MaxTestSubsetSize = mMethQty;

    if (MaxTestSubsetSize > MaxClosedTestMethodQty) {
        stringstream err;

        err << "The maximum allowed number of method combinations (" << MaxClosedTestMethodQty << ") was exceeded!";
        throw runtime_error(err.str());
    }

    if (mMethQty > GetMaxMethQty()) {
        stringstream err;

        err << "The maximum allowed number of methods (" << GetMaxMethQty() << ") was exceeded!";
        throw runtime_error(err.str());
    } 

    vector<vector<FloatType> >      vvPVals(mMethQty, vector<FloatType>(mMethQty));

    unsigned long long TestedSubsetQty = 0;
    for (MaskType MethMask = (MaskType(1) << mMethQty) - 1; MethMask; --MethMask) {
        if (((MaskType(1) << BaselineId) & MethMask) == 0) {
            continue; // In the scenario where we compare only against the baseline, 
                      // the baseline method should be present in every subset
        }
        unsigned SubsetMethQty = GetSubsetQty(MethMask);

        if (SubsetMethQty < 2) continue;
        if (SubsetMethQty > MaxTestSubsetSize) continue;

        if (mVerbosityLevel) {
            cout << hex << "Method mask: " << MethMask << dec << endl;
        }

        bool f = false;
        for (unsigned j = 0; j < vNonPromis.size(); ++j) {
            if ((vNonPromis[j] & MethMask) == MethMask) {
                f = true;
                break;
            }
        }
        if (f) {
            if (mVerbosityLevel) {
                cout << "Skipping combinations of methods (a containing hypothesis was not rejected)" << endl;
            }
            continue;
        }

        ++TestedSubsetQty;
        if (mVerbosityLevel) {
            cout << "# of tested method combinations: " << TestedSubsetQty << endl;
        }

        FloatType pval = ComputeCompleteNullMethSubset(PermQty, MethMask, StatObj, BaselineId);

        if (pval > MaxPValue && SubsetMethQty > 2) {
            vNonPromis.push_back(MethMask);
            if (mVerbosityLevel) {
                cout << "Containing hypothesis wasn't rejected, p-value=" << pval << " the total # of such hypotheses so far: " << vNonPromis.size() << endl;
            }
        }

        if (mVerbosityLevel) {
            cout << "p-val: " << setprecision(OutputPrecision) << pval << " for methods: " << endl;
        }

        UpdatePVals(vvPVals, MethMask, true, BaselineId, pval);
    }

    OutputPVals(vvPVals, true, OutFile);

    if (mVerbosityLevel) {
        cout << "# of tested method combination: " << TestedSubsetQty << endl;
    }
}

struct CStateBase {
    MaskType     mMask;
    FloatType    mPVal;
    unsigned     mMethQty;

    CStateBase(MaskType Mask, FloatType pVal) : mMask(Mask), mPVal(pVal), mMethQty(GetSubsetQty(Mask)) {}
};

template <class TCompare>
struct CState : CStateBase {
    TCompare     mComp;

    CState(MaskType Mask, FloatType pVal) : CStateBase(Mask, pVal) {}

    bool operator<(const CState& That) const {
        return mComp(*this, That);
    }
};

struct CCompFrontier {
    bool operator()(const CStateBase& e1, const CStateBase& e2) const {
#if 0
        // Adjust p-values by the number of methods, we need
        // to prioritize processing of subsets of hypotheses with few components.
        FloatType p1 = e1.mPVal/max(e1.mMethQty, 1U);
        FloatType p2 = e2.mPVal/max(e2.mMethQty, 1U);
#else
        // It appears that this approach slightly better
        unsigned m1 = e1.mMethQty;
        unsigned m2 = e2.mMethQty;
        if (m1 != m2) return m1 > m2;

        FloatType p1 = e1.mPVal;
        FloatType p2 = e2.mPVal;
#endif

        if (p1 != p2) return p1 < p2;
        return e1.mMask > e2.mMask;

    }
};

struct CCompPVal {
    bool operator()(const CStateBase& e1, const CStateBase& e2) const {
        FloatType p1 = e1.mPVal;
        FloatType p2 = e2.mPVal;

        if (p1 != p2) return p1 < p2;
        return e1.mMask > e2.mMask; 
    }
};

struct CCompPValRev {
    bool operator()(const CStateBase& e1, const CStateBase& e2) const {
        FloatType p1 = e1.mPVal;
        FloatType p2 = e2.mPVal;

        if (p1 != p2) return p1 > p2;
        return e1.mMask < e2.mMask; 
    }
};

typedef CState<CCompFrontier>   CStateFrontier;
typedef CState<CCompPVal>       CPValWrap;
typedef CState<CCompPValRev>    CPValRevWrap;

template <class TStatistic>
void 
CMultiplePerm::DoStatClosedTestFrontier(size_t PermQty, 
            unsigned MaxApproxTestSubsetSize, FloatType ExtraWorkCoeff,
            FloatType MaxPValue, FloatType ImputedPValue,
            ostream& OutFile, 
            unsigned BaselineId,
            const TStatistic& StatObj)
{

    if (mMethQty < 2) {
        throw runtime_error("There should be at least 2 rows!");
    }

    if (MaxApproxTestSubsetSize == 0) {
        throw runtime_error("The maxmeth should be > 0"); 
    }

    if (BaselineId >= mMethQty) {
        cerr << "Baseline id: " << BaselineId << " is too high!" << endl;
        throw runtime_error("Baseline ID");
    }

    if (mMethQty > GetMethQty()) {
        stringstream err;

        err << "The maximum allowed number of methods (" << GetMethQty() << ") was exceeded!";
        throw runtime_error(err.str());
    } 

    unsigned long long TestedSubsetQty = 0;

    FloatType   work = 0;

    for (unsigned i = 0; i <= MaxApproxTestSubsetSize; ++i) {
        FloatType r = 1;
        for (unsigned j = 0; j < i; ++j) {
            r *= FloatType(mMethQty - j) / (j + 1);
        }
        work += r;
    }

    unsigned long long MaxTestedSubsetQty = static_cast<unsigned long long>(ceil(work * ExtraWorkCoeff));

    if (mVerbosityLevel) {
        cout << "Testing at most " << MaxTestedSubsetQty << " combinations." << endl;
    }


    priority_queue<CStateFrontier>                  state;
    map<MaskType, FloatType>                        WasProcessed; 
    map<MaskType, CStateFrontier>                   Frontier;
    map<MaskType, CStateFrontier>::const_iterator   FrontierIt;
    vector<MaskType>                                vAcceptedFrontier;

    // Creating seed states
    for (MaskType MaskId = 0, MaskCurr = 1; MaskId < mMethQty; ++MaskId, MaskCurr <<= 1) {
        // Ensure that the baseline method is included into every intersection hypothesis,
        // if we are comparing only against the baseline.
        if (MaskId != BaselineId) continue;
        // Use the largest p-values to ensure that the seed states are processed
        // before any other states
        CStateFrontier  elem(MaskCurr, 1.0);
        state.push(elem); // Add to the stack, but not to the frontier!!!
    }
            

    vector<vector<FloatType> >      vvPVals(mMethQty, vector<FloatType>(mMethQty));
    
    while (!state.empty() && TestedSubsetQty < MaxTestedSubsetQty) {
        CStateFrontier  BaseState = state.top();
        state.pop();

        if (BaseState.mMethQty - 1 == MaxApproxTestSubsetSize) {
        // 1) Don't remove it from the frontier!!!
        // 2) Don't extend the frontier either
            continue;
        } else {
            Frontier.erase(BaseState.mMask);
        }

        if (mVerbosityLevel) {
            cout << "Base p-val: " << BaseState.mPVal << " Base subset qty: " << BaseState.mMethQty << endl;
        }

        for (MaskType MaskId = 0, MaskCurr = 1; MaskId < mMethQty; ++MaskId, MaskCurr <<= 1) {
            if (MaskCurr & BaseState.mMask) continue; // This method is already included in the set
            MaskType   MethMask = MaskCurr | BaseState.mMask;

            // Ignore already processed elements, they are all "under" the frontier
            if (WasProcessed.find(MethMask) != WasProcessed.end()) {
                continue; 
            }

            if (((MaskType(1) << BaselineId) & MethMask) == 0) {
                // In the scenario where we compare only against the baseline, 
                // the baseline method should be present in every subset,
                // because it is included through seeding.
                // If it is not here, something went wrong.
                throw runtime_error("Internal error: baseline id is not included!"); 
            }

            unsigned SubsetMethQty = GetSubsetQty(MethMask);

            if (SubsetMethQty < 2) {
                throw runtime_error("Internal error: should be at least 2 methods in the combination at this point!"); 
            }
            if (mVerbosityLevel) {
                cout << hex << "Method mask: " << MethMask << dec << " SubsetHypQty: " << (SubsetMethQty - 1) << endl;
            }


            bool f = false;
            for (unsigned j = 0; j < vAcceptedFrontier.size(); ++j) {
                if ((vAcceptedFrontier[j] & MethMask) == MethMask) {
                    f = true;
                    break;
                }
            }
            if (f) {
                if (mVerbosityLevel) {
                    cout << "Skipping combinations of methods (a containing hypothesis was not rejected)" << endl;
                }
                continue;
            }

            ++TestedSubsetQty;
            if (mVerbosityLevel) {
                cout << "Tested: " << TestedSubsetQty << " hypothesis combinations." << endl;
            }

            FloatType pval = ComputeCompleteNullMethSubset(PermQty, MethMask, StatObj, BaselineId);

            WasProcessed.insert(make_pair(MethMask, pval));

            if (pval > MaxPValue) {
                vAcceptedFrontier.push_back(MethMask);
                if (mVerbosityLevel) {
                    cout << "Containing hypothesis wasn't rejected, p-value=" << pval << " the total # of such hypotheses so far: " << vAcceptedFrontier.size() << endl;
                }
            } else {
                // Do not add an accepted hypothesis to the frontier or to the queue!!!
                CStateFrontier elem(MethMask, pval);

                state.push(elem);
                Frontier.insert(make_pair(elem.mMask, elem));
            }

            if (mVerbosityLevel) {
                cout << "p-val: " << setprecision(OutputPrecision) << pval << " for methods: " << endl;
            }

            UpdatePVals(vvPVals, MethMask, true, BaselineId, pval);
        }
    }

    unsigned MinFrontSubsetQty = mMethQty;
    vector<CPValRevWrap>   vTopPValues;
    
    for(FrontierIt = Frontier.begin(); FrontierIt != Frontier.end(); ++FrontierIt) {
        const CStateFrontier& elem = FrontierIt->second;

        if (elem.mPVal <= MaxPValue) {
            vTopPValues.push_back(CPValRevWrap(elem.mMask, elem.mPVal));
            MinFrontSubsetQty = min(elem.mMethQty, MinFrontSubsetQty);
        } else {
            throw runtime_error("Internal error: must not have accepted hypotheses here!");
        }
    }

    unsigned ActiveMethQty = mMethQty;

    // In the case of a baseline-only adjustement, a lot of methods would be inactive,
    // i.e., non-distinguishable from the baseline. We can use this fact to get
    // better upper bounds for methods such that hypotheses that they
    // are different from the baseline are not accepted yet. 
    //
    if (1) {
        // Baseline should be excluded from consideration!
        MinFrontSubsetQty--;
        ActiveMethQty--;
        for (unsigned i = 0; i < mMethQty; ++i) {
            if (i != BaselineId) {
                if (vvPVals[i][BaselineId] > MaxPValue) {
                    --ActiveMethQty; 
                }
            }
        }
    }

    // ceil of ActiveMethQty / MinFrontSubsetQty
    size_t TopK = (ActiveMethQty + MinFrontSubsetQty - 1)/ (MinFrontSubsetQty);

    // Now we have to proces frontier to estimate the maximum p-value of the all-true-null hypothesis combination.
    // This hypotheses can be represented by as an interesection of at most TopK hypothesis combinations from
    // the frontier. Using the Bonferroni method, one can estimate the p-value of the all-true-null hypothesis
    // by summing up TopK p-values of frontier combinations THAT ARE REJECTED. A slightly better approach
    // is to use the Holm-Bonferoni method.

    sort(vTopPValues.begin(), vTopPValues.end());

    if (mVerbosityLevel) {
        cout << "MinFrontSubsetQty: " << MinFrontSubsetQty << " TopK: " << TopK << " # of active methods: " << ActiveMethQty << endl;
        cout << "The # of elements on the frontier accepted already: " << vAcceptedFrontier.size() << endl;
        cout << "The # of elements on the frontier not yet accepted: " << vTopPValues.size() << endl;
    }

    size_t  skip = 0;


    if (mVerbosityLevel) {
        cout << "vTopPValues.size() = " << vTopPValues.size() << endl;
        cout << "Top p-values & covered methods: ";
    }
    for (; skip < vTopPValues.size(); ++skip) {
        size_t Qty = min(TopK, vTopPValues.size() - skip);

        if (mVerbosityLevel) {
            cout << vTopPValues[skip].mPVal << " ";
            cout << " ";
        }

#if 0
        long double sum = 0.0;
        // Bonferroni method
        for (size_t i = 0, Curr = skip + Qty - 1;  i < Qty; ++i, --Curr) {
            sum +=  vTopPValues[Curr].mPVal;
        }
        bool bRej = (sum <= MaxPValue);
        if (mVerbosityLevel) {
            cout << "Skip: " << skip << " Qty: " << Qty << " Sum: " << sum << endl;
        }
#else
        // Holm-Bonferroni method
        bool bRej = true;
        for (size_t i = 0, Curr = skip + Qty - 1;  i < Qty; ++i, --Curr) {
            if (vTopPValues[Curr].mPVal * (Qty - i) > MaxPValue) {
                bRej = false;
                break;
            }
        }
#endif

        if (bRej) { 
            if (mVerbosityLevel) {
                for (size_t i = skip; i < min(vTopPValues.size(), skip + TopK); ++i) {
                    cout << vTopPValues[i].mPVal << " ";
                    PrintSubset(vTopPValues[skip].mMask);
                }
            }
            break;
        }
    }

    if (mVerbosityLevel) {
        cout << endl;
        cout << "Skip: " << skip << endl;
    }


    // Forcing accept of hypotheses with too-large p-values
    for (unsigned i = 0; i < skip; ++i) {
        MaskType    MethMask = vTopPValues[i].mMask;
        UpdatePVals(vvPVals, MethMask, true, BaselineId, ImputedPValue);
    }

    if (mVerbosityLevel) {
        cout << "# of hypotheses accepted forcibly: " << skip << endl;
    }

    OutputPVals(vvPVals, true, OutFile);

    if (mVerbosityLevel) {
        cout << "# of tested method combinations: " << TestedSubsetQty << endl;
    }
}

// This one works apparently much better than the MAX-type aggregate stat. (see code below)
//#define SUM_AGGR_STAT 

template <class TStatistic>
FloatType
CMultiplePerm::ComputeCompleteNullMethSubset(size_t PermQty, MaskType MethMask, 
                                            const TStatistic& StatObj, unsigned ExternalBaselineId)
{
    unsigned SubsetMethQty = 0;
    unsigned BaselineId = 0;

    //ofstream OutStat("outstat.csv");

    if (mVerbosityLevel > 1) {
        cout << "Methods compared: ";
    }

    for (MaskType i = 0, Mask = 1; i < mMethQty; ++i, Mask <<= 1) {
        if (i == ExternalBaselineId) {
            BaselineId  = SubsetMethQty;
            if ((Mask & MethMask) == 0) {
                stringstream err;
                err << "Internal error: baseline (" << BaselineId << ") is not included into the mask!"; 
                throw runtime_error(err.str());
            }
        }
        if (Mask & MethMask) {
            if (mVerbosityLevel > 1) {
                cout << i << " ";
            }
            ++SubsetMethQty;
        }
    }

    if (mVerbosityLevel > 1) {
        cout << endl;
    }

    if (SubsetMethQty < 2) {
        throw runtime_error("There should be at least 2 rows!");
    }

    vector<vector<FloatType> >  vvSubsetSourceData(mObservQty, vector<FloatType>(SubsetMethQty));
    vector<vector<FloatType> >  vvSourceDataByRows(SubsetMethQty, vector<FloatType>(mObservQty));

    for (MaskType i = 0, Mask = 1, SubsetMethId = 0; i < mMethQty; ++i, Mask <<= 1) {
        if (Mask & MethMask) {
            for (unsigned j = 0; j < mObservQty; ++j) {
                vvSubsetSourceData[j][SubsetMethId] = 
                vvSourceDataByRows[SubsetMethId][j] = mvvSourceData[j][i];
            }
            ++SubsetMethId;
        }
    }

    vector<vector<FloatType> >      vvStatMatr(SubsetMethQty, vector<FloatType>(SubsetMethQty));

    if (mVerbosityLevel > 1) {
        cout << "Means:" << endl;
        for (unsigned i = 0; i < SubsetMethQty; ++i)  {
            cout << setw(5) << setprecision(3) << CalcMean(vvSourceDataByRows[i]) << "\t";
        }
        cout << endl;
    }

    if (mVerbosityLevel > 1) {
        cout << "STD dev:" << endl;

        for (unsigned i = 0; i < SubsetMethQty; ++i)  {
            cout << setw(5) << setprecision(3) << CalcSTD(vvSourceDataByRows[i])  << "\t";
        }
        cout << endl;
    }

    if (mVerbosityLevel > 1) {
        cout << "The statistics matrix:" << endl;
    }

    for (unsigned i = 0; i < SubsetMethQty; ++i)  {
        for (unsigned j = 0; j < SubsetMethQty; ++j)  {
            FloatType fStat = fabs(StatObj(vvSourceDataByRows[i], vvSourceDataByRows[j]));
            //FloatType fStat = StatObj(vvSourceDataByRows[i], vvSourceDataByRows[j]);

            vvStatMatr[i][j] = fStat;

            if (mVerbosityLevel > 1) {
                cout << setw(5) << setprecision(3) << fStat << "\t" ; 
            }
        }

        if (mVerbosityLevel > 1) {
            cout << endl;
        }
    }

    long double   fGlobAggrStat = 0;

    for (unsigned i = 0; i < SubsetMethQty; ++i)  {
        unsigned j = BaselineId; if (j != i) {
        //for (unsigned j = 0; j < i; ++j)  {
#ifdef SUM_AGGR_STAT
            fGlobAggrStat += (long double)vvStatMatr[i][j];
#else
            fGlobAggrStat = max(fGlobAggrStat, (long double)vvStatMatr[i][j]);
#endif
        }
    }

    fGlobAggrStat = fabs(fGlobAggrStat);

    //OutStat << fGlobAggrStat << endl;

    if (mVerbosityLevel > 1) {
        cout << "Aggregate statistic value: " << fGlobAggrStat << endl;
    }

    vector<vector<FloatType> >  vvCurrVals(SubsetMethQty, vector<FloatType>(mObservQty));
    vector<unsigned>            vIndex(SubsetMethQty);
    unsigned                    EqualOrGreater = 0;


    for (size_t nIter = 1; nIter <= PermQty; ++nIter) {
        PrintIter(nIter);

        for (unsigned i = 0; i < mObservQty; ++i) {
            // Compute a random matrix permutation (before computing statistic values)
            for (unsigned k = 0; k < SubsetMethQty; ++k)  {
                vIndex[k] = k;
            }
            const FloatType*    pData = &vvSubsetSourceData[i][0];
            for (unsigned k = SubsetMethQty; k; --k)  {
                unsigned kSwap = k > 1 ? GenRand() % k : 0;
                swap(vIndex[kSwap], vIndex[k - 1]);
                vvCurrVals[vIndex[k - 1]][i] = pData[k - 1];
            }
        }

        long double fAggrStat = 0;

        for (unsigned i = 0; i < SubsetMethQty; ++i)  {
            unsigned j = BaselineId; if (j != i) {
            //for (unsigned j = 0; j < i; ++j)  {
                FloatType Delta = fabs(StatObj(vvCurrVals[i], vvCurrVals[j]));
                //FloatType Delta = StatObj(vvCurrVals[i], vvCurrVals[j]);

#ifdef SUM_AGGR_STAT
                fAggrStat += (long double)Delta;
#else
                fAggrStat = max(fAggrStat, (long double)Delta);
#endif
            }
        }

        fAggrStat = fabs(fAggrStat);
        //OutStat << fAggrStat << endl;
        // Update the empirical distribution
        if (fAggrStat >= fGlobAggrStat) ++EqualOrGreater;
    }
    if (mVerbosityLevel > 1) {
        cout << endl;
    }


    FloatType fPValue = FloatType(EqualOrGreater) / PermQty;

    if (mVerbosityLevel > 1) {
        cout << "Aggregated p-value: " << setw(7) << setprecision(OutputPrecision) << fPValue << endl;
    }

    return fPValue;
}

template <class TStatistic>
void 
CMultiplePerm::DoStat(size_t PermQty, ostream& OutFile, 
            bool bDoAdjustment, bool bOnlyBaselineAdj, unsigned BaselineId,
            const TStatistic& StatObj)
{
    vector<vector<FloatType> >      vvStatMatr(mMethQty, vector<FloatType>(mMethQty));
    vector<vector<size_t> >         vvEqualOrGreaterQty(mMethQty, vector<size_t>(mMethQty));

    //ofstream OutStat("outstat.csv");

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

    vector<vector<FloatType> >      vvSourceDataByRows(mMethQty, vector<FloatType>(mObservQty));

    for (unsigned i = 0; i < mMethQty; ++i)  {
        for (unsigned j = 0; j < mObservQty; ++j) {
            vvSourceDataByRows[i][j] = mvvSourceData[j][i];
        }
    }

    if (mVerbosityLevel) {
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
            FloatType fStat = StatObj(vvSourceDataByRows[i], vvSourceDataByRows[j]);
            //if (j== 0 && i ==1) OutStat << fStat << endl;
            fStat = fabs(fStat);

            vvStatMatr[i][j] = fStat;


            if (mVerbosityLevel) {
                cout << setw(5) << setprecision(3) << fStat << "\t" ; 
            }

            if (j < i && bDoAdjustment) {
                if (!bOnlyBaselineAdj || BaselineId == i || BaselineId == j) {
                    vStat.push_back(CStatistics(i, j, fStat));
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

    vector<vector<FloatType> >  vvCurrVals(mMethQty, vector<FloatType>(mObservQty));
    vector<unsigned>            vIndex(mMethQty);


    for (size_t nIter = 1; nIter <= PermQty; ++nIter) {
        PrintIter(nIter);

        for (unsigned i = 0; i < mObservQty; ++i) {
            // Compute a random matrix permutation (before computing statistic values)
            for (unsigned k = 0; k < mMethQty; ++k)  {
                vIndex[k] = k;
            }
            const FloatType*    pData = &mvvSourceData[i][0];
            for (unsigned k = mMethQty; k; --k)  {
                unsigned kSwap = k > 1 ? GenRand() % k : 0;
                swap(vIndex[kSwap], vIndex[k - 1]);
                vvCurrVals[vIndex[k - 1]][i] = pData[k - 1];
            }
        }
        if (bDoAdjustment) {
            FloatType       Delta = 0;

            for (unsigned p = 0; p < vStat.size(); ++p) {
                unsigned i = vStat[p].mMethI;
                unsigned j = vStat[p].mMethJ;

                Delta = max(Delta, fabs(StatObj(vvCurrVals[i], vvCurrVals[j])));

                if (Delta >= vvStatMatr[i][j]) {
                    vvEqualOrGreaterQty[i][j]++;
                    vvEqualOrGreaterQty[j][i]++;
                }
            }
        } else {
            // Update the empirical distribution
            // These are unadjusted p-values!!!
            for (unsigned i = 0; i < mMethQty; ++i)  {
                for (unsigned j = 0; j < i; ++j)  {
                    FloatType Delta = StatObj(vvCurrVals[i], vvCurrVals[j]);
                    //if (j== 0 && i ==1) OutStat << Delta << endl;
                    Delta = fabs(Delta);

                    if (Delta >= vvStatMatr[i][j]) {
                        vvEqualOrGreaterQty[i][j]++;
                        vvEqualOrGreaterQty[j][i]++;
                    }
                }
            }
        }
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


#endif
