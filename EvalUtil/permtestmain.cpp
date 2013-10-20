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
#include <getopt.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <string>

#include "permtest.hpp"

#define ADJUST_NONE     "none"
#define ADJUST_PAIRWISE "pairwise"
#define ADJUST_BASELINE "baseline"


#define PROC_NONE           "none"
#define PROC_MAXT           "maxt"
#define PROC_CLOSED         "cls"
#define PROC_CLOSED_FRONT   "clsf"

#define STAT_MEAN     "mean"
#define STAT_TSTAT    "t"
#define STAT_FSCORE   "fscore"
#define STAT_ACC      "acc"

const size_t DefaultPermNum     = 100000;
const string DefaultAdjustType  = ADJUST_BASELINE;
const string DefaultProcType    = PROC_NONE;
const string DefaultStatType    = STAT_TSTAT;

size_t      PermNum             = DefaultPermNum;
string      AdjustType          = DefaultAdjustType;
string      ProcType            = DefaultProcType;
string      StatType            = DefaultStatType;
unsigned    VerbosityLevel      = 0;

void Usage(const char *p, const char *err) {
    if (err) {
        cerr << "ERROR: " << err << endl;
    }
    cerr << "Usage: " << p << " <input file> <output file> " << endl << 
                              "\t\t\t[-v --verbose: <verbosity level>] [-h --help]" << endl << 
                              "\t\t\t[-n --nperm <# of permutations> default:" << DefaultPermNum << "] " << endl << 
                              "\t\t\t[-a --atype <adjust type: " ADJUST_NONE ", "ADJUST_PAIRWISE ", "ADJUST_BASELINE "> default: " << DefaultAdjustType << "] " << endl << 
                              "\t\t\t[-b --baseid <baseline id>] " <<  endl << 
                              "\t\t\t\tIf atype = baseline, baseline id should be specified!" << endl << 
                              "\t\t\t[-r --aproc (" << PROC_MAXT << "|" << PROC_CLOSED << "|" << PROC_CLOSED_FRONT << ") ]"<<  endl << 
                              "\t\t\t\tA type of the adjustement procedure:" << endl <<
                              "\t\t\t\t\t" << PROC_MAXT << "\t\tMaxT" << endl <<
                              "\t\t\t\t\t" << PROC_CLOSED << "\t\tClosed testing" << endl <<
                              "\t\t\t\t\t" << PROC_CLOSED_FRONT "\t\tExperimental frontier-based closed testing (must specify maxmeth)" << endl <<
                              "\t\t\t\t\tNOTE: closed testing procedures provide only baseline adjustment" << endl <<
                              "\t\t\t[-m --maxmeth ]" << endl << 
                              "\t\t\t\tClosed testing: the maximum # of method combinations" << endl <<
                              "\t\t\t\tFrontier-based closed testing: an approximate upper bound to determine the stopping point" <<  endl << 
                              "\t\t\t[-t --threshpval <the minimal p-value for the adaptive approximation of the closed test> ] " <<  endl << 
                              "\t\t\t[-p --maxpval <the maximum p-value, specifying maxpval allows one to decrease run-time by 10-20\% sometimes> ] " <<  endl << 
                              "\t\t\t[-i --pvalimp <the imputed p-value for the frontier-based method, should be greater than maxpval> ] " <<  endl << 
                              "\t\t\t[-s --stat <statistic: " STAT_MEAN ", " STAT_TSTAT ", " STAT_FSCORE ", " STAT_ACC " > default: " << DefaultStatType << "] " <<  endl << 
                              "\t\t\t\tNOTE: " STAT_FSCORE " and " STAT_ACC  " should be used for binary classification only!" << endl <<
                              "\t\t\t\t      The first row always represents ground truth labels." << endl <<
                              "\t\t\t\t      Other rows contains labels selected by classifiers," << endl <<
                              "\t\t\t\t      (each row corresponds to an outcome of a single classifier)." << endl <<
                                endl;
    exit(1);
};

static struct option aLongOpts[] = {
    {"nperm",       1, 0, 'n'},
    {"atype",       1, 0, 'a'},
    {"aproc",       1, 0, 'r'},
    {"baseid",      1, 0, 'b'},
    {"stat",        1, 0, 's'},
    {"verbose",     1, 0, 'v'},
    {"help",        0, 0, 'h'},
    {"maxmeth",     1, 0, 'm'},
    {"maxpval",     1, 0, 'p'},
    {"pvalimp",     1, 0, 'i'},
    {"threshpval",  1, 0, 't'},
    {"extrworkcoeff",1, 0, 'e'},
    {0,             0, 0,  0 },
};

enum EAvailStatType {
    kSampleMean,
    kTStat,
    kFScore,
    kAccuracy
};

string GetStatName(EAvailStatType type) {
  switch (type) {
    case kSampleMean: return STAT_MEAN;
    case kTStat:      return STAT_TSTAT;
    case kFScore:     return STAT_FSCORE;
    case kAccuracy:   return STAT_ACC;
    default: break; 
  }
  cerr << "Unknown statistic type: " << type << " ... aborting" << endl;
  exit(1);
  return "";
}

enum EProcType {
    kProcNone,
    kProcMaxT,
    kProcClosed,
    kProcClosedFrontier,
};

int main(int argc, char * pArgv[]) {
    bool            bBaselineSpec = false; 
    unsigned        BaselineId = 0;
    int             OptIndex;

    EAvailStatType  eStatType = kTStat;
    EProcType       eProcType = kProcNone;

    unsigned        MaxTestSubsetSize = 0;
    FloatType       ImputedPValue   = 1.0;
    FloatType       MaxPValue       = 1.0;
    FloatType       CutoffPValue    = 0.0;
    FloatType       ExtraWorkCoeff  = 1.0;

    while (true) {
        // The colon *AFTER* the option name means that there should be an argument, e.g., -n 10000
        int c = getopt_long(argc, pArgv, "n:a:r:b:s:v:m:p:t:e:i:h", aLongOpts, &OptIndex);
        if (c == -1) break;
        switch (c) {
            case 'v':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> VerbosityLevel)) Usage(pArgv[0], ("Invalid verbosity level: " + string(optarg)).c_str());
                     }
                     
                     break;
            case 'h':Usage(pArgv[0], NULL);
                     break;
            case 'm':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> MaxTestSubsetSize)) Usage(pArgv[0], ("Invalid maxmeth: " + string(optarg)).c_str());
                     }
                     break;
            case 'p':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> MaxPValue)) Usage(pArgv[0], ("Invalid maxpval: " + string(optarg)).c_str());
                     }
                     break;
            case 'i':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> ImputedPValue)) Usage(pArgv[0], ("Invalid pvalimp: " + string(optarg)).c_str());
                     }
                     break;
            case 'e':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> ExtraWorkCoeff)) Usage(pArgv[0], ("Invalid extrworkcoeff: " + string(optarg)).c_str());
                     }
                     break;
            case 't':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> CutoffPValue)) Usage(pArgv[0], ("Invalid threshpval: " + string(optarg)).c_str());
                     }
                     break;
            case 'n':assert(optarg); {
                        stringstream  str(optarg);
                        if (!(str >> PermNum)) Usage(pArgv[0], ("Invalid # of permutations: " + string(optarg)).c_str());
                     }
                     break;
            case 'a':assert(optarg);
                    AdjustType = optarg;
                    TransformToLower(AdjustType);
                    break;
            case 'r':assert(optarg);
                    ProcType = optarg;
                    TransformToLower(ProcType);
                    break;
            case 's':assert(optarg);
                    StatType = optarg;
                    TransformToLower(StatType);
                    break;
            case 'b':assert(optarg); {
                        stringstream str(optarg);
                        if (!(str >> BaselineId)) Usage(pArgv[0], ("Invalid baseline id: " + string(optarg)).c_str());
                        bBaselineSpec = true;
                    }
                    break;
            default: Usage(pArgv[0], NULL);
        }

    }

    bool        bDoAdjustment       = false;
    bool        bOnlyBaselineAdj    = false;


    if (ADJUST_NONE == AdjustType) {
    } else if (ADJUST_PAIRWISE == AdjustType) {
        bDoAdjustment       = true;
        bOnlyBaselineAdj    = false;
    } else if (ADJUST_BASELINE == AdjustType) {
        bDoAdjustment       = true;
        bOnlyBaselineAdj    = true;
        if (!bBaselineSpec) {
            Usage(pArgv[0], "Missing baseline id");
        }
    } else {
        Usage(pArgv[0], ("Invalid adjustment type: " + AdjustType).c_str());
    }

    if (PROC_NONE == ProcType) {
    } else if (PROC_MAXT == ProcType) {
        eProcType = kProcMaxT;
    } else if (PROC_CLOSED == ProcType) {
        eProcType = kProcClosed;
    } else if (PROC_CLOSED_FRONT == ProcType) {
        eProcType = kProcClosedFrontier;
        if (!MaxTestSubsetSize) {
            Usage(pArgv[0], "Invalid maxmeth, should be > 0 for this type of procedure");
        }
    } else {
        Usage(pArgv[0], ("Invalid procedure type: " + ProcType).c_str());
    }

    if (STAT_MEAN == StatType) {
        eStatType = kSampleMean;
    } else if (STAT_TSTAT == StatType) {
        eStatType = kTStat;
    } else if (STAT_FSCORE == StatType) {
        eStatType = kFScore;
    } else if (STAT_ACC == StatType) {
        eStatType = kAccuracy;
    } else {
        Usage(pArgv[0], ("Invalid statistic type: " + StatType).c_str());
    }

    if (optind + 1 >= argc) {
        Usage(pArgv[0], "Both input and output files should be specified");
    }
    const char*     pInputFileName = pArgv[optind];
    const char*     pOutputFileName = pArgv[optind + 1];

    cout << "Input file     : '" << pInputFileName << "'" << endl;
    cout << "Output file    : '" << pOutputFileName << "'" << endl;
    cout << "Permutation #  : " << PermNum << endl;
    cout << "Adjust for mt? : " << (bDoAdjustment ? "yes":"no") << endl;
    if (bDoAdjustment) {
        cout << "Adjust type    : " << (bOnlyBaselineAdj ? "against baseline" : "pairwise") << endl;
        if (bOnlyBaselineAdj) {
            cout << "Baseline id    : " << BaselineId << endl;
        }
        cout << "Procedure      : " << ProcType << endl; 
        if (MaxPValue < 1) {
            cout << "Maximum p-value: " << MaxPValue << endl;
        }
        if (ImputedPValue < 1) {
            cout << "p-value for imputation (frontier-based method): " << ImputedPValue << endl;
        }
        if (MaxTestSubsetSize) {
            cout << "Maxmeth: " << MaxTestSubsetSize << endl;
        }
        if (CutoffPValue) {
            cout << "Threshold p-value: " << CutoffPValue << endl;
        }
        if (ExtraWorkCoeff != 1.0) {
            cout << "Extra work coeff.: " << ExtraWorkCoeff << endl;
        }
    }
    cout << "Statistic      : " << GetStatName(eStatType) << endl;
    cout << "Verbosity level : " << VerbosityLevel << endl;
    

    ifstream        InpFile(pInputFileName); 

    if (!InpFile) {
        cerr << "Cannot open file: '" << pInputFileName << "'" << endl;
        return 1;
    }

    InpFile.exceptions(ios::badbit);

    ofstream        OutFile(pOutputFileName, ios::trunc | ios::out);

    if (!OutFile) {
        cerr << "Cannot create output file: '" << pOutputFileName << "'" << endl;
        return 1;
    }

    OutFile.exceptions(ios::badbit | ios::failbit | ios::eofbit);

    CMultiplePerm       Perm(VerbosityLevel);


    try {
        string line;

        bool bFirst = true;

        vector<FloatType> labels;

        while (getline(InpFile, line)) {
            if ((eStatType == kAccuracy || eStatType == kFScore) && bFirst) {
              if (!ReadVec(line, labels)) {
                cerr << "Cannot read labels (from the first line)" << endl;
                return -1;
              }
            } else {
              Perm.AddObservRow(line);
            }
            bFirst = false;
        }

        if (Perm.GetMethQty() <= BaselineId && bOnlyBaselineAdj) {
            Usage(pArgv[0], "Invalid BaselineId (>= # of rows)");
        }

        if (bDoAdjustment && kProcNone == eProcType) {
            Usage(pArgv[0], "You should specify the adjustment procedure");
        }

        if (!bDoAdjustment) {
            eProcType = kProcNone;
        }

        if (kProcClosed == eProcType ||
            kProcClosedFrontier == eProcType) {
            if (!bOnlyBaselineAdj) {
                Usage(pArgv[0], "Closed testing procedures provide adjustments only against a baseline!");
            }
        }
        
        if (kProcClosed == eProcType) {
            if (kSampleMean == eStatType) {
                Perm.DoStatClosedTest(PermNum, MaxTestSubsetSize, MaxPValue, OutFile, BaselineId, CSampleMean());
            } else {
                Perm.DoStatClosedTest(PermNum, MaxTestSubsetSize, MaxPValue, OutFile, BaselineId, CTStat());
            }
        } else if (kProcClosedFrontier == eProcType) {
            if (ImputedPValue < 1 && ImputedPValue <= MaxPValue) {
                Usage(pArgv[0], "Invalid pvalimp, it should be > than pval)");
            }
            if (kSampleMean == eStatType) {
                Perm.DoStatClosedTestFrontier(PermNum, MaxTestSubsetSize, ExtraWorkCoeff, MaxPValue, ImputedPValue, OutFile, BaselineId, CSampleMean());
            } else {
                Perm.DoStatClosedTestFrontier(PermNum, MaxTestSubsetSize, ExtraWorkCoeff, MaxPValue, ImputedPValue, OutFile, BaselineId, CTStat());
            }
        } else if (kProcMaxT == eProcType || !bDoAdjustment) {
            if (kSampleMean == eStatType) {
                Perm.DoMeanSpeedOptim(PermNum, OutFile, bDoAdjustment, bOnlyBaselineAdj, BaselineId);
                //Perm.DoStat(PermNum, OutFile, bDoAdjustment, bOnlyBaselineAdj, BaselineId, CSampleMean());
            } else if (kTStat == eStatType) {
                Perm.DoStat(PermNum, OutFile, bDoAdjustment, bOnlyBaselineAdj, BaselineId, CTStat());
            } else if (kFScore == eStatType) {
cout << "here" << endl;
                Perm.DoStat(PermNum, OutFile, bDoAdjustment, bOnlyBaselineAdj, BaselineId, CFScoreStat(labels));
            } else if (kAccuracy == eStatType) {
cout << "here" << endl;
                Perm.DoStat(PermNum, OutFile, bDoAdjustment, bOnlyBaselineAdj, BaselineId, CAccuracyStat(labels));
            }
        } else {
            Usage(pArgv[0], "Invalid option combination.");
        }

    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        cerr << "Failed to process the file: '" << pInputFileName << "'" << endl;
        return 1;
    }

    return 0;
}
