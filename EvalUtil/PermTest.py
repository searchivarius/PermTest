import numpy as np
from numpy import genfromtxt
import sys
import argparse
import subprocess
import tempfile


class PermTest(object):

    def __init__(self, inData, argDict):
        self.inData = inData
        self.args = argDict
        self.scores = tempfile.NamedTemporaryFile()
        self.pvalues = tempfile.NamedTemporaryFile()
        self.optList = ["binpath","verbose","nperm","atype","baseid","aproc",
                        "maxmeth","threshpval","maxpval","pvalimp","stat"]

    def prepInput(self):
        np.savetxt(self.scores.name, self.inData, delimiter='\t')

    def readOutpus(self):
        res = genfromtxt(self.pvalues.name, delimiter='\t')
        return res[:,:-1]

    def getArgs(self):
        specs = [self.args.get(self.optList[0])]
        for a in range(1,len(self.optList)):
            if self.args.get(self.optList[a]) is not None:
                specs.append("--"+self.optList[a])
                if isinstance(self.args.get(self.optList[a]),(int,float)):
                    specs.append(str(self.args.get(self.optList[a])))
                else:
                    specs.append(self.args.get(self.optList[a]))
        return specs

    def run(self):
        self.prepInput()
        specs = self.getArgs()
        specs.extend([self.scores.name, self.pvalues.name])
        subprocess.call(specs)
        return self.readOutpus()


def main(argv):
    parser = argparse.ArgumentParser(description='PermTest runner')
    parser.add_argument('--binpath', type=str,
                        default="permtest",
                        help='Path to binary executable file')
    parser.add_argument('--infile', type=str,
                        default="testInput.txt",
                        help='Path to input file')
    parser.add_argument('--outfile', type=str,
                        default="testOutput.txt",
                        help='Path to output file')
    parser.add_argument('--verbose', type=int,
                        default=0,
                        help='Verbosity level')
    parser.add_argument('--nperm', type=int,
                        default=100000,
                        help='Number of permutations')
    parser.add_argument('--atype', type=str,
                        choices = ["none", "pairwise", "baseline"],
                        default="baseline",
                        help='Multiplicity adjustment type. If baseline one must'
                             'specify baseline system row id.')
    parser.add_argument('--baseid', type=int,
                        default=0,
                        help='Baseline system row number')
    parser.add_argument('--aproc', type=str,
                        choices = ["maxt","cls","clsf"],
                        default="maxt",
                        help='Type of adjustment procedure')
    parser.add_argument('--maxmeth', type=int,
                        default=100,
                        help='Maximum number of system combinations (for cls|clsf only)')
    parser.add_argument('--threshpval', type=float,
                        default=0.01,
                        help='The minimal p-value for the adaptive approximation of the closed test')
    parser.add_argument('--maxpval', type=float,
                        default=0.50,
                        help='The maximal p-value for the adaptive approximation of the closed test')
    parser.add_argument('--pvalimp', type=float,
                        default=0.60,
                        help='The imputed p-value for the frontier-based method, '
                             'should be greater than maxpval')
    parser.add_argument('--stat', type=str,
                        choices = ["mean", "t", "fscore", "acc"],
                        default="t",
                        help='The statistic. NOTE: fscore and acc should be used for binary classification only. '
                             'The first row of the data always represents ground truth labels')

    args = parser.parse_args(argv)
    print(args)

    inData = genfromtxt(args.infile, delimiter='\t')
    test = PermTest(inData,vars(args))
    results = test.run()
    print(results)
    np.savetxt(args.outfile, results, delimiter='\t')


if __name__ == '__main__':
    main(sys.argv[1:])
