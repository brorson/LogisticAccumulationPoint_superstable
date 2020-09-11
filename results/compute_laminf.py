#! /usr/bin/python3

import gmpy2
from gmpy2 import mpfr
import string

class Sequence():
    """
    This class reads in the lambda sequence contains methods
    used to verify the results and compute the accumulation 
    point.
    """

    def __init__(self):
        # Set number of bits to keep around
        gmpy2.get_context().precision=1500
        # These are values taken from the OEIS which I use when 
        # checking my results.
        self.lambdainf = mpfr("3.5699456718709449018420051513864989367638369115148323781079755299213628875001367775263210342163")
        self.delta = mpfr("4.669201609102990671853203820466201617258185577475768632745651343004134330211314")

    def read_vals(self, filename):
        """
        This reads in the lines in filename, and builds
        a dictionary of {order, lambda} pairs.
        """

        # Local dictionary to build.
        ldict = dict()

        # Read first line and discard it
        F = open(filename, "r")
        F.readline()

        # Now read in successive lines and parse them up.
        while True:
            line = F.readline()
            if not line: 
                break
            #print("Line: {}".format(line.strip()))

            ordstr, valstr = line.split(",")
            name, ord = ordstr.split("=")
            name, val = valstr.split("=")
            #print("ord = {}, val = {}".format(ord, val))

            ordf = int(ord.strip())
            valf = mpfr(val.strip())

            ldict.update({ordf: valf})
        return ldict

    def compare_results(self, dict1, dict2):
        """
        This looks for differences between the two computed
        lambda dicts.  
        """
        for key in sorted(dict1.keys()):
            # Only compare if keys exist in both dicts.
            if key in dict2.keys():
                v1 = dict1[key]
                v2 = dict2[key]
            
                diff = v2-v1
                print("order = {}, diff = ".format(key))
                print(diff)

    def extract_lambdas(self, ldict):
        x = []
        # Now add lambdas to vector x
        for key in sorted(ldict.keys()):
            x.append(ldict[key])
        return x

        
    def compute_lambdainf(self, x):
        """
        This uses Aitken's method to compute a high-res 
        value of the accumulation point lambda_inf.  This
        fcn applies the method once.  Call it in a loop
        to repeat the algo for better convergence.
        """
        M = len(x)
        ax = []
        for n in range(M-2):
            num = (x[n+2]-x[n+1])*(x[n+2]-x[n+1])
            denom = (x[n]-x[n+1]) - (x[n+1]-x[n+2])
            lam = x[n+2] - num/denom
            # print(lam)
            ax.append(lam)

        return ax

    def compare_lambdas(self, l1, l2):
        """
        This accepts two mpfr lambdas and compares them.
        """
        diff = l2-l1
        print("Difference = ")
        print(diff)

    def compute_converging_diffs(self, x, lambdainf):
        """
        This takes as input:
        x = sequence of lambda_n
        lambdainf = accumulation point lambda computed using 
        sequence acceleration.
        It returns the difference between each lambda and the accum point
        as a list.
        """
        y = []
        M = len(x)
        for i in range(M):
            y.append(x[i] - lambdainf)
        return y

    def compute_converging_deltas(self, x):
        """
        This takes a sequence of converging lambda diffs, and computes the
        ratio between successive diffs, which should converge to 
        Feigenbaum's delta.
        """
        y = []
        M = len(x)
        for i in range(M-1):
            y.append(x[i]/x[i+1]) 
        return y
        

#*************************************************************************
# 

if (__name__ == "__main__"):
    seq = Sequence()


    print("============================================")
    # First read in files to analyze and compare.  Edit this
    # by hand.
    filename1 = "results.txt"

    # This is index of highest element of sequence to use.
    idx = 26
    print("ax1 = "+filename1)


    # Look at last lambdas of each sequence to see how they match up.  
    # This verifies
    # both computations give the same result to some number 
    # of digits.  Also verify
    # last lambdas are less than lambda_inf.
    print("Checking last lambda from each file....")
    dict1 = seq.read_vals(filename1)
    x1 = seq.extract_lambdas(dict1)
    print("len(x1) = "+str(len(x1)))

    print("idx = "+str(idx)+", lambda[idx] from x1 = ")
    print(x1[idx])
    if (x1[idx] < seq.lambdainf):
        print("lambda["+str(idx)+"] is less than lambdainf -- good!")
    else:
        print("!!! lambda["+str(idx)+"] is greater than lambdainf -- bad!!!")

    # Truncate sequences to remove parts I don't want and to make sure
    # they are the same length.
    x1 = x1[0:idx]


    print("============================================")
    print("Checking lambdas computed via Aitken acceleration....")
    # Now run a loop and do repeated Aitken acceleration until I run
    # out of terms in my sequence.  The final result will be held in 
    # the last element of ax1 and ax2.
    print("len(x1) = "+str(len(x1)))
    ax1 = seq.compute_lambdainf(x1)
    M = len(ax1)
    while (M > 2):
        ax1 = seq.compute_lambdainf(ax1)
        M = len(ax1)
    # Now we're done.  Print out result.
    print("Lambda_inf computed from accelerated sequence ax1 = ")
    print(ax1[-1])   # Result after several accelerations

    print("============================================")
    # Check my results against the OEIS value.
    print("Compare my acclerated lambda_inf values with OEIS value....")
    print("lambdainf from OEIS = ")
    print(seq.lambdainf)
    print("compare lambda_inf computed from ax1 against seq.lambdainf:")
    seq.compare_lambdas(ax1[-1], seq.lambdainf)

    print("============================================")
    # Now use my accelerated sequence to compute Feigenbaum's
    # delta using delta = (lam_n+1 - lam_inf)/(lam_n -lam_inf)
    # Use my lam_inf computed using Aitken acceleration.
    print("Compute delta using x1 sequence and my lambda_inf")
    lambdainf = ax1[-1]
    y1 = seq.compute_converging_diffs(x1, lambdainf)
    M = len(y1)

    print("-------------------------------------------")
    deltas = seq.compute_converging_deltas(y1)
    
    print("-------------------------------------------")
    diff = seq.delta - deltas[-1]
    print("Difference of my computed delta from OEIS delta = ")
    print(diff)

    
    print("============================================")
    # Now use my accelerated sequence to compute Feigenbaum's
    # delta using delta = (lam_n+1 - lam_inf)/(lam_n -lam_inf)
    # Use the OEIS lam_inf this time.
    print("Compute delta using sequence x1 and the OEIS lambda_inf.")
    y1 = seq.compute_converging_diffs(x1, seq.lambdainf)
    M = len(y1)

    print("-------------------------------------------")
    deltas = seq.compute_converging_deltas(y1)
    
    print("-------------------------------------------")
    diff = seq.delta - deltas[-1]
    print("Difference of my computed delta from OEIS delta = ")
    print(diff)

