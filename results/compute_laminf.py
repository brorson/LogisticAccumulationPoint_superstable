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
        gmpy2.get_context().precision=2000
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

        
    def compute_lambdainf_aitken(self, x):
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

    def compute_lambdainf_theta(self, x):
        """
        Accelerates convergence using the theta algorithm.
        Reference:  https://core.ac.uk/download/pdf/82687053.pdf
        """
        M = len(x)
        ax = []
        for n in range(M-3):
            Snp1 = x[n+1]
            Snp2 = x[n+2]
            DSn = x[n+1]-x[n]
            DSnp1 = x[n+2]-x[n+1]
            DSnp2 = x[n+3]-x[n+2]
            D2Sn = (x[n+2]-x[n+1]) - (x[n+1]-x[n])
            D2Snp1 = (x[n+3]-x[n+2]) - (x[n+2]-x[n+1])

            num = Snp1*DSnp2*D2Sn - Snp2*DSn*D2Snp1
            denom = DSnp2*D2Sn - DSn*D2Snp1
            Tn = num/denom
            ax.append(Tn)

        return ax

    def compute_lambdainf_epsilon(self, x):
        """
        Accelerates convergence using the epsilon
        References:
          Convergence acceleration methods: The past decade -- Claude BREZINSKI
          The epsilon algorithm and related topics -- P.R. Graves-Morrisa, et al.
        """
        M = len(x)
        k = 0
        e = [ [ 0 for k1 in range(M) ] for k2 in range(M) ]

        for n in range(M):
            e[k][n] = x[n]
        #print(e)
        #print("---------------------------------------------------------------------------")

        k = 1
        for n in range(0, M-1):
            e[k][n] = 1.0/(e[k-1][n+1] - e[k-1][n])
        #print(e)
        #print("---------------------------------------------------------------------------")


        for k in range(2, M):
            for n in range(M-k):
                e[k][n] = e[k-2][n+1] + 1.0/(e[k-1][n+1] - e[k-1][n])
        #print("---------------------------------------------------------------------------")  
        #print(e)
        #print("---------------------------------------------------------------------------")  
        if ( (M % 2) == 1 ):
            # odd
            return e[M-1][0]
        else:
            # even
            return e[M-1][1]


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
#*************************************************************************
if (__name__ == "__main__"):
    seq = Sequence()


    print("============================================")
    # First read in files to analyze and compare.  Edit this
    # by hand.
    filename2 = "results_superstable_prec200.txt"


    # This is index of highest element of sequence to use.
    idx = -1

    # Look at last lambdas of each sequence to see how they match up.  
    # This verifies
    # both computations give the same result to some number 
    # of digits.  Also verify
    # last lambdas are less than lambda_inf.
    print("Checking last lambda from each file....")
    dict2 = seq.read_vals(filename2)
    x2 = seq.extract_lambdas(dict2)
    print("len(x2) = "+str(len(x2)))

    print("idx = "+str(idx)+", lambda[idx] from x2 = ")
    print(x2[idx])
    if (x2[idx] < seq.lambdainf):
        print("last lambda["+str(idx)+"] is less than lambdainf -- good!")
    else:
        print("!!! lambda["+str(idx)+"] is greater than lambdainf -- bad!!!")


    print("============================================")
    print("Checking lambdas computed via convergence acceleration....")
    # Now run a loop and do repeated convergence acceleration until I run
    # out of terms in my sequence.  The final result will be held in 
    # the last element of ax2.
    
    #-------------------------------------------------
    # Aitken accelerated convergence
    print("len(x2) = "+str(len(x2)))
    ax2 = seq.compute_lambdainf_theta(x2)
    M = len(ax2)
    while (M > 2):
        ax2 = seq.compute_lambdainf_aitken(ax2)
        M = len(ax2)
    # Now we're done.  Print out result.
    print("Lambda_inf computed from Aitken accelerated sequence ax2 = ")
    print(ax2[-1])   # Result after several accelerations
    
    #-------------------------------------------------
    # Theta accelerated convergence
    print("len(x2) = "+str(len(x2)))
    ax2 = seq.compute_lambdainf_theta(x2)
    M = len(ax2)
    while (M > 3):
        ax2 = seq.compute_lambdainf_theta(ax2)
        M = len(ax2)
    # Now we're done.  Print out result.
    print("Lambda_inf computed from theta accelerated sequence ax2 = ")
    print(ax2[-1])   # Result after several accelerations

    #-------------------------------------------------
    # epsilon accelerated convergence
    print("len(x2) = "+str(len(x2)))
    ax2 = seq.compute_lambdainf_epsilon(x2)
    print("Lambda_inf computed from epsilon accelerated sequence ax2 = ")
    print(ax2)

    print("============================================")
    # Check my results against the OEIS value.
    print("Compare my acclerated lambda_inf values with OEIS value....")
    print("lambdainf from OEIS = ")
    print(seq.lambdainf)
    print("compare lambda_inf computed from ax2 against seq.lambdainf:")
    seq.compare_lambdas(ax2, seq.lambdainf)

    print("============================================")
    # Now use my accelerated sequence to compute Feigenbaum's
    # delta using delta = (lam_n+1 - lam_inf)/(lam_n -lam_inf)
    # Use my lam_inf computed using Aitken acceleration.
    print("Compute delta using x2 sequence and my lambda_inf")
    lambdainf = ax2
    y2 = seq.compute_converging_diffs(x2, lambdainf)
    M = len(y2)
    #print("Converging diffs = ")
    #for i in range(M):
    #    print(y2[i])

    print("-------------------------------------------")
    deltas = seq.compute_converging_deltas(y2)
    #print("Computing deltas")
    #print("Converging deltas = ")
    #for i in range(M-1):
    #    print(deltas[i])
    
    print("-------------------------------------------")
    diff = seq.delta - deltas[-1]
    print("Difference of my computed delta from OEIS delta = ")
    print(diff)

    
    print("============================================")
    # Now use my accelerated sequence to compute Feigenbaum's
    # delta using delta = (lam_n+1 - lam_inf)/(lam_n -lam_inf)
    # Use the OEIS lam_inf this time.
    print("Compute delta using sequence x2 and the OEIS lambda_inf.")
    y2 = seq.compute_converging_diffs(x2, seq.lambdainf)
    M = len(y2)
    #print("Converging diffs = ")
    #for i in range(M):
    #    print(y2[i])

    print("-------------------------------------------")
    deltas = seq.compute_converging_deltas(y2)
    #print("Computing deltas")
    #print("Converging deltas = ")
    #for i in range(M-1):
    #    print(deltas[i])
    
    print("-------------------------------------------")
    diff = seq.delta - deltas[-1]
    print("Difference of my computed delta from OEIS delta = ")
    print(diff)

