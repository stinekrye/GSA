def read_SA_LCP(fasta):
    import numpy as np
    data = np.loadtxt(f"SA_LCP_{fasta}.txt", dtype = int, skiprows = 1)
    SA = data[:,0]
    LCP = data[:,1]
    return SA,LCP