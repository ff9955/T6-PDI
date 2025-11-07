import numpy as np
from xsh_analysis_functions import build_sim_H, get_eigen
from file_parsers import xyz_parser

pseudoH_lines = xyz_parser('../test-eigen-order/run-fssh-0/run-pseudo-hamilt-1.xyz', 3)
hamiltonian = build_sim_H(pseudoH_lines, 420)

eigenvals, eigenvecs = np.linalg.eigh(hamiltonian)
idx = eigenvals.argsort()

eigenvals = eigenvals[idx]
eigenvecs = eigenvecs[:, idx]

xsh_eigenvals = np.loadtxt('../test-eigen-order/run-fssh-0/eigenvalues.txt')
xsh_eigenvecs = np.loadtxt('../test-eigen-order/run-fssh-0/eigenvectors.txt')

for j in range(len(eigenvecs)):
#    print(abs(xsh_eigenvals[j] - eigenvals[j]))

    xsh_ct_pop = np.sum(xsh_eigenvecs[j][:400]**2)
    python_ct_pop = np.sum(eigenvecs[:,j][:400]**2)

    print(abs(xsh_ct_pop - python_ct_pop))

