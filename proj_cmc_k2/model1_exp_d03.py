# %%
import numpy as np
import os
import sys

cur_dir = os.path.dirname(__file__)
sys.path.append(cur_dir)

from RiemannianBCD import *
from PRW import ProjectedRobustWasserstein
from gen_data import *
from utils import InitialStiefel


n = 300
model = 'model1'
trans_type = 'exp'
m = 3
k = 2
p = 50
nb_exp = 10

a = (1./n) * np.ones(n)
b = (1./n) * np.ones(n)
eta = 1
tau = 0.01

# standard normal sequence
normal_standard = np.random.normal(0, 1, (n, m))

times = np.zeros((nb_exp, p))
proj_RBCD = np.zeros((p, m, m))
values_subspace_RBCD = np.zeros((nb_exp, p))

# Real optimal subspace
proj = np.zeros((m,m))  
proj[0, 0] = 1
proj[1, 1] = 1


for i in range(nb_exp):
    np.random.seed(i)
    params = {'n': n, 'p': p, 'm': m, 'k': k, 'model': model, 'trans_type': trans_type}
    X_oracle, X = generate_data(**params)
    X_transformed = np.zeros((n, m*p))

    U0 = InitialStiefel(m, k)

    for j in range(p):
        algo = RiemannianBlockCoordinateDescent(eta=eta, tau=None, max_iter=5000, threshold=0.005, verbose=True)
        PRW = ProjectedRobustWasserstein(X[:, (j*m):((j+1)*m)], normal_standard, a, b, algo, k)
        PRW.run('RBCD', tau, U0)

        values_subspace_RBCD[i, j] = np.linalg.norm(PRW.get_Omega() - proj)

        pi = PRW.get_pi()
        X_transformed[:, (j*m):((j+1)*m)] = pi.dot(normal_standard)*n

        proj_RBCD[j, :, :] += PRW.get_Omega() / nb_exp


U = np.zeros((p, m*k))
for j in range(p):
    eigen_values, eigen_vectors = np.linalg.eig(proj_RBCD[j, :, :] + proj_RBCD[j, :, :].T)

    idx = eigen_values.argsort()[::-1]
    eigen_values = eigen_values[idx]
    eigen_vectors = eigen_vectors[:, idx]

    U[j, :] = eigen_vectors[:, :k].flatten('F')


np.savetxt('proj_subspace_model1_exp_d03.csv', U, delimiter=',')

# %%
