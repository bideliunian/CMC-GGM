# %%
import numpy as np
from scipy.stats import norm
from scipy.stats import qmc

from Optimization.RiemannianBCD import RiemannianBlockCoordinateDescent
from Optimization.RiemannianGAS import RiemannianGradientAscentSinkhorn
from PRW import ProjectedRobustWasserstein
from gen_data import *
from scipy.stats import qmc

n = 200
n_halton = 200
model = 'model1'
trans_type = 'exp'
m = 10
k = 2
p = 30
nb_exp = 10

a = (1./n) * np.ones(n)
b = (1./n_halton) * np.ones(n_halton)
eta = 1
tau = 0.01

# halton sequence
halton_sampler = qmc.Halton(d=m, scramble=False)
halton_sample = halton_sampler.random(n=n_halton)
thresh = 1 / (4 * (np.power(n, 0.25) * np.sqrt(3.14 * np.log(n))))
halton_sample[halton_sample <= thresh] = thresh
halton_sample[halton_sample >= (1-thresh)] = 1 - thresh
normal_standard = norm.ppf(halton_sample)

times = np.zeros((nb_exp, p))
proj_RBCD = np.zeros((p, m, m))
proj_RGAS = np.zeros((p, m, m))
values_subspace_RBCD = np.zeros((nb_exp, p))
values_subspace_RGAS = np.zeros((nb_exp, p))

proj = np.zeros((m,m))   # Real optimal subspace
proj[0, 0] = 1
proj[1, 1] = 1


for i in range(nb_exp):
    np.random.seed(i)
    params = {'n': n, 'p': p, 'm': m, 'k': k, 'model': model, 'trans_type': trans_type, 'B_type':1}
    X_oracle, X = generate_data(**params)
    X_transformed = np.zeros((n, m*p))

    U0 = InitialStiefel(m, k)

    for j in range(p):
        algo = RiemannianBlockCoordinateDescent(eta=eta, tau=None, max_iter=3000, threshold=0.01, verbose=True)
        PRW = ProjectedRobustWasserstein(X[:, (j*m):((j+1)*m)], normal_standard, a, b, algo, k)
        PRW.run('RBCD', tau, U0)

        algo1 = RiemannianGradientAscentSinkhorn(eta=eta, tau=None, max_iter=3000, threshold=0.01,
                                                 sink_threshold=1e-4, verbose=True)
        PRW1 = ProjectedRobustWasserstein(X[:, (j*m):((j+1)*m)], normal_standard, a, b, algo1, k)
        PRW1.run('RGAS', tau/eta, U0)

        values_subspace_RBCD[i, j] = np.linalg.norm(PRW.get_Omega() - proj)
        values_subspace_RGAS[i, j] = np.linalg.norm(PRW1.get_Omega() - proj)

        pi = PRW.get_pi()
        X_transformed[:, (j*m):((j+1)*m)] = pi.dot(normal_standard)*n

        proj_RBCD[j, :, :] += PRW.get_Omega() / nb_exp
        proj_RGAS[j, :, :] += PRW1.get_Omega() / nb_exp


U = np.zeros((p, m*k))
for j in range(p):
    eigen_values, eigen_vectors = np.linalg.eig(proj_RGAS[j, :, :] + proj_RGAS[j, :, :].T)

    idx = eigen_values.argsort()[::-1]   
    eigen_values = eigen_values[idx]
    eigen_vectors = eigen_vectors[:, idx]

    U[j, :] = eigen_vectors[:, :k].flatten('F')


np.savetxt('proj_subspace_model1_exp_d10.csv', U, delimiter=',')

