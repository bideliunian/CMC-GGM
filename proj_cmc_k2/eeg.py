# %%
import numpy as np
from scipy.stats import norm
from scipy.stats import qmc

# from Optimization.RiemannianBCD import RiemannianBlockCoordinateDescent
from Optimization.RiemannianGAS import RiemannianGradientAscentSinkhorn
from PRW import ProjectedRobustWasserstein
from gen_data import *
from scipy.stats import qmc

n = 200
n_halton = 200
model = 'model1'
trans_type = 'exp'
m = 5
k = 2
p = 50
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
# proj_RBCD = np.zeros((p, m, m))
proj_RGAS = np.zeros((p, m, m))
# values_subspace_RBCD = np.zeros((nb_exp, p))
values_subspace_RGAS = np.zeros((nb_exp, p))

proj = np.zeros((m,m))   # Real optimal subspace
proj[0, 0] = 1
proj[1, 1] = 1


def InitialStiefel(d, k):
    U = np.random.randn(d, k)
    q, _ = np.linalg.qr(U)
    return q


for i in range(nb_exp):
    np.random.seed(i)
    params = {'n': n, 'p': p, 'm': m, 'k': k, 'model': model, 'trans_type': trans_type}
    X_oracle, X = generate_data(**params)
    X_transformed = np.zeros((n, m*p))

    U0 = InitialStiefel(m, k)

    for j in range(p):
        # params = {'X':X[:, (j*m):((j+1)*m)], 'Y':normal_standard, 'a':a, 'b':b,
        #                            'U0': U0, 'reg':eta, 'k':k, 'max_iter':2000,
        #                            'stopThr':0.01, 'verbose':1}
        # ot.dr.projection_robust_wasserstein(**params)
        # algo = RiemannianBlockCoordinateDescent(eta=eta, tau=None, max_iter=3000, threshold=0.01, verbose=True)
        # PRW = ProjectedRobustWasserstein(X[:, (j*m):((j+1)*m)], normal_standard, a, b, algo, k)
        # PRW.run('RBCD', tau, U0)

        algo1 = RiemannianGradientAscentSinkhorn(eta=eta, tau=None, max_iter=3000, threshold=0.01,
                                                 sink_threshold=1e-4, verbose=True)
        PRW1 = ProjectedRobustWasserstein(X[:, (j*m):((j+1)*m)], normal_standard, a, b, algo1, k)
        PRW1.run('RGAS', tau/eta, U0)

        # values_subspace_RBCD[i, j] = np.linalg.norm(PRW.get_Omega() - proj)
        values_subspace_RGAS[i, j] = np.linalg.norm(PRW1.get_Omega() - proj)

        pi = PRW1.get_pi()
        X_transformed[:, (j*m):((j+1)*m)] = pi.dot(normal_standard)*n

        # proj_RBCD[j, :, :] += PRW.get_Omega() / nb_exp
        proj_RGAS[j, :, :] += PRW1.get_Omega() / nb_exp


# U_RBCD = np.zeros((p, m, k))
U_RGAS = np.zeros((p, m*k))

for j in range(p):
    # _, v_RBCD = np.linalg.eig(proj_RBCD[j, :, :])
    eigen_values_RGAS, eigen_vectors_RGAS = np.linalg.eig(proj_RGAS[j, :, :] + proj_RGAS[j, :, :].T)

    idx = eigen_values_RGAS.argsort()[::-1]   
    eigen_values_RGAS = eigen_values_RGAS[idx]
    eigen_vectors_RGAS = eigen_vectors_RGAS[:, idx]

    # U_RBCD[j, :, :] = v_RBCD[:, :k]
    U_RGAS[j, :] = eigen_vectors_RGAS[:, :k].flatten('F')


np.savetxt('proj_subspace_model1_exp_d05.csv', U_RGAS, delimiter=',')
# %%
# import matplotlib.pyplot as plt
# X, Y = X_oracle, X_transformed

# a = (1./n) * np.ones(n)
# b = (1./n) * np.ones(n)

# plt.scatter(X[:,0], X[:,1], s=X.shape[0]*20*a, c='r', zorder=10, alpha=0.7)
# plt.scatter(Y[:,0], Y[:,1], s=Y.shape[0]*20*b, c='b', zorder=10, alpha=0.7)
# plt.axis('equal')
# %%
