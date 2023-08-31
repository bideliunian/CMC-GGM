# %%
import numpy as np
from scipy.stats import norm
from scipy.linalg import sqrtm


def rnorm_precision(n, mu, theta):
    # generate multivariate gaussian data with precision matrix
    assert n >= 1

    sigma = np.linalg.inv(theta)
    X = np.random.multivariate_normal(mu, sigma, size=n)

    return X


def transform(x, type):
    # 1-dimensional montone transformation, standardized to keep std of x
    if type == 'power':
        y = np.power(x, 3)
    elif type == 'logit':
        y = 1. / (1. + np.exp(-x))
    elif type == 'exp':
        y = np.exp(x)
    elif type == 'cdf':
        y = norm.cdf(x)

    return y


def generate_data(n, p, m, k, model, trans_type, B_type):
    '''
    x~N(mu, sigma)2)
    Para:
        n: sample size
        p: number of nodes
        m: length of vector on each node
        k: dimension of non-gaussianality
        model: model1 - AR(2) banded; model2 - random sparse;
            model3 -  hub connected.

    Return:
        n * (m * p) gaussian data(oracle)
        n * (m * p) transformed data(observed)
    '''
    assert m >= 2
    assert n >= 1
    assert m >= 1
    assert p >= 1
    assert k >= 1
    assert model in ['model1', 'model2']

    Omega = np.zeros((m*p, m*p))
    if model == 'model1':
        A = np.zeros((m, m))
        A[0, 0] = A[1, 1] = 1

        for j in range(p-1):
            Omega[(j*m):((j+1)*m), ((j+1)*m):((j+2)*m)] = 0.4 * A

            if j != (p-2):
                Omega[(j*m):((j+1)*m), ((j+2)*m):((j+3)*m)] = 0.2 * A

        Omega = Omega + Omega.T + np.diag(np.ones(m*p))

    if model == 'model2':
        delta = 0.3
        A = np.eye(m)

        for i in range(p/2-2):
            for j in range(i+1, p/2-1):
                if np.random.uniform < 0.1:
                    Omega[(i*m):((i+1)*m), (j*m):((j+1)*m)] = delta * A

        for i in range(p/2, p-2):
            for j in range(i+1, p-1):
                if np.random.uniform < 0.1:
                    Omega[(i*m):((i+1)*m), (j*m):((j+1)*m)] = delta * A

        Omega = Omega + Omega.T
        np.fill_diagonal(Omega, 0)
        ee = np.linalg.eigvalsh(Omega)[0]
        if ee < 0:
            np.fill_diagonal(Omega, -ee + 1.0)
        else:
            np.fill_diagonal(Omega, 1.0)

    Sigma = np.linalg.inv(Omega)
    diag_Sigma = np.zeros((m*p, m*p))
    for i in range(p-1):
        diag_Sigma[(i*m):((i+1)*m), (i*m):((i+1)*m)] = \
        np.linalg.inv(sqrtm(Sigma[(i*m):((i+1)*m), (i*m):((i+1)*m)]))

    cov_X = diag_Sigma.dot(Sigma).dot(diag_Sigma)
    X_oracle = np.random.multivariate_normal(mean=np.zeros(m*p), cov=cov_X, size=n)

    X = np.zeros((n, m*p))

    if (B_type == 1):
        B = np.zeros((m, k))
        B[0, 0] = B[1, 1] = B[0, 1] = 1 / np.sqrt(2)
        B[1, 0] = -1 / np.sqrt(2)
    else:
        B = np.zeros((m, k))
        B[0, 0] = B[1, 0] = B[m-1, 0] = 1 / np.sqrt(3)
        B[2, 1] = B[3, 1] = 1 / np.sqrt(2)

    Pv = np.eye(m) - B.dot(B.T)

    for j in range(p):
        X_A = X_oracle[:, (j*m):((j+1)*m)].dot(B)
        X_B = transform(X_A, trans_type)
        X_B = (X_B - np.mean(X_B, axis=0)) / np.std(X_B, axis=0) * np.std(X_A, axis=0)
        X[:, (j*m):((j+1)*m)] = X_B.dot(B.T) + X_oracle[:, (j*m):((j+1)*m)].dot(Pv)
        # X[:, (j*m):(j*m+k)] = X_B.dot(B.T[:, :k])
        # X[:, (j*m+k):((j+1)*m)] = X_oracle[:, (j*m+k):((j+1)*m)]

    return X_oracle, X

# %%

# import matplotlib.pyplot as plt
# n = 200
# params = {'n': n, 'p': 4, 'm': 5, 'k': 2, 'model': 'model1', 'trans_type': 'exp',
#           'B_type': 1}
# X, Y = generate_data(**params)

# a = (1./n) * np.ones(n)
# b = (1./n) * np.ones(n)

# plt.scatter(X[:,0], X[:,1], s=X.shape[0]*20*a, c='r', zorder=10, alpha=0.7)
# plt.scatter(Y[:,0], Y[:,1], s=Y.shape[0]*20*b, c='b', zorder=10, alpha=0.7)
# plt.axis('equal')

# %%
