import numpy as np
from unittest import TestCase
from GMM import *
import matplotlib.pyplot as plt

'''
Some code is adapted from 
https://github.com/ethen8181/machine-learning/blob/master/clustering/GMM/GMM.ipynb,  
which has MIT license: https://github.com/ethen8181/machine-learning/blob/master/LICENSE
'''

def generate_data(n_data, means, covariances, weights):
    """creates a list of data points"""
    n_clusters, n_features = means.shape
    
    data = np.zeros((n_data, n_features))
    for i in range(n_data):
        # pick a cluster id and create data from this cluster
        k = np.random.choice(n_clusters, size = 1, p = weights)[0]
        x = np.random.multivariate_normal(means[k], covariances[k])
        data[i] = x
   
    return data

def plot_contours(data, means, covs, title):
    """visualize the gaussian components over the data"""
    plt.figure()
    plt.plot(data[:, 0], data[:, 1], 'ko')

    delta = 0.025
    k = means.shape[0]
    x = np.arange(-2.0, 7.0, delta)
    y = np.arange(-2.0, 7.0, delta)
    x_grid, y_grid = np.meshgrid(x, y)
    coordinates = np.array([x_grid.ravel(), y_grid.ravel()]).T

    col = ['green', 'red', 'indigo', 'blue', 'brown', 'orange']
    for i in range(k):
        mean = means[i]
        cov = covs[i]
        z_grid = multivariate_normal(mean, cov).pdf(coordinates).reshape(x_grid.shape)
        plt.contour(x_grid, y_grid, z_grid, colors = col[i])

    plt.title(title)
    plt.tight_layout()

def testA():
    init_mu = np.array([
        [5, 0],
        [1, 1],
        [0, 5]
    ])

    init_sigma = np.array([
        [[.5, 0.], [0, .5]],
        [[.92, .38], [.38, .91]],
        [[.5, 0.], [0, .5]]
    ])

    init_phi = [1 / 4, 1 / 2, 1 / 4]

    seed = 1234
    np.random.seed(seed)

    # generate data
    X = generate_data(100, init_mu, init_sigma, init_phi)

    plt.plot(X[:, 0], X[:, 1], 'ko')
    plt.tight_layout()

    gmm = GMM(k=3, tol=1E-7, seed=seed)
    p_x = gmm.fit(X)

    gmm.print_params()
    plot_contours(X, gmm.mu, gmm.sigma, 'Initial clusters')

    phi = gmm.get_phi()
    mu = gmm.get_mu()
    sigma = gmm.get_sigma()

    # assert expected within a tolerance.  asserting mu.  TODO: add sigma and phi
    '''
    found_j = set()
    mu_tol = 0.25
    for i in range(init_mu.shape[0]):
        best_j = -1
        min_d = 1E11*np.ones((mu[i].shape[0]))
        for j in range(mu.shape[0]):
            if j in found_j:
                continue
            d = abs(init_mu[i] - mu[j])
            if all(d < min_d):
                min_d = d
                best_j = j
            np.testing.assert_array_equal((init_mu[i]-mu[best_j] < mu_tol), np.ones((mu[best_j].shape[0])), "not equal")
            found_j.add(best_j)
    '''

    print(f'aic={gmm.aic(p_x)}')
    print(f'bic={gmm.bic(p_x)}')

    pxs = []
    aics = []
    bics = []
    js = []
    for j in range(1, 5) :
        gmm = GMM(k=j, tol=1E-7, seed=seed)
        p_x = gmm.fit(X)
        plot_contours(X, gmm.mu, gmm.sigma, 'Initial clusters')
        aics.append(gmm.aic(p_x))
        bics.append(gmm.bic(p_x))
        js.append(j)

    print(f'aics={aics}')
    print(f'bics={bics}')
    print(f'pxs={pxs}')
    i = np.argmin(aics)
    i2 = np.argmin(bics)
    print(f'best k={js[i]}, {js[i2]}')


if __name__ == "__main__":
    testA()
