import numpy as np
from scipy.stats import multivariate_normal

'''
Implements a Gaussian Mixture Model solver.
The code is adapted from 
https://github.com/ethen8181/machine-learning/blob/master/clustering/GMM/GMM.ipynb,  
which has MIT license: https://github.com/ethen8181/machine-learning/blob/master/LICENSE

adaptations here:
- changed parameter names to match notation from Stanford CS229 Machine Learning lecture notes
- minor change to the log likelihood calculation
'''
class GMM:
    def __init__(self, k: int = 3, tol: float = 1E-4, seed=None):
        '''
        initialize the Gaussian Mixture Model
        Args:
            k: number of gaussians to use in model.  default is 3
            tol: amount at which the difference between log loss of successive
                 iterations is considered negligible, for convergence.  default is 1E-7.
            seed: random number seed if want reproducible results.  default is None.
        '''
        self.k = k
        self.tol = tol
        if seed is not None:
            np.random.seed(seed)

    def fit(self, X : np.array, max_iter : int=1000) :
        '''
        fit a Gaussian Mixture Model's parameters φ, μ and Σ given the hyperparameter k from initialization, and
        the data X using Expectation-Maximization.
        Args:
            X: matrix of data of size [n x p] where n is the number of training sets and p is the
               dimension of the datasets (i.e. data features).
            max_iter: maximum number of iterations to use in Expectation-Maximization.  
               default is 1000.
        Returns:
            True if converged, else False
        '''
        self._init_params(X)

        log_l_prev = np.log(1.E-323)
        ni = 0
        while ni < max_iter:

            # E-step: calc posterior w = p(z|x)
            self._calc_posterior(X)
            log_l = self._calc_log_likelihood()

            # M-step: estimate params
            self._calc_params(X)

            ni += 1

            if (log_l - log_l_prev) <= self.tol:
                print(f'converged in {ni} iterations')
                return True

            print(f'{ni} {log_l}')
            log_l_prev = log_l

        return False
        
    def _init_params(self, X : np.array) :
        self.n = X.shape[0] 
        self.p = X.shape[1]

        # init mu with random selection from X
        # k mu's
        # shape [k x p]
        self.mu = X[np.random.choice(self.n, self.k, replace = False)] 
        # shape [k x 1]
        self.phi = np.full(self.k, 1. / self.k)
        # np.cov(X, rowvar=False) is same as np.matmul(X_zc.T, X_zc)/( (1.-p)**2) where X_zc has zero-centered cols
        # k arrays of size [pxp] initialized with covariance of entire dataset X
        # shape [k x p x p]
        self.sigma = np.full((self.k, self.p, self.p), np.cov(X, rowvar=False))
        # shape [n x k].  p(z|x)
        self.w = np.zeros((self.n, self.k))
        
    def _calc_posterior(self, X : np.array) :
        '''
        Q_i(z^(i) == j) = w^(i)_j = p(z^(i) = j | x^(i); φ, μ, Σ) 
        from Bayes Rule = p(z|x) = p(x|z)*p(z)/p(x)
                                  
        w^(i)_j = N(x^(i) | μ_j, Σ_j) * φ_j / sum_over ell = 1 to_k( φ_ell * N(x^(i) | μ_ell, Σ_ell) )
        '''
        #shape [N x p x k]
        for j in range(self.k):
            # length p.  p(z)
            prior = self.phi[j]
            # shape [n x 1].  p(x|z)
            likelihood = multivariate_normal(self.mu[j], self.sigma[j]).pdf(X)
            # shape [n x 1] in column k
            self.w[:, j] = prior * likelihood

        # where p(x; θ) = sum over z of p(x|z)*p(z)
        self.px = np.sum(self.w, axis=1)

        #there should not be 0's in self.p(x)
        #these should sum to 1 for a single data point x^(i)
        self.w = (self.w.T / self.px).T

    # from https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/sklearn/mixture/_gaussian_mixture.py#L510
    # and https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/sklearn/mixture/_base.py
    #uses the BSD-3 Clause New or Revised license: https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/COPYING
    def bic(self):
        n_params = self.k * self.p * (self.p + 1) / 2.0
        return -2 * np.mean(np.log(self.px)) * self.n + n_params * np.log(self.n)
        return b

    # from https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/sklearn/mixture/_gaussian_mixture.py#L510
    # and https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/sklearn/mixture/_base.py
    #uses the BSD-3 Clause New or Revised license: https://github.com/scikit-learn/scikit-learn/blob/77aeb825b6494de1e3a2c1e7233b182e05d55ab0/COPYING
    def aic(self):
        n_params = self.k * self.p * (self.p + 1) / 2.0
        return -2 * np.mean(np.log(self.px)) * self.n + 2 * n_params

    def _calc_log_likelihood(self) :
        #log_l = l(θ) = sum over j=1,n ( log p(x; θ) ) 
        return np.sum(np.log(self.px))

    def _calc_params(self, X : np.array) :
        # phi(j) = sum over i=1,n of w^(i)_j
        wsum = np.sum(self.w, axis=0)

        self.phi = wsum/self.n

        self.mu = np.dot(self.w.T, X) / wsum.reshape(-1,1)

        for j in range(self.k):
            diff = (X - self.mu[j]).T
            weighted_sum = np.dot(self.w[:, j] * diff, diff.T)
            self.sigma[j] = weighted_sum / wsum[j]

    def print_params(self):
        print(f'mean=mu={self.mu}')
        print(f'phi={self.phi}')
        print(f'cov=sigma={self.sigma}')

    def get_mu(self):
        return self.mu.copy()

    def get_sigma(self):
        return self.sigma.copy()

    def get_phi(self):
        return self.phi.copy()

