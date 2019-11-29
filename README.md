# Graphical Component Analysis

This is an R implementation of the algorithm described in <a href="gca.pdf" download>Graphical component analysis via score matching</a>. A sample script to test the algorithm is below. 

## External dependencies

[orthopolynom](https://cran.r-project.org/web/packages/orthopolynom/index.html)

[igraph](https://igraph.org/r/)

[huge](https://cran.r-project.org/web/packages/huge/index.html)

## Example
```
from lhood import *
from sm import *

mu0 = -4.0
mu1 = 4.0
sigma_prior = 1.0
sigma = 1.0
x = 0.0
knots = np.arange(-5.0, 5.0, 1.0/20)
degree = 3
theta_grid = np.arange(-6, 6, .01) 

# score-matching
tmp = sm_spline(theta_grid, joint_gaussian_gmm, knots, degree)
post_sm_normalized = tmp.estimate_unbdd()

# plot the true posterior and the estimate
plt.plot(theta_grid, safe_log(post_sm_normalized), color = "orange", label = "score-matching spline")
plt.plot(theta_grid, safe_log(post_normalized), color = "black", label = "true")
plt.legend(loc = "lower right")
plt.suptitle("Log posterior estimates: Gaussian likelihood and Gaussian mixture prior on means")
plt.show()
```
