# Compare different quadrature rules for integration

There are two examples provided for calculating the weights and abscissas for gaussian quadrature rules, try:

```
make
./gqconstants
```

or

```
python gqconstants.py
```

You can also use the C++ example as a guide to build your own executable

There is no need to look at rules >~25 for Gaussian quadrature.  And you can also stop at ~ 1000 divisions for the trapezoidal and Simpson's rules.  If you run much longer you'll see the numerical errors bevome visible for the trapezoidal, but hyou'll need to think about how to code efficiently or the running time may be very long.


Landau question 4.) From the plots, the power law dependence of the error on the number of points for the trapezoidal method is alpha = -2.47 and alpha = -2.01 in the algorithmic and round-off error regimes, whereas for Simpson's rule it is alpha = -5.02 and alpha = -3.98 in the algorithmic and round-off error regimes.

The hard to integrate functions used to create bad error with these algorithms were mostly functions with either (or both) many discontinuities or points of non-differentiability. It is clear why discontinuities will create bad error for these algorithms because they all approximate the function via continuous functions. Furthermore, the points of non-differentiability also create bad error for the Simpson's rule and Gaussian Quadrature methods in particular because those methods approximate the integrand via smooth polynomials. 

This error could be improved by integrating over the connected components of the graph of the function and then summing the results - which would help resolve the first main cause of the error (in the example used, for instance, the domains integrated over would be (-1,-1/2), (-1/2,1/2), and (1/2,1) because |x| = 1/2 were the points of discontinuity of the function). The second way to decrease error in the regions of non-differentiability would be to perform a seperate integration in a neighborhood of the point with a larger number of evaluation points (or smaller step size) and add the result to the integral performed over the rest of the domain, as suggested in the last way of reducing error. Alternatively, the step size could be adaptively increased/decreased where appropriate in a single integration over the entire domain.


