The Discrete Laplace mechanism is similar to the Laplace mechanism except it operates only on the integers.
Rather than sampling from the Laplace distribution, it samples from the [Discrete Laplace](https://docs.scipy.org/doc/scipy/tutorial/stats/discrete_dlaplace.html) distribution:
```math
DLaplace_{a,\mu}(x) = tanh\left(\frac{a}{2}\right)e^{-a |x-\mu|}
```

Here is a python implementation. Like the normal Laplace mechanism it operates on $l_1$ bounded sensitivity vectors (denoted by $\Delta$).
```python
import scipy
def discrete_laplace_mech(v, sensitivity, epsilon):
  return scipy.stats.dlaplace.rvs(loc=v, a=epsilon/sensitivity)
```

The Discrete Laplace mechanism satisfies $\epsilon$-differential privacy. The proof is nearly identical to the [proof of the normal Laplace mechanism](laplace-proof.md).
The only difference is the $tanh\left(\frac{a}{2}\right)$ factor in the PMF, which is different from the $\frac{1}{2 b}$ factor in the Laplace distribution because
the discrete Laplace only operates on the integers. These factors drop out immediately in the proof.

Recall the definition of differential privacy:
```math
Pr[A(D) = x] \le e^\epsilon Pr[A(D’) = x]
```
For $D$ and $D’$ neighoring datasets (that differ on one row / user), a privacy mechanism $A$, and all outputs $x$.

```math
\begin{align*}
\frac{Pr[A(D) = x]}{Pr[A(D') = x]} &=
\frac{ \prod_i tanh\left(\frac{a}{2}\right)e^{-a |x-D_i|} } { \prod_i tanh\left(\frac{a}{2}\right)e^{-a |x-D'_i|} } \\
&= \frac{ \prod_i exp(-a |x-D_i|) } { \prod_i exp(-a |x-D'_i|) } \\
&= \prod_i exp(-a |x-D_i| + a|x - D'_i|) \\
&= \prod_i exp(a (|x - D'_i| - |x-D_i|)) \\
&\le \prod_i exp(a (|D_i - D'_i|)) &&\text{by triangle inequality}\\
&= exp\left(a \sum_i |D_i - D'_i|\right) \\
&\le exp\left(a \Delta\right) &&\text{by $l_1$ sensitivity}\\
&= exp\left(\frac{\epsilon}{\Delta} \Delta\right) \\
&= e^\epsilon \\
\end{align*}
```
