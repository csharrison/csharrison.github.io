## Laplace mechanism
The Laplace distribution is a continuous probability distribution with following PDF:

```math
Laplace_{b, \mu}(x) := \frac{1}{2 b} exp\left(-\frac{|x-\mu|}{b}\right)
```
For scale param $b$ and mean $\mu$.

The Laplace mechanism operates in the setting of vector aggregation.
The input is a vector $V \in \mathbb{R}^d$ with a fixed $l_1$ sensitivity $\Delta$: for neighboring datasets $D, D' \in \mathbb{R}^d$,
$||D-D'||_1 \le \Delta$. In other words, each user can only contribute $\Delta$ across all of the $d$ buckets in the vector.

Here is a buggy implementation of the Laplace mechanism in python:
```python
import scipy
def laplace_mech(v, sensitivity, epsilon):
  return scipy.stats.laplace.rvs(loc=np.array(v), scale=sensitivity/epsilon)
```
This is buggy due in part to the results in [Mironov 2012](https://www.microsoft.com/en-us/research/wp-content/uploads/2012/10/lsbs.pdf)
where the DP guarantees could be broken by the floating point encoding used in practice for these algorithms. See the doc
published by the Google DP library [Secure Noise Generation](https://github.com/google/differential-privacy/blob/main/common_docs/Secure_Noise_Generation.pdf)
for more information on how to fix these bugs.

For now we will assume generously that computers can operate on real numbers perfectly, and analyze the Laplace mechanism in this idealized scenario.

## Proof
The Laplace mechanism satisfies $\epsilon$-differential privacy. Here is a proof. Recall the definition of differential privacy:
```math
Pr[A(D) \in S] \le e^\epsilon Pr[A(D’) \in S]
```
For $D$ and $D’$ neighoring datasets (that differ on one row / user), a privacy mechanism $A$, and all subsets $S \subseteq image(A)$.

```math
\begin{align*}
\frac{Pr[A(D) \in S]}{Pr[A(D') \in S]} &\le
\frac{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D_i|}{b}\right)}{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D'_i|}{b}\right)} \\
&= \frac{\prod_{i} exp\left(-\frac{|x-D_i|}{b}\right)}{\prod_{i} exp\left(-\frac{|x-D'_i|}{b}\right)} \\
&= \prod_{i} exp\left(-\frac{|x-D_i|}{b}\right) exp\left(\frac{|x-D'_i|}{b}\right) \\
&= \prod_{i} exp\left(\frac{|x-D'_i| - |x - D_i|}{b}\right) \\
&\le \prod_{i} exp\left(\frac{|D_i-D'_i|}{b}\right) &&\text{by triangle inequality}\\
&= exp\left(\frac{1}{b} \sum_i |D_i-D'_i|\right) \\
&= exp\left(\frac{1}{b} ||D-D'||_1\right) \\
&\le exp\left(\frac{1}{b} \Delta \right) &&\text{by $l_1$ sensitivity}\\
&= exp\left(\frac{\epsilon}{\Delta} \Delta \right) \\
&= e^\epsilon
\end{align*}
```

Note that the triangle inequality states that the sum of any two sides of a triangle is greater than or equal to the third side.
You can write this as the following:
```math
\begin{align*}
|A - C| &\le |A - B| + |B - C| \\
|A - C| - |A - B| &\le |B - C|
\end{align*}
```
In other words, the distance between the points $A$ and $C$ will be less than (or equal to) the distance from $A$ to $B$, and then $B$ to $C$.
The proof follows by substituting $A = x$, $B = D_i$, and $C = D'_i$.
