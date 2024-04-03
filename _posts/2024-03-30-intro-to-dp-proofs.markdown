---
layout: post
title:  "Introduction to differential privacy proofs"
date:   2024-03-30
needsmath: true
---

Differential privacy is a robust, formal notion of privacy. For a friendly introduction, it's difficult to beat
this [blog post series](https://desfontain.es/privacy/friendly-intro-to-differential-privacy.html) by Damien Desfontaines.
In this post I'll assume a basic familiarity with the concept, and dive into a few examples of _proving_ that mechanisms are differentially private. I've found that getting a handle on a few simple proofs helps build a better intuition for how differential privacy works.

## The definition, and some helpful simplifications
This is the definition of differential privacy:

$$
\begin{align}
a_1& =b_1+c_1\\
a_2& =b_2+c_2-d_2+e_2
\end{align}
$$


$$
\begin{align*}
a_1& =b_1+c_1\\
a_2& =b_2+c_2-d_2+e_2
\end{align*}
$$

For all neighboring datasets $D$ and $D’$ (that differ on one row / user), a privacy mechanism $A$, and all (measureable) subsets $S \subseteq image(A)$.

For mechanisms that induce discrete probability distributions, this definition is often simplified to the following:

$$
Pr[A(D) = x] \le e^\epsilon Pr[A(D’) = x]
$$

for all possible outputs $x$. [This follows](https://cstheory.stackexchange.com/a/50824) from the fact that

$$
Pr[A(D) \in S] = \sum_{x \in S} Pr(A(D) = x) \le \sum_{x \in S} e^\epsilon Pr(A(D') = x) = e^\epsilon Pr[A(D’) \in S]
$$

Similarly, in the continuous case, $P(A(D) = x)$ will always equal zero at any given output $x$. However, if $p_{A(D)}$ is the PDF of the distribution induced by $A(D)$, we can say

$$
Pr[A(D) \in S] = \int_{x \in S} p_{A(D)}(x) dx \le \int_{x \in S} e^\epsilon p_{A(D')}(x) dx = e^\epsilon Pr[A(D’) \in S]
$$

So to show that $A$ satisfies differential privacy, it suffices to show that for all possible outputs $x$:

$$
p_{A(D)}(x) \le e^\epsilon p_{A(D')}(x)
$$

## Proofs

### K-ary randomized response
K-ary randomized response is a privacy mechanism with the following algorithm:
```python
import math
import random
def kary_rr(value: Any, possible_outputs: Set[Any], epsilon: float) -> Any:
  """Execute the k-rr mechanism on value

  Arguments:
  value: the value to randomize. Must be a member of `possible_outputs`
  possible_outputs: a set of possible outputs of the mechanism
  epsilon: the privacy loss parameter.
  """

  p = k / (k - 1 + math.exp(epsilon))
  assert value in possible_outputs
  if random.random() <= p:
    return random.choice(possible_outputs)
  return value
```

This algorithm satisfies $\epsilon$-differential privacy. Because k-RR is a local differentially private mechanism (i.e. it operates on only a single row), we can assume neighboring datasets $D$ and $D’$ differ only in that they are different elements of the `possible_outputs` set.

$$
\begin{align*}
\frac{Pr[A(D) = x]}{Pr[A(D') = x]} &\le \frac{Pr[A(x) = x]}{Pr[A(D') = x]} &&\text{ for $D' \ne x$}\\
&= \frac{(1 - p) + \frac{p}{k}}{\frac{p}{k}} \\
&= \frac{k}{p} - k + 1 \\
&= k \frac{k - 1 + e^\epsilon}{k} - k + 1 \\
&= e^\epsilon
\end{align*}
$$

The first inequality follows from the fact that the mechanism will always output the true input with higher probability than any other output.

### RAPPOR

Basic RAPPOR ([Erlingsson et al 2014](https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/42852.pdf))
is a privacy mechanism with the following algorithm:

```python
import math
import random
def rappor(value: int, domain_size: int, epsilon: float):
  """Execute the RAPPOR mechanism on value

  Arguments:
  value: the value to randomize. Must be an integer between 0 and domain_size - 1
  domain_size: a positive integer > 1
  epsilon: the privacy loss parameter.
  """
  # Generate a one-hot vector of length `histogram_length`
  # All zeros except for the `value`th entry.
  out = [0] * domain_size
  out[value] = 1

  # With probability f, set the bit randomly.
  # Equivalently, flip the bit with probability f/2
  f = 2 / (1 + math.exp(epsilon/2))
  for i in range(domain_size):
    if random.random() < f/2:
      out[i] ^= 1 # flip the bit
  return out
```

This algorithm satisfies $\epsilon$-differential privacy.

Because RAPPOR is also a local differentially private mechanism (i.e. it operates on only a single row), $D$ and $D’$ differ only in `value`. For any value of $x$ we choose, $A(D)$ will need to set two bits differently from how $A(D')$ sets them. In other words,
where $A(D)$ leaves the bit unchanged, $A(D')$ must flip the bit. Where $A(D)$ flips the bit, $A(D')$ must leave the bit unchanged.

As an example: take $D = 0$, $D' = 1$, and $x = [1, 0, 0, ...]$. It is clear that for $A(D') = x$ the first two bits need to be flipped, unlike in $A(D)$
where they are left unchanged.

Therefore the ratio of $\frac{Pr[A(D) = x]}{Pr[A(D') = x]}$ is maximized when the $D'$ flips those two bits with probability $f/2$, since $f/2 \le 1 - f/2$.

$$
\begin{align*}
\frac{Pr[A(D) = x]}{Pr[A(D') = x]} &\le
\frac{\left(1-\frac{f}{2}\right)^2 M}{\left(\frac{f}{2}\right)^2 M} &&\text{ for $M$ the probability all other bits are flipped to $x$}\\
&= \frac{\left(1-\frac{f}{2}\right)^2}{\left(\frac{f}{2}\right)^2} \\
&= \frac{\left(1-\frac{1}{1 + e^{\epsilon/2}}\right)^2}{\left(\frac{1}{1 + e^{\epsilon/2}}\right)^2} \\
&= \left(1-\frac{1}{1 + e^{\epsilon/2}}\right)^2 \left(1 + e^{\epsilon/2}\right)^2 \\
&= \left(\frac{e^{\epsilon/2}}{1 + e^{\epsilon/2}}\right)^2 \left(1 + e^{\epsilon/2}\right)^2 \\
&= \frac{e^{\epsilon}}{\left(1 + e^{\epsilon/2}\right)^2} \left(1 + e^{\epsilon/2}\right)^2 \\
&= e^\epsilon
\end{align*}
$$

### Laplace mechanism

The Laplace distribution is a continuous probability distribution with following PDF:

$$
Laplace_{b, \mu}(x) := \frac{1}{2 b} exp\left(-\frac{|x-\mu|}{b}\right)
$$
For scale param $b$ and mean $\mu$.

The Laplace mechanism operates in the setting of vector aggregation.
The input dataset is a vector $V \in \mathbb{R}^d$ with a fixed $l_1$ sensitivity $\Delta$: for neighboring datasets $D, D' \in \mathbb{R}^d$,
$||D-D'||_1 \le \Delta$. In other words, each user can only contribute $\Delta$ across all of the $d$ buckets in the vector.

The mechanism simply adds independent noise to each vector coordinate sampled from the Laplace distribution with $b = \Delta/\epsilon$ and $\mu = 0$.

Here is a (buggy) implementation of the Laplace mechanism in python:
```python
import scipy
def laplace_mech(v, sensitivity, epsilon):
  return scipy.stats.laplace.rvs(loc=np.array(v), scale=sensitivity/epsilon)
```
This is buggy due in part to the results in [Mironov (2012)](https://www.microsoft.com/en-us/research/wp-content/uploads/2012/10/lsbs.pdf)
where the DP guarantees could be broken by the floating point encoding used in practice for these algorithms. See the doc
published by the Google DP library [Secure Noise Generation](https://github.com/google/differential-privacy/blob/main/common_docs/Secure_Noise_Generation.pdf)
for more information on how to fix these bugs.

For now we will assume generously that computers can operate on real numbers perfectly, and analyze the Laplace mechanism in this idealized scenario.

The Laplace mechanism satisfies $\epsilon$-differential privacy.

$$
\begin{align*}
\frac{p_{A(D)}(x)}{p_{A(D')}(x)} &\le \frac{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D_i|}{b}\right)}{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D'_i|}{b}\right)} \\
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
$$

Note that the triangle inequality states (informally) that the sum of any two sides of a triangle is greater than or equal to the third side.
You can write this as the following:
$$
\begin{align*}
|A - C| &\le |A - B| + |B - C| \\
|A - C| - |A - B| &\le |B - C|
\end{align*}
$$

In other words, the distance between the points $A$ and $C$ will be less than (or equal to) the distance from $A$ to $B$, and then $B$ to $C$. This holds for scalars as well as vectors.


The proof follows by substituting $A = x$, $B = D_i$, and $C = D'_i$.

### Discrete Laplace mechanism
The Discrete Laplace mechanism is similar to the Laplace mechanism except it operates only on the integers.
Rather than sampling from the Laplace distribution, it samples from the [Discrete Laplace](https://docs.scipy.org/doc/scipy/tutorial/stats/discrete_dlaplace.html) distribution:
$$
DLaplace_{a,\mu}(x) = tanh\left(\frac{a}{2}\right)e^{-a |x-\mu|}
$$

Here is a python implementation. Like the normal Laplace mechanism it operates on $l_1$ bounded sensitivity vectors (denoted by $\Delta$).
```python
import scipy
def discrete_laplace_mech(v, sensitivity, epsilon):
  return scipy.stats.dlaplace.rvs(loc=v, a=epsilon/sensitivity)
```

The Discrete Laplace mechanism satisfies $\epsilon$-differential privacy. The proof is nearly identical to the proof of the normal Laplace mechanism.
The only difference is the $tanh\left(\frac{a}{2}\right)$ factor in the PMF, which is different from the $\frac{1}{2 b}$ factor in the Laplace distribution's PDF because
the discrete Laplace only operates on the integers. These factors drop out quickly in the proof.

$$
\begin{align*}
\frac{Pr[A(D) = x]}{Pr[A(D') = x]} &=
\frac{ \prod_i tanh\left(\frac{a}{2}\right)e^{-a |x-D_i|} } { \prod_i tanh\left(\frac{a}{2}\right)e^{-a |x-D'_i|} } \\
&= \prod_i exp(a (|x - D'_i| - |x-D_i|)) \\
&\le \prod_i exp(a (|D_i - D'_i|)) &&\text{by triangle inequality}\\
&= exp\left(a \sum_i |D_i - D'_i|\right) \\
&\le exp\left(a \Delta\right) &&\text{by $l_1$ sensitivity}\\
&= e^\epsilon \\
\end{align*}
$$