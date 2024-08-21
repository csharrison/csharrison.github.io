---
layout: post
title:  "Tight RDP and zCDP bounds for RAPPOR"
date:   2024-08-20
needsmath: true
---

We briefly covered the basic RAPPOR mechanism in a [previous post]({% post_url 2024-03-30-intro-to-dp-proofs %}#rappor). We'll copy the example code for reference:

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

Previously we showed that RAPPOR satisfies $\epsilon$ differential privacy. We'll attempt here to show
Renyi DP ([Mironov 2017](https://arxiv.org/pdf/1702.07476)) and zCDP bounds ([Bun & Steinke 2016](https://arxiv.org/pdf/1605.02065)). First, lets go over some definitions. These definitions are a bit simplified in that
we will assume discrete output spaces and distributions with full support over those spaces.

<a id="def1" href="#def1">Definition 1</a>: A mechanism $A$ satisfies $(\alpha, \hat{\epsilon})$ Renyi
Differential privacy (RDP) if, for all pairs of neighboring inputs $D$ and $D'$, and all possible outputs $x$:

$$
\frac{1}{\alpha-1}\log\left(\sum_{x} Pr[A(D) = x]^\alpha Pr[A(D')=x]^{1-\alpha}\right) \le \hat{\epsilon}
$$

<a id="def2" href="#def2">Definition 2</a>: A mechanism $A$ satisfies $\rho$ zero concentrated differential
privacy (zCDP) if $A$ satisfies $(\alpha, \rho \alpha)$-RDP for all $\alpha > 1$.


## RDP bounds for RAPPOR

Let's start by seeing what we can show about the RDP of RAPPOR. For RAPPOR, we can consider the input
and output as $d$-length binary vectors. So

$$
Pr[A(D) = x] = \prod_{i=1}^d
    \begin{cases}
    \frac{e^{\epsilon/2}}{e^{\epsilon/2} + 1} & \text{ if $D_i = x_i$}\\
    \frac{1}{e^{\epsilon/2} + 1} & \text{ if $D_i \neq x_i$}\\
    \end{cases}
$$

<a id='thm1' href='#thm1'>Theorem 1:</a> The RAPPOR algorithm where bits are kept with probability
$p=\frac{e^{\epsilon/2}}{1 + e^{\epsilon/2}}$ and flipped with probability $q = 1-p$ satisfies $(\alpha, \hat{\epsilon(\alpha)})$-RDP for all $\alpha > 1$ where:

$$
\hat{\epsilon}(\alpha) = 
\frac{2}{\alpha - 1}\log \left(
    \frac{e^{\alpha \epsilon / 2} + e^{(1 - \alpha) \epsilon/2}}{e^{\epsilon/2} + 1}
 \right)
$$

Without loss of generality (due to symmetry), consider $D$ and $D'$ that differ
only on the first and second index, i.e. $D = [1, 0, ...]$ and $D' = [0, 1, 0, ...]$.

There are 4 special outputs we need to consider:
- [0, 0, ...]
- [1, 1, ...]
- [1, 0, ...]
- [0, 1, ...]

Each of these outputs has the same possible $2^{d-2}$ suffixes, so

$$
\begin{align*}
e^{(\alpha - 1)D_\alpha\left(P||Q\right)} &= \sum_{x} Pr[A(D) = x]^\alpha Pr[A(D')=x]^{1-\alpha}\\
&= \sum_{x \in \{[0, 0, ...], [1, 1, ...]\}} p\ q\ Pr[A(D)[2:]=x[2:]]
    + \sum_{x \in [1, 0, ...]} p^{2\alpha} q^{2(1 - \alpha)} Pr[A(D)[2:]=x[2:]]
    + \sum_{x \in [0, 1, ...]} q^{2 \alpha} p^{2(1 - \alpha)} Pr[A(D)[2:]=x[2:]]
    \\\\
&= 2 p q + p^{2\alpha} q^{2(1 - \alpha)} +  q^{2\alpha} p^{2(1 - \alpha)}\\

&= 2 \frac{e^{\epsilon/2}}{\left(1 + e^{\epsilon/2}\right)^2}
    + \left(\frac{e^{\epsilon/2}}{1 + e^{\epsilon/2}}\right)^{2 \alpha} \left(\frac{1}{1 + e^{\epsilon/2}}\right)^{2 (1-\alpha)}
    + \left(\frac{e^{\epsilon/2}}{1 + e^{\epsilon/2}}\right)^{2 (1-\alpha)} \left(\frac{1}{1 + e^{\epsilon/2}}\right)^{2 \alpha}\\\\

&= 2 \frac{e^{\epsilon/2}}{\left(1 + e^{\epsilon/2}\right)^2}
    + \frac{e^{\alpha \epsilon}}{\left(1 + e^{\epsilon/2}\right)^{2}}
    + \frac{e^{(1 - \alpha) \epsilon}}{\left(1 + e^{\epsilon/2}\right)^{2}}\\\\

&= \left(\frac{e^{\alpha \epsilon / 2} + e^{(1-\alpha)\epsilon/2}}{1 + e^{\epsilon/2}}\right)^2
\end{align*}
$$

The theorem immediately follows.

## zCDP bounds for RAPPOR

We'll make use of recent results from [Steinke](https://differentialprivacy.org/pdp-to-zcdp/) which
derives tight zCDP bounds for the randomized response algorithm. Summarizing here:

<a id='thm2' href='#thm2'>Theorem 2</a>: The binary randomized response algorithm which flips the bit with
probability $q = \frac{1}{1 + e^\epsilon}$ satisfies $(\alpha, \hat{\epsilon}(\alpha))$ RDP for all $\alpha > 1$ where:

$$
\hat{\epsilon}(\alpha) = \frac{1}{\alpha - 1}\log \left(
    \frac{e^{\alpha \epsilon} + e^{(1 - \alpha) \epsilon}}{e^{\epsilon} + 1}
 \right)
$$

This is precisely the bounds we showed for RAPPOR, with two differences:
- We substitute $\epsilon/2$ for $\epsilon$
- We multiply the overall term by 2

That is, $RDP(RAPPOR(\epsilon)) = 2 RDP(RR(\epsilon/2))$.[^alt] It is equivalent (from the RDP perspective)
to invoking two instances of binary randomized response with half the "epsilon budget" each! With this perpsective,
we can leverage the additional Steinke result:

<a id='thm3' href='#thm3'>Theorem 3</a>: $\epsilon$-DP binary randomized response satifies
$\rho$-zCDP where $\rho = \frac{e^\epsilon - 1}{e^\epsilon + 1}\epsilon$.

In fact, the result in the blog post holds (as an upper bound) across all pure $\epsilon$-DP mechanisms, but it
is tight for binary randomized response. We can use this to immediately show the following:

<a id='thm4' href='#thm4'>Theorem 4</a>: $\epsilon$-DP RAPPOR satisfies (tight) $\rho$-zCDP where

$$
\rho = \frac{e^{\epsilon/2} - 1}{e^{\epsilon/2} + 1}\epsilon
$$

This follows immediately from the zCDP composition of two binary randomized response mechanisms,
which we know has the same Renyi divergence as RAPPOR. Tightness is witnessed by taking the limit
$\alpha \to 1$ for RAPPOR's Renyi divergence:

$$
\lim_{\alpha \to 1} \frac{2}{\alpha - 1}\log \left(
    \frac{e^{\alpha \epsilon / 2} + e^{(1 - \alpha) \epsilon/2}}{e^{\epsilon/2} + 1}
 \right) = \frac{e^{\epsilon/2} - 1}{e^{\epsilon/2} + 1}\epsilon
$$

Let's see what these bounds look like:

![Figure 1](/images/rappor-zcdp.svg)

This is a substantial improvement over the zCDP upper bounds.

## Leveraging zCDP for composition

Sure we can show improved zCDP bounds for RAPPOR, but what's the point? One big reason is that
zCDP has very nice advanced composition properties. That is, we can leverage zCDP to show better
privacy in the setting where we invoke RAPPOR many times over the dataset. In the plot below, we use [Bun & Steinke proposition 1.3](https://arxiv.org/pdf/1605.02065)
to translate zCDP bounds to $(\epsilon, \delta)$-DP.

![Figure 2](/images/rappor-zcdp-composition.svg)

With the improved zCDP analysis in hand, RAPPOR under k-fold composition can be shown to have much
better privacy than using simple composition results. Interestingly enough, the tight zCDP bounds
from <a href="#thm3">theorem 3</a> improve on the original advanced composition bound from Dwork & Roth
while maintaining relative simplicity.

## Aside: $k$-RAPPOR

We can consider the generalized RAPPOR algorithm where $k \le d$ bits are allowed to be 1, and the rest 0. We can
call this $k$-RAPPOR. This mechanism may be useful in the setting where users can contribute to multiple
histogram buckets at once.

<a id='thm5' href="#thm5">Theorem 5</a>: $k$-RAPPOR with flip probability $q = \frac{1}{1 + exp\left(\frac{\epsilon}{2 k}\right)}$ satisfies $\epsilon$-DP.

This is straightforward to show using a similar approach in our [previous post]({% post_url 2024-03-30-intro-to-dp-proofs %}#rappor). Briefly:

$$
\begin{align*}
\frac{Pr[A(D) = x]}{Pr[A(D') = x]} &\le \frac{(1 - q)^{2k}}{q ^ {2k}}\\
&= \frac{\left(\frac{exp\left(\frac{\epsilon}{2 k}\right)}{1 + exp\left(\frac{\epsilon}{2 k}\right)}\right)^{2k}}
        {\left(\frac{1}{1 + exp\left(\frac{\epsilon}{2 k}\right)}\right)^{2k}}\\
&= e^\epsilon
\end{align*}
$$

Similarly, we can show RDP / zCDP bounds, which we will present without proof[^sorry].

<a id='thm6' href='#thm6'>Theorem 6</a>: $k$-RAPPOR satifies $(\alpha, \hat{\epsilon}(\alpha))$-RDP for all
$\alpha > 1$ where 

$$
\hat{\epsilon}(\alpha) = 
\frac{2 k}{\alpha - 1}\log \left(
    \frac{e^\frac{\alpha \epsilon}{2 k} + e^\frac{(1 - \alpha) \epsilon}{2 k}}{e^\frac{\epsilon}{2 k} + 1}
 \right)
$$


It follows that $\epsilon$-DP $k$-RAPPOR satisfies $\rho$-zCDP for

$$
\rho = \frac{e^{\epsilon/(2k)} - 1}{e^{\epsilon/(2k)} + 1}\epsilon
$$


It's a tad unweildy to parameterize RAPPOR with $\epsilon$. Let's rewrite in terms of $q$:

$$
\begin{align*}
q &= \frac{1}{1 + e^{\epsilon / (2k)}}\\
\frac{\epsilon}{2k} &= \log \frac{1 - q}{q}\\
\epsilon &= 2 k \log\left(\frac{1-q}{q}\right)
\end{align*}
$$

So, 

$$
\begin{align*}
\rho &= \frac{e^{\log\left(\frac{1-q}{q}\right)} - 1}{e^{\log\left(\frac{1-q}{q}\right)} + 1}2 k \log\left(\frac{1-q}{q}\right)\\
&= \frac{\frac{1 - q}{q} - 1}{\frac{1 - q}{q} + 1} 2 k \log\left(\frac{1-q}{q}\right)\\
&= (1 -2q) 2 k \log \left(\frac{1-q}{q}\right)\\
\end{align*}
$$

For a fixed value of $q$, $\rho$ scales linearly with $k$. This behavior is similar to e.g. the Gaussian
mechanism, where $\rho$ scales proportional to $\Delta_2^2$, where $k$-hot binary vectors have $\Delta_2 = \sqrt{k}$. Unfortunately, some quick empirical results[^lazy] don't really show this beating the simplest baseline of the
(local) gaussian mechanism. Still, it's a somewhat interesting mechanism. Maybe there is some setting it could
be useful in :)

![Figure 3](/images/rappor-mae.svg)

[^alt]: Rather than a direct computation of the RAPPOR Renyi divergence for <a href="#thm1">Theorem 1</a>,
    I believe it should be possible to to prove it via something like RDP composition of $d$ instances of
    binary randomized response where only one coordinate can be 1. In this post I chose to do a direct 
    computation as practice, but maybe I will explore these other techniques in a future post.

[^sorry]: Sorry, it's a little tedious to do the direct computation.
[^lazy]: It should be pretty straightforward to compute e.g. per-bucket variance or some other
    utility metric for $k$-RAPPOR analytically. Maybe if it performed well empirically, I'd
    be more motivated!