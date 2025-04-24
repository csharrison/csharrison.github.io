---
layout: post
title:  "Privacy amplification by repeated noise addition"
date:   2025-04-23
needsmath: true
tag: privacy
---

Imagine we're trying to privately estimate some statistic by adding independent noise to it. Furthermore, as in
our [previous post]({% post_url 2024-08-05-distributed-noise %}), we'll assume noise is added across $n$ parties.
Unlike the previous post, we will explicitly _not_ take advantage of infinite divisibility. Rather, we will assume
each party adds iid noise satisfying a _local_ $(\alpha, \epsilon_0)$ Renyi DP bound ([Mironov 2017](https://arxiv.org/pdf/1702.07476)), and
use a neat amplification argument to show a privacy improvement. We'll begin by stating the main result.

<a id="thm1" href="#thm1">Theorem 1</a>: Let $P$ be a continuous probability distribution on $\mathbb{R}$, and let $P^{*n}$ denote
its $n$-fold convolution. Then

$$
D_\alpha(P^{*n} + 1 || P^{*n}) \le n \cdot D_\alpha(P + 1/n || P)
$$

Where $D_\alpha$ denotes the Renyi divergence of order $\alpha$.

### Amplification by iteration and the shift-reduction lemma

[Theorem 1](#thm1) follows from the work of [Feldman et al. 2018](https://arxiv.org/abs/1808.06651), which shows a generic
privacy amplification result under contractive iteration, which allowed them to better analyze the privacy of
stochastic gradient descent.

We'll consider a much simpler setting, but using a key lemma they proved called the _shift-reduction lemma_.
Unfortunately, understanding the lemma requires introducing a few new concepts.

#### $\infty$-Wasserstein Distance


We'll denote the $\infty$-Wasserstein distance between two distributions $U$ and $V$ as $W_\infty(U, V)$.
Given that we are working with distributions in a relatively simple space (the reals), I'll refrain from a formal measure-theoretic
definition of the Wasserstein distance. Instead we will simply say that the $W_\infty(U, V) \le s$ is equivalent to saying that
$\Pr[|U - V| \le s] = 1$, where we abuse notation a bit between a random variable and its distribution.

In our case we can think about the distance as measuring _how shifted_ $U$ is with respect to $V$. If $U = V + x$ for $x$ a constant, $W_\infty(U, V) = x$.

#### Shifted Renyi divergence

The $z$-shifted Renyi divergence between distributions $U$ and $V$ is defined as 

$$
D_{\alpha}^{(z)}(U || V) = \inf_{U' : W_\infty(U, U') \le z} D_\alpha(U' || V)
$$

In other words, we are minimizing the Renyi Divergence by considering all possible distributions $U'$ that are $z$-close to $U$
in the Wasserstein distance.

#### The shift-reduction lemma

Finally we can state a simplified version of the shift-reduction lemma:

<a id="thm1" href="#lem2">Lemma 2</a>: Let $U, V$ and $\zeta$ be distributions on $\mathbb{R}$, then for any $a \ge 0$,

$$
D_{\alpha}^{(z)}(U + \zeta|| V + \zeta) \le D_{\alpha}^{(z + a)}(U || V) + R_\alpha(\zeta, a)
$$

Where 
$$
R_\alpha(\zeta, a) = \sup_{x: |x| \le a} D_\alpha(\zeta + x || \zeta)
$$ measures how well $\zeta$ hides a constant shift.


### Proof of [Theorem 1](#thm1)

The proof is a fairly simple application of the shift-reduction lemma.

$$
\begin{align*}
D_{\alpha}(P^{*n} + 1 || P^{*n}) &\le  D_{\alpha}^{(a_1)}(P + 1 || P) + R_\alpha(P^{*n-1}, a_1)\\
&= D_{\alpha}(P + 1 - a_1 || P) + D_\alpha(P^{*n-1} + a_1 || P^{*n-1})\\
&\le D_{\alpha}(P + 1 - a_1 || P) +  D_{\alpha}(P + a_1 - a_2 || P) + D_\alpha(P^{*n-2} + a_2 || P^{*n-2})\\
&\cdots\\
&\le D_\alpha(P + a_{n-1} || P) + \sum_{i=1}^{n-1} D_\alpha(P + a_{i-1} - a_i || P) &\text{ where $a_0 = 1$}
\end{align*}
$$.

Now set $a$ values such that $a_{i-1} - a_i = 1/n$ i.e. $\\{a_0, a_1, \dots, a_{n-1}\\} = \\{1, \frac{n-1}{n}, \dots, 1/n\\}$.
Then this expression is simplified to 

$$
D_\alpha(P^{*n} + 1 || P^{*n}) \le n \cdot D_\alpha(P + \frac{1}{n} || P)
$$.

### Applications

#### Gaussian noise

Let $P \sim N(0, \sigma^2)$.
We know from Mironov 2017 that $D_\alpha(P + \Delta || P) = \frac{\alpha \Delta^2}{2 \sigma^2}$. From Theorem 1 we can say that

$$
D_\alpha(P^{*n} + 1 || P^{*n}) \le n \cdot D_\alpha(P + 1/n || P) = n \cdot \frac{\alpha (1/n)^2}{2 \sigma^2} = \frac{\alpha}{2 n \sigma^2}.
$$

We know this is tight since $P^{*n} \sim N(0, n \sigma^2)$.

#### Laplace noise

Again from Mironov 2017 we know that if $P \sim Lap(\lambda)$:

$$
D_\alpha(P + \Delta || P) =  \frac{1}{\alpha - 1}\log\left(
    \frac{\alpha}{2 \alpha -1} \exp\left(\frac{\Delta (\alpha-1)}{\lambda}\right) + 
    \frac{\alpha-1}{2\alpha - 1}\exp\left(-\frac{\Delta \alpha}{\lambda}\right)
\right).
$$

So from Theorem 1 we can say:
$$
\begin{align*}
D_\alpha(P^{*n} + 1 || P^{*n}) &\le n \cdot D_\alpha(P + 1/n || P)\\
&= \frac{n}{\alpha - 1}\log\left(
    \frac{\alpha}{2 \alpha -1} \exp\left(\frac{\alpha-1}{n \lambda}\right) + 
    \frac{\alpha-1}{2\alpha - 1}\exp\left(-\frac{\alpha}{n \lambda}\right)
\right)\\
&\sim \frac{\alpha}{2 n \lambda^2} \text{(as  $n \to \infty$).}
\end{align*}
$$

Since $$Var(P^{*n}) = 2 n \lambda^2$$
, we can say the RDP is bounded by $\frac{\alpha}{Var(P^{*n})}$ as $n \to \infty$, which is
a factor of 2 worse than what we get with Gaussian noise. My gut is that this is not tight, but
right now I am too lazy to run the numerical accounting tools on Laplace convolutions.
Still, it's pretty cool this technique gets within a constant factor of the Gaussian mechanism.

#### Discrete distributions

Unfortunately the technique does not extend naturally to discrete distributions on $\mathbb{Z}$.
One thing we can do is say that $D_\alpha(P + \Delta || P) \le \Delta D_\alpha(P + 1 || P)$,
but unless we are doing $\Delta$-summation anyway, it's not really as natural. For instance,
if we are doing binary summation we could intentionally crank up $\Delta$ and rescale, but at that
point we're pretty much in the "continuous approximation" regime anyways. 