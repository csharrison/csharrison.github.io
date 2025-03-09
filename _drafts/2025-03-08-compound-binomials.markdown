---
layout: post
title:  "Compound binomial distributions and binomial splitting"
date:   2025-03-08
needsmath: true
tag: probability
---

In this post we'll explore the binomial distribution, denoted $\Bin(n, p)$, which is a discrete
probability distribution on the non-negative integers. Its probability mass function (PMF) is

$$
f_{Bin(n, p)}(k) = \binom{n}{k} p^k (1-p)^{n-k}.
$$

The binomial distribution models the number of _successful_ trials out of $n$ total trials, where
a successful trial occurs with probability $p$. In other words, if $X_1, \cdots, X_n$ are i.i.d.
$\Bern(p)$ random variables, $\Bin(n, p) \sim \sum_{i=1}^n X_i$.

A _compound_ probability distribution (also sometimes called a scale mixture, or a parameter mixture distribution) is
formed when one of the _parameters_ of the distribution is another random variable. We'll consider compounding the
binomial dis

## Binomial thinning: randomizing the number of trials

In this section we'll consider $N$ trials where $N$ is some non-negative discrete random variable. The
process here is sometimes called _binomial thinning_, and $\Bin(N, p)$ the $p$-thinned version of $N$.
If we interpret $N$ as counting some number of events, the term "thinning" makes sense since $\Bin(N, p)$
takes $p$ fraction of the events, and discards the rest.

### Geometric number of trials

Let's assume $N \sim \Geo(q)$ is distributed according to the geometric distribution (with support on the non-negative integers) whose PMF is

$$
f_{\Geo(q)}(k) = q (1-q)^k.
$$

The geometric distribution models the number of failed trials until the first success, for
a sequence of events with success probability $q$. 
Now consider $\Bin(N, p)$, we can compute the PMF as follows:

$$
\begin{align*}
f_{\Bin(N, p)}(k) &= \sum_{n=k}^\infty f_{\Geo(q)}(n) \cdot f_{\Bin(n, p)}(k)\\
&= \sum_{i=0}^\infty f_{\Geo(q)}(k+i) \cdot f_{\Bin(k+i, p)}(k)\\
&= \sum_{i=0}^\infty \binom{k+i}{k} q(1-q)^{k+i} \cdot p^k (1-p)^i\\
&= q \cdot p^k \cdot \sum_{i=0}^\infty \binom{k+i}{k} (1-q)^{k+i} (1-p)^i\\
&= \frac{q \cdot p^k (1-q)^k}{(p + q -q \cdot p )^{k+1}}\\
&= \frac{q}{p + q - q \cdot p} \cdot \left( \frac{p (1-q)}{p + q -q \cdot p} \right)^k\\
&= f_{\Geo\left(\frac{q}{p + q - q \cdot p}\right)}(k)
\end{align*}
$$

Curiously, a $p$-thinned geometric distribution is still geometric!

### Negative binomial number of trials

Let $N \sim \NB(r, q)$, where $\NB$ refers to the negative binomial distribution (with support 
on the non-negative integers) whose PMF is

$$
f_{\NB(r, q)}(k) = \binom{k + r - 1}{k} (1-q)^k q^r.
$$

The negative binomial distribution can be seen as a generalization of the geometric distribution with $r = 1$,
as it models the number of failed trials until $r$ successes, for a sequence of events with success probability $q$.
If $X_1, \cdots, X_r$ are i.i.d. $\Geo(q)$, then $\NB(r, q) \sim \sum_{i=1}^r X_i$.

At least for integer $r$, this allows us to show that

$$
\begin{align*}
\Bin(N, p) &\sim \Bin\left(\sum_{i=1}^r \Geo(q), p\right)\\
&\sim \sum_{i=1}^r \Geo\left(\frac{q}{p + q - q \cdot p}\right)\\
&\sim \NB\left(r, \frac{q}{p + q - q \cdot p}\right).
\end{align*}
$$

The third step follows because for non-negative $x, y$, $\Bin(x + y, p) \sim \Bin(x, p) + \Bin(y, p)$.

### Poisson number of trials

Let $N \sim \Poi(\lambda)$, where $\Poi$ refers to the Poisson distribution (with support 
on the non-negative integers) whose PMF is

$$
f_{\Poi(\lambda)}(k) = \frac{\lambda^k e^{-\lambda}}{k!}.
$$

The Poisson distribution models the number of events occuring in an interval of time if
events occur independently from one another and at a constant mean rate.

Let's directly try to compute $f_{\Bin(N, p)}$.

$$
\begin{align*}
f_{\Bin(N, p)}(k) &= \sum_{n=k}^\infty f_{\Poi(\lambda)}(n) \cdot f_{\Bin(n, p)}(k)\\
&= \sum_{n=k}^\infty \frac{\lambda^n e^{-\lambda}}{n!} \cdot \binom{n}{k} p^k (1-p)^{n-k}\\
&= \frac{e^{-\lambda} p^k}{k!}\sum_{n=k}^\infty \frac{\lambda^n}{(n-k)!} (1-p)^{n-k}\\
&= \frac{e^{-\lambda} p^k}{k!}\sum_{i=0}^\infty \frac{\lambda^{k+i}}{i!} (1-p)^i\\
&= \frac{p^k \lambda^k e^{-p \lambda}}{k!}\\
&= f_{\Poi(p \lambda)}(k).
\end{align*}
$$

This fact shouldn't be too surprising given what the Poisson distribution models!
Since each event occurs independently from one another, we can consider "thinning" each
one independently in the Poisson model.

## Randomizing the success probability

In the previous section, we randomized the number of trials of a binomial distribution. Here we can
randomize the _success probability_ by first sampling from a random variable $P$ on $[0, 1]$,
then sampling from $\Bin(n, P)$ for some fixed $n$.

The classic example here is when $P \sim \Beta(\alpha, \beta)$ i.e. $P$ is sampled from the beta distribution.
The beta distribution is continuous and has support on $[0, 1]$. Its probability density function (PDF) is

$$
f_{\Beta(\alpha, \beta)}(x) = \frac{x^{\alpha - 1}(1-x)^{\beta - 1}}{\B(\alpha, \beta)},
$$

where $\B(\alpha, \beta) = \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha + \beta)}$ is the beta function.

## Binomial splitting