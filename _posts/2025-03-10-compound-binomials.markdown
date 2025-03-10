---
layout: post
title:  "Compound binomial distributions and binomial splitting"
date:   2025-03-10
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

A [_compound_ probability distribution](https://en.wikipedia.org/wiki/Compound_probability_distribution) (also sometimes called a scale mixture, or a [parameter mixture distribution](https://reference.wolfram.com/language/ref/ParameterMixtureDistribution.html)) is
formed when one of the _parameters_ of the distribution is another random variable. We'll consider compounding the
binomial distribution by randomizing both the number of trials $n$, and the success probability $p$.

## Binomial thinning: randomizing the number of trials

Let's start with randomizing the number of trials.
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
&= q \cdot p^k (1-q)^k \sum_{i=0}^\infty \binom{k+i}{k} ((1-q)(1-p))^{i}\\
\text{(negative binomial series)}&= \frac{q \cdot p^k (1-q)^k}{(p + q -q \cdot p )^{k+1}}\\
&= \frac{q}{p + q - q \cdot p} \cdot \left( \frac{p (1-q)}{p + q -q \cdot p} \right)^k\\
&= f_{\Geo\left(\frac{q}{p + q - q \cdot p}\right)}(k)
\end{align*}
$$

Curiously, a $p$-thinned $\Geo(q)$ distribution is just $\Geo\left(\frac{q}{p + q - q \cdot p}\right)$.

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
\text{(Maclaurin series for $e^x$)} &= \frac{p^k \lambda^k e^{-p \lambda}}{k!}\\
&= f_{\Poi(p \lambda)}(k).
\end{align*}
$$

This fact shouldn't be too surprising given what the Poisson distribution models!
Since each event occurs independently from one another, we can consider "thinning" each
one independently in the Poisson model.

### Binomial number of trials

Let's go meta, and consider $N \sim \Bin(n, q)$, i.e. we are $p$-thinning a binomial distribution itself.
Note that $q$ is the success probability of the input random variable, and $p$ the success probability of the
thinning process. In this case

$$
\begin{align*}
f_{\Bin(N, p)}(k) &= \sum_{x=k}^n f_{\Bin(n, q)}(x) \cdot f_{\Bin(x, p)}(k)\\
&= \sum_{x = k}^n \binom{n}{x} q^x (1-q)^{n-x} \cdot \binom{x}{k} p^k (1-p)^{x-k}\\
&= \sum_{i = 0}^{n-k} \binom{n}{k + i} q^{k + i} (1-q)^{n-k-i} \cdot \binom{k+i}{k} p^k (1-p)^{i}\\
&= \sum_{i = 0}^{n-k} \frac{n!}{(k+i)!(n-k-i)!} q^{k + i} (1-q)^{n-k-i} \cdot \frac{(k+i)!}{k!\cdot i!} p^k (1-p)^{i}\\
&= \sum_{i = 0}^{n-k} \frac{n!}{(n-k-i)!} q^{k + i} (1-q)^{n-k-i} \cdot \frac{1}{k!\cdot i!} p^k (1-p)^{i}\\
&= \frac{n!}{k!}(q \cdot p)^k(1-q)^{n-k} \sum_{i = 0}^{n-k} \frac{1}{(n-k-i)! \cdot i!} \left(\frac{q (1-p)}{1-q}\right)^i\\
\text{(binomial theorem)} &= \frac{n!}{k!}(q \cdot p)^k(1-q)^{n-k} \cdot \frac{1}{(n-k)!} \cdot \left(1+\frac{q (1-p)}{1-q}\right)^{n-k}\\
&= \frac{n!}{k! (n-k)!}(q \cdot p)^k(1- q \cdot p)^{n-k} \\
&= f_{\Bin(n, q \cdot p)}.
\end{align*}
$$

This should make sense. If the binomial process filters each event independently with probability $p$,
then stacking $S$ binomial processes on top of each other should be equivalent to filtering each event with success
probability equal to $\prod_{i=1}^S p_i$, since each event needs to pass through all filters to be accepted.

## Randomizing the success probability

In the previous section, we randomized the number of trials of a binomial distribution. Here we can
randomize the _success probability_ by first sampling from a random variable $P$ on $[0, 1]$,
then sampling from $\Bin(n, P)$ for some fixed $n$.

The classic example here is when $P \sim \Beta(\alpha, \beta)$ i.e. $P$ is sampled from the beta distribution.
The beta distribution is continuous and has support on $[0, 1]$. Its probability density function (PDF) is

$$
f_{\Beta(\alpha, \beta)}(x) = \frac{x^{\alpha - 1}(1-x)^{\beta - 1}}{\B(\alpha, \beta)},
$$

where $\B(\alpha, \beta) = \int_{0}^1 t^{\alpha-1} (1-t)^{\beta-1} \mathrm{d}t = \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha + \beta)}$ is the beta function.
Evaluating its PMF, we find that

$$
\begin{align*}
f_{\Bin(n, P)}(k) &= \int_{0}^1 f_{\Beta(\alpha, \beta)}(p) \cdot f_{\Bin(n, p)}(k) \mathrm{d}p\\
&= \frac{\binom{n}{k} }{\B(\alpha, \beta)}\int_{0}^1 p^{\alpha - 1}(1-p)^{\beta - 1} p^k (1-p)^{n-k} \mathrm{d}p\\
&= \frac{\binom{n}{k} }{\B(\alpha, \beta)}\int_{0}^1 p^{k + \alpha - 1}(1-p)^{n-k + \beta - 1}\mathrm{d}p\\
&= \binom{n}{k} \frac{\B(k + \alpha, n -k +\beta)}{\B(\alpha, \beta)} \\
&= f_{\BetaBin(n, \alpha, \beta)}(k).
\end{align*}
$$

The binomial distribution with beta distributed success probability is the _beta binomial_ distribution!
If we have $X_1 \sim \NB(r_1, p)$ and $X_2 \sim \NB(r_2, p)$, then $X_1 | X_1 + X_2 \sim \BetaBin(n, r_1, r_2)$
i.e. the beta binomial is the negative binomial, conditional on its sum with another negative binomial.

$$
\begin{align*}
\Pr[X_1 = x | X_1 + X_2 = n] &= \frac{\Pr[X_1 + X_2 = n | X_1 = x] Pr[X_1 = x]}{\Pr[X_1 + X_2 = n]}\\
&= \frac{\Pr[X_2 = n - x] Pr[X_1 = x]}{\Pr[X_1 + X_2 = n]}\\
&= \frac{f_{\NB(r_2, p)}(n-x) \cdot f_{\NB(r_1, p)}(x)}{f_{\NB(r_1 + r_2, p)}(n)}\\
&= \frac{\binom{n-x + r_2 - 1}{n-x} (1-q)^{n-x} q^{r_2} \cdot \binom{x + r_1 - 1}{x} (1-q)^x q^{r_1}}{\binom{n + r_1 + r_2 - 1}{n} (1-q)^n q^{r_1 + r_2}}\\
&= \frac{\binom{n-x + r_2 - 1}{n-x} \cdot \binom{x + r_1 - 1}{x}}{\binom{n + r_1 + r_2 - 1}{n}}\\
&= \frac{(n-x + r_2-1)!}{(n-x)! (r_2 - 1)!} \cdot \frac{(x + r_1 - 1)!}{x! (r_1 - 1)!} \cdot \frac{n! (r_1 + r_2 - 1)!}{(n + r_1 + r_2 - 1)!}\\
&= \binom{n}{x} \cdot \frac{\Beta(x+r_1, n-x+r_2)}{\Beta(r_1, r_2)} \\
&= f_{\BetaBin(n, r_1, r_2)}(x)
\end{align*}
$$

The beta binomial is the one-dimensional version of the [Dirichlet multinomial distribution](https://en.wikipedia.org/wiki/Dirichlet-multinomial_distribution) which is the distribution
of _all_ $X_1 + \cdots + X_n$ conditioned on their sum.


## Binomial splitting and the magic of the Poisson distribution

In this section we'll explore a related concept to binomial thinning, which we'll call _splitting_. We will consider the
two-dimensional random variable $(X, Y)$, where $X \sim \Bin(N, p)$ and $ Y \sim N - X$. In other words, we will have a
process which classifies each value in $N$ into the first index with probability $p$, and the second index otherwise.
Note that this is a special case of the (compound) [multinomial distribution](https://en.wikipedia.org/wiki/Multinomial_distribution).

For the most part, we can't say anything too meaningful about splitting a distribution in this way. In general, $X$ and $Y$
will be correlated in some way (i.e. learning something about $X$ tells you something about $Y$), making the split distribution
unweildy to analyze. However, there is one distribution where splitting leads to _independent_ $X$ and $Y$, the Poisson distribution!

We know from above that if $N \sim \Poi(\lambda)$, then individually $X \sim \Poi(\lambda p)$ and $Y \sim \Poi(\lambda (1-p))$, but how can
we show independence? Let's consider their joint probability distribution.

$$
\begin{align*}
\Pr[(X, Y) = (x, y)] &= \Pr[N = x+y] \Pr[\Bin(x+y, p) = x]\\
&= \frac{\lambda^{x+y} e^{-\lambda}}{(x + y)!} \cdot \binom{x + y}{x} p^x (1-p)^{y}\\
&= \lambda^{x+y} e^{-\lambda} \cdot \frac{p^x (1-p)^y}{x! y!}\\
&= e^{-\lambda} \cdot \frac{(p \lambda)^{x}}{x!} \cdot \frac{((1- p) \lambda)^y}{y!} \\
&= e^{-\lambda(p + (1-p))} \cdot \frac{(p \lambda)^{x}}{x!} \cdot \frac{((1- p) \lambda)^y}{y!} \\
&= \frac{(p \lambda)^{x} e^{-p \lambda}}{x!} \cdot \frac{((1- p) \lambda)^y e^{-(1-p)\lambda}}{y!} \\
&= f_{\Poi(p \lambda)}(x) \cdot f_{\Poi((1-p)\lambda)}(y)
\end{align*}
$$

Thus, $X$ and $Y$ are independent, as we can write their joint probability distribution as the product of
$f_X$ and $f_Y$. A natural question is whether any other binomially split distributions end up with independence
on each split portion. Let's try to answer this question following the approach in [these lecture notes](https://soumendu041.gitlab.io/notes/Poisson.pdf) from Soumendu Mukherjee.

We'll now assume $N$ is _any_ distribution on the non-negative integers where $\Pr[N = 0] < 1$, with $X \sim \Bin(N, p)$ and $Y \sim N - X$.
Let $G_N(z) = \mathbb{E}(z^N)$ be the [probability generating function](https://en.wikipedia.org/wiki/Probability-generating_function) (PGF) of $N$. Before we move on we will just recall a few facts about PGFs:

1. $G_N(0) = \Pr[N = 0]$
1. $G_N(1) = 1$

Computing the joint PGF of $X$ and $Y$ yields

$$
\begin{align*}
G_{X, Y}(z_1, z_2) &= \mathbb{E}(z_1^X \cdot z_2^Y)\\
&= \sum_{x_1, x_2} f_{X, Y}(x_1, x_2) z_1^{x_1} z_2^{x_2}\\
&= \sum_{x_1, x_2} f_{N}(x_1 + x_2) f_{\Bin(x_1+x_2)}(x_1) z_1^{x_1} z_2^{x_2}\\
&= \sum_{x_1, x_2} f_{N}(x_1 + x_2) \binom{x_1+x_2}{x_1}p^{x_1}(1-p)^{x_2} z_1^{x_1} z_2^{x_2}\\
&= \sum_{x_1, x_2} f_{N}(x_1 + x_2) \binom{x_1+x_2}{x_1} (p \cdot z_1)^{x_1} ((1-p) z_2)^{x_2}\\
&= \sum_{n=0}^\infty f_{N}(n) \sum_{x_1 = 0}^n \binom{n}{x_1} (p \cdot z_1)^{x_1} ((1-p) z_2)^{n - x_2}\\
&= \sum_{n=0}^\infty f_N(n) (p\cdot z_1 + (1-p)z_2)^n\\
&= G_N(p\cdot z_1 + (1-p)z_2).
\end{align*}
$$

If $X$ and $Y$ are independent, it _must_ be the case that

$$
G_{X,Y}(z_1, z_2) = \mathbb{E}[z_1^X z_2^Y] = \mathbb{E}[z_1^X]\cdot\mathbb{E}[z_2^Y] = G_{X,Y}(z_1, 1)\cdot G_{X,Y}(1, z_2).
$$

Or, via the equivalence above,

$$
G_N(p\cdot z_1 + (1-p)z_2) = G_N(p\cdot z_1 + (1-p)) \cdot G_N(p + (1-p)z_2).
$$

Nice, we have a handle on a non-trivial functional relationship for $G_N$. Let's try to an additive relationship by analyzing the related function $g(z) = \log(G_N(z+1))$.

$$
\begin{align*}
g(z_1 + z_2) &= \log(G_N(z_1 + z_2 + 1))\\
&= \log(G_N(p \cdot (x_1+1) + (1-p)(x_2+1)) ) & \text{where } z_1 = p \cdot x_1 \text{ and } z_2 = (1-p)  x_2 \\
&= \log(G_N(p\cdot (x_1 +1) + (1-p))) + \log(G_N(p + (1-p)(x_2 + 1)))\\
&= \log(G_N(p\cdot x_1 + 1)) + \log(G_N((1-p)x_2 + 1))\\
&= g(z_1) + g(z_2)\\
\end{align*}
$$

This satisfies [Cauchy additive functional equation](https://en.wikipedia.org/wiki/Cauchy's_functional_equation),
which is known to be _linear_ under very weak conditions, e.g. $g(z) \ge 0$ whenever $z \ge 0$, which is clearly the case.
Thus we can say $g(z) = \lambda z$ for some positive constant $\lambda$. Why does $\lambda$ need to be positive? Well,
it certainly can't be 0, otherwise, $\log(G_N(z+1)) = 0$ and therefore $G_N(z) = 1$, but $G_N(0) = \Pr[N=0]$ which by construction is less than 1. Similarly, $\lambda$ cannot be negative, as then $g(0) < g(-1)$ and $1 = G_N(1) < G_N(0)  = \Pr[N=0] < 1$ which is also a contradiction.

So where does that leave us? Well if $g(z) = \lambda z$, then $e^{g(z)} = G_N(z+1) = e^{\lambda z}$, so $G_N(z) = e^{\lambda (z-1)}$, which indeed is the PGF for the Poisson distribution! This tells us that the _independence_ of the two random variables
from binomial splitting uniquely characterizes the Poisson distribution. There aren't any other distributions with this
property!