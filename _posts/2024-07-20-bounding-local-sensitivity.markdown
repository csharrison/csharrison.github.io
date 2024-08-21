---
layout: post
title:  "Privately bounding local sensitivity"
date:   2024-07-20
needsmath: true
tag: privacy
---

Typically in differentially private mechanisms, we scale noise proportional to the _global sensitivity_ of the measurement:

$$
GS(f) = \max_{x, x':dist(x, x') \le 1} |f(x)-f(x')|
$$

Where $x$ and $x'$ are all possible neighboring datasets. However, this is often overkill, especially in cases where only
"edge case" datasets have high sensitivity. In many cases, like median estimation, the measurement is quite robust to moving to
a local neighbor, and we would prefer to scale noise proportional to this "local" sensitivity.

It turns out this is possible in some cases, if we relax our notion of privacy to approximate differential privacy.

## Scaling noise to local sensitivity

The _local sensitivity_ of a function $f$ with a given input $x$ is defined as:

$$
LS(f, x) = \max_{x':dist(x, x') \le 1} |f(x)-f(x')|
$$

i.e. we only look at neighbors to a _specific_ input database, not all possible databases. The biggest challenge
with scaling noise to local sensitivity is any measurement that reveals $LS(f, x)$ can reveal something about $x$.
The general approach we will take is to _privately_ measure an upper bound on $LS(f, x)$, $\hat{\Delta}$ and scale
noise according to $\hat{\Delta}$.

<a href="#thm1" id='thm1'>Theorem 1:</a> Let $Pr[g(f, x) \le LS(f, x)] \ge 1-\delta$ for some $\epsilon_1$-differentially private function $g$
that estimates the upper bound of the local sensitivity. The Laplace mechanism over $f$ with scale
parameter $g(f, x) / \epsilon_2$ satisfies $(\epsilon_1 + \epsilon_2, \delta)$-differential privacy.

**Proof:** Recall the typical proof for the Laplace mechanism:

$$
\begin{align*}
\frac{p_{A(D)}(x)}{p_{A(D')}(x)} &\le \frac{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D_i|}{b}\right)}{\prod_{i} \frac{1}{2 b} exp\left(-\frac{|x-D'_i|}{b}\right)} \\
&= \prod_{i} exp\left(\frac{|x-D'_i| - |x - D_i|}{b}\right) \\
&\le \prod_{i} exp\left(\frac{|D_i-D'_i|}{b}\right) &&\text{by triangle inequality}\\
&\le exp\left(\frac{1}{b} LS(f, D) \right)\\
&\le exp\left(\frac{1}{b} g(f, D) \right) && \text{conditioned on the event that $LS(f, D) \le g(f, D)$}\\
&= e^{\epsilon_2}
\end{align*}
$$

And so, conditioning on the event $E$ that $LS(f, D) \le g(f, D)$:

$$
\begin{align*}
Pr[A(D) \in S] &= Pr[A(D) \in S | E] Pr(E) + Pr[A(D) \in S | \neg E] Pr(\neg E)\\
&\le Pr[A(D) \in S | E] Pr(E) + \delta\\
&\le e^{\epsilon_2} Pr[A(D') \in S | E] Pr(E) + \delta\\
&\le e^{\epsilon_2} Pr[A(D') \in S] + \delta
\end{align*}
$$

The theorem follows from basic composition between $\epsilon_1$ and $(\epsilon_2, \delta)$.


## Private ratios

Let's use this techniquee to privately compute the ratio $a/b$, when users can either contribute a value in {0, 1} to $b$, or $a$ and $b$.
The simplest approach[^1] in this case would be to compute $a/b$ by privately querying both $a$ and $b$ separately with the Laplace mechanism[^2].

```python
def private_ratio_naive(a, b, epsilon):
    noisy_a = scipy.stats.laplace.rvs(loc=a, scale=2/epsilon)
    noisy_b = scipy.stats.laplace.rvs(loc=b, scale=2/epsilon)
    return noisy_a / noisy_b
```

By basic composition, this mechanism satisfies $\epsilon$ differential privacy. Can we do better?
We'll try scaling noise to local sensitivity. Fix a given dataset $x$ which has the given ratio $a/b$. It has four possible neighbors:

- $x'_1 = \frac{a+1}{b+1}$
- $x'_2 = \frac{a}{b+1}$
- $x'_3 = \frac{a-1}{b-1}$
- $x'_4 = \frac{a}{b-1}$

If we assume $b > 1$, we can make some claims about the local sensitivity of $A(x)=a/b$.

$$
\begin{align*}
LS(f, x) &= \max_{x':dist(x, x') \le 1} |A(x)-A(x')|\\
&= max\left(\left|\frac{a}{b} - \frac{a+1}{b+1}\right|, \left|\frac{a}{b} - \frac{a}{b+1}\right|, \left|\frac{a}{b} - \frac{a-1}{b-1}\right|, \left|\frac{a}{b} - \frac{a}{b-1}\right|\right)\\
&= max\left(\left|\frac{a}{b} - \frac{a-1}{b-1}\right|, \left|\frac{a}{b} - \frac{a}{b-1}\right|\right) &\text{these terms dominate}\\
&= max\left(\left|\frac{a-b}{b-b^2}\right|, \left|\frac{a}{b-b^2}\right|\right)\\
&= max\left(\frac{b-a}{b^2 - b}, \frac{a}{b^2-b}\right)\\
\end{align*}
$$

The best we can hope to do is scale noise proportional to this quantity, but we'll proceed by computing an upper bound[^3].

## An upper bound on local sensitivity

We know $LS(f, x) \le max\left(\frac{b-a}{b^2 - b}, \frac{a}{b^2-b}\right)$ as long as $b>1$. This points to a possible way to upper bound $LS$:

1. Privately query $a$ and $b$ separately, which are both sensitivity 1 queries.
1. For the second term, we know this is maximized by taking high probability upper bounds of $a$ ($a_u$) and lower bounds of $b$ ($b_l$).
1. For the first term, we know this is maximized by taking lower bounds of $a$ ($a_l$). However, choosing bounds for $b$ is a little tricky since the term is non-linear.
    However, we should be able to upper bound this quantity by solving for the maximum within the valid range. 

$$
\begin{align*}
\frac{d}{db} \frac{b-a}{b^2 - b} &= \frac{a(2b - 1)-b^2}{(b-1)^2 b^2} = 0\\
a(2b - 1) - b^2 &= 0\\
b &= \pm \sqrt{a^2-a} + a
\end{align*}
$$

Note that because $b \ge a$, the solution $-\sqrt{a^2-a} + a$ is invalid.
So we can just take integer solution around $b_x = \sqrt{a^2-a} + a$. If
this value is not within bounds, it suffices to just check the endpoints of the bounds of $b$.

$$
\max_{b \in [b_l,b_u]} \frac{b-a_u}{b^2 - b}  = \left\{\begin{array}{lr}
    \frac{b_x-a_u}{b_x^2 - b_x}, & \text{if } b_x \in [b_l, b_u]\\
    \max(\frac{b_l-a_u}{b_l^2 - b_l}, \frac{b_u-a_u}{b_u^2 - b_u}) & \text{otherwise}
    \end{array}\right\}
$$

This points to a possible mechanism:

```python
def private_ratio_local_sensitivity(a, b, epsilon, delta, epsilon_1_fraction=0.1):
    epsilon_1 = epsilon * epsilon_1_fraction
    epsilon_2 = epsilon - epsilon_1    

    # Use half of epsilon_1 for a and b querying.
    noisy_a = scipy.stats.dlaplace.rvs(loc=a, a=epsilon_1/2)
    noisy_b = scipy.stats.dlaplace.rvs(loc=b, a=epsilon_1/2)

    a_upper = ceil(noisy_a + 2 * np.log(2/delta)/epsilon_1)
    a_lower = max(0, floor(noisy_a - 2 * np.log(2/delta)/epsilon_1))
    b_upper = ceil(noisy_b + 2 * np.log(2/delta)/epsilon_1)
    b_lower = floor(noisy_b - 2 * np.log(2/delta)/epsilon_1)

    # Bail out since our approach fails for b <= 1
    if b_lower <= 1:
        return private_ratio_naive(dataset, epsilon_2)

    bound_1 = a_lower / (b_lower**2 - b_lower)

    bx = sqrt(a_lower**2 - a_lower) + a_lower
    bound_2 = max((bi-a_upper) / (bi**2 - bi) 
        for bi in [b_upper, b_lower, ceil(bx), floor(bx)]
        if bi <= b_upper and bi >= b_lower)

    local_sensitivity_upper_bound = max(bound_1, bound_2)
    return scipy.stats.laplace.rvs(loc=a/b, scale=local_sensitivity_upper_bound / epsilon_2)
```

<a id="thm2" href="#thm2">Theorem 2:</a> `private_ratio_local_sensitivity` is $(\epsilon, \delta)$ differentially private.

Proof: We are doing basic composition across 3 queries:

1. The query for $a$, yielding $\hat{a}$
1. The query for $b$, yielding $\hat{b}$
1. The query for $a/b$

The first two are $\frac{\epsilon_1}{2}$ differentially private via basic usage of the Discrete Laplace mechanism (each query is sensitivity 1).
It suffices to prove privacy for the last one. From <a href="#thm1">Theorem 1</a>, it suffices to show that `local_sensitivity_upper_bound` is a
valid upper bound for $LS(f, D)$ with probability at least $1-\delta$.

First, recall that $Pr(|x| \ge T) = e^{-T a}$ for $x \sim DLaplace(a)$,
and therefore 

$$
\begin{align*}
Pr(|\hat{a} - a| \ge 2 log(2/\delta)/\epsilon_1) &= exp(- 2 a\ log(2/\delta)/\epsilon_1)\\
&= exp(-log(2/\delta))\\
&= \frac{\delta}{2}\\
\end{align*}
$$

By the same logic, $Pr(|\hat{b} - b| \ge 2 log(2/\delta)/\epsilon_1) = \frac{\delta}{2}$.
Therefore,

$$
Pr(|\hat{a} - a| \ge 2 log(2/\delta)/\epsilon_1 \land |\hat{b} - b| \ge 2 log(2/\delta)/\epsilon_1) \le \delta
$$

due to independence. This means we know $a \in [a_l, a_u]$ and $b \in [b_l, b_u]$ with probability $1-\delta$.
From the bounds derived above, with $g(f, D) = max\left(\frac{a_u}{b_l^2 - b_l}, \max_{b \in [b_l, b_u]} \frac{b-a_l}{b^2 - b}\right)$:
$Pr[g(f, D) \le LS(f, D)] \le \delta$.

<a href="#thm1">Theorem 1</a> tells us that the final query for $a/b$ satisfies $(\epsilon_2, \delta)$-DP. Basic composition across
the three queries completes the proof[^4].


## Some empirical results

We first plot[^5] the mean absolute error of the naive approach vs. our local sensitivity attempt, across 100,000 trials.

![Figure 1](/images/mae-1.svg)

As is shown in the plot, we can "afford" to spend less privacy budget on the query privately bounding the local sensitivity when the
sensitivity is very "stable". Even still, we handily outperform the naive mechanism and come pretty close to the non-private baseline.

How does our attempt stack up to some state-of-the-art results? We'll compare to the simple Algorithm 3 from [Kulesza, Suresh, and Wang 2024](https://openreview.net/pdf?id=cwIhvoTzuK).
This genius algorithm exploits the fact that in the add-remove adjacency relationship, a single two-dimensional query of the form $(a, b-a) $ has sensitivity 1.
The algorithm looks something like this (slightly tweaked for our simplified setting):

```python
def private_ratio_KSW(a, b, epsilon):
    noisy_ones = scipy.stats.laplace.rvs(loc=a, scale=1/epsilon)
    noisy_zeros = scipy.stats.laplace.rvs(loc=b-a, scale=1/epsilon)
    return noisy_ones / (noisy_ones + noisy_zeros)
```

![Figure 3](/images/mae-3.svg)

![Figure 2](/images/mae-2.svg)

The KSW algorithm has some clear advantages:

- It outperforms our algorithm in many regmies, especially when $b$ is small and our upper bounds on the local sensitivity become very loose.
- Its performance is predictable.
- It is simple, and doesn't require so many finicky bits like we have, e.g. our "bail out" strategy when our lower bounds on $b$ are too small.
- It does not require any hyperparameter tuning, unlike our algorithm which requires figuring out the optimal budget split for bounding the sensitivity.
- It satisfies pure DP, rather than our algorithm which relaxes to approximate DP.

Still, it is nice to see the approach to privately bound the local sensitivity do well in some regimes!
Our algorithm also does not have a degradation in the substitution relation unlike the KSW algorithm which requires querying
with $\lambda = 2/\epsilon$, as the sensitivity of their two-dimensional query becomes 2.

## Conclusion

I hope this helped show privately bounding local sensitivity can be a powerful tool for desigining differentially private
mechanisms! Special thanks to Gautam Kamath's [public lectures and lecture notes](http://www.gautamkamath.com/CS860-fa2020.html)
on this topic for inspiration.

[^1]: Note here and elsewhere (in the KSW algorithm) we could use some other noise mechanism (e.g. Discrete Laplace, staircase) for better performance, but we use the regular
    Laplace mechanism for simplicity and for easier comparisons. In particular, the KSW paper proposes a new hourglass mechanism.

[^2]: Like some of our previous posts, these are not production-grade implementations and are just for reference.
    For example, using `scipy.stats.laplace` naively will have privacy bugs outlined in [Mironov 2012](https://www.microsoft.com/en-us/research/wp-content/uploads/2012/10/lsbs.pdf)

[^3]: Note that the quantity here is upper bounded by $1/(b-1)$, which should allow us to generate an upper bound
    Without querying two values. However, it makes it impossible for us to optimize for the reduced sensitivity as $a$ approaches
    $b/2$.

[^4]: We skip some details for the "bailing out" approach when $b_l \le 1$, but it should be clear this path is also $(\epsilon_1 + \epsilon_2)$-DP from basic composition.

[^5]: Plots generated in [this colab notebook](https://colab.research.google.com/drive/14z-rLkL1pTt9CYQbCoOvwDYdJM0GnDvH?usp=sharing)