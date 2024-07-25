---
layout: post
title:  "Private continual counting and lower bounds"
date:   2024-07-25
needsmath: true
---

In this post, we'll discuss the continual counting problem in differential privacy and outline the proof of the
lower bound found in [Dwork & Roth 2014](https://www.cis.upenn.edu/~aaroth/Papers/privacybook.pdf) (Theorem 12.6).
This problem is quite rich, but we'll discuss a simplified version: there is a continual stream of binary data
coming in, and for every new input, we issue a query returning the number of 1's seen so far.

Let $X = (x_1, ... , x_T)$ be the inputs and $Q = (q_1, ... , q_T)$ are the query outputs.
Here is an example, with no privacy:

$X = (0, 0, 1, 0, 1, 1, 1, 0 ,1)$

$Q = (0, 0, 1, 1, 2, 3, 4, 4, 5)$

There are a few simple ways to do this privately:

- Use the Laplace mechanism independently for each $q_i$. Each query $q_i$ would have sensitivity 1, but we need to pay for
    $T$-fold composition, since we would be querying $x_1$, $T$ times. In this case, the variance of each $q_i$ would be
    $O(T^2 / \epsilon^2)$, since the Laplace distribution with scale parameter $\lambda$ has variance $O(\lambda^2)$.
- Use the Laplace mechanism (or binary randomized response, etc) to privitize $X$, and emit each $q_i$ as post-processing the privitized
    vector of inputs. In this case each mechanism has sensitivity 1 and we pay no composition, but emitting each $q_i$
    requires summing $i$ iid random variables, so variance will grow linearly with $i$. For the Laplace mechanism the
    variance will be $O(T / \epsilon^2)$ in the worst case for $q_T$. This is similar to how to do basic summation
    in the _local_ privacy model, since we are adding noise to each individual data point.

It turns out we can do better than these ideas, e.g. with a binary tree datastructure.
The Dwork & Roth book linked above has a description of the basic binary mechanism, which we may
explore in some future post. It enjoys worst case query variance of $O(\log^3 T)$.

## A more sophisticated notion of utility

Above we discussed per-query variance as a means to compare algorithms. In the literature for continual counting,
a different metric is often used: the $\ell_\infty$ error, or _worst case_ error. That is, for noisy outputs
$\hat{Q} = (\hat{q_1}, ..., \hat{q_T})$, we measure
$err_\infty = \left\lVert Q-\hat{Q} \right\rVert_\infty = max_{i \in [0, T]} \vert q_i - \hat{q_i} \vert $. Furthermore, we
want to show that we only exceed this error with some probability $\beta$.

Let's use this utility metric in an example. Take the naive continual counting algorithm outlined above of
outputing each $q_i$ with the Laplace mechanism directly (the one with $T$-fold composition). In this case,
each query output $\hat{q_i} \sim Laplace\left(q_i, \frac{T}{\epsilon}\right)$. Recall also that for $x \sim Laplace(\mu, \lambda)$:

$$
Pr[|x - \mu| \ge t] = e^{-t / \lambda}
$$

Plugging this in for our mechanism, we see that:

$$
Pr[|\hat{q_i} - q_i| \ge t] = e^{-t\ \epsilon / T}
$$

By a union bound, we can say that:

$$
\beta = Pr\left[\bigcup_{i=1}^T |\hat{q_i} - q_i| \ge t\right] \le T e^{-t\ \epsilon / T}
$$

Solving for $err_\infty = t$:

$$
err_{\ell_\infty} \le \frac{T\ \log\left(\frac{T}{\beta}\right)}{\epsilon} = O_{\epsilon, \beta}(T\ \log(T))
$$

A similar argument for the binary mechanism shows that it enjoys $err_\infty \le O_{\epsilon, \beta}(\log^{2.5} T)$.

## Lower bounds on error

For the remainder of this post we will focus on _lower bounds_, i.e. what kind of error is unavoidable
for this problem if we want to satisfy differential privacy? As far as I can tell, for pure differential
privacy the best bound we know is originally from
[Dwork et al 2010](https://guyrothblum.wordpress.com/wp-content/uploads/2014/11/dnpr10.pdf)
(theorem 4.2) where they show at least $\Omega(\log{T})$ error on this problem.
The Dwork & Roth book's proof is a bit clearer. We'll cover that one.

<a id="thm1" href="#thm1">Theorem 1</a>: There is no $\epsilon=1$ algorithm for continual counting
with $err_\infty \le \frac{\log T}{4} = k$ that holds with probability at least $1 - \beta = 2/3$.

**Proof**: We proceed with a proof by contradiction. Assume we _do_ have an algorithm $A$ as described above.
We'll construct a set of $\frac{T}{k}$ inputs $S = (x_1, ... , x_{T/k})$ by dividing the $T$ time periods into $T/k$ consecutive phases of length $k$. Each input $x_i$ will have the $i$th phase set to all 1s, with the rest set to zeros:

$x_i = 0^{k i} \circ 1^k \circ 0^{k ((T/k) - (i + 1))}$

An output $x$ _matches_ $i$ if, just before the $i$th phase, the output is less than $k/2$,
and at the end of the $ith$ phase, the output is at least $k/2$. Formally, for an output $\hat{Q}$:

$$
matches(Q, i) = \hat{q}_{k i - 1} \lt k/2 \land \hat{q}_{k (i + 1) - 1} \ge k/2
$$

Let $O_i = \left\\{Q : matches(Q, i)\right\\}$. From our utility assumptions:

$$
\begin{align*}
Pr[A(x_i) \in O_i] \ge 1 - \beta
\end{align*}
$$

From $\epsilon$ differential privacy (and group privacy), we know for all $j \neq i$, and set of outputs $O_j$ that _match_ $x_j$:

$$
\begin{align*}
Pr[A(x_j) \in O_j] &\le e^{2 k \epsilon} Pr[A(x_i) \in O_j] \\
Pr[A(x_i) \in O_j] &\ge e^{-2 k \epsilon} Pr[A(x_j) \in O_j] \\
&\ge e^{-2 k \epsilon}(1 - \beta)
\end{align*}
$$

Furthermore for any $i \neq j$, $O_j$ and $O_i$ are disjoint. i.e. any output _matching_ $i$ cannot _match_ $j$. This means that $\sum_{j} Pr[A(x_i) \in O_j] \le 1$.

$$
\begin{align*}
\sum_{j} Pr[A(x_i) \in O_j] &\ge (1-\beta) + \sum_{j \neq i} Pr[A(x_i) \in O_j]\\
&\ge (1-\beta) + \left(\frac{T}{k}-1\right) e^{-2 k \epsilon} (1 - \beta)\\
&\ge \frac{T}{k} e^{-2 k \epsilon} (1 - \beta)\\
\end{align*}
$$

Setting[^1] $\beta=\frac{1}{3}$, $k = \frac{\log{T}}{4}$ and $\epsilon = 1$:

$$
\begin{align*}
\sum_{j} Pr[A(x_i) \in O_j] &\ge \frac{2}{3} \frac{T}{k}\ exp\left(- 2 k \epsilon\right)\\
&= \frac{2}{3} \frac{4 T}{\log(T)}\ exp\left(- \log(T) / 2\right)\\
&= \frac{2}{3} \frac{4 T}{\log(T)} \frac{1}{\sqrt{T}}\\
&= \frac{8 \sqrt{T}}{3 \log(T)}\\
&\gt 1 && \text{ for $T > 1$}
\end{align*}
$$

So $1 \ge \sum_{j} Pr[A(x_i) \in O_j] \gt 1$, leading to a contradiction.

#### Poking at the bound

It's instructive to try to see how this proof changes with tweaks to $k$. For example,
linearly scaling $k = \log(T) / 2$ gives us 
$\sum_{j} Pr[A(x_i) \in O_j] \ge \frac{4}{3 \log(T)}$ which does not lead to a contradiction for $T>e^{4/3}$.
This indicates that the $\log$ lower bound is probably the best we can show with
this construction.

Other notes related to this lower bound:
- [Dwork et al 2015](https://www.iacr.org/archive/asiacrypt2015/94520352/94520352.pdf) showed an algorithm
    which achieves $O(\log{T} + \log^2{n})$ error where $n$ is the number of 1's in the dataset.
- [Cohen et al 2024](https://arxiv.org/pdf/2403.00028) have an improved lower bound for approximate DP, showing $\Omega(\min(n, \log{T}))$ error. 

[^1]: I believe there is a typo of setting $k = \frac{\log_2{T}}{4}$ in the D&R book, when they
    should be using the natural log.