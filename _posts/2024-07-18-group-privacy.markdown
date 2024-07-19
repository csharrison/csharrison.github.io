---
layout: post
title:  "Group privacy"
date:   2024-07-18
needsmath: true
---

Within the context of differential privacy, group privacy is the notion of protecting _groups_ of users, rather than the "standard" notion of privacy which protects a single user. I've always considered this aspect of differential privacy somewhat confusing, especially as it relates to DP composition:

- Composition is running multiple mechanisms over the same dataset
- Group privacy is protecting multiple users / rows in the same dataset

Intuitively, both composition and achieving group privacy seems like they ought to amplify the privacy parameter (e.g. $\epsilon$). However, they differ in subtle respects. In this post we will outline how group privacy works under a few different DP definitions, and give some intution about some of the properties of group privacy (across privacy defintions and vs. composition).

## Pure differential privacy

"Protecting multiple users" in differential privacy relates to the notion of adjacent datasets. Typically we want to show, for two neighboring datasets $D$ and $D'$ that differ on one row:

$$
Pr[A(D) \in S] \le e^\epsilon Pr[A(D’) \in S]
$$

Under $k$ group privacy, we simply extend then notion of a neighbor to be a database that differs on $k$ rows. It is easy to see that a mechanism satisfying $\epsilon$-DP also satisfies $k \epsilon$ group privacy.

Consider a database $D^k$ that differs on $k$ records from $D$. Let $D^1, D^2, ... D^k$ represent the series of intermediate 
databases, each differing on only a single record from the next, needed to get to $D^k$ from $D$.

$$
Pr[A(D) \in S] \le e^\epsilon Pr[A(D^1) \in S] \le e^{2 \epsilon} Pr[A(D^2) \in S] ... \le e^{k \epsilon} Pr[A(D^k) \in S]
$$

## Approximate differential privacy

The story is largely similar with approximate differential privacy, except here we have an additional term $\delta$ to worry about.
A mechanism $A$ satisfies $(\epsilon, \delta)$ approximate differential privacy if:

$$
Pr[A(D) \in S] \le e^\epsilon Pr[A(D') \in S] + \delta
$$

Group privacy bounds follow directly from the definition:

$$
\begin{align*}
Pr[A(D) \in S] &\le e^\epsilon Pr[A(D^1) \in S] + \delta \\
&\le e^\epsilon\left(e^\epsilon Pr[A(D^2) \in S] + \delta\right) + \delta \\
&\le e^\epsilon\left(e^\epsilon \left( e^\epsilon Pr[A(D^3) \in S] +\delta \right) + \delta\right) + \delta \\
&...\\
&\le e^{k \epsilon} Pr[A(D^k) \in S] +  \delta (1 + e^\epsilon + e^{2 \epsilon} + ... + e^{(k-1)\epsilon}) \\
&= e^{k \epsilon} Pr[A(D^k) \in S] + \delta \frac{e^{k \epsilon} - 1}{e^\epsilon - 1}
\end{align*}
$$

Therefore any $(\epsilon, \delta)$-DP mechanism should give $(k \epsilon, \delta \frac{e^{k \epsilon} - 1}{e^\epsilon - 1})$ group privacy. See Lemma 2.2 of the [Complexity of Differential Privacy](https://privacytools.seas.harvard.edu/files/complexityprivacy_1.pdf) for a similar proof.

## zero-concentrated differential privacy

A mechanism $A$ satisfies $\rho$-zCDP if, for all $\alpha \in (1, \infty)$ and all neighboring datasets $D$ and $D'$:
$$
D_\alpha(A(D) || A(D')) \le \rho \alpha
$$

Where $D_\alpha$ is the Renyi Divergence of order $\alpha$. The zCDP definition also protects groups. To prove this we will use a technical lemma (Lemma 5.2 in [Bun & Steinke 2016](https://arxiv.org/abs/1605.02065)). Let $P, Q, R$ be probability distributions, then:

$$
\begin{equation}
D_\alpha(P||Q) \le \frac{k \alpha}{k \alpha - 1}D_{\frac{k \alpha -1 }{k - 1}}(P||R) + D_{k \alpha}(R||Q)\\
\end{equation}
$$

for all $k, \alpha \in (1, \infty)$. This is a "triangle-like" inequality for Renyi Divergence.

We can use this lemma to iteratively show group privacy like we did for pure and approximate dp.
We will show that a $\rho$-zCDP mechanism satifies $k^2 \rho$ zCDP for groups of size $k$. Like Bun & Steinke we will proceed via induction on $k$, although unlike their proof we will only concern ourself with $\rho$-zCDP mechanisms rather than $(\xi, \rho)$-zCDP mechanisms, making the proof a bit simpler.

The base case is trivial, zCDP trivially protects groups of size $k=1$. We will assume the induction hypothesis, that:

$$
D_\alpha(A(D)||A(D^{k-1})) \le (k-1)^2 \rho \alpha
$$

It remains to show that the if groups of size $k-1$ are protected with $(k-1)^2 \rho$-zCDP, then groups of size $k$ are protected with $k^2 \rho$-zCDP:

$$
\begin{align*}
D_\alpha(A(D)||A(D^k)) &\le \frac{k \alpha}{k \alpha - 1}D_{\frac{k \alpha -1 }{k - 1}}(A(D)||A(D^{k-1})) + D_{k \alpha}(A(D^{k-1})||A(D^k)) && \text{By the triangle-like inequality}\\
&\le \frac{k \alpha}{k \alpha - 1}D_{\frac{k \alpha -1 }{k - 1}}(A(D)||A(D^{k-1})) + k \rho \alpha && \text{By zCDP}\\
&\le \frac{k \alpha}{k \alpha - 1} (k-1)^2 \rho \frac{k \alpha -1 }{k - 1} + k \rho \alpha && \text{By the induction hypothesis}\\
&= \rho \left(\frac{k \alpha}{k \alpha - 1} (k-1)^2 \frac{k \alpha -1 }{k - 1} + k \alpha\right)\\
&= \rho \left(\frac{k \alpha}{k - 1} (k-1)^2 + k \alpha\right)\\
&= \rho \left(k \alpha (k-1) + k \alpha\right)\\
&= \rho \left(k^2 \alpha - k\alpha + k \alpha\right)\\
&= k^2 \rho \alpha\\
\end{align*}
$$

## Why do non-pure DP definitions have "worse" group privacy properties?

To protect groups of size $k$, pure DP enjoys a linear increase in $\epsilon$, whereas approximate DP and zCDP have worse blowups in their privacy parameters. Why is that?

One way to think about this is via _sensitivity_. Most private mechanisms tune noise to the sensitivity of the measurement i.e.
how much any one person can contribute. For instance, the Laplace mechanism with standard deviation $\sigma$ satisfies $\frac{\sqrt{2} \Delta_1}{\sigma}$ differential privacy, where $\Delta_1$ is the $\ell_1$ sensitivity of the measurement. If we think of group privacy as
allowing one person to contribute to the measurement with the inputs of $k$ people, the effective sensitivity is $k \Delta_1$, and we see a linear blowup satisfying $\frac{\sqrt{2} k \Delta_1}{\sigma}$ differential privacy.

On the other hand, the Gaussian mechanism with standard deviation $\sigma$ satisfies $\rho = \frac{\Delta_2^2}{2 \sigma^2}$ zCDP, Where $\Delta_2$ is the $\ell_2$ sensitivity. If one person is allowed to contribute with the inputs of $k$ other people, the $\ell_2$ sensitivity is increased linearly by $k$. In that case, $\rho_k = \frac{(k \Delta_2)^2}{2 \sigma^2} = k^2 \rho$.

## Group privacy vs. composition

It is tempting to think of group privacy "as if" we amplify a person's contributions by querying them $k$ times and sum up each output, but this is _not_ the case. We can use the Gaussian mechanism again to showcase the difference, where each person inputs a vector with sensitivity $\Delta_2$.

$k$-fold composition and summing is similar to running the Gaussian mechanism with scale parameter $\sqrt{k} \sigma$ on a vector whose $\ell_2$ sensitivity is $k \Delta_2$. This is because we are adding $k$ independent samples of the Gaussian for each element, and the result follows from the fact that variances add up. It makes sense in that case that $k$-fold composition will result in $k \rho$-zCDP, since $\rho' = \frac{(k \Delta_2)^2}{(\sqrt{k} \sigma)^2} = \frac{k \Delta_2^2}{\sigma^2}= k\rho$.

On the other hand, using group privacy we only add noise a single time, with scale parameter $\sigma$, even though the underlying $\ell_2$ sensitivity is $k \Delta_2$. This will result in $k^2 \rho$-zCDP as shown above.

## Group privacy superpower: bounding multiple privacy units

A useful property of group privacy is that it allows us to more easily protect coarser grained "privacy units". For instance, if there is a mechanism
that satisfies $\epsilon$-DP for individual records in a database, but users can contribute $k$ records, group privacy automatically gives us a user-level privacy bound. Similarly if a mechanism protects a person's contribution for a day, but a new data release occurs every day, group privacy can give bounds on the privacy of all a person's contributions across multiple days / weeks. 