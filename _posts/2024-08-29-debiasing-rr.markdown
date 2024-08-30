---
layout: post
title:  "Debiasing randomized response"
date:   2024-08-29
needsmath: true
tag: privacy
---

In this post I wanted to derive the formula to debias the randomized response / RAPPOR algorithms.
I often find myself reaching for these formulas so it's a good exercise to re-derive them and have a reference.
We'll start with $k$-ary randomized response, which generalizes to the others. Here's a reference implementation
for convenience:

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
  k = len(possible_outputs)
  p = k / (k - 1 + math.exp(epsilon))
  assert value in possible_outputs
  if random.random() <= p:
    return random.choice(possible_outputs)
  return value
```

We showed in a [previous post]({% post_url 2024-03-30-intro-to-dp-proofs %}#k-ary-randomized-response).
I'll now show that it is _biased_. Consider the function $A$ which invokes $k$-rr and encodes its output
as a $k$-length one-hot vector (suitable for e.g. aggregation). Let $x_i$ be the $i$th possible input value.

$$
\begin{align*}
E[A(x_i)_j] &= 
    \begin{cases}
    (1-p) + \frac{p}{k} & i = j\\
    \frac{p}{k} & i \neq j
    \end{cases}\\

&= \begin{cases}
\frac{e^\epsilon}{e^\epsilon + k - 1} & i = j\\
\frac{1}{e^\epsilon + k - 1} & i \neq j
\end{cases}\\
\end{align*}
$$

Hence this is biased since we don't recover the original (one-hot) vector in expectation.
Fixing this is straightforward. Simply post-process $A(x_i)$ by via a translation / scaling function $F$.
We want to:
- Subtract $\frac{p}{k} = \frac{1}{e^\epsilon + k - 1}$ so that the expectation of wrong items is correctly 0
- Rescale by a factor of $c$ so that the expectation of correct items is 1

Let's derive $c$.

$$
\begin{align*}
c \left((1-p) + \frac{p}{k} - \frac{p}{k}\right) &= 1\\\\
c &= \frac{1}{1-p}
\end{align*}
$$

So, with an abuse of notation ($F$ is applied to each element independently):

$$
\begin{align*}
F(x) &= \frac{1}{1-p} \left(x - \frac{p}{k}\right)\\
\end{align*}
$$

Now, as we hoped, we have an unbiased output:

$$
\begin{align*}
E[F(A(x_i))] = \begin{cases}
    1 & i = j\\
    0 & i \neq j
\end{cases}
\end{align*}
$$

Let's see how this looks in practice. We'll generate some sample of inputs, put them through $k$-rr with $k=10$, then
debias them and check the error. We'll use $\epsilon = \log 3$ and $N = 100000$ samples. Here is an example for a 
single element:
```python
print(", ".join([f'{x:.2f}' for x in debias_krr(1, k=10, epsilon=math.log(3))]))
> -0.50, 5.50, -0.50, -0.50, -0.50, -0.50, -0.50, -0.50, -0.50, -0.50
```

and here's the whole lot of 'em:

![Debiased k-rr](/images/debiased-rr.svg)

As you can see, the raw noisy output is nearly unusable, but the debiased output closely follows the distribution true
distribution.

## Extending to RAPPOR

It is clear that binary randomized response can be debiased by taking the formulas above and substituting $k=2$, i.e.
$F(x) = \frac{x - p/2}{1-p}$.
This gives us an approach to debias an individual RAPPOR report, which is just element-wise randomized response across
a vector. Just keep in mind that for RAPPOR, $p$ is computed a little differently than in binary randomized response
(it has a factor of $e^{\epsilon  / 2}$ rather than $e^\epsilon$).

## Debiasing aggregates

Since the debiasing functions are linear, we can apply them to individual outputs or aggregate outputs.
Let $F_A$ be the aggregate debiasing function over $N$ items. Let's look at $F$ applied to each element individually:

$$
\begin{align*}
\sum_{x \in X} F(x) &= \sum_{x\in X} \frac{1}{1-p} \left(x - \frac{p}{k}\right)\\
&= \sum_{x \in X} \frac{x}{1-p} - \sum_{x \in X} \frac{p/k}{1-p}\\
&= \sum_{x \in X} \frac{x}{1-p} - N \frac{p/k}{1-p}
\end{align*}
$$

So we can debias one bucket of the aggregate histogram via:

$$
F_A(Y) = \frac{Y - N p / k}{1-p}
$$

## Conclusion

This post was fairly basic, but actually dealing with the output of differentially private mechanisms like randomized
response is not straightforward and requires a bit of care. It makes me question whether a system even ought to
output biased output if it can be post-processed by that system before any emission takes place. There are
downsides like increased communication cost of this approach, but it does eliminate a possible footgun of
someone blindly using the biased output.