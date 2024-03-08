Basic RAPPOR ([Erlingsson et al 2014](https://static.googleusercontent.com/media/research.google.com/en//pubs/archive/42852.pdf))
is a privacy mechanism with the following algorithm:

```python
import math
import random
def rappor(value: int, domain_size: int, epsilon: float):
  """Execute the RAPPOR mechanism on value

  Arguments:
  value -- the value to randomize. Must be an integer between 0 and domain_size - 1
  domain_size -- a positive integer > 1
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

This algorithm satisfies $\epsilon$-differential privacy. Here is a proof. Recall the definition of differential privacy:
```math
Pr[A(D) = x] \le e^\epsilon Pr[A(D’) = x]
```
For $D$ and $D’$ neighoring datasets (that differ on one row / user), a privacy mechanism $A$, and any output $x$.

Because RAPPOR is a local differentially private mechanism (i.e. it operates on only a single row), we can assume $D$ and $D’$ differ
only in `value`. For any value of $x$ we choose, $A(D)$ will need to set two bits differently from how $A(D')$ sets them. In other words,
where $A(D)$ leaves the bit unchanged, $A(D')$ must flip the bit. Where $A(D)$ flips teh bit, $A(D')$ must leave the bit unchanged.

As an example: take $D = 0$, $D' = 1$, and $x = [1, 0, 0, ...]$. It is clear that for $A(D') = x$ the first two bits need to be flipped, unlike in $A(D)$
where they are left unchanged.

Therefore the ratio of $\frac{Pr[A(D) = x]}{Pr[A(D') = x]}$ is maximized when the $D'$ flips those two bits with probability $f/2$, since $f/2 \le 1 - f/2$.

```math
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
```
