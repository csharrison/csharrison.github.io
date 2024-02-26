# K-ary randomized response
K-ary randomized response is a privacy mechanism with the following algorithm:
```python
import math
import random
def kary_rr(epsilon: float, input: Any, possible_outputs: Set[Any]) -> Any:
  p = k / (k - 1 + math.exp(epsilon))
  assert input in possible_outputs
  If random.random() <= p:
    return random.choice(possible_outputs)
  return input
```

This algorithm satisfies $\epsilon$-differential privacy. In this post I will present a simple proof. Recall the definition of differential privacy:
```math
Pr[A(D) = x] \le e^\epsilon Pr[A(D’) = x]
```
For $D$ and $D’$ datasets that differ on one row / user, a privacy mechanism $A$, and any output $x$. Because k-RR is a local differentially private mechanism (i.e. it operates on only a single row), we can assume $D$ and $D’$ differ only in that they are different elements of the `possible_outputs` set. Call that set $O$, and the differing inputs $o$ and $o’$.

It is obvious that if $o \ne x$ and $o’ \ne x$, the inequality holds. It suffices to show then that the inequality holds when either $o = x$ or $o’ = x$.

```math
\begin{align*}
\frac{Pr[A(o) = x]}{Pr[A(o') = x]} &\le \frac{Pr[A(x) = x]}{Pr[A(o') = x]} &&\text{ for $o' \ne x$}\\
&= \frac{(1 - p) + \frac{p}{k}}{\frac{p}{k}} \\
&= \frac{k}{p} - k + 1 \\
&= k \frac{k - 1 + e^\epsilon}{k} - k + 1 \\
&= e^\epsilon
\end{align*}
```

The first inequality follows from the simple observation that the mechanism will always output the true input with higher probability than any other output.
