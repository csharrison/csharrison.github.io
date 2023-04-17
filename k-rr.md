# K-ary randomized response
K-ary randomized response is a privacy mechanism with the following algorithm:
```python
def kary_rr(p: float, input: Any, possible_outputs: Set[Any]) -> Any:
  assert input in possible_outputs
  If random.random() <= p:
    return random.choice(possible_outputs)
  return input
```

This algorithm satisfies $\epsilon$-differential privacy. In this post I will present a proof. Recall the definition of differential privacy:
```math
Pr[A(D) = x] \le e^\epsilon Pr[A(D’) = x]
```
For $D$ and $D’$ datasets that differ on one row / user, a privacy mechanism $A$, and any output $x$. Because k-RR is a local differentially private mechanism (i.e. it operates on only a single row), we can assume $D$ and $D’$ differ only in that they are different elements of the `possible_outputs` set. Call that set $O$, and the differing inputs $o$ and $o’$.

It is obvious that if $o \ne x$ and $o’ \ne x$, the inequality holds, and in fact without the $e^\epsilon$ term it’s an equality.

It suffices to show that the inequality holds when either $o = x$ or $o’ = x$.

### $o = x$ and $o’ \ne x$
```math
\begin{align*}
Pr[A(o) = x] &= (1 - p) + \frac{p}{k} \\
Pr[A(o') = x] &= \frac{p}{k} \\
\end{align*}
```
So
```math
\begin{align*}
(1 - p) + \frac{p}{k} &\le e^{\epsilon} \frac{p}{k} \\ 
\frac{(k + p - k p)}{p} &\le e^{\epsilon} \\
\end{align*}
```
### $o \ne x$ and $o’ = x$
```math
\begin{align*}
Pr[A(o) = x] &= \frac{p}{k} \\
Pr[A(o’) = x] &= (1 - p) + \frac{p}{k} \\
\end{align*}
```
So
```math
\begin{align*}
\frac{p}{k} &\le e^{\epsilon} ((1 - p) + \frac{p}{k}) \\ 
\frac{p}{(k + p - k p)} &\le e^{\epsilon} \\
\end{align*}
```
### Putting it together

$\forall k > 0, p \in [0, 1]:  k + p - kp \ge p$, so
```math
\begin{align*}
e^\epsilon &\ge \frac{(k + p - kp)}{p} \\
\epsilon &\ge ln(\frac{(k + p - kp)}{p}) \text{ so}\\
\end{align*}
```

Therefore $p \le \frac{k}{k - 1 + e^\epsilon}$ satisfies $\epsilon$-differential privacy.

