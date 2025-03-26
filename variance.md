Consider a gene whose expression fraction follows some distribution with mean $q$ and variance $v$ (non log scale), and counts $K_i$ for this gene for cell $i$, with total count $s_i$:

$$\begin{align}
 K_i|Q_i &\sim \text{Pois}(s_i Q_i) \\
 \operatorname{E}(Q_i) &= q \\
 \operatorname{Var}(Q_i) &= v
\end{align}$$ The total expectation and variance of $K_i$ is
$$\begin{align}
\text{E}(K_i) &= \text{E}(\text{E}(K_i|Q_i))=E(s_i Q_i)=s_iq \\
\text{Var}(K_i) &= \text{E}\left(\text{Var}(K_i|Q_i)\right) + \text{Var}\left(\text{E}(K_i|Q_i)\right) = \text{E}(s_iQ_i) + \text{Var}(s_iQ_i) =
s_iq+s_i^2v
\end{align}$$
The obvious estimator for $q$ is the mean of the $K_i/s_i$:
$$ \hat q = \frac{1}{n}\sum_i\frac{K_i}{s_i}. $$
We check that it has expectation $q$:
$$ \text{E}(\hat q) = \frac{1}{n}\sum_i\frac{\text{E}(K_i)}{s_i} = \frac{1}{n}\sum_i \frac{s_iq}{s_i} = q $$
For later use, we find 
$$ \text{E}(K_i^2)=[\text{E}(K_i)]^2+\text{Var}(K_i) = 
s_i^2 q^2  + s_iq + s_i^2v$$
For $i \neq j$, we have $$E(K_i K_j)=\text{E}(K_i)\text{E}(K_j)=s_is_jq^2,$$because they are independent.

To estimate the variance, we study the following estimator, which is the formula for the sample variance of $K_i/s_i$:
$$\hat w := \frac{1}{n-1}\sum_i\left(\frac{K_i}{s_i}-\hat q\right)^2 $$

Its expectation is (see after the calculation for explanations):

$$\begin{align}
\text{E}(\hat w) 

&= \text{E}\left(\frac{1}{n-1}\sum_i\left(\frac{K_i}{s_i}-\frac{1}
{n}\sum_j\frac{K_j}{s_j}\right)^2\right) \\


&=\frac{1}{n-1}\sum_i\text{E}\left(\left(\frac{1}{n}\sum_j\left(\frac{K_i}{s_i}-\frac{K_j}{s_j}\right)\right)^2\right)\\

&=\frac{1}{n-1}\sum_i\text{E}\left(\frac{1}{n^2}\sum_j\sum_l\left(\frac{K_i}{s_i}-\frac{K_j}{s_j}\right)\left(\frac{K_i}{s_i}-\frac{K_l}{s_j}\right)\right)\\

&=\frac{1}{n^2(n-1)}\sum_i\sum_j\sum_l\left( \frac{\text{E}(K_i^2)}{s_i^2} -
\frac{\text{E}(K_iK_l)}{s_is_l} - \frac{\text{E}(K_iK_j)}{s_is_j} + \frac{\text{E}(K_jK_l)}{s_js_l}
\right)\\

&\stackrel{(a)}{=}\frac{1}{n^2(n-1)}\left(
n^2\sum_i\frac{\text{E}(K_i^2)}{s_i^2}-n\sum_i\sum_j\frac{\text{E}(K_iK_j)}{s_is_j}
\right)\\

&\stackrel{(b)}{=}\frac{1}{n^2(n-1)}\left(
(n^2-n)\sum_i\frac{\text{E}(K_i^2)}{s_i^2}-n\sum_i\sum_{\substack{j\neq i}}\frac{\text{E}(K_iK_j)}{s_is_j}
\right)\\

&\stackrel{}{=}\frac{1}{n^2(n-1)}\left(
(n^2-n)\sum_i\frac{s_i^2 q^2  + s_iq + s_i^2v}{s_i^2}-n\sum_i\sum_{\substack{j\neq i}}\frac{s_is_jq^2}{s_is_j}
\right)\\

&\stackrel{}{=}\frac{1}{n^2(n-1)}\left(
(n^2-n)\sum_i\left( q^2 + \frac{q}{s_i} + v\right) - n\sum_i\sum_{\substack{j\neq i}}q^2
\right)\\

&\stackrel{}{=}\frac{1}{n^2(n-1)}\left(
n(n-1)(nq^2+q\sum_i\frac{1}{s_i}+nv) - n^2(n-1)q^2
\right)\\

&\stackrel{}{=}
q^2+q\frac{1}{n}\sum_i\frac{1}{s_i}+v - q^2
\\

&=q\Xi+v

\end{align} $$

In (a), we made use of the fact that the three mixed terms are all the same, with permuted indices, and hence two of them cancel. In (b), we moved the diagonal part ($i=j$) from the second term to the first. In the last line, we have introduced $\Xi$, the reciprocal harmonic mean of the $s_i$:

$$ \Xi:=\frac{1}{n}\sum_i\frac{1}{s_i}$$

We see that the sample variance of the $K_i/s_i$ is the sum of two components: $q\Xi$ captures the expected Poisson part of the variance and $v$ is the expected extra-Poisson part.
