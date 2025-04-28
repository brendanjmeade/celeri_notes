#import "@preview/cetz:0.2.0"

= Technical Note: Formulations for Constrained Least Squares Optimization

This note explores different mathematical formulations for solving the
constrained least squares optimization problem:

$min_x |C x - b|_2 quad "subject to" quad A x = d, quad G x >= f$

where $C in bb(R)^(m times n)$, $b in bb(R)^m$, $A in bb(R)^(p times n)$, $d in
bb(R)^p$, $G in bb(R)^(q times n)$, $f in bb(R)^q$, and $m >> n$.

== 1. Direct Quadratic Programming Approach (direct)

The 2-norm can be squared to obtain a standard quadratic programming formulation:

$min_x |C x - b|_2 = min_x |C x - b|_2^2 = min_x (C x - b)^T (C x - b)$

Expanding gives
$(C x - b)^T (C x - b) = x^T C^T C x - 2b^T C x + b^T b.$

The constant term does not affect the minimization, so we can drop it. We get

$min_x x^T C^T C x - 2b^T C x quad "subject to" quad A x = d, quad G x >= f.$

This can be written in standard form as:

$min_x 1/2x^T P x + q^T x quad "subject to" quad A x = d, quad G x >= f$

where $P = 2C^T C$ and $q = -2C^T b$.

*Note on numerical issues*: While dropping the constant term $b^T b$ is
mathematically valid, it can lead to numerical issues when $|b|$ is large. The
objective function value is offset by a large constant, which may cause
convergence criteria to behave poorly.

*KKT System Dimensionality*: The KKT system for this approach involves:
- $n$ variables from the original problem
- $p$ dual variables for equality constraints
- Up to $q$ dual variables for inequality constraints.

The resulting KKT matrix has dimensions $(n+p+q) times (n+p+q)$.

== 2. Variable Splitting Method (splitting)

We introduce auxiliary variables $r in bb(R)^m$ where $r = C x - b$, and get

$min_(x,r) |r|_2^2 quad "subject to" quad r = C x - b, quad A x = d, quad G x >= f.$

This formulation uses additional variables to avoid the numerical issues of the
previous approach.

*KKT System Dimensionality*: The KKT system for this approach involves:
- $n + m$ primal variables ($n$ from $x$ and $m$ from $r$)
- $m$ dual variables for the constraint $r = C x - b$
- $p$ dual variables for equality constraints $A x = d$
- Up to $q$ dual variables for inequality constraints $G x >= f$

The resulting KKT matrix has dimensions $(n+2m+p+q) times (n+2m+p+q)$.

== 3. Variable splitting with Dimension Reduction via QR (reduced-splitting-qr)

In this approach, we use an economic QR decomposition of $C$ to reduce the
dimensionality of the problem:

$C = Q R$

where $Q in bb(R)^(m times n)$ has orthonormal columns and $R in bb(R)^(n times
n)$ is upper triangular.

We introduce new variables $hat(r) in bb(R)^n$ defined as:

$hat(r) = R x - Q^T b$

Our problem becomes:

$min_(x,hat(r)) |hat(r)|_2 quad "subject to" quad hat(r) = R x - Q^T b, quad A x
= d, quad G x >= f$.

=== Proof of Equivalence

We need to show that minimizing $|hat(r)|_2$ with the constraint $hat(r) = R x -
Q^T b$ leads to the same minimizer as the original problem $min_x |C x - b|_2$.

Starting from the original objective $|C x - b|_2 = |Q R x - b|_2$.

Since $Q$ has orthonormal columns, we can decompose $b$ into a component in the
column space of $Q$ and a component orthogonal to it:

$b = Q Q^T b + (I - Q Q^T)b$

This gives

$|Q R x - b|_2 = |Q R x - Q Q^T b - (b - Q Q^T b)|_2$.

Since $Q R x - Q Q^T b$ is in the column space of $Q$ and $b - Q Q^T b$ is orthogonal to it we get
$|Q R x - b|_2^2 = |Q R x - Q Q^T b|_2^2 + |b - Q Q^T b|_2^2.$

The second term $|b - Q Q^T b|_2^2$ is independent of $x$, so it doesn't affect
the minimizer. But the first term is $|Q R x - Q Q^T b|_2^2 = |Q hat(r)|_2^2 =
hat(r)^T Q^T Q hat(r) = |hat(r)|_2, $ because the columns of $Q$ are
orthonormal.

*KKT System Dimensionality*: The KKT system for this approach involves:
- $n + n$ primal variables ($n$ from $x$ and $n$ from $hat(r)$)
- $n$ dual variables for the constraint $hat(r) = R x - Q^T b$
- $p$ dual variables for equality constraints $A x = d$
- Up to $q$ dual variables for inequality constraints $G x >= f$

The resulting KKT matrix has dimensions $(3n+p+q) times (3n+p+q)$.

== 4. Second-Order Cone Programming (socp)

We can reformulate the original problem as a second-order cone program using the epigraph form:

$min_(x,t) t quad "subject to" quad |C x - b|_2 <= t, quad A x = d, quad G x >= f$.

This is equivalent to:

$min_(x,t) t quad "subject to" quad (C x - b, t) in cal(Q)^(m+1), quad A x = d,
quad G x >= f$

where $cal(Q)^(m+1)$ is the second-order cone in $bb(R)^{m+1}$.

This formulation allows the use of SOCP solvers.

*KKT System Dimensionality*: The KKT system for the SOCP formulation involves:
- $n+1$ primal variables ($n$ from $x$ and $1$ from $t$)
- $m+1$ dual variables for the second-order cone constraint
- $p$ dual variables for equality constraints
- Up to $q$ dual variables for inequality constraints

The resulting KKT matrix has dimensions $(n+m+p+q+2) times (n+m+p+q+2)$.
