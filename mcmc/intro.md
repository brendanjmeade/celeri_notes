
The geophysical forward model in celeri provides maps from geodynamic parameters
(block rotations, elastic slip, Mogi sources, strain rates, etc.) to predicted
ground velocities at GPS stations.

Mathematically, if $\theta$ denotes the vector of parameters, the predicted
velocity field at all stations is

$$
v_\text{total}(\theta) \in \mathbb{R}^{2N_\text{stations}},
$$

where each station contributes two observed components (east and north).

The state vector $\theta$ is assembled from multiple physical contributions.
Each process contributes an additive term to the total expected velocity
$v_\text{total}(\theta)$. The components are build from linear maps that encode
elastic properties, that don't depend on the parameters $\theta$ and are
precomputed by celeri.

This decomposes $v_\text{total}(\theta)$ into a sum

$$
v_\text{total}(\theta) = v_\text{strain} + v_\text{rot} + v_\text{rot,okada} + v_\text{mogi} + v_\text{elastic}.
$$

The PyMC model turns this deterministic forward map into a probabilistic model
for inference.

## Observations
At $N_\text{stations}$ GPS sites we observe noisy velocities

$$
y_i \in \mathbb{R}^2, \quad i=1,\dots,N_\text{stations}.
$$

We assume a Gaussian likelihood with shared variance $\sigma^2$:

$$
y_i \sim \mathcal{N}\left(v_\text{total, $i$}(\theta), \sigma^2 I_2\right).
$$

Here $v_\text{total, $i$}(\theta)$ is the forward model prediction for station
$i$, restricted to the horizontal $(x,y)$ directions.

## Geophysical Model Components

### Block strain rates
Each deformable block has a homogeneous 2D strain tensor with 3 parameters. We
use a standard normal prior on the strain parameters after rescaling the
operator columns. The contribution to station velocities is linear:

$$
v_\text{strain} = O_\text{strain} \theta_\text{strain}.
$$

### Block rotations
Each block rotates about the Earth’s center (3 parameters), again with a
standard normal prior after rescaling. The contribution has two parts:
- Station velocities from rigid rotation: $v_\text{rot} = O_\text{rot} \, \theta_\text{rot}$.
- Slip along faults induced by rotation, mapped through Okada dislocation
  operators $v_\text{rot,okada} = O_\text{rot,okada} \, \theta_\text{rot}$.

### Mogi sources
Point sources of pressure at fixed locations, with a scalar intensity parameter
each. Again, we use standard normal prior after rescaling. The contribution to
the velocity is again linear:

$$
v_\text{mogi} = O_\text{mogi} \theta_\text{mogi}.
$$

### Elastic slip on faults (TDEs)

On meshed fault segments, we want to model how much of the kinematically
expected slip rate (from rigid block rotations) is actually accommodated
elastically at the fault.

#### Coupling field as a Gaussian process

Define the *coupling ratio field*

$$
c(s) = \frac{v_\text{elastic,fault}(s)}{v_\text{kinematic,fault}(s)}, \quad s \in \text{fault surface},
$$

where $v_\text{kinematic,fault}(s)$ is the relative slip rate implied by block
rotations and $v_\text{elastic,fault}(s)$ is the actual elastic slip rate at the
fault. The kinematic velocities are computed as a linear function of the
rotation parameters $v_\text{kinematic,fault} = O_\text{kinematic,fault} \theta_\text{rot}$.

$c(s)$ is given a Gaussian process prior, encoding smoothness along the fault mesh:

$$
c(s) \sim \mathcal{GP}(0, K(s, s')),
$$

with a squared-exponential kernel and constant variance and length scale
hyperparameters. Those hyperparameters are currently fixed, and configured the
same way as in the SQP solver.

We use separate GPs for each fault. At each fault, we compute one GP for the
strike-slip component and one for the dip-slip component, and assume
independence between them.

#### Eigen-expansion (Karhunen–Loève truncation)

In practice, the GP covariance matrix $K \in \mathbb{R}^{M \times M}$ is
computed on the mesh nodes ($M$ = number of fault DOFs). But to avoid very
high-dimensional inference, we use a low-rank approximation. We perform an
eigen-decomposition

$$
K = U \Lambda U^\top, \quad \Lambda = \text{diag}(\lambda_1, \dots, \lambda_M),
$$

and truncate after $n$ leading eigenmodes (largest $\lambda_i$).
Thus,

$$
c(s) \approx \sum_{i=1}^n \alpha_i u_i(s),
$$

with coefficients $\alpha_i \sim \mathcal{N}(0, \lambda_i)$.

This reduces the high-dimensional GP to a low-dimensional expansion that still
captures smooth variation.

#### Bounded transform for constraints

Physically, coupling ratios should fall within plausible limits, e.g.

$$
0 \leq c(s) \leq 1,
$$

or some more general $(\ell, u)$ bounds.

We use smooth transforms of $c(s)$ to enforce these bounds. Logistic (sigmoid)
if both lower and upper bounds are finite, softplus if only a single bound is
given. This ensures the posterior remains differentiable for HMC. So for instance

$$
c_\text{bounded}(s) = \text{expit}(c(s)) (u - \ell) + \ell.
$$

#### Resulting elastic slip and station velocities

The elastic slip rate field on the fault surface is then

$$
v_\text{elastic,fault}(s) = c(s) v_\text{kinematic,fault}(s).
$$

This field is then mapped through a linear operator to yield predicted
velocities at GPS stations:

$$
v_\text{elastic} = \sum_{\text{faults}}\sum_{\text{dip slip, strike slip}} O_\text{tde} v_\text{elastic,fault}.
$$

#### Alternative direct elastic formulation

Depending on the model configuration, we also support modeling the elastic slip
field as a GP, instead of a coupling ratio. In this case, the elastic slip field
itself** is given a GP prior and eigen-expansion:

$$
v_\text{elastic,fault}(s) \approx \sum_{i=1}^n \beta_i \, u_i(s),
\quad \beta_i \sim \mathcal{N}(0, \lambda_i).
$$

Bounds (if any) are again imposed via softplus or sigmoid transforms, and the
slip field is mapped linearly to station velocities.

Which formulation is used is configured in the model setup, and can differ by
fault or dip/strike component if desired.

## Likelihood and Posterior
The complete likelihood is

$$
\begin{gathered}
\sigma \sim \text{HalfNormal}(\text{scale}=2) \\
y \sim \mathcal{N}\left(v_\text{strain} + v_\text{rot} + v_\text{rot,okada} + v_\text{mogi} + v_\text{elastic}, \sigma^2 I \right).
\end{gathered}
$$

Since all inequality constraints (e.g. on slip rates or coupling ratios) are
enforced by bounding transforms rather than hard truncation, the posterior
density is smooth and differentiable, so we can use gradient-based MCMC.

The posterior is sampled with **HMC/NUTS** (via PyMC + nutpie). Since there are
highly correlated parameters (e.g. slip on adjacent fault segments), we use the
diagonal plus low-rank mass matrix adaptation scheme from nutpie.


## Output
The posterior predictive mean $\hatv_\text{total}(\theta)$ is projected back into
celerí’s state vector format for downstream analysis.
Individual draws can be accessed as `estimation.mcmc_draw(chain_idx, draw_idx)`.
