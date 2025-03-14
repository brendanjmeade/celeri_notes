# Notes on Celeri

## Adrian, 2025-01-22

### **Block**
- A block represents a tectonic plate or a part of a tectonic plate.
- Blocks can:
  1. Rotate globally around the center of the earth.
  2. Deform globally as an elastic body.
  3. Deform due to forces at their boundaries. (the tde things, I don't really get those yet though...)
  4. Deformation due to point sources (mogi). Locations of these sources
     are configured in a csv file.

Data about blocks is stored in CSV files, but their exact outlines seem to be
inferred from their boundaries (the segments).

The blocks themselves do not have discretizations, they are assumed to be
homogeneous.

I am still a little confused about how I can think about those blocks. They are
3D objects, but have only partial boundaries? For instance, they don't have a
boundary below? And only some segments have associated meshes, so how does the
interaction between blocks work between blocks that don't intersect at a mesh? I
think some of this can probably be explained because we don't really look at how
the blocks change over time, but only the derivative, so only infinitesimal
changes? That seems to make it possible to keep some properties of the blocks
more vague?

---

### **Faults and Segments**

- Segments represent the boundaries between blocks.
- Each segment is a single line between two points, and each segment has
  properties defined in a csv file.

Some subsets of segments have associated 2D meshes of boundaries between blocks.
Other segments are only defined on the surface? Not sure what that implies for
the structure of the blocks...

Is it correct to call each of those meshes a fault?

---

### **Stations**

- Stations are the locations where GPS data measures ground motion.
- For a given run, each station has a single 3D velocity vector that the model
  tries to fit.

---

## **Model Parameters**

The model tries to find the following parameters:

### Block Rotations

- For each block, there is a global rotation.
- The rotation is defined around the Earth's center with three angles.
- In the Japan example, this corresponds to `3 * n_blocks` parameters.
- The velocities are computed using `operators["rotation_to_velocities"] @
  rotation_parameters`, so the velocities are a linear function of the parameters.

There are 954 stations, so 3 * 954 = 2862 velocity values.
`operators["rotation_to_velocities"]` has shape (2862, 63), so we can think of
this as a reshaped (954, 3, 21, 3), or `("station", "direction", "block",
"rotation_parameter")` tensor.

There are some linear constraints about the rotation as well:

`block_motion_constraints` I guess can give bounds on how much the block can
move? In the Japan example, these are empty. What would typically go here?

`slip_rate_constraints`: I guess they control how much relative movement of the
blocks there can be at the faults? Empty in the Japan example.

`rotation_to_slip_rate`: This matrix computes the slip rate (relative velocity)
at each fault segment? It contains pretty gigantic values. What units are these
numbers in?

### Mogi

Mogi refers to point-like pressure sources below the block. The locations of
those point sources are defined in a csv file. In the japan example there are in
total 3 of those point sources. We have one parameter per point source (what is
that parameter?), so in total 3 mogi parameters.

The station velocities are given by `operators["mogi_to_velocities"] @ mogi_parameters`.

### Block strain rate

We assume that every block is being homogeneously deformed (as a 2D object?).
In 2D, the strain rate is a symmetric 2 by 2 tensor, so 3 degrees of freedom
per block. Again, this gives us a linear mapping from strain rate parameters
in each block to station velocities `operators["block_strean_rate_to_velocities"]`.

Only blocks with the `strain_rate_flag` set can be deformed.

### TDE slip

This is I think the hardest to understand, I'm still struggling with it a bit.

We have parameters that describe what happens on the boundaries between the
blocks (only the meshes or also the non-meshed segments?), and infer from that,
how the interior of the block will behave.

From what I get (or guess...) so far:
- We build basis functions on the meshes, something like taking the eigenvalues
  of the kernel of a GP or so? I'm not sure what those functions are exactly.
  The slip rate on the boundaries?
- We then say that the actual function on the boundary is a linear span of those
  basis functions, and use the coefficients as free parameters.
- There is then a linear mapping from those coefficients to the velocities at
  the station locations.

I think the operators that results from this is hidden a bit in `get_full_dense_operator_eigen`.

There are also a *lot* of constraints about these functions. I don't know yet
what those constraints actually do. Most are pretty harmless simple linear
constraints, but there is also a non-convex constraint...

`cutde` seems to be a related software package.

Nice reference for TDEs: https://doi.org/10.1093/gji/ggv035
I haven't read it yet though.


## Fitting

We can combine the linear operators to get a linear function from all our
parameters to the expected velocity at each station, and a couple of linear
constraints on our parameters.

In addition, we also have a non-convex constraint, that specifies that the ratio
of two velocities at the faults has to be between 0 and 1.
(this involves some smoothing function as well...)

`celeri` combines the individual operators into a larger (dense) matrix $C$. In
the japan example this has shape `(1948, 299)` with dimensions ("station",
"free_parameter").

Note: Since we only need $P = C^TC$ in the optimization, I don't think we really
have to build this full matrix? We could just build $P$ directly from the individual
operators.

Note: The resulting matrix $P$ seems to have some scaling issue. If I throw them
into `cvxpy` unchanged, it isn't too happy. Maybe it would be better to change
some units to prevent precision loss? We can also reparameterize our free
parameters using something like `scale = np.sqrt((C_val ** 2).sum(0))`. This
still leaves us with some pretty small eigenvalues in $P$ though... When we
rescale like this clarabel seems fine optimizing the linear problem, ignoring
the non-convex constraints and takes with sparse matrices about as long as
cvxopt with (I think) dense matrices:

```python
import cvxpy as cp

params = cp.Variable(name="params", shape=operators.eigen.shape[1])

C_val = operators.eigen * np.sqrt(weighting_vector_eigen[:, None])
C = C_val

d_val = data_vector_eigen * np.sqrt(weighting_vector_eigen)
d = d_val

A_val = qp_inequality_constraints_matrix
A = A_val

b_val = qp_inequality_constraints_data_vector
b = b_val

scale = np.sqrt((C_val ** 2).sum(0))

params_scaled = params / scale

C_hat = C / scale
P = C_hat.T @ C_hat

# Maybe reparameterize based on the eigenvalue decomposition of P?
# Something like (but not exactly) this?
#vals, vecs = linalg.eigh(P)
#U = vecs * np.sqrt(vals)
#P_hat = np.eye(len(vals))
#params_trafo = (vecs.T / np.sqrt(vals)) @ params_scaled

b_scale = np.abs(b)

objective = cp.Minimize(cp.quad_form(params, 0.5 * P, True) - (d.T @ C_hat) @ params)
constraint = (sparse.csr_array((A / b_scale[:, None]) / scale)) @ params <= b / b_scale

kinematic = operators.rotation_to_tri_slip_rate[0] @ params[0 : 3 * len(block)]
estimated = operators.eigenvectors_to_tde_slip[0] @ params[
    index.start_col_eigen[0]:index.end_col_eigen[0]
]

problem = cp.Problem(objective, [constraint])

problem.solve(verbose=True, solver="CLARABEL")
```

This takes the solver about 2s (roughly the same as the cvxopt solver in the
code)

A second matrix $A$ describes the linear constraints of the problem. It has
dimensions ("constraint", "free_parameter") with shape (19320, 299), and about
85% zeros.

### Nonconvex constraints

There are also some non-convex constraints though, but I don't really understand
what they represent.

For some reason it is required that the ratio of the slip rate between blocks
due to rotation $v_k$ of the whole blocks and the estimated slip $v_e$ is
between 0 and 1. Both those slip rates are linear functions of the parameters.
I guess this tries to make sure the movement at the boundary is not more or
of the opposite sign as we would expect from the overall movement of the blocks?

The code solves the optimization problem once without those constraints, and
then iteratively adds linearized constraints and resolves.

Note: This looks somewhat similar to https://en.wikipedia.org/wiki/Sequential_quadratic_programming

Note: We can reformulate $0 < v_e/v_k < 1$ as $v_e * v_e - v_e * v_k > 0$. Maybe
there are other ways to rewrite this? The space of valid $v_e, v_k$ is a double
cone, so not convex. If we knew the sign of one of $v_e$ or $v_k$, it would be
convex though... It also reminds me of McCormick envelopes. Maybe they could be
used? (https://optimization.cbe.cornell.edu/index.php?title=McCormick_envelopes)







# Comments

Brendan:

Wow, @aseyboldt this is a very good summary and analysis!  I'm glad the existing code was tolerably interpretable.

The method I made up is definitely very close to sequential quadratic programming and I termed it such for a while!  The only reason I don't call it that is because in the optimization community that phrase seems to be reserved for a method with some sort of proper gradient in the non-linear mapping and I just have some sort of ad hoc reduction in the count of some parameters rather than a formal gradient.

Let me add some notes:

1. The matrix does indeed have scaling issues.  Many times I column normalize and that really improves the condition number.  I recommend it.  Something like:
```python
def normalize_columns(A):
    # Compute column norms
    col_norms = np.sqrt(np.sum(A**2, axis=0))

    # Avoid division by zero
    col_norms[col_norms == 0] = 1.0

    # Normalize columns
    A_normalized = A / col_norms[None, :]

    return A_normalized, col_norms
```
Then do the solve and project the estimated vector, `x_est`, back into the unormalized space,

```python
x_est = x_est / col_norms
```

2. It's great that the solve time is comparable to the QP solver.  Are the estimated solutions similar?  Hopefully!

3. As you note, adding constraints is the hard part of this problem.  There are generally two types of constraints that we really want.
- Linear inequality constraints $\mathrm{C} \mathrm{x} \leq \mathrm{d}$.  We use these in the notebook I shared, both directly and as part of the iterative algorithm.
- Non-linear non-convex constraints where we want $a\leq f({\mathbf{x}}) \leq b$.  In the current version I do this via an iterative approach that does a non-linear mapping between the constraint that we want and the linear inequality constraints that be applied to quadratic programming problems.

4.  The non-linear non-convex constraints.  A challenge to cleverness here is that we never know the signs of $v_\mathrm{e}$ or $v_\mathrm{k}$ a priori.  This is challenging.  If we know the signs, we could actually write this as a set of two linear inequality constraints (one each for the upper and lower bounds)!  The reason we can't do that is illustrated by considering this upper bound example,

$$
\frac{v_\mathrm{e}}{v_\mathrm{k}} \leq b
$$

The obvious thing to do is to multiply both sides by $v_\mathrm{k}$ and then rearrange.  However, not knowing the sign of $v_\mathrm{k}$ is a problem because if it were negative, and we multiplied both sides by that it would flip the sign of the inequality constraint!  This is why the reformulation you suggested in your last note isn't a way forward: Because not knowing the sign of $v_\mathrm{k}$ means that we can't know whether the inequality should be $\leq$ or $\geq$.  Does this make sense?  I'd love to be able to write the problem as a constraint like this but the sign ambiguity seems to prevent it.

To me, this puts us in a great place to ask the big question of whether you think it's worth trying and MCMC-type approach to solving the problem with constraints, including the non-linear non-convex ones.

Thoughts?

I'd certainly be willing to invest in this if you think there's some potentially positive prospect worth exploring here.  I respect that there may not be, but I'm hoping that there might be!  I'll follow up via email



Adrian:

I replied to the mail, but about the constraint:
Totally possible that I'm missing something, but I think $x^2 - xy \leq 0$ iff ($\frac{x}{y} \geq 0$ and $\frac{x}{y} \leq 1$), no matter what signs of $x$ and $y$ are...

```python
x, y = np.mgrid[-3:3:500j, -3:3:500j]
plt.pcolor(x, y, (x / y > 0) & (x / y < 1), vmin=0, vmax=1)
plt.pcolor(x, y, (x * x - x * y <= 0), vmin=0, vmax=1)
```

Maybe it is work checking what gurobi is doing to solve those?
https://docs.gurobi.com/projects/optimizer/en/current/concepts/modeling/constraints.html#quadratic-constraints



Brendan:


It took me a while but I understand this now!  The inequality constraints:

$$
a < \frac{x}/{y} < b
$$

is equivalent to the quadratic form:

$$
(x - ay) (x-by) < 0
$$

Following your code snippet, here are a couple of examples.

<img width="833" alt="Screenshot 2025-01-28 at 10 26 10 PM" src="https://gist.github.com/user-attachments/assets/cc06f670-8cbe-47cc-984f-70a98dfbfd2f" />

This is super clever, and I never would have recognized this!

Now, let's consider my specific problem.  I have:

$$
a < v_\mathrm{e} / v_\mathrm{k} < b
$$

In general, $v_\mathrm{e}$ is a linear function of the estimated state vector, $\mathbf{Ex}$.  In theory $v_\mathrm{k}$ is a linear function of the state vector as well, $\mathbf{Kx}$.  However, for pragmatic reasons, I'm currently writing this as a non-linear relationship $v_\mathrm{k} = k(\mathrm{x})$.  There are two reasons for this:

1. $v_\mathrm{e}$ is implicitly smooth by assumption but $\mathbf{Kx}$ so that the "raw" $v_\mathrm{k}$ Can have very short wavelength variations that would be impossible to ever match with a handful of smooth functions.  To mitigate this, I smooth $v_\mathrm{k}$.  This part is currently done in a non-linear way, but there are linear smoothing operators that could be used instead!
2.  $v_\mathrm{k}$ changes sign over meshes, and this means that there are zero crossings.  When there are zero crossings that $v_\mathrm{e}/v_\mathrm{k}$ term can blow up and make it so that $v_\mathrm{e}$ must effectively be very close to zero to satisfy the inequality constraint.  So in the current algorithm, I limit how small the magnitude of $v_\mathrm{k}$ can get as a type of regularization.  This step is non-liner and proved to help the convergence of the current approach.  However, it may not be necessary with your newly proposed idea because there's no division by a small number, at least in the apparent formulation of the problem!

So, let's assume for a moment that we can write the classic inequality constraint as:

$$
a < \frac{\mathbf{Ex}}{\mathbf{Kx}} < b
$$

The quadratic form would then be:

$$
(\mathbf{Ex} - a \mathbf{Kx}) (\mathbf{Ex} - b \mathbf{Kx}) < 0
$$

After some algebra I think (?) this would reduce to something like:

$$
\mathbf{x}^\mathrm{T} \left[ \mathbf{E}^\mathrm{T} \mathbf{E} - (a+b) \mathbf{E}^\mathrm{T} \mathbf{K} + ab \mathbf{K}^\mathrm{T}\mathbf{Kx} \right] \mathbf{x}  < 0
$$

Which can be rewritten as,

$$
\mathbf{x}^\mathrm{T} \mathbf{Q} \mathbf{x}  < 0
$$

And that's one of the forms that seems compatible with the [Gurobi link](https://docs.gurobi.com/projects/optimizer/en/current/concepts/modeling/constraints.html#quadratic-constraints) that you shared!

This sounds amazing and is probably worth trying!  I have three caveats (because I'm a worrier):

1. Gurobi is closed-source and commercial.  It's free for some academics (Harvard has a license).  I'm happy to consider and get this working with Gurobi if there's also an open source path.  And I think there might be calling [SCS](https://github.com/cvxgrp/scs) from [CVXPY](https://github.com/cvxpy/cvxpy), which is similar to CVXOPT, which I'm currently using. Probably worth trying!
2. It's not clear whether the regularized form of the problem where $v_\mathrm{k}$ is not a linear function of the state vector, $\mathbf{x}$ immediately compatible with this solution method.  Probably worth trying!
3. This would be awesome if it works and is fast (who knows).  I'd still be constructing large constraint matrices that consume a lot more memory than the core linear operator!  This is one of the reasons that I'm interested in finding out whether there's a sampling approach to solving the problem with this constraint where I don't have to build and store these large constraint matrices.  Who knows!

In summary, I think your insight is exceptionally clever and may make a big difference!  Let's figure out how to move forward!


Adrian:

That's a nice summary, thanks.
I saw the smoothing in the code, and decided to just ignore it for the moment, I hope that doesn't get in the way later.
I'd also be quite disappointed if we only call gurobi in the end, after working on open source software for so long, that would feel like a defeat. :-)
Sounds like in general this quadratically constrained problem is NP hard unfortunately. But this here looks like it might be useful: https://web.stanford.edu/~boyd/papers/pdf/qcqp.pdf
The software looks pretty much abandoned however.

I've also been thinking about another approach:
Right now, the parameters are the velocities on the mesh (right?). We then have a nice linear function from those velocities to the velocities at the stations. But the mapping to the ratio v_e / v_k is nonlinear, so we have nonlinear constraints.

Could we switch that up a bit? Could we use a scalar field v_e / v_k (or two, for the components of the velocity) on the mesh as parameter, and then compute the velocities based on that? This would give us nicer constraints, at the cost of a more complicated (non-linear) function from parameters to station velocities.
That might not be a good idea for optimization, but maybe it allows us to run MCMC much more easily. Disconnected parameter spaces and MCMC really don't mix well, but nonlinearity in the logp function really isn't a problem. We might then run into multi modality though...


Brendan:


> That's a nice summary, thanks.

You're the one with the good ideas!

> I saw the smoothing in the code, and decided to just ignore it for the moment, I hope that doesn't get in the way later.

I don't think it will be an issue because we can write smoothing as a linear operator.  With the help of my good friend ChatGPT it was easy to find a linear operator to do smoothing,

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

# Generate irregularly spaced 2D points
num_points = 500
x = np.random.uniform(-5, 5, num_points)
y = np.random.uniform(-5, 5, num_points)
points = np.vstack((x, y)).T

# Generate noisy function values
z_noisy = np.sin(x) * np.cos(y) + 0.2 * np.random.randn(num_points)

# Compute pairwise Euclidean distance matrix
D = distance_matrix(points, points)

# Define Gaussian weight function
sigma = 1.0  # Controls the spread of influence
W = np.exp(-(D**2) / (2 * sigma**2))  # Gaussian kernel

# Normalize rows so each row sums to 1 (diffusion-like smoothing)
W /= W.sum(axis=1, keepdims=True)

# Direct smoothing operation (L @ z_noisy)
z_smooth = W @ z_noisy  # No matrix inversion needed!

# Plot results
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc1 = axes[0].scatter(x, y, c=z_noisy, cmap="viridis", edgecolor="k", s=30)
axes[0].set_title("Noisy Data")
fig.colorbar(sc1, ax=axes[0])

sc2 = axes[1].scatter(x, y, c=z_smooth, cmap="viridis", edgecolor="k", s=30)
axes[1].set_title("Smoothed Data (Dense Operator)")
fig.colorbar(sc2, ax=axes[1])

plt.show()
```
<img width="995" alt="Screenshot 2025-01-31 at 1 51 35 PM" src="https://gist.github.com/user-attachments/assets/32055321-c1ee-45b6-96df-1b79663c5a0e" />


> I'd also be quite disappointed if we only call gurobi in the end, after working on open source software for so long, that would feel like a defeat. :-)

Strong agree.  However, it might be an ok place to start because it has clear support.


> Sounds like in general this quadratically constrained problem is NP hard unfortunately. But this here looks like it might be useful: https://web.stanford.edu/~boyd/papers/pdf/qcqp.pdf

I've skimmed this paper now and if the constraints matrix is PSD then it's polynomial time so there's reason for hope!

> The software looks pretty much abandoned however.

Drat.

>I've also been thinking about another approach:
Right now, the parameters are the velocities on the mesh (right?). We then have a nice linear function from those velocities to the velocities at the stations. But the mapping to the ratio v_e / v_k is nonlinear, so we have nonlinear constraints.
>
> Could we switch that up a bit? Could we use a scalar field v_e / v_k (or two, for the components of the velocity) on the mesh as parameter, and then compute the velocities based on that? This would give us nicer constraints, at the cost of a more complicated (non-linear) function from parameters to station velocities.

   I'll think about this.  My initial guess is that this is unlikely (maybe?) because both $v_e$ and $v_k$ both have strong spatial self-correlations.  Niether field is allowed to vary freely in space rather $v_e$ results from summing a set of eigenmodes and $v_k$ is a super complicated, albeit linear, function of the geometry of the the entire mesh.  Again I'm not sure and will think more about this.
