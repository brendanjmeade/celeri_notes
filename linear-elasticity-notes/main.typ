#import "@preview/cetz:0.2.0"
#import "@preview/plotst:0.2.0": *

= Linearized displacement fields

Consider a smooth deformation of some material.
Given a point with coordinates $(x, y, z)^T$, the position of the point after deformation is given by

$ mat(x'; y'; z') = mat( x; y; z ) + mat(u_x (x, y, z); u_y (x, y, z); u_z (x, y, z)). $

Taylor expanding $bold(u)$ at the origin,

$ 
mat(u_x (x, y, z); u_y (x, y, z); u_z (x, y, z)) = 

mat(u_x (0, 0, 0); u_y (0, 0, 0); u_z (0, 0, 0))
 + mat(
    diff_x u_x (0, 0, 0), diff_y u_x (0, 0, 0), diff_z u_x (0, 0, 0);
    diff_x u_y (0, 0, 0), diff_y u_y (0, 0, 0), diff_z u_y (0, 0, 0);
    diff_x u_z (0, 0, 0), diff_y u_z (0, 0, 0), diff_z u_z (0, 0, 0)
 )
 mat(x; y; z) \
 + O(x^2 + y^2 + z^2).
$

More compactly, we write

$ x_i ' = x_i + (u^0_i + (diff_j u_i) x_j), $

neglecting the higher order terms. Here $u^0_i$ is the displacement vector at the origin, and $diff_j u_i$ is the displacement gradient tensor.
The displacement gradient tensor determines a linear vector field. A linear vector field is determined by a matrix with real entries. (The corresponding linear map converts a position to the vector at that position.) For example, if $u = (2y, 0, 0)$, then
$ diff_j u_i = mat(0, 2, 0; 0, 0, 0; 0, 0, 0), $
and the vector field in the $x$-$y$ plane is as shown in @fig-linear-vector-field:

#figure(
    caption: [The vector field of the displacement $u = (2y, 0, 0) = mat(0, 2, 0; 0, 0, 0; 0, 0, 0) mat(x; y; z)$ in the $x$-$y$ plane.],
    box(image("figures/vector_field_traditional_shear.svg"), width: 50%)
) <fig-linear-vector-field>


= Scalar-Vector-Tensor decomposition

In order to analyze the displacement gradient tensor and its connection to stress and strain, we decompose it in terms of irreducible representations of the rotation group $"SO"(3)$. These are subspaces such that first projecting and then rotating is equivalent to first rotating and then projecting.

There is a well-developed theory for how to decompose any representation and construct the projections. A general rank-2 tensor in three dimensions decomposes into a scalar, vector, and tensor part. The irreducible representations of $"SO"(3)$ are uniquely characterized by their dimensions, and the decomposition can be expressed as:

$ bold(3) times.circle bold(3) = bold(1) plus.circle bold(3) plus.circle bold(5). $

In terms of the corresponding 3×3 matrices, these are the diagonal, antisymmetric, and symmetric-traceless parts. The corresponding orthogonal projection operators $P^s$, $P^v$, and $P^t$ are:

$
P^s_(i j k l) &= 1/3 delta_(i j) delta_(k l) \
P^v_(i j k l) &= 1/2 delta_(i k) delta_(j l) - 1/2 delta_(i l) delta_(j k) \
P^t_(i j k l) &= 1/2 delta_(i k) delta_(j l) + 1/2 delta_(i l) delta_(j k) - 1/3 delta_(i j) delta_(k l)
$

For a rank-2 tensor $T_(i j)$ and an operator $P$, the applied tensor $P T$ is given by
$ (P T)_(i j) = P_(i j k l) T_(k l). $

Note that the sum of projections $P^s + P^v + P^t$ is the identity operator $delta_(i k) delta_(j l)$.

The scalar part represents isotropic compression or expansion.
The vector part represents rotation.
The tensor part represents pure shear.

#figure(
    caption: [
      Scalar, vector, and tensor parts of a general rank-2 tensor.
      The matrices are
      $mat(-1, 0, 0; 0, -1, 0; 0, 0, 0)$,
      $mat(0, 1, 0; -1, 0, 0; 0, 0, 0)$, and
      $mat(0, 1, 0; 1, 0, 0; 0, 0, 0)$, respectively,
      restricted to the $x$-$y$ plane.
      ],
    grid(
        columns: 3,
        box(image("figures/vector_field_compression.svg"), width: 90%),
        box(image("figures/vector_field_rotation.svg"), width: 90%),
        box(image("figures/vector_field_pure_shear.svg"), width: 90%)
    )
) <fig-tensor-decomp>

As an example, consider the shear transformation given by the matrix

$ mat( 0, 2, 0; 0, 0, 0; 0, 0, 0 ) $

This can be decomposed into an antisymmetric rotational part and a symmetric pure shear part as in @fig-shear-decomp:

$ mat( 0, 2, 0; 0, 0, 0; 0, 0, 0 ) = mat( 0, 1, 0; -1, 0, 0; 0, 0, 0 ) + mat( 0, 1, 0; 1, 0, 0; 0, 0, 0 ). $

#figure(
    caption: [Traditional shear transformation decomposed as the sum of rotational and pure shear components],
    grid(
        columns: 3,
        box(image("figures/vector_field_traditional_shear.svg"), width: 90%),
        box(image("figures/vector_field_rotation.svg"), width: 90%),
        box(image("figures/vector_field_pure_shear.svg"), width: 90%)
    )
) <fig-shear-decomp>

In terms of elasticity, the overall translational part $u^0$ and the rotational/antisymmetric part of $diff_j u_i$ combine to form a rigid motion, and thus don't contribute to the stress. The remaining symmetric part of $diff_j u_i$ is the strain tensor

$ epsilon_(i j) = 1/2 (diff_j u_i + diff_i u_j) $

and determines the stresses in the material.

= Elastic Moduli

The response of a material to a particular strain can be described by the rank-2 stress tensor $sigma_(i j)$ giving the pressure vector components in a direction $i$.

A consequence of the following representation theory is that any isotropic material is characterized by two elastic moduli: the bulk modulus $K$ and the shear modulus $mu$.
The bulk modulus $K$ describes the incompressibility of a material, while the shear modulus $mu$ describes its resistance to shearing.

The stress tensor is, like the strain tensor, also symmetric. This is typically argued in terms of conservation of angular momentum, but it also follows as a consequence of isotropy from purely representation theoretic considerations.

Both stress and strain transform under the action of the rotation group $"SO"(3)$ as the symmetric square of the vector representation.  and their irreducible subrepresentations are the diagonal (1-D) and traceless (5-D) representations.
The elasticity tensor, which intertwines the stress and strain tensors, must act as a scalar on each irreducible representation.
Schur's lemma implies that there is a block decomposition along these irreducible representations, and blocks between isomorphic irreducible subrepresentations are diagonal, and all other blocks are zero.

The stress-strain tensor giving the stress as a function of strain is given explicitly by

$
C &= 3K P^s + 2mu P^t \
C_(i j k l) &= 3K (1/3 delta_(i j) delta_(k l)) + 2mu (1/2 delta_(i k) delta_(j l) + 1/2 delta_(i l) delta_(j k) - 1/3 delta_(i j) delta_(k l)).
$

so that

$ sigma_(i j) = C_(i j k l) epsilon_(k l) $

The corresponding Schur multipliers are $3K$ on the diagonal and $2mu$ on the traceless.
(The factors of $3$ and $2$ are due to repetitions when the tensors are contracted.)
These are the eigenvalues of the mass matrix of the coupled harmonic oscillators, with multiplicities 1 and 5, respectively.
The mass matrix must be positive definite, and this happens precisely when $K > 0$ and $mu > 0$.


The inverse is given by the compliance tensor for strain as a function of stress is given by

$ S &= 1/(3K) P^s + 1/(2mu) P^t \
S_(i j k l) &= 1/(3K) (1/3 delta_(i j) delta_(k l)) + 1/(2mu)(1/2 delta_(i k) delta_(j l) + 1/2 delta_(i l) delta_(j k) - 1/3 delta_(i j) delta_(k l)), $

so that

$ epsilon_(i j) = S_(i j k l) sigma_(k l). $

= Poisson's Ratio

Poisson's ratio is defined as the ratio of transverse strain to axial strain under uniaxial stress.

If $1$ is the axial direction and $2$ is a transverse direction, then unit axial stress is given by $sigma_(k l)=delta_(k 1) delta_(l 1)$. Axial strain is given by $epsilon_(11) = S_(11 k l) sigma_(k l) = S_(11 11)$. Transverse strain is typically negative, and is thus counted with a minus sign by $-epsilon_(22) = S_(2211)$. Poisson's ratio is therefore given by

$ nu = (-epsilon_(22)) / epsilon_(11) = (-S_(22 11)) / S_(11 11) = (1/(6mu) - 1/(9K)) / (1/(9K) + 1/(3mu)) = (3K - 2mu)/(6K + 2mu) = (3 - 2(mu/K))/(2(3 + (mu/K))). $

The Poisson's ratio depends only on the ratio of the shear modulus to the bulk modulus, with the inverse given by

$ mu / K = (3(1-2nu))/(2(1+nu)). $

The requirement that both $K$ and $mu$ are positive restricts the possible range of Poisson's ratio to:

$ -1 < nu < 0.5. $

Normal materials have $nu > 0$, and correspondingly $K > 2/3 mu$. Curiously, there do exist exotic materials such as auxetic foams that have $nu < 0$. 

In the plot below we show these rays in the positive $K$–$mu$ plane. In particular, the rays represent:
- $nu = 0$ (typical for cork),
- $nu = 0.5$ (as in rubber), and
- $nu = -1$ (an idealized limit as used in the Hoberman Sphere mechanism).

#figure(
  caption: [Contour plot of Poisson's ratio in the $mu$-$K$ plane.],
  box(
    image("poisson-ratio-plot.svg"),
    width: 80%
  )
)

= Navier-Cauchy Equation

To be in equilibrium, the stress tensor must be divergence-free:

$ diff_i sigma_(i j) = 0. $

This is known as the Navier-Cauchy equation.

To verify:

$
mu med (nabla dot nabla) thin u + (K + mu/3) nabla (nabla dot u) = 0.
$

$
(1 - 2 nu) (nabla dot nabla) thin u + nabla (nabla dot u) = 0.
$

$
F_i = sigma_(i j) u_j
$

$
0 &= 2 "div"(F) = 2 diff_i (sigma_(i j) u_j) = 2 u_j diff_i sigma_(i j) + 2 sigma_(i j) diff_i u_j \
  &= u_j diff_i (diff_j u_i + diff_i u_j) + (diff_j u_i + diff_i u_j) diff_i u_j \

$



