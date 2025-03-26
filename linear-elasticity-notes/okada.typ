= Okada's rectangular dislocations

The Python function `dc3dwrapper` computes the displacement and stress fields due to dislocation along a rectangular of a very particular type.

== Cartesian coordinates in half-space

We use a right-handed Cartesian coordinate system, with the $z$-axis denoting height. The elastic material is assumed to occupy the half-space $z <= 0$. The observation points will typically be on the boundary $z = 0$ because GPS stations are on the surface, but mathematically there's no reason to restrict the observation points to the surface.

Side remark: The boundary condition used along $z = 0$ for the Navier-Cauchy equation is that for a free surface, known as "traction-free", mathematically given by $sigma_(i j) n_j = 0$ where $n$ is a normal vector and $sigma$ is the stress tensor. In terms of strain tensor $epsilon_(i j)$, this implies that along the boundary $z = 0$, the strain tensor satisfies $epsilon_(x z) = 0 = epsilon_(y z)$, and $epsilon_(z z) = -nu/(1-nu) (epsilon_(x x) + epsilon_(y y))$.

== Geometry of the rectangle

In the conventions of `dc3dwrapper`, there are several implicit constraints on the orientation of the rectangle relative to the coordinate system.

Two of the edges of the rectangle (bottom and top) are assumed to be parallel to the $x$-axis, and the other two (left and right) are perpendicular to the $x$-axis. The shared angle of the left and right pair of edges when projected onto the $y$-$z$ plane is given by `dip` degrees above the $y$-axis, towards the $z$-axis. Thus when `dip` is $0$, the left and right edges are parallel to the $y$-axis, and when `dip` is $90$ degrees, the left and right edges are parallel to the $z$-axis.

== `dc3dwrapper` arguments

The `dc3dwrapper` function accepts the following arguments:

- `alpha`: this is defined in terms of Lame parameters as $alpha = (lambda + mu) / (lambda + 2 mu)$. Note that Poisson's ratio $nu = lambda / (2 (lambda + mu))$ is equivalent, with $nu = 1 - 1/(2 alpha)$ and $alpha = 1/(2(1 - nu))$.
- `xo`: a triple $(x, y, z)$ giving the Cartesian coordinates of the observation point.
- `depth`: the (constant!) $z$-coordinate of the bottom edge of the rectangle. (In the current data, `dip` is always between $0$ and $180$ degrees, but if this were not the case, then the corresponding edge would actually be the top edge.)
- `dip`: the angle of the perpendicular pair of edges, in degrees, above the $y$-axis, towards the $z$-axis.
- `strike_width`: a pair containing the minimum and maximum $x$-coordinates of the rectangle.
- `dip_width`: the minimum and maximum coordinate values in the "dip" direction. The dip coordinate preserves length along the rectangle, so the difference between these two values gives the length of the left and right edges. The dip coordinate is zero when $y=0$ and $z=-mono("depth")$.
- `dislocation`: the slip vector of the dislocation in components of strike, dip, and tensile/opening. The strike component is parallel to the $x$-axis, the dip component is along the dip direction, and the tensile/opening component is normal to the rectangle.

== Precise mathematical description of the dislocation

The unit vector in the dip direction is given by

$ (0, cos(mono("dip")), sin(mono("dip"))), $

and the unit normal to the plane of the rectangle is given by

$ hat(n) = (0, -sin(mono("dip")), cos(mono("dip"))). $

The slip vector $b$ is the (constant) value of the jump discontinuity in the displacement vector field across the rectangle. Specifically, it's the value of displacement infinitesimally offset in the $hat(n)$ direction minus the value in the $-hat(n)$ direction.

The plane of the rectangle intersects the point $(0, 0, -mono("depth"))$. Thus in strike-dip coordinates, the Cartesian coordinates for a point with strike-dip coordinates $(s, d)$ are given by

$ (s, d cos(mono("dip")), d sin(mono("dip")) - mono("depth")). $

The $z$-coordinate of the rectangle below a point on the surface with coordinates $(x, y, 0)$ is given by

$ z = y tan(mono("dip")) - mono("depth"). $

The vertices of the rectangle in Cartesian coordinates in order of bottom-left, bottom-right, top-right, top-left are:

$
(mono("strike_width[0]"), mono("dip_width[0]") cos(mono("dip")), mono("dip_width[0]") sin(mono("dip")) - mono("depth")), \
(mono("strike_width[1]"), mono("dip_width[0]") cos(mono("dip")), mono("dip_width[0]") sin(mono("dip")) - mono("depth")), \
(mono("strike_width[1]"), mono("dip_width[1]") cos(mono("dip")), mono("dip_width[1]") sin(mono("dip")) - mono("depth")), \
(mono("strike_width[0]"), mono("dip_width[1]") cos(mono("dip")), mono("dip_width[1]") sin(mono("dip")) - mono("depth")).
$

Given a slip vector specified by strike-dip-opening components $s = mono("dislocation[0]")$,  $d = mono("dislocation[1]")$, and $o=mono("dislocation[2]")$, the corresponding Cartesian components are given by

$
mat(b_x; b_y; b_z)
= mat(
    s;
    d cos(mono("dip")) - o sin(mono("dip"));
    d sin(mono("dip")) + o cos(mono("dip"))
).
$

Conversely, given a Cartesian slip vector $(b_x, b_y, b_z)$, the strike-dip-opening components are given by

$
mat(s; d; o)
= mat(b_x; b_y cos(mono("dip")) + b_z sin(mono("dip")); -b_y sin(mono("dip")) + b_z cos(mono("dip"))).
$
