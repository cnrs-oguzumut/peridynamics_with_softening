# Peridynamics With Softening: Model and Code Details

This document explains the active model implemented in this repository: the
equations, discretization, boundary conditions, damage laws, numerical solver,
outputs, and the role of each source file.

The current executable is a one-dimensional tensile peridynamic bar. By default
the left grip is fixed, the right grip is pulled in displacement control, and
the virtual substrate coupling is disabled. Many older two-dimensional, thermal,
and substrate-grid ideas are still present as commented code, but the active
calculation is one-dimensional in `x`.

## Source File Map

| File | Purpose |
| --- | --- |
| `main.cpp` | Defines the simulation parameters, allocates fields, builds geometry and correction factors, calls `mechanical`, then frees memory. |
| `Geometry.h` | Builds the 1D nodal coordinates and the peridynamic family list for each node. |
| `surface_correction_factors.h` | Computes boundary/horizon correction factors used in the bond integration weight. |
| `matrices.h` | Allocates, initializes, precomputes, and frees all vector and pointer-based arrays. |
| `only_mechanical.h` | Contains the main quasi-static mechanical solve: load increments, boundary conditions, bond forces, softening/damage, substrate coupling, ADR integration, energy and output writing. |
| `DefectUtils.h` | Builds a node mask for a possible defect region. In the current force loop, the defect critical stretch is not applied because the corresponding code is commented out. |
| `DiagnosticsUtils.h` | Helpers for selecting diagnostic nodes and writing optional diagnostic traces. |
| `time_settings.h` | Defines the mechanical and thermal time-step constants. Only the mechanical values are active. |
| `variable_initialization.h` | Initializes scalar scratch variables used by the included-style C++ code. |
| `write_ovito_files.h` | Legacy OVITO writer. It is not included by the active `main.cpp`. |
| `plot_strain.py` | Post-processing script for `result/conf*.dat` files. It can plot or animate strain, displacement, and damage. |
| `Makefile` | Builds the executable `peridynamics_with_softening`. |

## Build and Run

Build:

```bash
make
```

Run with defaults:

```bash
./peridynamics_with_softening
```

Run with command-line control:

```bash
./peridynamics_with_softening [increments_to_run] [load_increments] [damage_model]
```

The first argument controls how many outer increments are executed. The second
argument controls the denominator used for the strain increment,
`target_total_strain / load_increments`.

The third argument selects the damage law:

| `damage_model` | Law |
| --- | --- |
| `-1` | Elastic, no damage |
| `0` | Brittle failure |
| `1` | Linear softening, default |
| `2` | Bilinear softening |
| `3` | Cornelissen softening |
| `4` | Exponential softening |

For example:

```bash
./peridynamics_with_softening 40 80 3
```

runs 40 increments while using the load step size for an 80-increment ramp and
the Cornelissen softening law.

## Active Geometry

The active geometry is a 1D bar/film:

```text
x in [0, L]
y = 0
```

Current values from `main.cpp`:

| Quantity | Code variable | Current value |
| --- | --- | --- |
| Length | `length` | `1` |
| Number of material points | `ndivx`, `totnode`, `totint` | `300` |
| Cross-sectional area | `area_bar` | `1` |
| Thickness | `thickness` | `1`, currently not used in the active force law |
| Maximum family storage | `maxfam` | `20` |

The grid spacing is

```math
\Delta x = \frac{L}{N-1}.
```

The nodal coordinates are

```math
x_i = i\Delta x,\qquad y_i = 0,\qquad i = 0,\ldots,N-1.
```

The code names these nondimensional coordinates `coord_bar[i][0]` and
`coord_bar[i][1]`.

Each node has volume

```math
V = A\Delta x.
```

With the active constants:

```math
A=1,\qquad V=\Delta x.
```

## Peridynamic Horizon and Families

The horizon is

```math
\delta = 3\Delta x.
```

In code:

```cpp
double delta_bar = 3.0 * dx_bar;
double m = delta_bar / dx_bar;
```

Thus the horizon ratio is

```math
m = \frac{\delta}{\Delta x} = 3.
```

Node `j` belongs to the family of node `i` if

```math
i \ne j,\qquad |x_j-x_i| \le \delta + 10^{-12}.
```

Family data are stored in a compact, one-based style inherited from the older
code:

```text
numfam[i]      = number of neighbors of node i
pointfam[i]    = one-based offset into nodefam for node i
nodefam[...]   = one-based neighbor node index
```

When the code retrieves a neighbor, it converts back to a zero-based C++ index:

```cpp
const int cnode = nodefam[pointfam[i] + j - 1] - 1;
```

## Material Parameters

Current active values:

| Quantity | Code variable | Current value |
| --- | --- | --- |
| Poisson ratio | `pratio` | `1.0 / 3.0` |
| Film density | `rho_film` | `1` |
| Substrate density | `rho_substrate` | `1`, not active in the current virtual-substrate force |
| Film Young's modulus | `emod_film` | `1` |
| Substrate Young's modulus | `emod_substrate` | `1`, mostly inactive |
| Film plane-stress modulus | `emod_bar_film` | `E_f / (1 - nu^2)` |
| Film shear modulus | `mu_film` | `E_f / (2(1+nu))`, currently not used by active force law |

The active bond micromodulus is the 1D expression

```math
c_0 = \frac{2E}{A\delta^2}.
```

A discrete correction is then applied:

```math
c = c_0\frac{m}{m+1}
  = \frac{2E}{A\delta^2}\frac{m}{m+1}.
```

In code:

```cpp
double discrete_correction = m / (m + 1.0);
double bc_mechanical_film_bar =
    2 * emod_film / (area_bar * pow(delta_bar, 2)) * discrete_correction;
```

The code uses `bc_mechanical_film_bar` as the active peridynamic bond stiffness.

## Surface and Volume Correction

Near the two ends of the bar, a node has an incomplete horizon. The code applies
a geometric correction based on the visible horizon length.

For node `i`,

```math
\delta_i^- = \min(\delta,x_i),
```

```math
\delta_i^+ = \min(\delta,L-x_i),
```

and

```math
G_i = \frac{2\delta}{\delta_i^-+\delta_i^+}.
```

Interior nodes have `G_i = 1`. Boundary-near nodes usually have `G_i > 1`.

However, the first and last three displacement-prescribed grip nodes are forced
back to

```math
G_i = 1.
```

For a bond `(i,j)`, the correction is averaged:

```math
G_{ij} = \frac{G_i+G_j}{2}.
```

The code stores this as `scr_mechanical[i][j]`. The separate volume-correction
array `fac[i][j]` is currently always `1.0`.

The effective bond integration weight is

```math
w_{ij} = V\,G_{ij}\,\mathrm{fac}_{ij}.
```

In code this is precomputed as

```cpp
bond_weight[i][j] = vol_bar * scr_mechanical[i][j] * fac[i][j];
```

## Kinematics

The displacement of node `i` is

```math
u_i.
```

The current position is

```math
y_i = x_i + u_i.
```

For bond `(i,j)`, the initial separation is

```math
\xi_{ij} = x_j - x_i,
```

with length

```math
r_{ij}^0 = |\xi_{ij}|.
```

The current separation is

```math
\eta_{ij} = (x_j+u_j) - (x_i+u_i),
```

with length

```math
r_{ij} = |\eta_{ij}|.
```

The bond stretch is

```math
s_{ij} = \frac{r_{ij}-r_{ij}^0}{r_{ij}^0}.
```

The 1D direction is

```math
n_{ij} = \frac{\eta_{ij}}{r_{ij}}.
```

The current implementation computes the bond length using the `x` component
only:

```cpp
const double nlen2 = nx * nx;
const double nlength = std::sqrt(nlen2);
```

This matches the active 1D geometry, where all active `y` coordinates and
`y` displacements are zero.

For damage, the code stores the historical maximum stretch:

```math
s_{ij}^{max}(t) = \max_{\tau \le t} s_{ij}(\tau).
```

In code this is `bond_s_history[i][j]`.

## Loading and Boundary Conditions

The simulation is quasi-static. The applied strain starts at

```math
\varepsilon_0 = 0.
```

The target total strain is

```math
\varepsilon_{target} = 0.004.
```

The default number of load increments is

```math
N_{load} = 80.
```

The strain increment is

```math
\Delta\varepsilon =
\frac{\varepsilon_{target}}{N_{load}}.
```

With defaults:

```math
\Delta\varepsilon = \frac{0.004}{80} = 5\times 10^{-5}.
```

At every outer load increment:

```math
\varepsilon_k = \varepsilon_{k-1} + \Delta\varepsilon.
```

### Predictor Before Relaxation

Before the inner relaxation loop starts, the code applies a displacement
predictor to every node:

```math
u_i \leftarrow u_i + 0.75\,\Delta\varepsilon\,x_i.
```

This is controlled by

```cpp
const double load_predictor_fraction = 0.75;
```

The previous converged displacement field is therefore retained and nudged in
the direction of the next affine strain state before dynamic relaxation begins.

### Prescribed Grip Nodes

The first and last three nodes are displacement-controlled:

```text
i = 0, 1, 2
i = N-3, N-2, N-1
```

During every inner relaxation step, these nodes are reset to

```math
u_i = 0,
\qquad i=0,1,2,
```

and

```math
u_i = \varepsilon_k L,
\qquad i=N-3,N-2,N-1.
```

Their velocities are reset to

```math
v_i = 0.
```

Thus the left end satisfies

```math
u(0)=0,
```

and the right end satisfies

```math
u(L)=\varepsilon_k L.
```

The two neighboring grip nodes on each side are part of the grip blocks: the
left grip block is fixed, and the right grip block is pulled by the same
prescribed displacement.

The code also skips force, stress, strain-energy, substrate-energy, and damage
calculation for those prescribed force nodes:

```cpp
if (prescribed_force_node) {
    dmg[i] = 0.0;
    strain_energy[i] = 0.0;
    substrate_energy[i] = 0.0;
    stress[i][0] = 0.0;
    continue;
}
```

### Body Forces

The active body-force array is initialized to zero and no body force is applied:

```math
b_i = 0.
```

The commented traction-style boundary-force code is inactive.

## Virtual Substrate Coupling

The virtual substrate code path is optional and disabled by default:

```math
\texttt{enable_substrate_coupling}=\texttt{false},
\qquad
\texttt{c_int_ratio}=0.
```

Therefore, in the simpler tensile-bar scenario,

```math
p_i^{sub}=0,
\qquad
W_i^{sub}=0.
```

If the flag is re-enabled in the source, the substrate is represented by a
virtual layer below the bar with the same horizontal coordinates and a vertical
offset

```math
h = \Delta x.
```

In code:

```cpp
double h_interface_bar = 1.0 * dx_bar;
```

The virtual substrate is prescribed to follow the affine displacement field:

```math
u_j^{sub} = \varepsilon_k x_j.
```

For a film/bar node `i` and virtual substrate point `j`, the initial distance is

```math
r_{ij}^{0,sub} =
\sqrt{(x_j-x_i)^2+h^2}.
```

Only substrate points inside the horizon are included:

```math
r_{ij}^{0,sub} \le \delta + 10^{-12}.
```

The current horizontal separation is

```math
d_{ij}^{sub} =
(x_j+u_j^{sub}) - (x_i+u_i).
```

The current distance is

```math
r_{ij}^{sub} =
\sqrt{(d_{ij}^{sub})^2+h^2}.
```

The interface stretch is

```math
s_{ij}^{sub} =
\frac{r_{ij}^{sub}-r_{ij}^{0,sub}}{r_{ij}^{0,sub}}.
```

The horizontal direction cosine is

```math
n_{ij}^{sub} =
\frac{d_{ij}^{sub}}{r_{ij}^{sub}}.
```

The interface stiffness is

```math
c_{int} = \alpha_{int}c,
```

where

```math
\alpha_{int}=1.
```

When enabled, the substrate force contribution added to node `i` is

```math
p_i^{sub} =
\sum_j c_{int}s_{ij}^{sub}n_{ij}^{sub}V.
```

The optional implementation uses `vol_bar` directly for the virtual-substrate sum. It
does not use the bar-bar surface correction `scr_mechanical` in this interface
term.

## Bar-Bar Bond Force

For an intact elastic bond, the active bond force contribution is

```math
f_{ij} =
c\,s_{ij}\,n_{ij}\,w_{ij}.
```

The internal force density at node `i` is accumulated as

```math
p_i^{bar} =
\sum_{j\in\mathcal{H}_i} f_{ij}.
```

The total force used for the relaxation update is

```math
p_i = p_i^{bar} + p_i^{sub}.
```

## Damage and Degradation Models

The active model is selected by

```cpp
int degradation_model = 1;
```

The options implemented in `only_mechanical.h` are:

| Value | Meaning |
| --- | --- |
| `-1` | Pure elastic, no damage |
| `0` | Brittle failure |
| `1` | Linear softening |
| `2` | Bilinear softening, also called trilinear in comments |
| `3` | Cornelissen nonlinear softening |
| `4` | Exponential softening |
| other | Fallback to brittle-style failure |

### Code-to-Equation Mapping Inside the Bond Loop

For each directed bond from node `i` to neighbor `j`, the code computes:

```math
r_{ij}^0 = \texttt{idist},
```

```math
w_{ij} = \texttt{w_const}
= V\,G_{ij}\,\mathrm{fac}_{ij},
```

```math
n_{ij} = \texttt{dirx}
= \frac{\eta_{ij}}{|\eta_{ij}|},
```

```math
s_{ij} = \texttt{stretch}
= \frac{|\eta_{ij}|-r_{ij}^0}{r_{ij}^0}.
```

Before the damage branch is evaluated, the history variable
`bond_s_history[i][j]` is updated:

```math
\hat{s}_{ij}^{n+1}
=
\max\left(\hat{s}_{ij}^{n},s_{ij}^{n+1}\right).
```

Below, `\hat{s}_{ij}` means this updated history stretch. The code then checks
whether the bond is on the softening envelope:

```math
\mathrm{onEnvelope}_{ij}
=
\left(|s_{ij}-\hat{s}_{ij}|<10^{-14}\right).
```

The scalar force contribution added to node `i` is

```math
\Delta p_{ij} = \texttt{dforce1_mechanical}.
```

It enters the nodal force balance as

```math
p_i^{bar} \leftarrow p_i^{bar}+\Delta p_{ij}.
```

It also enters the stress moment:

```math
M_i \leftarrow M_i+\Delta p_{ij}\xi_{ij},
\qquad
\sigma_i=\frac{1}{2}M_i.
```

For nodal damage, the code accumulates two local sums:

```math
A_i \leftarrow A_i + q_{ij}V\,\mathrm{fac}_{ij},
```

```math
B_i \leftarrow B_i + V\,\mathrm{fac}_{ij},
```

where

```math
q_{ij} =
\begin{cases}
1, & \text{elastic/intact bond},\\
\mu_{ij}, & \text{softened bond},\\
0, & \text{fully failed bond}.
\end{cases}
```

The nodal damage is then

```math
D_i=1-\frac{A_i}{B_i}.
```

Because the family table stores directed bonds, the code also checks the
reciprocal entry `(j,i)`. If either `fail[i][j]` or the reciprocal entry has
failed, the bond is treated as failed and carries no force.

If `fail[i][j] == 0` when the loop reaches the bond, the bond is already failed
and the code uses

```math
\Delta p_{ij}=0,\qquad
W_{ij}=0,\qquad
\Delta G_{ij}^{fracture}=0,\qquad
q_{ij}=0.
```

Already failed bonds therefore keep contributing to the denominator of the
damage ratio, but not to force, strain energy, or the damage numerator.

### Common Elastic Update

Whenever a branch treats the bond as elastic, the code uses:

```math
\Delta p_{ij}
=
c\,s_{ij}\,w_{ij}\,n_{ij},
```

```math
W_{ij}
=
\frac{1}{4}c\,s_{ij}^{2}\,r_{ij}^0\,w_{ij},
```

```math
\Delta G_{ij}^{fracture}=0,
\qquad
q_{ij}=1.
```

### Common Softening Update

For the linear, bilinear, Cornelissen, and exponential softening branches, the
irreversibility variable is the history stretch

```math
\kappa_{ij}(t)=\max_{\tau\le t}s_{ij}^{+}(\tau),
```

stored in the code as `bond_s_history[i][j]`. The prescribed softening envelope
is written as

```math
f_{env,ij}(\kappa)=c\,p(\kappa),
```

where, for the implemented linear, bilinear, Cornelissen, and exponential laws,
the code has `p(kappa) = s0 * mu_history` after the elastic limit. The degraded
unloading/reloading stiffness is

```math
g(\kappa)=\frac{p(\kappa)}{\kappa},
\qquad
k_{unload,ij}=c\,g(\kappa)
=\frac{f_{env,ij}}{\kappa}.
```

If the bond is on the envelope:

```math
\Delta p_{ij}
=
f_{env,ij}\,w_{ij}\,n_{ij}.
```

If the bond is unloading or reloading:

```math
\Delta p_{ij}
=
k_{unload,ij}\,s_{ij}\,w_{ij}\,n_{ij}.
```

In both cases, the strain energy stored by the code is the directed-bond
version of the free energy:

```math
W_{ij}
=
\frac{1}{4}
c\,g(\kappa_{ij})\,s_{ij}^{2}
r_{ij}^0
w_{ij}.
```

The factor `1/4` appears because this code accumulates directed bond
contributions at nodes; it is the pair-split counterpart of the usual
`\frac12 c r g(\kappa)s^2` bond free energy.

The damage numerator uses

```math
q_{ij}=\mu_{ij}.
```

### Common Complete-Failure Update

For the brittle, linear, bilinear, Cornelissen, and fallback branches, complete
failure sets

```math
\mathrm{fail}_{ij}=0,
```

and uses

```math
\Delta p_{ij}=0,\qquad
W_{ij}=0,\qquad
q_{ij}=0.
```

For the brittle, bilinear, Cornelissen, exponential, and fallback branches, the
current code adds the branch-local complete-failure value stored in
`d_fracture_energy`. For `degradation_model == 1`, the code instead uses the
incremental history-energy update described in the linear-softening subsection.

For linear softening at complete failure, the limiting accumulated value is

```math
G_{ij}^{fracture}
=
\frac14 c\,s_0s_c\,r_{ij}^0w_{ij}.
```

### Critical Stretch and Fracture Parameters

The critical elastic stretch is

```math
s_0 =
\sqrt{\frac{3G_c}{E\delta}}.
```

In code this is `s_initial`.

The current fracture-energy parameter is

```math
G_c = 10^{-8}.
```

The softening energy parameter is

```math
R_d = G_c.
```

The complete-failure stretch `s_c` depends on the selected model:

Linear softening:

```math
s_c =
1.5\left(s_0 + \frac{2R_d}{cs_0}\right).
```

Bilinear softening:

```math
s_c =
s_0 + \frac{2R_d}{cs_0}.
```

Cornelissen softening:

```math
s_c =
s_0 + \frac{5.1361R_d}{cs_0}.
```

Exponential softening:

```math
s_c = 10s_0.
```

For the exponential law, this `s_c` is only a numerical cutoff value in
`main.cpp`; the active exponential branch in `only_mechanical.h` now uses it as
the complete-failure cutoff.

The bilinear transition stretch is

```math
s_k = s_0 + 0.15(s_c-s_0).
```

The bilinear residual ratio is

```math
\beta = 0.25.
```

### Elastic Model, `degradation_model < 0`

The bond never damages or fails:

```math
f_{ij} = c\,s_{ij}\,n_{ij}\,w_{ij}.
```

Damage remains zero because the intact weight and total weight are the same.

### Brittle Model, `degradation_model == 0`

If

```math
s_{ij} \le s_0,
```

the bond is elastic.

If

```math
s_{ij} > s_0,
```

the bond is marked failed:

```math
\mathrm{fail}_{ij}=0,
```

and its force becomes zero.

The fracture-energy accumulator receives

```math
G_{ij}^{fracture} =
\frac{1}{4}c\,s_0^2\,r_{ij}^0\,w_{ij}.
```

### Linear Softening, `degradation_model == 1`

This is the current default.

If the history stretch is still elastic,

```math
s_{ij}^{max} \le s_0,
```

the bond uses the elastic law.

If

```math
s_0 < s_{ij}^{max} < s_c,
```

the degradation factor is

```math
\mu_{ij} =
\frac{s_c-s_{ij}^{max}}{s_c-s_0}.
```

The softening-envelope force magnitude is

```math
f_{env} = cs_0\mu_{ij}.
```

If the current stretch is on the envelope, the force is

```math
f_{ij} = f_{env}n_{ij}w_{ij}.
```

If the bond is unloading or reloading below its historical maximum, the code
uses a degraded secant stiffness

```math
k_{unload} = \frac{f_{env}}{s_{ij}^{max}},
```

and

```math
f_{ij} = k_{unload}s_{ij}n_{ij}w_{ij}.
```

If

```math
s_{ij}^{max} \ge s_c,
```

the bond is fully failed:

```math
f_{ij}=0,\qquad \mathrm{fail}_{ij}=0.
```

For `degradation_model == 1`, fracture energy is stored incrementally with a
per-bond history variable `bond_fracture_energy_history[i][j]`. During partial
softening the code computes

```math
W_{input}
=
\frac12 c s_0^2
+
\frac12 c s_0(1+\mu_{ij})(s_{ij}^{max}-s_0),
```

```math
W_{recoverable}
=
\frac12 f_h s_{ij}^{max},
\qquad
f_h=cs_0\mu_{ij},
```

and

```math
G_{now}
=
\frac12
\left(
W_{input}-W_{recoverable}
\right)
r_{ij}^0w_{ij}.
```

Only the positive increment

```math
\Delta G
=
G_{now}
-
G_{history}
```

is added to `fracture_energy[i]`, and then `G_history` is updated. At complete
failure, the total target value is

```math
G_{total}
=
\frac14 c s_0s_c r_{ij}^0w_{ij}.
```

### Bilinear Softening, `degradation_model == 2`

For

```math
s_0 < s_{ij}^{max} < s_k,
```

the degradation factor is

```math
\mu_{ij} =
1 - (1-\beta)
\frac{s_{ij}^{max}-s_0}{s_k-s_0}.
```

For

```math
s_k \le s_{ij}^{max} < s_c,
```

the degradation factor is

```math
\mu_{ij} =
\beta\,\frac{s_c-s_{ij}^{max}}{s_c-s_k}.
```

Then the same envelope/unloading logic is used:

```math
f_{env}=cs_0\mu_{ij},\qquad
k_{unload}=\frac{f_{env}}{s_{ij}^{max}}.
```

### Cornelissen Softening, `degradation_model == 3`

The normalized history variable is

```math
\eta =
\frac{s_{ij}^{max}-s_0}{s_c-s_0}.
```

The degradation factor is

```math
\mu_{ij} =
\left(1 + (3\eta)^3\right)\exp(-6.93\eta)
- 28\eta\exp(-6.93).
```

The code clamps this at zero:

```math
\mu_{ij} = \max(0,\mu_{ij}).
```

Then it uses the same envelope/unloading logic as the other softening laws.

### Exponential Softening, `degradation_model == 4`

The exponential branch uses

```math
\alpha = \frac{cs_0}{R_d},
```

and

```math
\mu_{ij} =
\exp\left[-\alpha(s_{ij}^{max}-s_0)\right].
```

Then

```math
f_{env}=cs_0\mu_{ij},\qquad
k_{unload}=\frac{f_{env}}{s_{ij}^{max}}.
```

If

```math
s_{ij}^{max} \ge s_c,
```

the bond is fully failed:

```math
f_{ij}=0,\qquad \mathrm{fail}_{ij}=0.
```

## Nodal Damage Variable

The code computes a scalar nodal damage-like variable:

```math
D_i =
1 -
\frac{
\sum_{j\in\mathcal{H}_i}\mu_{ij}V\,\mathrm{fac}_{ij}
}{
\sum_{j\in\mathcal{H}_i}V\,\mathrm{fac}_{ij}
}.
```

For elastic or intact bonds,

```math
\mu_{ij}=1.
```

For softened bonds, `mu_history` is used. For fully failed bonds, the numerator
does not receive that bond's contribution.

Thus:

```math
D_i = 0
```

means all counted bonds around node `i` are intact, while larger values mean
local stiffness has degraded or bonds have failed.

The grip nodes are assigned

```math
D_i=0
```

because their force loop is skipped.

## Stress Measure

The code computes a stress-like virial measure from the bar-bar bond forces:

```math
\sigma_i =
C_\sigma
\sum_{j\in\mathcal{H}_i} f_{ij}\xi_{ij}.
```

The coefficient is

```math
C_\sigma = \frac{1}{2}.
```

In code:

```cpp
sum_ac_of_stress[i] += dforce1_mechanical * dix;
stress[i][0] = stress_coeff * sum_ac_of_stress[i];
```

The substrate/interface force is added to the nodal force balance, but it is not
included in this stress-like moment sum.

## Output Strain

The stored `strain[i]` is a post-processed finite-difference strain, not the
bond stretch.

For the left endpoint:

```math
\varepsilon_0 =
\frac{u_1-u_0}{\Delta x}.
```

For interior nodes:

```math
\varepsilon_i =
\frac{1}{2}
\left[
\frac{u_i-u_{i-1}}{\Delta x}
+
\frac{u_{i+1}-u_i}{\Delta x}
\right].
```

For the right endpoint:

```math
\varepsilon_{N-1} =
\frac{u_{N-1}-u_{N-2}}{\Delta x}.
```

## Energy Quantities

### Elastic Strain Energy

For an elastic intact bond, the strain-energy contribution accumulated at node
`i` is

```math
W_{ij} =
\frac{1}{4}
c\,s_{ij}^2\,r_{ij}^0\,w_{ij}.
```

For softened unloading/reloading, the code uses

```math
W_{ij} =
\frac{1}{4}
k_{unload}\,s_{ij}^2\,r_{ij}^0\,w_{ij}.
```

### Fracture Energy

For `degradation_model == 1`, fracture energy is accumulated incrementally from
the bond history. The code stores a previous per-bond value
`bond_fracture_energy_history[i][j]`, computes the current history value
`G_now`, and adds only

```math
\Delta G = G_{now}-G_{history}
```

when this increment is positive. This prevents the same history energy from
being added repeatedly during ADR iterations. In the current code, the other
damage branches keep their branch-local complete-failure increments.

### Substrate Energy

If substrate coupling is enabled, the virtual substrate energy density at node
`i` is

```math
W_i^{sub} =
\sum_j
\frac{1}{4}
c_{int}
\left(s_{ij}^{sub}\right)^2
r_{ij}^{0,sub}
V.
```

For the current default tensile-bar scenario, this term is identically zero.

### Kinetic Energy

The kinetic energy stored per node is

```math
K_i =
\frac{1}{2}m_i v_i^2.
```

The current lumped mass-like value is

```math
m_i = 1.
```

### Damping Dissipation

Damping dissipation is accumulated as a history variable:

```math
\mathcal{D}_i^{n+1}
=
\mathcal{D}_i^n
+ m_i c_d v_i^2\Delta t.
```

### Total Energy Density

The total stored nodal energy density is

```math
E_i^{tot} =
W_i + W_i^{sub} + K_i + \mathcal{D}_i + G_i^{fracture}.
```

System-level totals are accumulated by multiplying nodal values by `vol_bar`.

## Dynamic Relaxation Solver

Each load increment is equilibrated by an explicit damped dynamic relaxation
loop.

The active equation is

```math
m_i\ddot{u}_i + c_d\dot{u}_i = p_i + b_i.
```

Since `b_i = 0`,

```math
m_i\ddot{u}_i + c_d\dot{u}_i = p_i.
```

The acceleration used by the code is

```math
a_i^n =
\frac{p_i^n + b_i - c_d v_i^n}{m_i}.
```

The velocity update is

```math
v_i^{n+1} =
v_i^n + a_i^n\Delta t.
```

The displacement update is

```math
u_i^{n+1} =
u_i^n + v_i^{n+1}\Delta t.
```

The time step used by the mechanical solver is

```math
\Delta t = \mathrm{dt\_mechanical\_bar}.
```

`dt_mechanical_bar` is computed from

```math
\tau = \frac{L}{\sqrt{E/\rho}},
```

and

```math
\mathrm{dt\_mechanical\_bar}
=
\frac{\mathrm{dt\_mechanical}}{\tau}.
```

With the current default `L = E = rho = 1`, `tau = 1`, so

```math
\mathrm{dt\_mechanical\_bar} = 5\times 10^{-4}.
```

The code also computes `dt_pd`, `dt_lambda`, and `dt_bar`, but the active
mechanical update uses `dt_mechanical_bar`, not `dt_bar`.

## Damping

The default damping coefficient is

```math
c_d = 5.
```

In code:

```cpp
double cn = 5.0;
bool fixed_damping_coefficient = true;
```

If `fixed_damping_coefficient` is set to `false`, the code estimates an
adaptive damping value from changes in force acceleration:

```math
c_{opt} =
2\sqrt{\frac{\mathrm{num}}{\mathrm{den}}},
```

where the implementation accumulates

```math
\mathrm{num}
=
-\sum_i
\frac{u_i^2}{v_i}
\frac{
(p_i^n/m_i)-(p_i^{n-1}/m_i)
}{\Delta t},
```

and

```math
\mathrm{den} = \sum_i u_i^2.
```

The result is capped and smoothed:

```math
c_d \leftarrow 0.7c_d + 0.3c_{opt}.
```

This adaptive path is inactive by default.

## Convergence

The inner relaxation loop stops when the maximum nodal speed is below the
tolerance:

```math
\max_i |v_i| \le 10^{-6}.
```

In code:

```cpp
if (!fixed_nt_mechanical && maxval <= tolerance) {
    break;
}
```

Current defaults:

| Quantity | Code variable | Current value |
| --- | --- | --- |
| Tolerance | `tolerance` | `1.0e-6` |
| Maximum ADR iterations | `max_iteration` | `100000` |
| Fixed ADR iteration flag | `fixed_nt_mechanical` | `false` |
| Fixed ADR count if enabled | `nt_mechanical` | `10000` |

## Algorithm Summary

Initialization:

1. Read constants from `main.cpp`, `variable_initialization.h`, and `time_settings.h`.
2. Allocate arrays with `matrices`.
3. Build the 1D grid and family list with `build_Geometry`.
4. Compute surface correction factors with `surface_correction_factors`.
5. Precompute bond geometry and integration weights with `precompute_bond_invariants`.
6. Initialize displacement and velocity fields to zero.
7. Set the lumped mass-like coefficient to `1` for every node.
8. Call `mechanical`.

For each outer load increment:

1. Increase the applied strain by `del_epsilonn`.
2. Add the 0.75 affine predictor to all displacements.
3. Reset velocities and old force state for dynamic relaxation.
4. Run the ADR loop until convergence or maximum iteration count.
5. Write final fields for the increment.
6. Optionally detect first crack and early-stop after five more increments.

For each ADR iteration:

1. Reimpose displacement and zero velocity on the first and last three nodes.
2. For every unconstrained node, compute bar-bar bond forces.
3. Update each bond's historical maximum stretch.
4. Apply the selected degradation law.
5. Add virtual-substrate/interface force only if substrate coupling is enabled.
6. Compute nodal damage, stress, strain energy, and substrate energy.
7. Update damping if adaptive damping is enabled.
8. Update velocity and displacement explicitly.
9. Compute kinetic energy, damping dissipation, total energy, and convergence speed.
10. Compute finite-difference output strain.
11. Write optional diagnostics or snapshot buffers.

## Output Files

The simulation writes several text files in the repository root plus two output
directories.

### Per-Increment Configuration Files

For every outer increment, the code writes:

```text
result/confK.dat
```

where `K = loop_counter + 1`.

Each file begins with:

```text
# position displacement strain damage
```

Each data row contains:

```text
x_i  u_i  strain_i  damage_i
```

These files are the main input for `plot_strain.py`.

### OVITO Files

For every outer increment, the code writes:

```text
ovito_files/for_ovito_K.xyz
```

Each file contains:

```text
x y ux uy x_def y_def damage stress strain_energy strain
```

The active model is 1D, so `y`, `uy`, and `y_def` are written as `0`.

### Final ADR Field Matrices

The final converged values for each outer increment are appended to:

| File | Contents |
| --- | --- |
| `Final_ADR_stress.txt` | Final `stress[i][0]` for each node and increment |
| `Final_ADR_strain.txt` | Final finite-difference strain |
| `Final_ADR_strain_energy.txt` | Final strain energy |
| `Final_ADR_kinetic_energy.txt` | Final kinetic energy |
| `Final_ADR_fracture_energy.txt` | Final fracture energy |
| `Final_ADR_damping_dissipation_energy.txt` | Final damping dissipation |
| `Final_ADR_substrate_energy.txt` | Final virtual-substrate energy |

Rows correspond to outer increments. Columns correspond to nodes.

### Energy Histories During ADR

For every outer increment, the code writes a row of ADR-step system totals to:

| File | Contents |
| --- | --- |
| `total_energy_density_matrix.txt` | Total system energy |
| `total_kinetic_energy_density_matrix.txt` | Kinetic energy |
| `total_damping_dissipation_energy_density_matrix.txt` | Damping dissipation |
| `total_fracture_energy_density_matrix.txt` | Fracture energy |
| `total_strain_energy_density_matrix.txt` | Strain energy |
| `total_substrate_energy_density_matrix.txt` | Substrate/interface energy |

Each row is one outer increment. Each column is one ADR iteration.

### Stress Convergence

If

```cpp
bool track_stress_convergence = true;
```

the code writes:

```text
convergence_of_stress.txt
```

It tracks `stress[target_node][0]` at the node nearest

```cpp
target_x_coord_bar = 0.3;
```

Each row is one outer increment. Each value in the row is one ADR iteration.

### Snapshot Files

The snapshot mechanism can write all ADR-step fields for selected outer
increments:

```text
stressK.txt
strainK.txt
strain_energyK.txt
kinetic_energyK.txt
fracture_energyK.txt
damping_dissipation_energyK.txt
substrate_energyK.txt
total_energy_densityK.txt
```

By default:

```cpp
std::vector<int> write_snaps = {100000};
```

so no normal run reaches the requested snapshot increment. If
`stop_after_first_crack` is enabled, the cracking increment can be flushed.

The assignment to `write_snaps` after `mechanical` returns does not retroactively
write a snapshot for an already completed run.

### Average Stress File

The file

```text
Average_stress_of_interface.txt
```

is opened and given a header, but the current average-stress write loop is
commented out. In the active code, it does not receive per-increment data beyond
the header.

## Optional Diagnostics

If

```cpp
bool track_numerical_noise = true;
```

the code writes diagnostic traces:

```text
diagnostic_stress.txt
diagnostic_pforce.txt
diagnostic_disp_bar.txt
diagnostic_max_stretch.txt
```

For snapshot increments it can also write all-ADR-step versions:

```text
diagnostic_stress_all_ADR_steps_at_loop=K.txt
diagnostic_pforce_all_ADR_steps_at_loop=K.txt
diagnostic_max_stretch_all_ADR_steps_at_loop=K.txt
diagnostic_disp_bar_all_ADR_steps_at_loop=K.txt
```

The diagnostic nodes are the node nearest `target_x_coord_bar` and its immediate
neighbors, when available.

## First-Crack Early Stop

The code has an early-stop mode:

```cpp
bool stop_after_first_crack = false;
```

When enabled, the outer loop is allowed to run up to `max_extended_loop = 10000`
instead of the normal requested load count. The first crack is detected when any
node has

```math
D_i > 0.
```

After first crack detection, the simulation stops five increments later.

This mode is disabled by default.

## Defect Handling

`DefectUtils.h` can mark nodes inside a circular defect region:

```math
(x_i-x_c)^2 + (y_i-y_c)^2 \le r_{defect}^2.
```

It also computes a reduced critical stretch:

```math
s_0^{defect} = s_0 \times \mathrm{defect\_strength\_factor}.
```

Current defect defaults:

| Quantity | Code variable | Current value |
| --- | --- | --- |
| Enabled | `introduce_defect` | `false` |
| Center | `defect_x_center`, `defect_y_center` | `(0.0, 0.0)` |
| Radius | `defect_radius` | `dx_bar` |
| Strength factor | `defect_strength_factor` | `0.6` |

Important implementation detail: the force loop currently comments out the code
that would replace `s_initial` with the defect critical stretch. Therefore,
turning `introduce_defect` on currently builds and prints the defect mask, but it
does not change the active bond failure or softening law unless the commented
critical-stretch logic is restored.

## Post-Processing Script

`plot_strain.py` reads `result/conf*.dat`.

Static plot examples:

```bash
python3 plot_strain.py --steps 1 20 80
python3 plot_strain.py --steps 1 20 80 --displacement --damage
python3 plot_strain.py --steps 1 20 80 --strain-zoom --zoom-x-min 0.05 --zoom-x-max 0.95
```

Save a figure:

```bash
python3 plot_strain.py --steps 1 20 80 --damage --save strain_profiles.png
```

Create a movie:

```bash
python3 plot_strain.py --movie peridynamics.mp4 --damage --strain-zoom --fps 12
```

The script expects four columns:

```text
position displacement strain damage
```

If an older file has only three columns, the script fills damage with zeros.

## Active Versus Dormant Code

Several pieces of older model infrastructure remain in the source but are not
active in the current run:

| Dormant item | Current status |
| --- | --- |
| 2D film/substrate grid | Commented out; active model uses `y=0` for all nodes. |
| Thermal solve | `time_settings.h` has thermal values, but no active thermal loop is called. |
| Explicit traction boundary conditions | Present as comments; active loading is fixed-left and pulled-right displacement control. |
| Physical separate substrate material points | Disabled by default; optional virtual affine substrate coupling remains in the solver. |
| `write_ovito_files.h` | Not included by active `main.cpp`. |
| Defect critical stretch | Defect mask exists, but the local critical stretch substitution is commented out. |
| `dt_bar` stability estimate | Computed for reference, but active ADR uses `dt_mechanical_bar`. |
| Half-step velocity ADR variant | Older half-step update code remains commented out; active update is direct explicit Euler-style. |
| `heaviside` helper | Defined in `main.cpp`, but unused. |

## Key Implementation Notes

1. The active model is nondimensionalized through variables with the `_bar`
   suffix, but the code still uses familiar names such as `emod_film`, `rho_film`,
   and `Gc`.
2. The force and damage data are directional: the pair `(i,j)` and `(j,i)` have
   separate entries in `fail` and `bond_s_history`.
3. The stress measure includes bar-bar bond force moments only, not the optional
   virtual substrate force.
4. Grip nodes are skipped by the force loop; damage can occur in the free bar
   bonds once the selected law reaches its damage threshold.
5. Grip nodes are excluded from the force loop and report zero damage/stress
   from the active code path.
6. `maxfam` must be large enough to store every node's family. With
   `m = 3` in 1D, `maxfam = 20` is comfortably larger than needed.
7. The generated `result/`, `ovito_files/`, energy text files, and plot/movie
   outputs are run products, not source code.
