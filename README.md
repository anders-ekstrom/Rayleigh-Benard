## Governing Equations

The flow is governed by the incompressible Boussinesq equations in nondimensional form.

### Continuity
$\frac{\partial u_j}{\partial x_j} = 0$

### Momentum
\[
\frac{\partial u_i}{\partial t}
+ \frac{\partial}{\partial x_j}(u_i u_j)
= -\frac{\partial p}{\partial x_i}
+ \left(\frac{\mathrm{Pr}}{\mathrm{Ra}}\right)^{1/2}
\frac{\partial^2 u_i}{\partial x_j \partial x_j}
+ \theta \delta_{i2}
\]

### Temperature Transport
\[
\frac{\partial \theta}{\partial t}
+ \frac{\partial}{\partial x_j}(u_j \theta)
= \left(\frac{1}{\mathrm{Pr}\,\mathrm{Ra}}\right)^{1/2}
\frac{\partial^2 \theta}{\partial x_j \partial x_j}
\]

## Boundary Conditions

The boundary conditions are:

- **Periodic** in the horizontal direction: \(x=0\) and \(x=L\)
- **No-slip walls** at \(y=0\) and \(y=H\):
  \[
  u = v = 0
  \]
- **Constant temperature** at \(y=0\) and \(y=H\):
  \[
  \theta = \theta_{\text{bottom}}, \quad \theta = \theta_{\text{top}}
  \]
