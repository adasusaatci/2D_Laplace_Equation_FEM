# Numerical Solution of the 2D Laplace Equation Using the Finite Element Method (FEM)

## Project Overview

This project focuses on the numerical solution of the two-dimensional Laplace equation around a circular cylinder using the Finite Element Method (FEM) in MATLAB.

The study investigates the application of finite element discretization to partial differential equations (PDEs), with particular emphasis on mathematical formulation, mesh generation, numerical accuracy, and convergence analysis.

The numerical solution is compared with the analytical solution to evaluate the accuracy of the implemented method.

## Mathematical Formulation

The governing equation is the two-dimensional Laplace equation expressed in polar coordinates:

$$
\nabla^2 u =
\frac{1}{r}\frac{\partial}{\partial r}
\left(r\frac{\partial u}{\partial r}\right)
+\frac{1}{r^2}\frac{\partial^2 u}{\partial\theta^2}=0
$$

The computational domain is defined as:

- Radial direction: r ∈ [1, 2]
- Angular direction: θ ∈ [0, 2π]

Dirichlet boundary conditions are imposed on the inner and outer boundaries of the computational domain.

## Numerical Method

The Finite Element Method is employed to transform the governing partial differential equation into a discrete system of algebraic equations.

The numerical methodology includes:

- Weak formulation of the governing equation
- Computational domain discretization
- Finite element shape functions
- Element stiffness matrix formulation
- Global system assembly
- Implementation of Dirichlet boundary conditions
- Numerical solution of the resulting algebraic system

## Mesh Discretization

Three different mesh resolutions are considered to investigate the influence of spatial discretization on numerical accuracy:

| Mesh | Resolution |
|---|---|
| Coarse | 81 × 41 |
| Medium | 161 × 81 |
| Fine | 321 × 161 |

The mesh refinement study provides a basis for evaluating numerical convergence and discretization error.

## Error and Convergence Analysis

The numerical results are compared with the analytical solution using a normalized L2 error measure.

The analysis focuses on:

- Numerical solution accuracy
- Spatial discretization error
- Mesh dependency
- Numerical convergence
- Analytical and numerical solution comparison

## Results and Visualization

The project investigates the numerical solution distribution within the computational domain and the relationship between mesh resolution and solution accuracy.

## Tools and Technologies

- MATLAB
- Finite Element Method (FEM)
- Partial Differential Equations
- Numerical Linear Algebra
- Scientific Computing
- Computational Fluid Dynamics
- Numerical Error Analysis

## Project Documentation

The accompanying technical report presents the mathematical background, finite element formulation, numerical implementation, simulation results, and error analysis.

## Related Project

A complementary study using the Finite Difference Method (FDM) is available in the following repository:

[2D Laplace Equation — Finite Difference Method](https://github.com/adasusaatci/2D_Laplace_Equation_FDM)

Both projects investigate the same governing equation and computational domain using different numerical discretization techniques.

## Author

**Ada Su Saatci**

Aeronautical Engineer | Modeling & Simulation | Control Systems

[LinkedIn](https://www.linkedin.com/in/adasusaatci/) | [GitHub](https://github.com/adasusaatci)

---

This project was developed as part of graduate-level engineering coursework at Istanbul Technical University.
