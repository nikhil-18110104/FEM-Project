# Finite Element Analysis Projects


Project 1: Bar Problem Problem statement

Figure 1 shows an elastic bar under traction load and constrained at the ends and . Develop a generic finite element code to get the approximate solution to the resulting governing differential equation for the bar shown in Figure 1.
The code should have the following capabilities:
1. Boundary conditions/End Constraints: Both ends can be constrained by specifying (a) primary variable (Dirichlet/Displacement/Essential),
(b) secondary variable or force (Force/Neumann/Natural) and
(c) springs (Mixed/Robin)
2. The variables ( ), ( ) and ( ) can vary from a constant to a quadratic function.
3. The length and the number of elements will be input values. Discretize the domain into given number of elements with equal lengths.
4. There should be a provision to put at least one concentrated load at any given location (excluding the ends).
5. Use of either Lagrange interpolation or hierarchic shape functions up to quartic order should be possible.
6. Postprocessing must be able to represent the primary, secondary and other variables over the domain either continuously or discretely as required.


Project 2: One-Dimensional Code for Beam Bending Problem

Write a one-dimensional finite element code using Hermite cubic shape functions with the following details for the beam bending problem.
1) Uniform cross-section: 1 cm X 1 cm
2) Length of the beam: 10 cm
3) E = 200 GPa
4) The code should be capable of handling the transverse loads of the type
a. Concentrated/point load
b. Uniformly distributed load
c. Point moments at the center of the beam length only
5) Further, it should be capable of applying the appropriate combination of boundary conditions at either of the ends as:
a. Specified transverse displacement
b. Specified slope of the transverse displacement c. Shear force
d. Bending moment
Now, take the appropriate values of loads as mentioned in Point # 4 above and perform the following finite element analysis using your code for 1, 4, 10, 50, and 100 elements.
1) Give continuous variation of transverse displacement and its slope
2) Give continuous variation of shear force and bending moment
3) Bending stress on the topmost line of the beam along its entire length.


