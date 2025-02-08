---
project: section_soot
author: Joshua Christopher
---

This is a Fortran implementation of a discrete sectional soot model.
We solve the General Dynamic Equations (GDE), also sometimes called the populaation balance equation.
For this study, a spatially homogeneous aerosol is assumed.
No variation in x- or y- coordinates.
All of the particles composed of molecules (monomers) are treated as of the same kind.
A particle containing $k$ monomers is named a $k$-mer.

We track the mass concentration variations from monomers up to $ND$-mer with 'discrete bins', and larger particles are divided up into $NS$ 'sectional bins', or 'sections'.
In total we will solve $NT=ND+NS$ eqautions.
Picture a 2D-axis.
On the x-axis is the particle mass expressed in number of molecules contained in the particles.
On the y-axis is the particle mass distributino function, $q_m$.

Discrete bins:
1. Have particle mass concentrated on discrete points
2. Particle mass is concentrated at discrete points of the x-axis

Sectional bins:
1. Particle mass is uniformly distributed within each section.
2. Within th ebin, $q_m dm$ gives the total particulate mass within the mass range of $m$ to $m+dm$ (where both $m$ and $q_m$ are continuous variables).
3. The uniformity of $q_m$ within bins is an assumption. 

So we use mass to split up the bins.
We could have also split up the bins by number, surface area, or mass concentration instead.
The CERFACS team split the soot up by volume.
With a uniform $q_m$, when converted to number distribution function using $n_m = q_m/m$ will give an apparently unphysical number distribution shape.
This doesn't make much of a difference to simulation results as long as you have enough bins.

## Assumptions
I just need to list all the assumptions so I don't forget about them.
1. Spatiallly homogenous aerosol.
2. All partiles composed of molecules are treated as of the same kind.
3. Uniformity of $q_m$ within a sectional bin.
4. Monomers do not coagulate

## Coagulation
The collision of two particles.
Particles from any two bins (or the same bin) can coagulate.
After coagulation, they will likely be put into a separate bin, or bins.
This is divided into two steps:
1. Calculate the mass of the particles that are colliding.
2. Decide which bins the coagulation products go into and how to distribute the particle mass among these bins.

