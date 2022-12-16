# Floquet-PDE
Solves for Floquet modes and exponents of a 1-dimensional PDE instead of decoupled ODEs

Uses Fourier pseudospectral discretization scheme to solve for the Floquet modes and Floquet exponents of linear fluctuations
around background solutions with spatial structure.  Specifically, the equation is of the form

$$\frac{\partial^2\delta\phi}{\partial t^2} - \frac{\partial^2\delta\phi}{\partial x^2}\delta\phi + \left(k^2 + V''(\bar{\phi}_{\rm bg})\right)\delta\phi = 0$$

and we want to decompose $\delta\phi$ into an appropriate set of basis functions.
