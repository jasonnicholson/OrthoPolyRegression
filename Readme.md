# Orthogonal Polynomial Regression

The monomial basis is exponential ill conditioned. A better approach is using an orthogonalization procedure that produces a basis that is orthogonal on the regression points. This has the effect of diagonlizing the least squares regression problem. This can also be done in multiple variables.

The plan here is to develop the code to implement Forsythe's paper from 1957. Then, extend the code using Weisfeld's paper for multiple variables.

Foresythe's algorithm is working and mostly documented. Give it a shot.

TODO:
- Document polyvalOrtho
- Write a packaging script.
- Implemented the multidimensiona version.
- Consider adding scaling to improve robustness.
- Consider converting to an object oriented library.
- 
