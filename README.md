# Numerical Accuracy of the Derivative-Expansion-Based Functional Renormalization Group

Program analyzing the precision of a numerical implementation of the derivative-expansion-based functional renormalization group in the three-dimensional O(N) models. Copyright (C) 2024 Andrzej Chlebicki

This program serves as an appendix to the article under the same name. 
* Article status: Upublished
* ArXiv DOI: Unassigned
* Article DOI: Unassigned

The program takes the form of a Wolfram Language (.wl) package file. As such, it can be run via the terminal command
```console
foo@bar:fRGNumericalPrecision$ /path/to/Mathematica/WolframKernel -script fRGNumericalPrecision.wl
```
executed from the project's main directory. In the command, /path/to/Mathematica/ has to be replaced with the path to Mathematica executable files, for the default installation directory see [Wolfram System File Organization](https://reference.wolfram.com/language/tutorial/WolframSystemFileOrganization.html). The program does not accept any command-line arguments.

Alternatively, the program can be run via the Mathematica front end. This allows one to run only specific analyses but one has to be careful to initialize all the libraries at the top before performing any calculations.

<h2>Main analyses</h2>
The main part of the program is responsible for performing six analyses discussed in the article:

1. Tracking the compactification error (&rho;<sub>Max</sub> dependence)
2. Propagation of the discretization error (optimal &rho;<sub>Max</sub>=3.5&rho;<sub>0</sub>)
3. Propagation of the discretization error (low &rho;<sub>Max</sub>=2&rho;<sub>0</sub>)
4. Comparison of the discretization error between subsequent eigenvalues
5. Propagation of the loop-integral error
6. Error of the finite-difference stability-matrix approximation

Except for analyses 2 and 3, each analysis is defined by a separate pair of functions; the first function performs the analysis at the LPA level and the second at the order O(&part;<sup>2</sup>). The analyses 2 and 3 are performed by the same functions with an additional parameter &rho;MaxTo&rho;0 set to either 2 or 3.5 to distinguish between the two.

Each analysis automatically generates several plots presenting its results and inserts them into a corresponding directory. When calculations are performed for N=2, the plots used in the article are also copied into the "Article figures" directory.

The analyses are self-contained and should not have any side effects. Therefore, it should be possible to run them separately and in any order.

<h2>Extras</h2>
At the end of the program, we included two extra sections. The first one is devoted to testing the accuracy of the Gauss-Legendre integrals. It calculates the fixed-point solutions for N={1,2,3} and calculates the errors in the integration of the &beta; functions with different Gauss-Legendre parameter sets. This code was used to calculate the errors of the GL quadrature used in this study:

1. Low precision, q &isin;[0, 5], 20 integrand-evaluation points, maximal error: 10<sup>-5</sup>; 
2. Standard precision, q &isin;[0, 5], 35 integrand-evaluation points, maximal error: 10<sup>-10</sup>; 
3. Reference precision, q &isin;[0, 5.5], 50 integrand-evaluation points, maximal error: 10<sup>-13</sup>.

The second extra section presents the finite derivative operators included in the appendix of the paper.
