# P.D.E.-Poisson-OpenMP
## Program Objective

Solve a 2D Poisson partial differential equation by finite differences and using OpenMP.

## Solution

The implemented program is based on the pseudocode of the book (Burden and Faires, 2010). The solution for the PDE represented below:

![SerialResult.JPG](https://github.com/M-MSilva/PDE-OpenMP/blob/main/SerialResult.JPG)

where the subtitles numerico ( is the numerical solution), and analitico ( is the analytical solution). The average uncertainty of the resolution is in the ninth decimal place.

## Motivation

This program serves both for learning the gauss seidel and finite difference methods, as well as for testing any category of P.D.E Elliptic with the appropriate boundary condition. Boundary conditions and initial function are shown at the beginning of the code and can be changed as required. 

Any questions my e-mail to contact is: marcosmatheusdepaivasilva@gmail.com.

## Remarks about the program

The code using the shared memory API has a chunk size = 80, just for safety and to avoid possible future errors, since later we define the chunk size = multiple of the number of threads.
