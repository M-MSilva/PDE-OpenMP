# P.D.E.-Poisson-OpenMP
## Program Objective

Solve a 2D Poisson partial differential equation by finite differences and using OpenMP.

## Solution

The implemented program is based on the pseudocode of the book (Burden and Faires, 2010). The solution for the PDE represented below:

![SerialResult.JPG](https://github.com/M-MSilva/PDE-OpenMP/blob/main/SerialResult.JPG)

where the subtitles numerico ( is the numerical solution), and analitico ( is the analytical solution). The average uncertainty of the resolution is in the ninth decimal place.

## Conceptual Model

![SerialResult.JPG](https://github.com/M-MSilva/PDE-OpenMP/blob/main/Fluxogram.JPG)

The code in c solves the Poisson elliptical P.D.E by calculating the matrix that approximates the solution to the P.D.E. The steps to calculate the matrix are illustrated above.

## Motivation

This program serves to learn the methods of gauss seidel and finite differences, as well as to solve any elliptic P.D.E category with the appropriate boundary condition. Boundary conditions and start function are shown at the beginning of the code and can be changed as needed. :smiley:

## Remarks about the program

The code using the shared memory API has a chunk size = 80, just for safety and to avoid possible future errors, since later we define the chunk size = multiple of the number of threads.

## Instructions for compiling and executing the code

### Initial requirements

To compile the code you need any compiler that has the OpenMP shared memory API (compilers list: https://www.openmp.org/resources/openmp-compilers-tools/), however it is recommended to use GNU compiler, because some specific optimization flags will be introduced later.

### Compiling and running the code

To compile the code in serial using the GNU compiler do:

```bash
gcc PDE_MM_SerialCode.c -lm
```

and run:

```bash
./a.out
```

To compile the code in parallel using the GNU compiler do:

```bash
gcc -fopenmp -o prg.x PDE_MM_OpenMPCode.c -lm
```
export the number of threads:

```bash
export OMP_NUM_THREADS=number of threads
```

and run the program at the end:

```bash
./prg.x
```
## Contributing 

Reviews and suggestions feel free to send me:

e-mail: marcosmatheusdepaivasilva@gmail.com

LinkedIn: linkedin.com/in/marcos-silva-089699b3 :hugs:

## Author

Marocs Matheus de Paiva Silva

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
