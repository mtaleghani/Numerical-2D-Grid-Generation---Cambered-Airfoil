# Elliptic O-Type Grid Generation with Control Functions for Various Geometries

### By
Seyed MohammadAmin Taleghani


## Table of Contents
1. [Description of the Developed Code](#description-of-the-developed-code)
   1. [Vinokur’s Stretching Function](#vinokurs-stretching-function)
   2. [NACA 4-Digit Symmetrical Airfoil](#naca-4-digit-symmetrical-airfoil)
   3. [NACA 4-Digit Cambered Airfoil](#naca-4-digit-cambered-airfoil)
   4. [Circular Cylinder](#circular-cylinder)
   5. [Domain Initializer](#domain-initializer)
   6. [Elliptic Solver](#elliptic-solver)
2. [Results](#results)
   1. [NACA 0012 Algebraic Grid (51x31)](#naca-0012-algebraic-grid-51x31)
   2. [NACA 4412 Algebraic Grid (51x31)](#naca-4412-algebraic-grid-51x31)
   3. [NACA 4412 Elliptic Grid (auto control function)](#naca-4412-elliptic-grid-auto-control-function)

## Description of the Developed Code

The coding was performed using C++ language. Inside the `main` function, various functions are called to perform each task.

### Vinokur’s Stretching Function

This function returns a custom clustered set of points. In this example, it was determined that the upper surface of the airfoil should have `imaxu` points, and the function provides `imaxu` points that are `dsa` and `dsb` units apart at the beginning and the end of the set.

![image](https://github.com/user-attachments/assets/f3592a4b-6774-4e85-a58f-f5d91ade5a37)


### NACA 4-Digit Symmetrical Airfoil

This function discretizes the surface of a four-digit symmetrical NACA airfoil with a thickness `t` relative to the unit chord and `x_vink` points on its chord line. The input data above discretize a NACA 0012 airfoil.

![image](https://github.com/user-attachments/assets/420aefd2-cb56-4236-8467-0dc1e11666f7)


### NACA 4-Digit Cambered Airfoil

Uses the formulation below to find the upper and lower surface points of a cambered NACA airfoil, where `m` is the maximum camber (100 * m is the first of the four digits), and `p` is the location of maximum camber (10 * p is the second digit in the NACA xxxx description).

![image](https://github.com/user-attachments/assets/914e28b2-330f-479a-97ec-20203c09d00d)

![image](https://github.com/user-attachments/assets/45f87646-d77c-4db6-a079-ecb2e9dad494)

![image](https://github.com/user-attachments/assets/3cb8d0c3-6804-411d-9a0c-80a084f14fdf)


The figure below shows the function input parameters set to generate NACA 2412 airfoil

![image](https://github.com/user-attachments/assets/90875252-2464-4341-b402-c7db8e1fff2b)


### Circular Cylinder

Provides the discretization of a circular cylinder with radius `r`.

### Domain Initializer

After the inner boundary has been chosen, `Initialize` creates an O-type outer boundary, with a predetermined radius inside the function. Then, if `algebric` was set to true, it will use linear interpolation to generate an initial algebraic grid. If `algebric` was set to false, the function would use an initial guess of (0,0) for every point inside the domain. Although elliptic solvers are able to handle this kind of initial guess, it is desirable to use a better initial guess for better convergence rate.

![image](https://github.com/user-attachments/assets/1dcdcecf-8b81-4019-970b-e6e101aa807a)


### Elliptic Solver

Uses Point Successive Over-Relaxation (SOR) to solve elliptic grid generation PDE. The SOR parameters for the control functions and grid points are `omega` and `omegaPQ`, respectively. The user can manually set the Q1 and P1 coefficients of the control functions to a desirable constant number for the entire boundary. Conversely, by setting `Auto = true`, the algorithm will determine the proper distribution of P1 and Q1 to satisfy orthogonality and the determined distance of the initial grid.

![image](https://github.com/user-attachments/assets/9a3149a5-e9e0-4bf9-9cfb-8691962901cd)

## Results

### NACA 0012 Algebraic Grid (51x31)

![image](https://github.com/user-attachments/assets/f1774d2f-130f-419b-9a9b-5f799405e05b)

![image](https://github.com/user-attachments/assets/116a841a-90a2-4de6-8325-2d0f2fdbe294)


### NACA 4412 Algebraic Grid (51x31)

![image](https://github.com/user-attachments/assets/ca07eafc-2a7c-4bd7-bbb7-f0f29509c7a1)

![image](https://github.com/user-attachments/assets/14ea20cd-0a3a-4fc0-a61d-c44d835d8b3b)

### NACA 4412 Elliptic Grid (auto control function)

![image](https://github.com/user-attachments/assets/e7a1bd1f-99e2-4e98-a8cf-c4b6e3a36e23)

![image](https://github.com/user-attachments/assets/42e3f5b0-5fbc-412a-b43d-5e4852c666f3)
