# Matrix Models with Quantum Computation
## Introduction

This project is the first part of my master thesis where we try to simulate Matrix Models which corresponds to 2D Quantum Gravity Models. In the following, you can find a small introduction to the project. If you want to learn further, please look up to my research presentation pdf "2024_05_16_Varol_Burak". 
  
  This project will be further developed and optimised, in the hope of demonstrating a quantum advantage. Even if it does not, it gives a way of modelling some important integrals that is commonly encountered in physics.

## Theory

Our Goal is to simulate (hermitian) Random Multi Matrix Models of the kind:

![equation](https://latex.codecogs.com/svg.image?Z=\int\prod_{\alpha=1...\nu-1;i=1...N}d\lambda_i^{(\alpha)}\Delta(\Lambda^{(1)})e^{-\sum_i&space;V(\lambda_i^{(\alpha)})}\Delta(\Lambda^{(\nu-1)}).)

with 
![equation](https://latex.codecogs.com/svg.image?V(M^{(\alpha)})=\sum_{\alpha=1}^{\nu-1}V_{\alpha}M^{(\alpha)}-\sum_{\alpha=1}^{\nu-2}c_{\alpha}M^{(\alpha)}M^{(\alpha&plus;1)})

which are shown to be dual to matter coupled with 2D Quantum Gravity. Our end goal is to model and calculate meaningful expectation values without assuming coupling strength "g" is small. We therefore want to explore the landscape of different couplings and search for interesting phenomena. (eg. a phase transition, indicating a possible Conformal Field Theory ).

Because of the large dimensions of the NxN random matrices, providing a speedup with quantum computing can result into a possible investigation of Large N Limit (which is known analytically) and "middle" N limit of the theories.

## Model 

Our latest results are shown for the One Matrix Model:
![equation](https://latex.codecogs.com/svg.image?&space;Z=\int&space;dM\exp{-tr(V(M))})

with ![equation](https://latex.codecogs.com/svg.image?V(M)=M^2&plus;\sum_{k\geq&space;3}\alpha_k&space;M^k)

which would describe pure gravity and tessellation of the space. Integral is  parameterised  with eigenvalues and unitary matrices and reduce the integral to consecutive real integrals with the Vandermonde determinant: 
![equation](https://latex.codecogs.com/svg.image?Z=\int\prod_i^N&space;d\lambda_i\Delta(\Lambda)^2\exp{(-\sum_i&space;V(\lambda_i))})

We model this with operators

![equation](https://latex.codecogs.com/svg.image?\begin{align*}\lambda&=\sigma_z^{\pm}\sum_{m=0}^\infty\lambda_m&space;2^{m}\quad\text{with}\quad\lambda_m:=\frac{1-\sigma_z^{(m)}}{2}\\hat{x}_\lambda&=a\left(\frac{1}{2}\sigma_{z}^{(\pm)}&plus;\lambda\right)\end{align*})

discretises one eigenvalue with (ket is the state in bracket notation)

![equation](https://latex.codecogs.com/svg.image?ket({\pm\lambda})=ket({\lambda_\pm,\lambda_0,\lambda_1,\lambda_2,\dots}))

For multiple eigenvalues we have toi tensor product this state with similar states that describes the discretisation of other eigenvalues. For an NxN Random Matrix we need DxN qubits, where D is the number of discretisation points and N the dimension of Matrix.

Performing the Imaginary Time Evolution of the Hamiltonian:
![equation](https://latex.codecogs.com/svg.image?H=\sum_{i=0}^{N-1}1^{\otimes&space;i}\otimes&space;V(\hat{x}_\lambda)\otimes&space;1^{\otimes((N-1)-i)})

creates us the distributions we want to evaluate. 

Because of the current state of the quantum hardware in NISQ era, we use Variational Quantum Imaginary Time Evolution (varQITE) algorithm to perform the simulations. The Ansatz can be found in the repository 








## Usage

There are several scripts that are designed to increase given parameter and holds the others constant. Mainly, you must create a .dipf file in the example template:

    'depth': int,
    'D': int,
    "qubits_per_dim": int,
    'beta': float,
     "a": float,
     "g": float,
     "power": int,
     "power_interaction": int,
     "num_timesteps":int,
and then use (for example):

    python3 -m one_matrix_main $path_to_dipf_file

the code then automatically saves the result of the evolution. After doing this one time, the optimal parameter values are found and saved, so every distribution can be immediately created by putting the values to the Ansatz. Then one can use the observables module to calculate expectation values. 
