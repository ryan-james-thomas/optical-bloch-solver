# Optical Bloch Solver

This suite of MATLAB code solves the optical Bloch equations (OBEs) for either the steady state or for evolution for a given time.  It is mainly designed for solving problems involving optical transitions in alkali metal atoms, although it could be extended to other systems.  It is *not* well-designed for solving problems involving time-varying parameters such as frequency sweeps or pulsed lasers - if you need that, I suggest meta-programming as the solution for these kinds of problems.

The main class for solving the time-evolution of the density matrix for a given quantum system is the `densityMatrix` class.  The main idea behind this class is that the user provides a "bare" Hamiltonian which is the Hamiltonian in the absence of coupling between levels due to external fields; a "coupling" Hamiltonian which describes the interaction with external fields; and a matrix of "decay" rates which indicate the rate at which populations and/or coherences decay between states.  With these supplied, the class can then calculate the Lindblad term representing the loss of coherence and transfer of populations due to the decay rates, and then combine that with the rate of change of the density matrix due to the unitary evolution, to get a total time derivative of the density matrix.  This time derivative is represented as a super-operator that is applied to the density matrix.  So if we flatten the density matrix from NxN to N^2x1, where N is the number of eigenstates, then the super-operator is now a matrix of size N^2xN^2.  This means that for constant parameter problems, such as constant laser fields with constant detunings, the solution for any time is then the matrix exponential of the super-operator multiplied by the time step, left-multiplied with the initial density matrix (flattened).
