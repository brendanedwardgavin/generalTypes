# generalTypes

This module should prepare the basic objects that might be used in view of Portable Algorithmic Design with Julia (A suggestion for the Future name of the Library might be jPAD)

## Present considerations and blueprint for future actions
At present, in view of MPI compliancy at the lower level, the basic operations on the vectors have to be implememented by using the abstract type "myVec", which has been defined as a placeholder in the same way of "generalOp".
Such type should replace the "generalVec" as most as possible in the vector-specific operations, such as we have the list of the operations that might have to be implemented in the "userVec" type (which like "generalVec" with inherit from "myVec").

At the same time there is the need to provide some minimal "ergonomy" to the library at high level, that make the user aware of the quality of the design of the algorithm, either in prototyping and in the HPC-ready immplementation.

My suggestion for the next implementatation steps are therefore the following:

### Enable the insertion of other low-level components
 * Identify the list of myVec abstract operations that have to be implemented to preserve the functionalities for the code;
 * Use the operations defined in this way to introduce a "userVec" vector that allows flexibility in the distribution of the components;
 * Implement the MPI parallelisation of the components in the matmul_futile library, as a demonstrator;
 * Connect the BigDFT code with the wavefunctions as userVec and the Hamiltonian as userOp.


### Add high-level components to validate and visualize outputs

 * Define a set of non-regression tests for the two algorithm presented here, such as to validate future modifications in the library;
 * Define the non-regression also at the level of the external MPI layer in Julia;
 * Normalize the logfile such as to be inspected by the user to understand what is ongoing (a suggestion might be a yaml format for the logfile);
 * Introduce a set of "visualisation tools" that is able to provide at a glance the status of the convefgence with FEAST algorithm.
