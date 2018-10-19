# MD_Code_One_Sided_Communication_MPI
converting the two-sided communication in the code to one-sided communication and compare the performance of the modified code with the original code.

## Description
You are given a MPI code (md_mpi.c) which does molecular dynamics simulations of inter-atomic interactions. You have to convert the two-sided communication in the code to one-sided communication and compare the performance of the modified code with the original code.

(MD is a computer simulation technique where the time evolution of a set of
interacting atoms is followed by integrating their equations of motion. Simply,
we can follow the laws of classical mechanics, and most notably Newton’s
law: Fi = miai for each atom i in a system constituted by N atoms. Here, mi
is the atom mass, ai = d2ri/dt2 its acceleration, and Fi the force acting upon
it, due to the interactions with other atoms. Interactions can be divided
into categories as intramolecular interaction and intermolecular interaction.
The former is modeled using harmonic approximations to describe the bonds
and bends between atoms. The latter is modeled using Lennard-Jones 12-
6 potential and coulombic interactions to describe the interaction of point
charges.)

