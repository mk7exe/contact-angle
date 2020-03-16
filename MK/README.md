# Method of Khalkhali et al.
*Khalkhali M, Kazemi N, Zhang H, Liu Q, The Journal of Chemical Physics 146, 114704 (2017)*

All algorithms were coded in __MATLAB Version 8.6.0.267264__ and designed to read __LAMMPS dump files__ as inputs.
Functions:

..* read_LAMMPS_traj: This function reads a LAMMPS trajectory file and returns the coordinates of the centre of masses of water molecules.
..* hit_and_count: This Function applies the hit-and-count algorithm to remove outliers from a data set.
..* fine_percision: This function applies the fine precision droplet identification process to remove gas atoms which are very close to the liquid droplet surface.
..* contact_angle: This function calculates the contact angle of a liquid droplet of a solid surface through convex hull triangulation.
..* weighted_distribution: This function constructs a histogram of weighted data.
..* run.m: This is a sample run file to perform contact angle calculation.
