# Method of Khalkhali et al.
*Khalkhali M, Kazemi N, Zhang H, Liu Q, The Journal of Chemical Physics 146, 114704 (2017)*

All algorithms were coded in __MATLAB Version 8.6.0.267264__ and designed to read __LAMMPS dump files__ as inputs.
Functions:

* __read_LAMMPS_traj__: This function reads a LAMMPS trajectory file and returns the coordinates of the centre of masses of water molecules.

* __hit_and_count__: This Function applies the hit-and-count algorithm to remove outliers from a data set.

* __fine_percision__: This function applies the fine precision droplet identification process to remove gas atoms which are very close to the liquid droplet surface.

* __contact_angle__: This function calculates the contact angle of a liquid droplet of a solid surface through convex hull triangulation.

* __weighted_distribution__: This function constructs a histogram of weighted data.

* __run.m__: This is a sample run file to perform contact angle calculation.
