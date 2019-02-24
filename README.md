*********************************INTRODUCTION**********************************

Note the code is only for academic research use. Any direct use or modification of the code in commercial activity should contact the authors.

This code is for recovering complex field at focus (phase and amplitude) from a coherent defocus stack. The defocus stack includes multiple intensity images, which are modulus of thorough-focus coherent complex fields. It could be measured on a microscope by changing focus distance. 

The reconstruction method in the code obtains an estimated complex field at focus by solving a nonlinear optimization problem. The optimization problem minimizes the mean square error between the measured defocus stack and the predicted measurement, which is obtained by propagating the complex field at focus. The nonlinear optimization problem is solved with a modified Gauss-Newton method in the code

In the code, a Fresnel propagation kernel is used, since the NA of the microscope objective lens is small. For high NA case (e. g. NA=0.9), the propagation kernel code could be modified. But this modification is not included.

The derivation of the method in the code could be found in our paper: 
J. Zhong, L. Tian, P. Varma and L. Waller, "Nonlinear Optimization Algorithm for Partially Coherent Phase Retrieval and Source Recovery," in IEEE Transactions on Computational Imaging, vol. 2, no. 3, pp. 310-322, Sept. 2016.
doi: 10.1109/TCI.2016.2571669

*****************************HOW TO USE THECODE**************************

Example data sets:
The downloaded folder should include two defocus stack data sets: 
•	a simulation data set, named ‘SimulationCoherentDefocusStack.mat’; 
•	an experimental data, named ‘ExperimentCoherentDefocusStack.mat’. 
•	
The datasets contain:
Ividmeas: a defocused intensity stack
z: the positions of the measured intensity images [m]
ps: pixel size [m]
lambda: wavelength [m]
nfocus: the index of images which is focused

Open CoherentGradientMain.m and run it in Matlab. Just load the data set you would like to try, and click run.

Parameters to set:

Improve speed by running on GPU by setting:
IsGPUprogram=1; %choose whether run on gpu. 1 = yes; 0 = no.

Parameters affect the convergence speed:
MaxIter=50; % at main function; the iteration number the nonlinear optimization problem (outer loop)

maxiterCG=50; % at main function;  the iteration number for solving search direction with conjugate gradient method (inner loop)






