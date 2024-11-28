# flux_vanilla
Code for running the flux emergence simulations provided by David MacTaggart. That original fortran code has been modified to output files in the netcdf format, to make interaction with python much easier. Many additional files are included, which I shall summarise shortly...

The main fortran code will only run on Archie-West because it's old. The makefile for that isn't included here as it baffles my local machine. However, the field line tracing code in the folder 'fltrace' works fine everywhere and the makefile determines the correct machine before beginning. 

The emergence code isn't set up for multiple runs as it takes a while. All data is saved in "./Data/" unless this is changed in oliver_output.F90. This is as normal netcdfs and saves the magnetic field, internal energy and density, with staggering as appropriate. Unlike most of my code, there is no "run.py" with the parameters -- these are set in "input.deck" instead. 

The bash script to run on Archie is:

module purge
module load intel/intel-2020.4
module load netcdf-fortran/intel-2020.4/4.5.4
module load anaconda/python-3.8.8/2021.05
make

rm -rf ./Data/*.nc
mpirun ./bin/lare3d 0

There are various options for plotting:

"plotslice.py" crudely does some slices through the field. Dimensions etc. will need to be put in manually for now. The field line tracer this talks to appears to have gone missing, but it wasn't very good anyway. To do such things, look in the folder 'fltrace'.

"lorentz_emergence.py" is a new one which plots slices of the lorentz force (either normalised or not) at various altitudes. Compares two fields with specified data directories (which must be changed manually).

"emergence_differences.py" uses Chris' topological code to plots the differences in FLH, Twist and Winding etc., using Daining's fancy fourier transform method. This has now been fixed.

The fltrace folder contains my fortran field line tracer, which is speedy and uses the staggered grid proper like. It is currently set up to obtain data from two different runs, calculate the differences between the FLH on the surface and plot the field lines which originate mainly from that area. Outputs are pyvista screenshots saved to "plots".

All of these python scripts deal with a single snapshot, which is determined by typing "plotslice.py 300" etc., for plot number 300. 

Carry on.


