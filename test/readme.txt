scons -j2 co=1 b=GccOpt ts=projects/ozzy/test/VascularRemodelling/TestRunFlowInPipe.hpp
scons -j2 co=1 b=GccOpt ts=projects/ozzy/test/VascularRemodelling/TestSetupFlowInPipe.hpp


To setup a simulation

./projects/ozzy/build/optimised/VascularRemodelling/TestSetupFlowInPipeRunner 

To run it for the next timestep.

./projects/ozzy/build/optimised/VascularRemodelling/TestRunFlowInPipeRunner -start_time 1.0 -duration 1.0 -traction_file /home/ozzy/workspace_miguel/Chaste/projects/ozzy/test/data/straight_vessel_traction_it1.dat 



