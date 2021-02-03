#!/bin/bash
# chmod 700 RunHemeLB
mpirun -np 15 hemelb -in /data/vascrem/testoutput/TestFlowThroughNonSymetricCollapse/2/HemeLBFluid/config.xml -out /data/vascrem/testoutput/TestFlowThroughNonSymetricCollapse/2/HemeLBFluid/results/ >/data/vascrem/testoutput/TestFlowThroughNonSymetricCollapse/2/HemeLBFluid/HemeLBTerminalOutput.txt
echo 'HemeLB has finished' > /data/vascrem/testoutput/TestFlowThroughNonSymetricCollapse/2/HemeLBFluid/WaitFile.txt
echo 'HemeLB simulation complete' 