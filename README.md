# Pi0EventDisplay
The LArSoft module and Python plotting scripts found here allow users to view simulated pi0 events in the DUNE far detector.  Useage instructions are as follows:

# Use analysis module to create output ROOT tree
The analysis module [Pi0EventDisplay_module.cc](https://github.com/lkashur/Pi0EventDisplay/blob/main/Pi0EventDisplay_module.cc) is used on the dunegpvms and creates an output ROOT tree titled FDPiZeroShower.root.  Credit to Jacob Larkin for developing the pi0 event/shower selection framework seen within the module.

The module is ran using the FHiCL files [runPi0EventDisplay.fcl](https://github.com/lkashur/Pi0EventDisplay/blob/main/runPi0EventDisplay.fcl) and [Pi0EventDisplay.fcl](https://github.com/lkashur/Pi0EventDisplay/blob/main/Pi0EventDisplay.fcl).

```
lar -c runPi0EventDisplay.fcl -S <text file of input files>
```

# View specific event using Python scripts
Running EventStigator.py will generate the plot for the selected event.  The scripts Pi0.py and EventPlotting.py will also need to be present in the working directory.

```
python3 EventStigator.py -id <event number:z-vertex> Pi0EventDisplay.root
```
