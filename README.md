# Pi0EventDisplay
The LArSoft module and Python plotting scripts found here allow users to view simulated pi0 events in the DUNE far detector.  Useage instructions are as follows:

# Use analysis module to create output ROOT tree
The analysis module [Pi0EventDisplay_module.cc](https://github.com/lkashur/Pi0EventDisplay/blob/main/Pi0EventDisplay_module.cc) is used on the dunegpvms and creates an output ROOT file Pi0EventDisplay.root.  Credit to Jacob Larkin for developing the pi0 event/shower selection framework seen within the module.

The module is ran using the FHiCL files [runPi0EventDisplay.fcl](https://github.com/lkashur/Pi0EventDisplay/blob/main/runPi0EventDisplay.fcl) and [Pi0EventDisplay.fcl](https://github.com/lkashur/Pi0EventDisplay/blob/main/Pi0EventDisplay.fcl).

```
lar -c runPi0EventDisplay.fcl -S <text file of input files>
```

I use dunetpc v09_18_00 and my working area can be found at /dune/app/users/lkashur/pi0pi0StudiesSpring2021/srcs/dunetpc/dune/pi0Analysis.


# View specific event using Python scripts
To view a specific event, the event number and z-coordinate of the true interaction vertex is needed.  The script [getEventList.py](https://github.com/lkashur/Pi0EventDisplay/blob/main/getEventList.py) will create a text file containing event numbers and z-coordinates given the output root file from the previous step.

```
python3 getEventList.py Pi0EventDisplay.root
```

Running EventStigator.py will then match true photons to reconstructed showers and generate the plot for an event selected from the text file.  The scripts Pi0.py and EventPlotting.py will also need to be present in the working directory.

```
python3 EventStigator.py -id <event number:z-vertex> Pi0EventDisplay.root
```

You will then be prompted to press "Enter" to exit the display, or press "s" to save the display (and then exit).


For any questions, please email me at lkashur@colostate.edu or message me on the DUNE Slack workspace.
