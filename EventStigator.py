# IMPORTS {{
#
import ROOT
from contextlib import contextmanager
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

from Pi0 import Pi0
from Colors import Colors
import eventPlotting as ep
from pprint import pprint

#
# }}


# SESSION FUNCTIONS FOR MAIN {{
#


def print_header():
    print(Colors.CLEAR)
    print("If you wish to move forward one choice,"
          " please enter 'n' or 'next'."
          "\nTo load 10 more choices enter 'c' or 'continue'."
          "\nTo save your last choice, enter 's' or 'save'."
          "\nIf you wish to quit, please enter 'q' or 'quit'\n\n")
    print(f"{Colors.BOLD}{Colors.RED}{'CHOICE':<18}"
          f"{Colors.MAGENTA}{'EVENT NUM':<16}"
          f"{Colors.YELLOW}{'RECO MASS':<16}"
          f"{Colors.GREEN}{'TRUE MASS':<16}"
          f"{Colors.ENDC}")
    print(70 * '=')


def print_list(pion_list):
    for i, pion in enumerate(pion_list):
        print(
            f"{Colors.RED}[+]    {i:<12}"
            f"{Colors.MAGENTA}{pion.evNum:<16}"
            f"{Colors.YELLOW}{pion.reco_pion_mass:<16}"
            f"{Colors.GREEN}{pion.true_pion_mass:<12}"
            f"{Colors.ENDC}"
        )
        print(70 * '-')


def run_interactive_session(pion_list):

    last_pion = None
    print_list(pion_list)

    while True:
        choice = input(
            f"\n{Colors.BLUE}Choose a photon ID to plot, "
            f"or choose to continue or quit: {Colors.ENDC}"
        )
        if choice in ['q', 'quit']:
            exit(0)
        elif choice in ['n', 'next']:
            pion_list.pop(0)
            return None
        elif choice in ['c', 'continue']:
            del pion_list[:]
            return None
        elif choice in ['s', 'save']:
            if last_pion is None:
                print(
                    f"{Colors.RED}"
                    f"Cannot save pion, you haven't chosen one yet!"
                    f"{Colors.ENDC}"
                )
            else:
                print(
                    f"{Colors.BOLD}{Colors.BLUE}"
                    f"Plotting pion event {last_pion.evNum}..."
                    f"{Colors.ENDC}"
                )
                ep.plot_pion(last_pion, save_fig=True)
                print_header()
        elif choice not in list(map(str, range(10))):
            print(f"{Colors.BOLD}Please choose a valid option...{Colors.ENDC}")
            # print_header()
        else:
            ep.plot_pion(last_pion := pion_list[int(choice)])
            print_header()
            print_list(pion_list)

#################################################
##### Adds info from selected pi0 event to Pi0.py
#################################################
def create_pion(event, dic):

    ## Grab relevant variables

    # sim::IDEs from true photons
    true_showerX = [dic[0]['true_showerX'], dic[1]['true_showerX']]
    true_showerY = [dic[0]['true_showerY'], dic[1]['true_showerY']]
    true_showerZ = [dic[0]['true_showerZ'], dic[1]['true_showerZ']]
    
    # simIDEs from all recob::Showers in event (found with backtracker)
    trur_showerX = event.ShowerSimIDEX
    trur_showerXs = [x for x in trur_showerX]
    trur_showerY = event.ShowerSimIDEY
    trur_showerYs = [y for y in trur_showerY]
    trur_showerZ = event.ShowerSimIDEZ
    trur_showerZs = [z for z in trur_showerZ]
            
    # all SpacePoints from event
    reco_hitsX = event.EventSpacePointX
    reco_hitsY = event.EventSpacePointY
    reco_hitsZ = event.EventSpacePointZ
    
    # all sim::IDEs in event (found from backtracking all recob::Hits in event)
    trueIDEX = event.TrueIDEX
    trueIDEY = event.TrueIDEY
    trueIDEZ = event.TrueIDEZ

    # SpacePoints from all recob::Showers in event
    reco_showerX = event.ShowerSpacePointX
    reco_showerXs = [x for x in reco_showerX]
    reco_showerY = event.ShowerSpacePointY
    reco_showerZ = event.ShowerSpacePointZ
    
    showerIDs = [x for x in event.PFShowerBestMCParentID]

    # The following recob::Hit info is for the pion.plot_planes function (not currently in use)
    '''
    # recob::Hits from entire event
    hit_peaktime0 = event.HitPeakTime0
    hit_wire0 = event.HitWire0
    hit_plane0 = event.HitPlane0
    hit_peaktime1 = event.HitPeakTime1
    hit_wire1 = event.HitWire1
    hit_plane1 = event.HitPlane1
    hit_peaktime2 = event.HitPeakTime2
    hit_wire2 = event.HitWire2
    hit_plane2 = event.HitPlane2

    # recob::Hit info from all recob::Showers in event
    showerHit_wire0 = event.ShowerHitWire0
    showerHit_wire0s = [w for w in showerHit_wire0]
    showerHit_time0 = event.ShowerHitPeakTime0
    showerHit_wire1 = event.ShowerHitWire1
    showerHit_time1 = event.ShowerHitPeakTime1
    showerHit_wire2 = event.ShowerHitWire2
    showerHit_time2 = event.ShowerHitPeakTime2

    # all recob::Hits
    thisPi.add_wires_and_time_all_hits(hit_wire0, hit_peaktime0, hit_wire1, hit_peaktime1, hit_wire2, hit_peaktime2)

    # recob::Hits from selected recob::Showers
    for i in range(len(dic)):
        thisPi.add_wires_and_time_sel_shower_hits(
            dic[i]['showerHit_wire0'],
            dic[i]['showerHit_peaktime0'],
            dic[i]['showerHit_wire1'],
            dic[i]['showerHit_peaktime1'],
            dic[i]['showerHit_wire2'],
            dic[i]['showerHit_peaktime2']
        )
     
    # recob::Hits from non-selected recob::Showers
    for w in showerHit_wire0s:
        i = showerHit_wire0s.index(w)
        if (showerHit_wire0[i] != dic[0]['showerHit_wire0']) & (showerHit_wire0[i] != dic[1]['showerHit_wire0']):
            thisPi.add_wires_and_time_nosel_shower_hits(
                showerHit_wire0[i],
                showerHit_time0[i],
                showerHit_wire1[i],
                showerHit_time1[i],
                showerHit_wire2[i],
                showerHit_time2[i]
            )
    '''


    # to be filled later
    direction_vectors = []
    muon_vectors = []
    recoShower_energies = []    
    
    # Initialize event
    thisPi = Pi0(dic[0]['ID'], dic[0]['Event'])

    # Loop through photons/selected recob::Showers
    for i in range(len(dic)):

        # Add sim::IDEs from true photons
        thisPi.add_shower(
            dic[i]['true_showerX'],
            dic[i]['true_showerY'],
            dic[i]['true_showerZ'],
            "true"
        )

        # Add sim::IDEs from selected recob::Showers (backtracked)
        thisPi.add_shower(
            [x for x in dic[i]['trur_showerX'] if x in true_showerX[i]], 
            [y for y in dic[i]['trur_showerY'] if y in true_showerY[i]], 
            [z for z in dic[i]['trur_showerZ'] if z in true_showerZ[i]], 
            'trur'
        )

        # Add recob::SpacePoints from selected recob:Showers
        thisPi.add_shower(
            dic[i]['reco_showerX'],
            dic[i]['reco_showerY'],
            dic[i]['reco_showerZ'],
            "reco"
        )        

        # Add info about direction vectors, energy, reco mass, etc.
        direction_vectors.append([
            dic[i]['showerDirX'],
            dic[i]['showerDirY'],
            dic[i]['showerDirZ']
        ])
        muon_vectors.append([
            dic[i]['muonDirX'],
            dic[i]['muonDirY'],
            dic[i]['muonDirZ']
        ])
        recoShower_energies.append(
            1000*dic[i]['showerEnergy']
        )

    thisPi.add_reco_angle_and_mass(
        *direction_vectors, recoShower_energies)
    thisPi.add_muon_angle_and_mass(
        *muon_vectors, recoShower_energies)

    photon1_vec = [
        [dic[0]['photonStartX'], dic[0]['photonStartY'], dic[0]['photonStartZ']],
        [dic[0]['photonEndX'], dic[0]['photonEndY'], dic[0]['photonEndZ']]
    ]

    photon2_vec = [
        [dic[1]['photonStartX'], dic[1]['photonStartY'], dic[1]['photonStartZ']],
        [dic[1]['photonEndX'], dic[1]['photonEndY'], dic[1]['photonEndZ']]
    ]

    thisPi.add_photons(photon1_vec, photon2_vec)
    
    thisPi.add_pion(
        dic[0]['pion_endX'],
        dic[0]['pion_endY'],
        dic[0]['pion_endZ'],
        [dic[0]['photonE'], dic[1]['photonE']]
    )
        
    thisPi.add_vertex(dic[0]['vtxX'], dic[0]['vtxY'], dic[0]['vtxZ'])
    
    # Add SpacePoints from non-selected recob::Showers
    for x in reco_showerXs:
        i = reco_showerXs.index(x)
        # if recob::Shower from event doesn't match with recob::Shower from dictionary, categorize as reco1
        if (reco_showerX[i] != dic[0]['reco_showerX']) & (reco_showerX[i] != dic[1]['reco_showerX']):
            thisPi.add_shower(
                reco_showerX[i],
                reco_showerY[i],
                reco_showerZ[i],
                "reco1"
            )
            
    # sim::IDEs that were backtracked from all recob::Hits in event (but now in list form)
    trueIDEXlist = []
    trueIDEYlist = []
    trueIDEZlist = []
    for x in trueIDEX:
        for i in x:
            trueIDEXlist.append(i)
    for y in trueIDEY:
        for i in y:
            trueIDEYlist.append(i)
    for z in trueIDEZ:
        for i in z:
            trueIDEZlist.append(i)

    # Grab showers that were not selected
    trurXlist_nosel = []
    trurYlist_nosel = []
    trurZlist_nosel = []

    for i in range(len(trur_showerX)):
        if (trur_showerX[i] != dic[0]['trur_showerX']) & (trur_showerX[i] != dic[1]['trur_showerX']):
            for x in trur_showerX[i]:
                trurXlist_nosel.append(x)
    for i in range(len(trur_showerY)):
        if (trur_showerY[i] != dic[0]['trur_showerY']) & (trur_showerY[i] != dic[1]['trur_showerY']):
            for y in trur_showerY[i]:
                trurYlist_nosel.append(y)
    for i in range(len(trur_showerZ)):
        if (trur_showerZ[i] != dic[0]['trur_showerZ']) & (trur_showerZ[i] != dic[1]['trur_showerZ']):
            for z in trur_showerZ[i]:
                trurZlist_nosel.append(z)
                
    # Grab sim::IDEs that match to non-selected recob::Showers and other recob::Hits in event
    for i in range(len(dic)):
        trurXlist_ns = []
        trurYlist_ns = []
        trurZlist_ns = []
        trurXlist_notSh = []
        trurYlist_notSh = []
        trurZlist_notSh = []
        trurXlist_notSh1 = []
        trurYlist_notSh1 = []
        trurZlist_notSh1 = []
        # Grab sim::IDEs that match to non-selected recob::Showers and other recob::Hits in event
        for x in true_showerX[i]:
            # if true photon sim::IDE matches to sim::IDE from non-selected recob::Shower, categorize as 'trur1'
            if x in trurXlist_nosel: 
                trurXlist_ns.append(x)
        for y in true_showerY[i]:
            if y in trurYlist_nosel:
                trurYlist_ns.append(y)
        for z in true_showerZ[i]:
            if z in trurZlist_nosel:
                trurZlist_ns.append(z)
            
        trurXlist_notSh = [x for x in trueIDEXlist]
        trurYlist_notSh = [y for y in trueIDEYlist]
        trurZlist_notSh = [z for z in trueIDEZlist]

        
        trueTups = list(zip(trueIDEXlist, trueIDEYlist, trueIDEZlist))
        trueShTups = list(zip(true_showerX[i], true_showerY[i], true_showerZ[i]))
        trurShTups = list(zip(dic[i]['trur_showerX'], dic[i]['trur_showerY'], dic[i]['trur_showerY']))

        for tup in trueTups:
            if tup in trueShTups and tup not in trurShTups:
                trurXlist_notSh1.append(tup[0])
                trurYlist_notSh1.append(tup[1])
                trurZlist_notSh1.append(tup[2])
        
        '''
        for x,y,z in zip(trurXlist_notSh, trurYlist_notSh, trurZlist_notSh):
            if x in true_showerX[i] and y in true_showerY[i] and z in true_showerZ[i] and x not in dic[i]['trur_showerX'] and y not in dic[i]['trur_showerY'] and z not in dic[i]['trur_showerZ']:
                trurXlist_notSh1.append(x)
                trurYlist_notSh1.append(y)
                trurZlist_notSh1.append(z)
        '''
        
        thisPi.add_shower(
            trurXlist_ns,
            trurYlist_ns,
            trurZlist_ns,
            'trur1'
        )
        thisPi.add_shower(
            trurXlist_notSh1, 
            trurYlist_notSh1,
            trurZlist_notSh1,
            'trur2'
        )

    # Grab SpacePoints from non-shower objects
    thisPi.add_shower(
        reco_hitsX,
        reco_hitsY,
        reco_hitsZ,
        "reco2"
    )
            
    return thisPi


@contextmanager
def open_rootFile(fileName):
    try:
        print(
            f"{Colors.BOLD}{Colors.RED}---\n  {Colors.BLUE}->{Colors.GREEN}"
            f"Loading root file {fileName} ..."
            f"{Colors.RED} \n---\n {Colors.ENDC}"
        )
        rootFile = ROOT.TFile.Open(fileName)
        evTree = rootFile.Get("trkUtil/PiZeroShower")
        yield evTree
    finally:
        rootFile.Close()

def search_for_pion(eventTree, a_evNum, a_vertex):

    next = 0
    evNum_multiple = a_evNum - 1

    for i, event in enumerate(eventTree):
        print(f"I have searched {i} events so far...", end="\r")

        # event info
        Event = event.event
        vtxX = event.vtx_x
        vtxY = event.vtx_y
        vtxZ = event.vtx_z

        # photon info
        photonID = event.PhotonID
        photonStartX = event.PhotonStartX
        photonStartY = event.PhotonStartY
        photonStartZ = event.PhotonStartZ
        photonEndX = event.PhotonEndX
        photonEndY = event.PhotonEndY
        photonEndZ = event.PhotonEndZ
        photonE = event.PhotonE
        photonIndex = 0

        # shower info
        showerID = event.PFShowerBestMCParentID
        showerWeight = event.PFShowerBestMCParentWeight
        showerEnergyCollection = event.PFShowerEnergyCollection
        showerDirectionX = event.PFShowerDirectionX
        showerDirectionY = event.PFShowerDirectionY
        showerDirectionZ = event.PFShowerDirectionZ
        showerStartX = event.PFShowerStartX
        showerStartY = event.PFShowerStartY
        showerStartZ = event.PFShowerStartZ
        showerIndex = 0

        # spatial info
        true_showerX = event.TruePhotonSimIDEX
        true_showerY = event.TruePhotonSimIDEY 
        true_showerZ = event.TruePhotonSimIDEZ 
        reco_showerX = event.ShowerSpacePointX
        reco_showerY = event.ShowerSpacePointY
        reco_showerZ = event.ShowerSpacePointZ
        trur_showerX = event.ShowerSimIDEX
        trur_showerY = event.ShowerSimIDEY
        trur_showerZ = event.ShowerSimIDEZ
        
        #while next < 100:
        #    next = next + 1
        #    continue

        #if (i % 100) != (evNum_multiple):
        #    continue
        
        if not event.isCC:
            continue

        # Check that there is only one pion and then check that it has
        #     2 daughters (both photons) and then check if its the right event
        #     and plot it if it is
        num_pis = event.Pi0CandidateID.size()
        if num_pis == 1:
            num_daughters = event.Pi0CandidateNumPhotonDaughters[0]
            if num_daughters == 2:
                
                vertex = round(event.vtx_z)

                #true photon 1 info
                truePh1ID = photonID[photonIndex]
                photon1E = photonE[photonIndex]
                photon1StartX = photonStartX[photonIndex]
                photon1StartY = photonStartY[photonIndex]
                photon1StartZ = photonStartZ[photonIndex]
                photon1EndX = photonEndX[photonIndex]
                photon1EndY = photonEndY[photonIndex]
                photon1EndZ = photonEndZ[photonIndex]
                photon1DirX = photonEndX[photonIndex] - photonStartX[photonIndex]
                photon1DirY = photonEndY[photonIndex] - photonStartY[photonIndex]
                photon1DirZ = photonEndZ[photonIndex] - photonStartZ[photonIndex]
                photon1DirList = [photonEndX[photonIndex] - photonStartX[photonIndex], photonEndY[photonIndex] - photonStartY[photonIndex], photonEndZ[photonIndex] - photonStartZ[photonIndex]]
                true_shower1X = true_showerX[photonIndex]
                true_shower1Y = true_showerY[photonIndex]
                true_shower1Z = true_showerZ[photonIndex]
                photonIndex += 1

                #true photon 2 info
                truePh2ID = photonID[photonIndex]
                photon2E = photonE[photonIndex]
                photon2StartX = photonStartX[photonIndex]
                photon2StartY = photonStartY[photonIndex]
                photon2StartZ = photonStartZ[photonIndex]
                photon2EndX = photonEndX[photonIndex]
                photon2EndY = photonEndY[photonIndex]
                photon2EndZ = photonEndZ[photonIndex]
                photon2DirX = photonEndX[photonIndex] - photonStartX[photonIndex]
                photon2DirY = photonEndY[photonIndex] - photonStartY[photonIndex]
                photon2DirZ = photonEndZ[photonIndex] - photonStartZ[photonIndex]
                photon2DirList = [photonEndX[photonIndex] - photonStartX[photonIndex], photonEndY[photonIndex] - photonStartY[photonIndex], photonEndZ[photonIndex] - photonStartZ[photonIndex]]
                photon2DirVector = np.array(photon2DirList)
                true_shower2X = true_showerX[photonIndex]
                true_shower2Y = true_showerY[photonIndex]
                true_shower2Z = true_showerZ[photonIndex]
                photonIndex += 1

                #reco shower info
                showerIDs = list(showerID)
                if len(showerIDs) < 2:
                    showerIDs = []
                else:
                    showerIDs = showerIDs
                numShowers = len(showerIDs)
                showerWeights = list(showerWeight)
                showerEnergies = list(showerEnergyCollection)
                showrWxEs = [a * b for a,b in zip(showerWeights, showerEnergies)]
                showerDirectionsX = list(showerDirectionX)
                showerDirectionsY = list(showerDirectionY)
                showerDirectionsZ = list(showerDirectionZ)

                # Ensure each true photon has a reco shower match
                if((sum(1 for shID in showerIDs if((shID == truePh1ID) or (shID == truePh2ID))) >= 2)):
                    shList = []
                    for x in range(numShowers):
                        shDict = {}
                        shDict['ID'] = showerID[showerIndex]
                        shDict['showerEnergy'] = showerEnergyCollection[showerIndex]
                        shDict['showerWeight'] = showerWeight[showerIndex]
                        shDict['showerWxE'] = showerWeight[showerIndex] * showerEnergyCollection[showerIndex]
                        shDict['showerDirX'] = showerDirectionX[showerIndex]
                        shDict['showerDirY'] = showerDirectionY[showerIndex]
                        shDict['showerDirZ'] = showerDirectionZ[showerIndex]
                        shDict['muonDirX'] = showerStartX[showerIndex] - vtxX
                        shDict['muonDirY'] = showerStartY[showerIndex] - vtxY
                        shDict['muonDirZ'] = showerStartZ[showerIndex] - vtxZ
                        shDict['reco_showerX'] = reco_showerX[showerIndex]
                        shDict['reco_showerY'] = reco_showerY[showerIndex]
                        shDict['reco_showerZ'] = reco_showerZ[showerIndex]
                        shDict['trur_showerX'] = trur_showerX[showerIndex]
                        shDict['trur_showerY'] = trur_showerY[showerIndex]
                        shDict['trur_showerZ'] = trur_showerZ[showerIndex]
                        shDict['pion_startX'] = event.Pi0CandidateStartX
                        shDict['pion_startY'] = event.Pi0CandidateStartY
                        shDict['pion_startZ'] = event.Pi0CandidateStartZ
                        shDict['pion_endX'] = event.Pi0CandidateEndX
                        shDict['pion_endY'] = event.Pi0CandidateEndY
                        shDict['pion_endZ'] = event.Pi0CandidateEndZ
                        shDict['vtxX'] = vtxX
                        shDict['vtxY'] = vtxY
                        shDict['vtxZ'] = vtxZ
                        shDict['Event'] = Event
                        # add this block back in if using pion.plot_planes
                        '''
                        #shDict['trueShowerPointsX'] = trueShowerPointsX[showerIndex]
                        #shDict['trueShowerPointsY'] = trueShowerPointsY[showerIndex]
                        #shDict['trueShowerPointsZ'] = trueShowerPointsZ[showerIndex]
                        #shDict['showerHit_peaktime0'] = showerHit_peaktime0[showerIndex]
                        #shDict['showerHit_peaktime1'] = showerHit_peaktime1[showerIndex]
                        #shDict['showerHit_peaktime2'] = showerHit_peaktime2[showerIndex]
                        #shDict['showerHit_wire0'] = showerHit_wire0[showerIndex]
                        #shDict['showerHit_wire1'] = showerHit_wire1[showerIndex]
                        #shDict['showerHit_wire2'] = showerHit_wire2[showerIndex]
                        '''
                        shList.append(shDict)
                        showerIndex += 1
                    
                    shList = [shDict for shDict in shList if((shDict['ID'] == truePh1ID) | (shDict['ID'] == truePh2ID))] # get rid of non-matches

                    # I use pandas here to finish off the rest of the matching algorithm...will do this in a cleaner way in the future
                    df1 = pd.DataFrame(shList)
                    df1 = df1.sort_values('showerWxE', ascending=False)
                    df1['dup'] = df1.duplicated(subset = 'ID')
                    df1 = df1[['ID', 'showerEnergy', 'showerWeight', 'showerWxE', 'showerDirX', 'showerDirY', 'showerDirZ', 'muonDirX', 'muonDirY', 'muonDirZ', 'reco_showerX', 'reco_showerY', 'reco_showerZ', 'trur_showerX', 'trur_showerY', 'trur_showerZ', 'pion_startX', 'pion_startY', 'pion_startZ', 'pion_endX', 'pion_endY', 'pion_endZ', 'vtxX', 'vtxY', 'vtxZ', 'Event']]
                    if (len(df1) > 1) & ((df1['showerWeight'] > 0.7).all()):
                        shList = df1.to_dict('records')
                        phList = []
                        ph1Dict = {'ID':truePh1ID, 'photonE':photon1E, 'photonStartX':photon1StartX, 'photonStartY':photon1StartY, 'photonStartZ':photon1StartZ, 'photonEndX':photon1EndX, 'photonEndY':photon1EndY, 'photonEndZ':photon1EndZ, 'photonDirX':photon1DirX, 'photonDirY': photon1DirY, 'photonDirZ': photon1DirZ, 'true_showerX': true_shower1X, 'true_showerY': true_shower1Y, 'true_showerZ': true_shower1Z}
                        phList.append(ph1Dict)
                        ph2Dict = {'ID':truePh2ID, 'photonE':photon2E, 'photonStartX':photon2StartX, 'photonStartY':photon2StartY, 'photonStartZ':photon2StartZ, 'photonEndX':photon2EndX, 'photonEndY':photon2EndY, 'photonEndZ':photon2EndZ, 'photonDirX':photon2DirX, 'photonDirY': photon2DirY, 'photonDirZ': photon2DirZ, 'true_showerX': true_shower2X, 'true_showerY': true_shower2Y, 'true_showerZ': true_shower2Z}
                        phList.append(ph2Dict)
                        d = defaultdict(dict)
                        for l in (phList, shList):
                            for elem in l:
                                d[elem['ID']].update(elem)
                        eventDict = dict(d)
                        df2 = pd.DataFrame(eventDict).T
                        df2.reset_index(inplace = True)
                        eventDict = df2.to_dict('records')
                        if (event.event == a_evNum) and (vertex == a_vertex):
                            return create_pion(event,eventDict)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fileName",
        type=str,
        help="Name of root file to be processed"
    )
    parser.add_argument(
        "-id",
        help="""Unique identifier for an event consisting of event number
                and pion true vertex. Should be called as
                -ev [evNum]:[vertex_z].
                (Default behavior is interactive plotting.)"""
    )

    args = parser.parse_args()

    with open_rootFile(args.fileName) as evTree:
        a_evNum, a_vertex = map(int, args.id.split(':'))
        print(
            f"{Colors.CYAN}"
            f"Searching for event {a_evNum} "
            f"with vertex (int-cast) {a_vertex}"
            f"{Colors.ENDC}\n"
        )
        pion = search_for_pion(evTree, a_evNum, a_vertex)
        if pion is None:
            print(
                f"{Colors.RED}{Colors.BOLD}"
                "Could not find the event given... exiting program"
                f"{Colors.ENDC}\n"
            )
            exit(2)
        ep.plot_pion(pion)
        #ep.plot_planes(pion)
        choice = input(
            f"Plotting event {a_evNum} "
            f"with vertex (int-cast) {a_vertex}\n"
            "Please press 's' to save or "
            "press <Enter> to close..."
        )
        if choice == 's':
            ep.plot_pion(pion, save_fig=True)
    

