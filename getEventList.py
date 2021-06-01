###########
### imports
###########
import sys
import ROOT
from ROOT import TH1F, TH2F, gStyle, gROOT, TStyle
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
import math
from collections import Counter, defaultdict
import pprint
import itertools as it
from itertools import chain
import operator as op
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pprint


##############################################
### Analyze pi0s from simulated DUNE FD events
##############################################
def main(*args):
    ###########################
    ###Load ROOT file and TTree
    ###########################
    filename = args[0]
    rootFile = ROOT.TFile.Open(filename, "READ")
    evTree = rootFile.Get('trkUtil/PiZeroShower')
    
    ##########################################################################################
    ### Create output file for event number and vertex info of specific events to be stored in
    ##########################################################################################
    print('Writing to file pi0Events.txt...')
    outF = open('pi0Events.txt', 'w')
    
    
    shListofLists = []
    phListofLists = []
    dfs = []
    
    ######################
    ### Loop through TTree
    ######################
    for i, event in enumerate(evTree):
        ####################
        ### Define variables
        ####################
        # event info
        Event = event.event
        vtxX = event.vtx_x
        vtxY = event.vtx_y
        vtxZ = event.vtx_z
        isCC = event.isCC
        pi0ID = event.Pi0CandidateID
        pi0StartX = event.Pi0CandidateStartX
        pi0StartY = event.Pi0CandidateStartY
        pi0StartZ = event.Pi0CandidateStartZ
        pi0EndX = event.Pi0CandidateEndX
        pi0EndY = event.Pi0CandidateEndY
        pi0EndZ = event.Pi0CandidateEndZ
        
        # photon info
        numPhotonDaughters = event.Pi0CandidateNumPhotonDaughters
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

        
        #####################
        ### Loop through pi0s
        #####################
        numPi = pi0ID.size()
        for j in range(len(pi0ID)):
            thisPi0ID = pi0ID[j]
            numDaughterPhotons = numPhotonDaughters[j]
            numDaughterShowers = len(showerID)

            if((isCC == 1) and (numPi == 1) and numDaughterPhotons == 2 and (vtxX > -310 and vtxX < 310) and (vtxY > -550 and vtxY < 550) and (vtxZ > 50 and vtxZ < 1244)):
                # True photon 1 info
                truePh1ID = photonID[photonIndex]
                photon1E = photonE[photonIndex]
                photon1DirX = photonEndX[photonIndex] - photonStartX[photonIndex]
                photon1DirY = photonEndY[photonIndex] - photonStartY[photonIndex]
                photon1DirZ = photonEndZ[photonIndex] - photonStartZ[photonIndex]
                photon1DirList = [photonEndX[photonIndex] - photonStartX[photonIndex], photonEndY[photonIndex] - photonStartY[photonIndex], photonEndZ[photonIndex] - photonStartZ[photonIndex]]
                photon1DirVector = np.array(photon1DirList)
                photonIndex += 1

                # True photon 2 info
                truePh2ID = photonID[photonIndex]
                photon2E = photonE[photonIndex]
                photon2DirX = photonEndX[photonIndex] - photonStartX[photonIndex]
                photon2DirY = photonEndY[photonIndex] - photonStartY[photonIndex]
                photon2DirZ = photonEndZ[photonIndex] - photonStartZ[photonIndex]
                photon2DirList = [photonEndX[photonIndex] - photonStartX[photonIndex], photonEndY[photonIndex] - photonStartY[photonIndex], photonEndZ[photonIndex] - photonStartZ[photonIndex]]
                photon2DirVector = np.array(photon2DirList)
                photonIndex += 1
                
                # Reco shower info
                showerIDs = list(showerID)
                if len(showerIDs) < 2: # would not be able to match to two true photons with < 2 showers
                    showerIDs = []
                else:
                    showerIDs = showerIDs
                numShowers = len(showerIDs)
                showerWeights = list(showerWeight)
                showerEnergies = list(showerEnergyCollection)
                showerWxEs = [a * b for a,b in zip(showerWeights, showerEnergies)]
                showerDirectionsX = list(showerDirectionX)
                showerDirectionsY = list(showerDirectionY)
                showerDirectionsZ = list(showerDirectionZ)
                
                ######################################################################################################
                ### Create a list of dictionaries (each dict is a shower) for each shower that matches a true photonID
                ######################################################################################################
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
                        shDict['showerStartX'] = showerStartX[showerIndex]
                        shDict['showerStartY'] = showerStartY[showerIndex]
                        shDict['showerStartZ'] = showerStartZ[showerIndex]
                        shDict['muonDirX'] = showerStartX[showerIndex] - vtxX
                        shDict['muonDirY'] = showerStartY[showerIndex] - vtxY
                        shDict['muonDirZ'] = showerStartZ[showerIndex] - vtxZ
                        #shDict['reco_showerX'] = reco_showerX[showerIndex]
                        #shDict['reco_showerY'] = reco_showerY[showerIndex]
                        #shDict['reco_showerZ'] = reco_showerZ[showerIndex]
                        #shDict['trur_showerX'] = trur_showerX[showerIndex]
                        #shDict['trur_showerY'] = trur_showerY[showerIndex]
                        #shDict['trur_showerZ'] = trur_showerZ[showerIndex]
                        shList.append(shDict)
                        showerIndex += 1

                    shList = [shDict for shDict in shList if((shDict['ID'] == truePh1ID) | (shDict['ID'] == truePh2ID))] # only matches
                    
                    ############################################################################################################################
                    ### Transfer shower data to pandas dataframe in order to choose showers with highest WxE if there are duplicate matching IDs
                    ############################################################################################################################
                    df1 = pd.DataFrame(shList)
                    df1 = df1.sort_values('showerWxE', ascending=False)
                    df1['dup'] = df1.duplicated(subset = 'ID')
                    df1 = df1[df1['dup'] != True]
                    df1 = df1[['ID', 'showerEnergy', 'showerWeight', 'showerWxE', 'showerDirX', 'showerDirY', 'showerDirZ', 'showerStartX', 'showerStartY', 'showerStartZ', 'muonDirX', 'muonDirY', 'muonDirZ']]
                    #print(df1)
                    #print(len(df1))
                    if (len(df1) > 1) & ((df1['showerWeight'] > 0.7).all()): # This is to enforce weight cut and eliminate pi0s with shower IDs like this: [4,4]
                        dfs.append(df1)
                        #print(df1)
                        shList = df1.to_dict('records')
                        shListofLists.append(shList)
                        # Create list of photon dicts in similar manner to showers
                        phList  = []
                        ph1Dict = {'ID':truePh1ID, 'photonE':photon1E, 'photonDirX':photon1DirX, 'photonDirY': photon1DirY, 'photonDirZ': photon1DirZ}
                        phList.append(ph1Dict)
                        ph2Dict = {'ID':truePh2ID, 'photonE':photon2E, 'photonDirX':photon2DirX, 'photonDirY': photon2DirY, 'photonDirZ': photon2DirZ}
                        phList.append(ph2Dict)
                        phListofLists.append(phList)
                        # Combine photon and shower lists of dicts
                        d = defaultdict(dict)
                        for l in (phList, shList):
                            for elem in l:
                                d[elem['ID']].update(elem)
                        comboList = dict(d)
                        df2 = pd.DataFrame(comboList).T
                        df2.reset_index(inplace = True) # As there are only two showers/photons per pi0, index them by 0 and 1
                        comboList = df2.to_dict('records')
                        #pprint.pprint(comboList)


                        ##################################################################################################
                        #### Write event number and z-vertex of relevant events to file (change filename at top of main())
                        ##################################################################################################
                        outF.write('Event Number: ')
                        outF.write(str(Event))
                        outF.write(' | ')
                        outF.write('Vertex-Z: ')
                        outF.write(str(round(vtxZ)))
                        outF.write('\n')
                        
                    
            else:
                for k in range(numDaughterPhotons):
                    photonIndex += 1 # This corresponds to the if statement specifying 1 pi0, 2photons, CC, and fiducial volume cut


    df = pd.concat(dfs) # returned by main...may do something worthwhile with this at some point
    
    ########################################
    #### Close ROOT file and output text file
    ########################################
    rootFile.Close()
    outF.close()
    return df

if __name__ == '__main__':
    df = main(*sys.argv[1:])
