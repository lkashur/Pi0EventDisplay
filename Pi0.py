import math
import itertools as it
import operator as op


# Rounding Function
def round_to_N(n):
    def rounder(x):
        return round(x, n)
    return rounder


# CLASSES FOR Pi0 AND FRIENDS {{
#


class Photon:
    def __init__(self):
        self.xs = []
        self.ys = []
        self.zs = []


class Shower:
    def __init__(self):
        self.xs = []
        self.ys = []
        self.zs = []
        self.sc = []

class HitShower:
    def __init__(self):
        self.w0 = []
        self.t0 = []
        self.w1 = []
        self.t1 = []
        self.w2 = []
        self.t2 = []

# Class which holds all information about the Pi0 and its daughters
class Pi0:

    def __init__(self, phID, n):

        self.pion_xs = []
        self.pion_ys = []
        self.pion_zs = []
        self.true_pion_mass = None
        self.reco_pion_mass = 'N/A'
        self.muon_pion_mass = 'N/A'

        '''
        self.wire0_all_hits = 'N/A'
        self.time0_all_hits = 'N/A'
        self.wire1_all_hits = 'N/A'
        self.time1_all_hits = 'N/A'
        self.wire2_all_hits = 'N/A'
        self.time2_all_hits = 'N/A'

        self.wire0_shower_hits = 'N/A'
        self.time0_shower_hits = 'N/A'
        self.wire1_shower_hits = 'N/A'
        self.time1_shower_hits = 'N/A'
        self.wire2_shower_hits = 'N/A'
        self.time2_shower_hits = 'N/A'

        self.sel_hit_showers = []
        self.nosel_hit_showers = []
        '''
        
        self.vertex_X = []
        self.vertex_Y = []
        self.vertex_Z = []

        self.photons = []
        self.true_photonE = []
        self.reco_photonE = []
        self.true_cosa = None
        self.reco_cosa = None
        self.muon_cosa = None

        self.true_showers = []
        self.reco_showers = []
        self.reco1_showers = []
        self.reco2_showers = []
        self.true2_showers = []
        self.trur_showers = []
        self.trur1_showers = []
        self.trur2_showers = []
        
        self.photonID = phID
        self.evNum = n

    def __mag(self, vec):
        return math.sqrt(
            sum(
                map(lambda x: x*x, vec)))

    def __cosEuclidianAngle(self, vec1, vec2):
        dotProd = sum(
            it.starmap(op.mul,
                       zip(vec1, vec2)))
        return dotProd / (self.__mag(vec1) * self.__mag(vec2))

    def add_photons(self, pointVec1, pointVec2):

        # Arrange the vectors for the photons into loopable data
        start_xs = [pointVec1[0][0], pointVec2[0][0]]
        start_ys = [pointVec1[0][1], pointVec2[0][1]]
        start_zs = [pointVec1[0][2], pointVec2[0][2]]
        end_xs = [pointVec1[1][0], pointVec2[1][0]]
        end_ys = [pointVec1[1][1], pointVec2[1][1]]
        end_zs = [pointVec1[1][2], pointVec2[1][2]]

        photonVec1 = list(
            it.starmap(op.sub,
                       zip(pointVec1[0],
                           pointVec1[1])))
        photonVec2 = list(
            it.starmap(op.sub,
                       zip(pointVec2[0],
                           pointVec2[1])))
        self.true_cosa = self.__cosEuclidianAngle(photonVec1, photonVec2)

        for i in range(2):
            thisPhoton = Photon()
            thisPhoton.xs.append(start_xs[i])
            thisPhoton.ys.append(start_ys[i])
            thisPhoton.zs.append(start_zs[i])
            thisPhoton.xs.append(end_xs[i])
            thisPhoton.ys.append(end_ys[i])
            thisPhoton.zs.append(end_zs[i])
            self.photons.append(thisPhoton)

    def add_reco_angle_and_mass(self, vec1, vec2, recoE):
        round5 = round_to_N(5)
        for rE in recoE:
            self.reco_photonE.append(round5(rE))
        self.reco_cosa = self.__cosEuclidianAngle(vec1, vec2)
        self.reco_pion_mass = round5(self.__pionMass(recoE[0],
                                                     recoE[1],
                                                     self.reco_cosa))
    '''
    def add_wires_and_time_all_hits(self, wire0, peaktime0, wire1, peaktime1, wire2, peaktime2):
        self.wire0_all_hits = wire0[0]
        self.time0_all_hits = peaktime0[0]
        self.wire1_all_hits = wire1[0]
        self.time1_all_hits = peaktime1[0]
        self.wire2_all_hits = wire2[0]
        self.time2_all_hits = peaktime2[0]
        
    def add_wires_and_time_sel_shower_hits(self, wire0, peaktime0, wire1, peaktime1, wire2, peaktime2):
        thisHitShower = HitShower()
        for a,b,c,d,e,f in zip(wire0, peaktime0, wire1, peaktime1, wire2, peaktime2):
            thisHitShower.w0.append(a)
            thisHitShower.t0.append(b)
            thisHitShower.w1.append(c)
            thisHitShower.t1.append(d)
            thisHitShower.w2.append(e)
            thisHitShower.t2.append(f)

        self.sel_hit_showers.append(thisHitShower)


    def add_wires_and_time_nosel_shower_hits(self, wire0, peaktime0, wire1, peaktime1, wire2, peaktime2):
        thisHitShower = HitShower()
        for a,b,c,d,e,f in zip(wire0, peaktime0, wire1, peaktime1, wire2, peaktime2):
            thisHitShower.w0.append(a)
            thisHitShower.t0.append(b)
            thisHitShower.w1.append(c)
            thisHitShower.t1.append(d)
            thisHitShower.w2.append(e)
            thisHitShower.t2.append(f)

        self.nosel_hit_showers.append(thisHitShower)

    '''
    
    def add_muon_angle_and_mass(self, vec1, vec2, recoE):
        round5 = round_to_N(5)
        for rE in recoE:
            #self.reco_photonE.append(round5(rE))
            self.muon_cosa = self.__cosEuclidianAngle(vec1, vec2)
            self.muon_pion_mass = round5(self.__pionMass(recoE[0], recoE[1], self.muon_cosa))

    def __pionMass(self, E1, E2, cosa):
        return math.sqrt(2 * E1 * E2 * (1-cosa))

    def add_pion(self, pion_x, pion_y, pion_z, trueE):
        self.pion_xs.append(pion_x[0])
        self.pion_ys.append(pion_y[0])
        self.pion_zs.append(pion_z[0])
        round5 = round_to_N(5)

        for tE in trueE:
            self.true_photonE.append(round5(tE))

        # Multiply by 1000 to convert GeV -> MeV
        self.true_pion_mass = round5(1000 * self.__pionMass(trueE[0],
                                                            trueE[1],
                                                            self.true_cosa))
            
    def add_vertex(self, vertex_x, vertex_y, vertex_z):
        self.vertex_X.append(vertex_x)
        self.vertex_Y.append(vertex_y)
        self.vertex_Z.append(vertex_z)
        round5 = round_to_N(5)
        


    def parse_shower_type(self, sh_type):
        if sh_type == "true":
            return self.true_showers
        elif sh_type == "reco":
            return self.reco_showers
        elif sh_type == "reco1":
            return self.reco1_showers
        elif sh_type == "reco2":
            return self.reco2_showers
        elif sh_type == "true2":
            return self.true2_showers
        elif sh_type == "trur":
            return self.trur_showers
        elif sh_type == "trur1":
            return self.trur1_showers
        elif sh_type == "trur2":
            return self.trur2_showers
        else:
            print(
                "In add_shower please include 'true', 'reco', or 'trur'"
                " depending on which shower you want to add"
            )
            return None

    def add_shower(self, X, Y, Z, sh_type, isSC=it.repeat(False)):
        if (showers := self.parse_shower_type(sh_type)) is None:
            return None

        thisShower = Shower()
        for x, y, z, sc in zip(X, Y, Z, isSC):
            
            thisShower.xs.append(x)
            thisShower.ys.append(y)
            thisShower.zs.append(z)
            thisShower.sc.append(sc)

        showers.append(thisShower)

#
# }}
