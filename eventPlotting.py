import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import itertools as it
import math
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

def make_label(thisLabel, i=0):
    if thisLabel is None:
        return None
    elif i == 0:
        return thisLabel
    else:
        return thisLabel + ' ' + str(i)

# sim::IDEs from true photons
def plot_trueShowers(ax, showers):
    
    shower_colors = it.cycle('rc')
    for i, shower in enumerate(showers):
        thisColor = next(shower_colors)
        ax[0].scatter(
            shower.xs,
            shower.ys,
            c=thisColor, 
            s=1,
            label=f"True Photon {i+1} IDEs that were not reconstructed" # this label is used because these points will be plotted over (any remaining points are not reconstucted)
        )
        ax[1].scatter(
            shower.zs,
            shower.xs,
            c=thisColor,
            s=1
        )
        ax[2].scatter(
            shower.zs,
            shower.ys,
            c=thisColor,
            s=1
        )
# recob::SpacePoints from selected recob::Showers
def plot_recoShowers(ax, showers):

    shower_label = 'Selected Reco Shower'
    shower_colors = it.cycle(mpl.colors.TABLEAU_COLORS)

    for i, shower in enumerate(showers):
        thisColor = next(shower_colors)
        ax[0].scatter(
            shower.xs,
            shower.ys,
            c=thisColor, #colors
            s=1, #sizes
            label = shower_label + ' ' + str(i+1) + ' SpacePoints'
        )
        ax[1].scatter(
            shower.zs,
            shower.xs,
            c=thisColor,
            s=1
        )
        ax[2].scatter(
            shower.zs,
            shower.ys,
            c=thisColor,
            s=1
        )

# Plot photons, vertices
def plot_axis(ax, sh_type, pion, want_label=True):

    photon_colors = it.cycle('gm')

    if want_label:
        photon_label = 'Photon'
        pion_label = 'pi0 Decay Vertex'
        vertex_label = 'Neutrino Interaction Vertex'
    else:
        photon_label = None
        pion_label = None
        vertex_label = None

    showers = pion.parse_shower_type(sh_type)
    if showers is None:
        return None

        
    #Plot interaction vertex in all plane views
    ax[0].scatter(
        pion.vertex_X,
        pion.vertex_Y,
        s=200,
        marker='*',
        c='r',
        label = vertex_label
    )
    ax[1].scatter(
        pion.vertex_Z,
        pion.vertex_X,
        s=200,
        marker='*',
        c='r'
    )
    ax[2].scatter(
        pion.vertex_Z,
        pion.vertex_Y,
        s=200,
        marker='*',
        c='r'
    )
           
    # Plot the pion end point in all views, which serves as the vertex for the photons
    ax[0].scatter(
        pion.pion_xs,
        pion.pion_ys,
        s=70,
        marker='*',
        c='k',
        label=pion_label
    )
    ax[1].scatter(
        pion.pion_zs,
        pion.pion_xs,
        s=70,
        marker='*',
        c='k'
    )
    ax[2].scatter(
        pion.pion_zs,
        pion.pion_ys,
        s=70,
        marker='*',
        c='k'
    )
    
    # Plot the photons
    for i, photon in enumerate(pion.photons):
        thisColor = next(photon_colors)
        ax[0].plot(
            photon.xs,
            photon.ys,
            linestyle="--",
            linewidth=1.5,
            c=thisColor,
            label=make_label(photon_label, i+1)
        )
        ax[1].plot(
            photon.zs,
            photon.xs,
            linestyle="--",
            linewidth=1.5,
            c=thisColor
        )
        ax[2].plot(
            photon.zs,
            photon.ys,
            linestyle="--",
            linewidth=1.5,
            c=thisColor
        )
    
    # Plot the showers
    if sh_type == 'true':
        plot_trueShowers(ax, showers)

    else:
        plot_recoShowers(ax, showers)

    ax[0].set_xlabel('X [cm]', labelpad=15, fontsize = 16)
    ax[0].set_ylabel('Y [cm]', labelpad=15, fontsize = 16)
    ax[1].set_xlabel('Z [cm]', labelpad=15, fontsize = 16)
    ax[1].set_ylabel('X [cm]', labelpad=15, fontsize = 16)
    ax[2].set_xlabel('Z [cm]', labelpad=15, fontsize = 16)
    ax[2].set_ylabel('Y [cm]', labelpad=15, fontsize = 16)
    
'''
def plot_planes(pion, save_fig=False):
    fig, (ax0, ax1, ax2) = plt.subplots(1,3)
    ax0.scatter(pion.wire0_all_hits, pion.time0_all_hits, marker=',', c='gold', s=0.75)
    ax1.scatter(pion.wire1_all_hits, pion.time1_all_hits, marker=',', c='gold', s=0.75)
    ax2.scatter(pion.wire2_all_hits, pion.time2_all_hits, marker=',', c='gold', s=0.75)

    for i,shower in enumerate(pion.nosel_hit_showers):
        ax0.scatter(shower.w0, shower.t0, marker=',', c='lightpink', s=0.75)
        ax1.scatter(shower.w1, shower.t1, marker=',', c='lightpink', s=0.75)
        ax2.scatter(shower.w2, shower.t2, marker=',', c='lightpink', s=0.75)

    
    shower_colors = it.cycle(mpl.colors.TABLEAU_COLORS)
    for i,shower in enumerate(pion.sel_hit_showers):
        thisColor = next(shower_colors)
        ax0.scatter(shower.w0, shower.t0, marker=',', c=thisColor, s=0.75)
        ax1.scatter(shower.w1, shower.t1, marker=',', c=thisColor, s=0.75)
        ax2.scatter(shower.w2, shower.t2, marker=',', c=thisColor, s=0.75)
    
    ax0.set_title('Induction Plane 1')
    ax1.set_title('Induction Plane 2')
    ax2.set_title('Collection Plane')
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel('Wire Number', labelpad = 20, size = 15)
    plt.ylabel('Peak Time', labelpad = 20, size = 15)
    plt.title('pi0 Event Display: recob::Hits', pad = 35, size = 20)
    plt.grid(False)
    plt.show()
'''

def plot_pion(pion, save_fig=False):

    # Set up some style
    mpl.style.use('seaborn-darkgrid')
    fig, ax = plt.subplots(2, 3,
                           figsize=(19, 10),
                           sharex='col',
                           sharey='col',
                           gridspec_kw={
                               'hspace': 0.05,
                               'wspace': 0.25
                           })

    fig.subplots_adjust(left=0.05)
    fig.subplots_adjust(right=0.75) #0.80

    # Plot sim::IDEs from true photons
    plot_axis(ax[0], "true", pion)
    
    # This sets the autoscale to focus on the true photon IDEs
    limitDict = {}
    for i, shower in enumerate(pion.true_showers):
        limitDict['minX{0}'.format(i)] = min(shower.xs)
        limitDict['maxX{0}'.format(i)] = max(shower.xs)
        limitDict['minY{0}'.format(i)] = min(shower.ys)
        limitDict['maxY{0}'.format(i)] = max(shower.ys)
        limitDict['minZ{0}'.format(i)] = min(shower.zs)
        limitDict['maxZ{0}'.format(i)] = max(shower.zs)

    minX = min(limitDict['minX0'], limitDict['minX1'], min(pion.pion_xs), min(pion.vertex_X))
    maxX = max(limitDict['maxX0'], limitDict['maxX1'], max(pion.pion_xs), max(pion.vertex_X))
    minY = min(limitDict['minY0'], limitDict['minY1'], min(pion.pion_ys), min(pion.vertex_Y))
    maxY = max(limitDict['maxY0'], limitDict['maxY1'], max(pion.pion_ys), max(pion.vertex_Y))
    minZ = min(limitDict['minZ0'], limitDict['minZ1'], min(pion.pion_zs), min(pion.vertex_Z))
    maxZ = max(limitDict['maxZ0'], limitDict['maxZ1'], max(pion.pion_zs), max(pion.vertex_Z))

    limitList = []
    limitList.append(minX)
    limitList.append(maxX)
    limitList.append(minY)
    limitList.append(maxY)
    limitList.append(minZ)
    limitList.append(maxZ)    

    
    # Plot true photon IDEs that come from any recob::Hit in event
    trur2_colors = iter(['lime','hotpink']) 
    for i, shower in enumerate(pion.trur2_showers):
        thisColor = next(trur2_colors)
        ax[0][0].scatter(
            shower.xs,
            shower.ys,
            c=thisColor, 
            s=1,
            label=f'True Photon {i +1} IDEs from Non-Shower Hits'
        )
        ax[0][1].scatter(
            shower.zs,
            shower.xs,
            c=thisColor,
            s=1
        )
        ax[0][2].scatter(
            shower.zs,
            shower.ys,
            c=thisColor,
            s=1
        )
    

    # Plot all recob::SpacePoints (some will be covered by shower SpacePoints)
    for i, shower in enumerate(pion.reco2_showers):
        ax[1][0].scatter(
            shower.xs,
            shower.ys,
            c='gold', 
            s=1,
            label=f"Non-Shower SpacePoints"
        )
        ax[1][1].scatter(
            shower.zs,
            shower.xs,
            c='gold',
            s=1
        )
        ax[1][2].scatter(
            shower.zs,
            shower.ys,
            c='gold',
            s=1
        )
    
    # Plot the trur1 showers (true photon sim::IDEs originating from non-selected recob::Showers in event
    trur1_colors = iter(['orange','mediumblue']) 
    for i, shower in enumerate(pion.trur1_showers):
        thisColor = next(trur1_colors)
        ax[0][0].scatter(
            shower.xs,
            shower.ys,
            c=thisColor,
            s=1,
            label=f"True Photon {i+1} IDEs from Non-Selected Shower Hits"
        )
        ax[0][1].scatter(
            shower.zs,
            shower.xs,
            c=thisColor,
            s=1
        )
        ax[0][2].scatter(
            shower.zs,
            shower.ys,
            c=thisColor,
            s=1
        )
    
    
    # Plot the trur showers (true photon sim::IDEs originating from selected recob::Showers)
    colors = it.cycle('gm')
    for i, shower in enumerate(pion.trur_showers):
        thisColor = next(colors)
        ax[0][0].scatter(
            shower.xs,
            shower.ys,
            c=thisColor,
            s=1,
            label=f"True Photon {i+1} IDEs from Selected Shower Hits"
        )
        ax[0][1].scatter(
            shower.zs,
            shower.xs,
            c=thisColor,
            s=1
        )
        ax[0][2].scatter(
            shower.zs,
            shower.ys,
            c=thisColor,
            s=1
        )
    
    
    # Plot recob::SpacePoints from non-selected recob::Showers
    for i, shower in enumerate(pion.reco1_showers):
        #thisColor = next(colors)
        ax[1][0].scatter(
            shower.xs,
            shower.ys,
            s = 1,
            c = 'lightpink',
            label = f"Non-Selected Reco Shower SpacePoints"
        )
        ax[1][1].scatter(
            shower.zs,
            shower.xs,
            c = 'lightpink',
            s=1
        )
        ax[1][2].scatter(
            shower.zs,
            shower.ys,
            c = 'lightpink',
            s=1
        )


    # Plot selected recob::SpacePoints from recob::Showers
    plot_axis(ax[1], "reco", pion, want_label=False)

    ax[0][0].set_title("YX-Plane View", fontsize = 16)
    ax[0][1].set_title("XZ-Plane View", fontsize = 16)
    ax[0][2].set_title("YZ-Plane View", fontsize = 16)

    ax[1][0].set_xlim((minX - 10, maxX + 10))
    ax[1][0].set_ylim((minY - 10, maxY + 10))
    ax[1][1].set_xlim((minZ - 10, maxZ + 10))
    ax[1][1].set_ylim((minX - 10, maxX + 10))
    ax[1][2].set_xlim((minZ - 10, maxZ + 10))
    ax[1][2].set_ylim((minY - 10, maxY + 10))
    
    if pion.true_cosa is None:
        tru_angle = "N/A"
    else:
        tru_angle = round(math.acos(pion.true_cosa) * 180/math.pi, 5)

    if pion.reco_cosa is None:
        reco_angle = "N/A"
    else:
        reco_angle = round(math.acos(pion.reco_cosa) * 180/math.pi, 5)
    if pion.muon_cosa is None:
        muon_angle = "N/A"
    else:
        muon_angle = round(math.acos(pion.muon_cosa) * 180/math.pi, 5)

    vertex = round(pion.vertex_Z[0])
    
    props = dict(boxstyle='round', alpha=0.25, facecolor = 'none', edgecolor = 'blue')
    fig.text(
        0.76, #0.005
        0.18, #0.88
        f"True Opening Angle: {round(tru_angle, 1)} [deg]\n"
        f"Reco Opening Angle: {round(reco_angle, 1)} [deg]\n"
        f"Muon Opening Angle: {round(muon_angle, 1)} [deg]\n"
        f"True Photon Energies: {' and '.join(str(round(1000 * e, 1)) for e in pion.true_photonE)} [MeV]\n"
        f"Reco Shower Energies: {' and '.join(str(round(e, 1)) for e in pion.reco_photonE)} [MeV]\n"
        f"True pi0 Mass: {round(pion.true_pion_mass, 1)} [MeV]\n"
        f"Reco pi0 Mass: {round(pion.reco_pion_mass, 1)} [MeV]\n"
        f"Muon pi0 Mass: {round(pion.muon_pion_mass, 1)} [MeV]",
        transform=fig.transFigure,
        fontsize=12,
        bbox=props
    )

    fig.text(0.76,
             0.82,
             f"Event Number: {pion.evNum}\n"
             f"Vertex-Z: {vertex} [cm]",
             fontsize=12
    )

    fig.suptitle(
        f"pi0 Event Display",
        fontsize=15
    )
    
    legend_elements = [Line2D([0], [0], marker='*', color='w', label='Neutrino Interaction Vertex', markerfacecolor='r', markersize=15),
                       Line2D([0], [0], marker='*', color='w', label='pi0 Decay Vertex', markerfacecolor='black', markersize=15),
                       Line2D([0], [0], linestyle='--', color='g', label='Photon 1', markersize=15),
                       Line2D([0], [0], linestyle='--', color='m', label='Photon 2', markersize=15),
                       Line2D([0], [0], color='w', label=' ', markersize=15),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 1 IDEs from Selected Shower Hits', markerfacecolor='g', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 2 IDEs from Selected Shower Hits', markerfacecolor='m', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 1 IDEs from Non-Selected Shower Hits', markerfacecolor='orange', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 2 IDEs from Non-Selected Shower Hits', markerfacecolor='mediumblue', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 1 IDEs from Non-Shower Hits', markerfacecolor='lime', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 2 IDEs from Non-Shower Hits', markerfacecolor='hotpink', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 1 IDEs that were not reconstructed', markerfacecolor='r', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='True Photon 2 IDEs that were not reconstructed', markerfacecolor='c', markersize=7),
                       Line2D([0], [0], color='w', label=' ', markersize=15),
                       Line2D([0], [0], marker='o', color='w', label='Selected Reco Shower 1 SpacePoints', markerfacecolor='tab:blue', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='Selected Reco Shower 2 SpacePoints', markerfacecolor='tab:orange', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='Non-Selected Reco Shower SpacePoints', markerfacecolor='lightpink', markersize=7),
                       Line2D([0], [0], marker='o', color='w', label='Non-Shower SpacePoints', markerfacecolor='gold', markersize=7)
    ]

    fig.legend(handles=legend_elements, loc={0.75, 0.35}, fontsize = 12)
    
    if save_fig:
        fig.savefig(
            f"pion_ev_{pion.evNum}_vtx_{vertex}"
        )
    else:
        fig.show()

