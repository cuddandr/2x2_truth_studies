import numpy as np
import ROOT as RT

def is_point_contained(pos):
    if abs(pos[0]) > 670: return False
    if abs(pos[1] - 430) > 670: return False
    if abs(pos[2]) > 670: return False
    return True

def is_track_contained(trk):
    end_point = trk.Points[-1].GetPosition().Vect()
    return is_point_contained(end_point)

def get_neutron_and_daughter_ids(event):

    neutrons  = set()
    daughters = set()

    for traj in event.Trajectories:

        if traj.GetPDGCode() == 2112:
            neutrons.add(traj.GetTrackId())
            continue
        par_id = traj.GetParentId()
        if par_id in neutrons or par_id in daughters:
            daughters.add(traj.GetTrackId())

    return neutrons.union(daughters)

def get_low_energy_ids(event):
    return set(x.GetTrackId() for x in event.Trajectories if x.GetInitialMomentum().E() < 10)

def get_traj_ids_for_pdg(particles, pdgs):
    return tuple(x.GetTrackId() for x in particles if x.GetPDGCode() in pdgs)

def is_hadronic_contained(event):

    neutron_ids = get_neutron_and_daughter_ids(event)
    low_energy_ids = get_low_energy_ids(event)
    muon_ids = get_traj_ids_for_pdg(event.Primaries[0].Particles, [13, -13])

    for seg in event.SegmentDetectors:

        nChunks = len(seg[1])
        for n in range(nChunks):

            key_contrib = seg[1][n].GetContributors()[0]
            par_contrib = seg[1][n].GetPrimaryId()

            if par_contrib in muon_ids:
                continue

            if key_contrib in neutron_ids:
                continue

            if key_contrib in low_energy_ids:
                continue

            pos = seg[1][n].GetStop()

            if not is_point_contained(pos):
                return False

    return True

def get_nu_vec(genie_tree):

    genie_evt = genie_tree
    for p in range(genie_evt.StdHepN):
        if genie_evt.StdHepStatus[p] != 0:
            continue

        if np.abs(genie_evt.StdHepPdg[p]) not in [12, 14, 16]:
            continue

        nu_vec = RT.TLorentzVector(genie_evt.StdHepP4[p*4 + 0]*1000,
                                   genie_evt.StdHepP4[p*4 + 1]*1000,
                                   genie_evt.StdHepP4[p*4 + 2]*1000,
                                   genie_evt.StdHepP4[p*4 + 3]*1000)
        nu_pdg = genie_evt.StdHepPdg[p]

        return (nu_vec, nu_pdg)

def calc_distance(trk):
    dist = 0
    n_pt = len(trk.Points)
    for i in range(1, n_pt):
        v0 = trk.Points[i-1].GetPosition().Vect()
        v1 = trk.Points[i].GetPosition().Vect()
        dist += (v1 - v0).Mag()

    return dist

def energy_by_range(range):
    r = 0
    T = 5
    dr = range / 100.0;
    while(r < range):
        T += dr * calc_energy_loss_cm(T)
        r += dr

    return T

def calc_energy_loss_cm(T):
    T *= (105.66 / 139.75)
    x = np.log10(T)
    c_lo = np.array([0.363907, -3.99702, 16.8216, -31.8385, 24.2120])
    c_hi = np.array([0.120316, -1.64161, 8.36222, -18.4671, 16.3644])
    density = 1.4
    c = c_lo if x < 3.0 else c_hi
    mev_per_cm = density * (c[0]*x*x*x*x + c[1]*x*x*x + c[2]*x*x + c[3]*x + c[4])
    return mev_per_cm

def energy_deposit_trk(segment_det, trk_id):
    reco_energy = 0

    # print("Edep for {}".format(trk_id))
    for k,v in segment_det:
        for edep in v:
            prim_id = edep.GetPrimaryId()
            contrib = edep.GetContributors()
            # print(prim_id)
            if prim_id != trk_id:
                continue
            # if trk_id not in contrib:
                # continue
            # if prim_id == trk_id or trk_id in contrib:
                # reco_energy += edep.GetEnergyDeposit()
            reco_energy += edep.GetEnergyDeposit()

    # print(reco_energy)
    return reco_energy

def edep_plus_children(event, trk_id):
    reco_energy = 0

    # list_particles = [trk_id]
    # list_particles = []
    for trk in event.Trajectories:
        if trk.GetParentId() == trk_id or trk.GetTrackId() == trk_id:
            temp_energy = energy_deposit_trk(event.SegmentDetectors, trk.GetTrackId())
            temp_energy = temp_energy / trk.GetInitialMomentum().Gamma()
            reco_energy += temp_energy
            # list_particles.append(trk.GetTrackId())

    # for pid in list_particles:
        # reco_energy += energy_deposit_trk(event.SegmentDetectors, pid)

    return reco_energy
