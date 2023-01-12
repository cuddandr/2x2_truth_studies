import ROOT as RT
import numpy as np
import sys

import lar_functions as lar

#ROOT.gSystem.Load("/opt/generators/edep-sim/install/lib/libedepsim_io.so")

edep_tree = RT.TChain("EDepSimEvents")
grtk_tree = RT.TChain("gRooTracker")
# grtk_tree = RT.TChain("DetSimPassThru/gRooTracker")

filelist = [sys.argv[x] for x in range(1, len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)

beam_angle = RT.TVector3(0, 0.05836, 1.0) # 3.343 degrees in the y-plane
pion_mass = 139.57
kaon_mass = 497.61
nevt = edep_tree.GetEntries()

h_kaon_pcos = RT.TH2D("k0_pcos", "k0_pcos;#theta; True KE (MeV)", 45, 0, 90.0, 50, 0, 10000)
h_muon_pcos = RT.TH2D("mu_pcos", "mu_pcos;#theta; True KE (MeV)", 45, 0, 90.0, 50, 0, 10000)
h_kaon_mass = RT.TH1D("k0_mass", "k0_mass;Mass (MeV); N", 50, 0, 1000)
h_pion_kint = RT.TH1D("pion_T", "pion_T;T (MeV); N", 100, 0, 5000)
h_evt_q2    = RT.TH1D("h_q2", "h_q2", 50, 0, 5.0)
h_vtx_dist  = RT.TH1D("vtx_dist", "vtx_dist;d (cm); N", 100, 0, 20)

print("Reading {} events...".format(nevt))
for evt in range(nevt):

    if evt % (int(nevt/10)) == 0:
        print("Processed event: ", evt)

    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)

    if not lar.is_hadronic_contained(edep_tree.Event):
        continue

    vtx = edep_tree.Event.Primaries[0]
    num_vtx = len(edep_tree.Event.Primaries)
    primary_pdg = [x.GetPDGCode() for x in vtx.Particles]

    if not np.any(np.isin(np.abs(primary_pdg), [13])):
        continue

    if not np.any(np.isin(np.abs(primary_pdg), [130, 310, 311, 321])):
        continue

    K0s_tid = -1
    traj = edep_tree.Event.Trajectories
    for trk in traj:
        if np.abs(trk.GetPDGCode()) == 310:
            K0s_tid = trk.GetTrackId()
            break

    if K0s_tid == -1:
        continue

    mu_tid = -1
    for trk in traj:
        if np.abs(trk.GetPDGCode()) == 13:
            mu_tid = trk.GetTrackId()
            break

    K0s_decay = []
    for trk in traj:
        if trk.GetParentId() == K0s_tid:
            K0s_decay.append(trk)

    if [x.GetPDGCode() for x in K0s_decay] != [-211, 211]:
        continue

    # pions_contained = True
    # for trk in K0s_decay:
        # pions_contained = pions_contained & lar.is_track_contained(trk)

    # if not pions_contained:
        # continue

    nu_vtx_pos = vtx.GetPosition().Vect()
    nu_vec, nu_pdg = lar.get_nu_vec(grtk_tree)

    k0_vec = traj[K0s_tid].GetInitialMomentum()
    k0_vtx_pos = traj[K0s_tid].Points[-1].GetPosition().Vect()

    k0_angle = k0_vec.Vect().Angle(beam_angle) * 180.0 / np.pi
    k0_KE = k0_vec.E() - k0_vec.M()
    h_kaon_pcos.Fill(k0_angle, k0_KE)

    if mu_tid != -1:
        mu_vec = traj[mu_tid].GetInitialMomentum()
        mu_angle = mu_vec.Vect().Angle(beam_angle) * 180.0 / np.pi
        mu_KE = mu_vec.E() - mu_vec.M()
        h_muon_pcos.Fill(mu_angle, mu_KE)

    if mu_tid != -1:
        q2 = -1 * (mu_vec - nu_vec).Mag2() / 1.0E6
        h_evt_q2.Fill(q2)

    vtx_dist = (nu_vtx_pos - k0_vtx_pos).Mag() / 10.0
    h_vtx_dist.Fill(vtx_dist)

    temp_vec = RT.TLorentzVector()
    for trk in K0s_decay:
        temp_vec += trk.GetInitialMomentum()

    reco_mass = 0
    # for trk in K0s_decay:
        # range = lar.calc_distance(trk) / 10.0
        # T = lar.energy_by_range(range)
        # T_true = trk.GetInitialMomentum().E() - trk.GetInitialMomentum().M()
        # print("Evt {}: T {} vs R {}".format(evt, T_true, T))
        # reco_mass += T + pion_mass

    for trk in K0s_decay:
        # T = lar.energy_deposit_trk(edep_tree.Event.SegmentDetectors, trk.GetTrackId())
        T = lar.edep_plus_children(edep_tree.Event, trk.GetTrackId())
        T_true = trk.GetInitialMomentum().E() - trk.GetInitialMomentum().M()
        # print("Evt {}: T {} vs R {}".format(evt, T_true, T))
        reco_mass += T + pion_mass

    for trk in K0s_decay:
        if np.abs(trk.GetPDGCode()) == 211:
            h_pion_kint.Fill(trk.GetInitialMomentum().E() - trk.GetInitialMomentum().M())

    # gamma = (K0s_decay[0].GetInitialMomentum() + K0s_decay[1].GetInitialMomentum()).Gamma()
    # reco_mass *= (1.0 / gamma)

    # K0s_mass = temp_vec.M()
    # print("Evt {}: {:3f} vs {:3f}".format(evt, K0s_mass, reco_mass))
    h_kaon_mass.Fill(reco_mass)

# can = RT.TCanvas("can", "can", 1000, 800)
# can.cd()
# h_kaon_mass.Draw()
# can.SaveAs("kaon_mass.png")
# h_pion_kint.Draw()
# can.SaveAs("pion_energy.png")
# h_vtx_dist.Draw()
# can.SaveAs("vtx_dist.png")

output_file = RT.TFile("kaon_output.root", "RECREATE")
h_kaon_pcos.Write()
h_muon_pcos.Write()
h_kaon_mass.Write()
h_pion_kint.Write()
h_evt_q2.Write()
h_vtx_dist.Write()
