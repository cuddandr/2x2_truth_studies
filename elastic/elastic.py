import ROOT as RT
import numpy as np
import sys
import time

import lar_functions as lar

#ROOT.gSystem.Load("/opt/generators/edep-sim/install/lib/libedepsim_io.so")

edep_tree = RT.TChain("EDepSimEvents")
grtk_tree = RT.TChain("DetSimPassThru/gRooTracker")

h_proton_ke = RT.TH1D("pr_ke", "pr_ke;True KE (MeV); N", 100, 0, 2500)
h_proton_tcos = RT.TH2D("pr_tcos", "pr_tcos;#theta; True KE (MeV)", 45, 0, 90, 100, 0, 2500)
h_pr_smearing = RT.TH2D("ke_smearing", "ke_smearing; Reco KE (MeV), True KE (MeV)", 100, 0, 2500, 100, 0, 2500)

filelist = [sys.argv[x] for x in range(1, len(sys.argv))]
for file in filelist:
    edep_tree.Add(file)
    grtk_tree.Add(file)

beam_angle = RT.TVector3(0, 0.05836, 1.0) # 3.343 degrees in the y-plane
nevt = edep_tree.GetEntries()
num_nc1p = 0
num_cont = 0

st = time.time()
print("Reading {} events...".format(nevt))
for evt in range(nevt):

    if evt % (int(nevt/10)) == 0:
        print("Processed event: ", evt)

    edep_tree.GetEntry(evt)
    grtk_tree.GetEntry(evt)

    vtx = edep_tree.Event.Primaries[0]
    num_vtx = len(edep_tree.Event.Primaries)
    # print("Num. verticies {}".format(num_vtx))

    primary_pdg = [x.GetPDGCode() for x in vtx.Particles]

    if primary_pdg != [14, 2212]:
        continue

    print("-------------------------")
    print("Event: ", evt, primary_pdg)
    num_nc1p += 1

    proton_tid = -1
    for p in vtx.Particles:
        if p.GetPDGCode() == 2212:
            proton_tid = p.GetTrackId()

    if proton_tid < 0:
        print("Something went wrong with the track ID...")
        continue

    traj = edep_tree.Event.Trajectories
    proton_track = traj[proton_tid]

    if not lar.is_track_contained(proton_track):
        print("Proton track not contained...")
        continue

    num_cont += 1
    edep_energy = 0.0
    for k,v in edep_tree.Event.SegmentDetectors:
        for edep in v:
            prim_id = edep.GetPrimaryId()
            contrib = edep.GetContributors()

            # if prim_id == proton_tid:
                # edep_energy += edep.GetEnergyDeposit()
            if contrib[0] == proton_tid:
                edep_energy += edep.GetEnergyDeposit()

    proton_mass = 938.272
    traj_energy = 0.0
    prev_momentum = proton_track.Points[0].GetMomentum().Mag()
    prev_energy = np.sqrt(prev_momentum**2 + proton_mass**2) - proton_mass
    for pt in proton_track.Points:
        curr_momentum = pt.GetMomentum().Mag()
        curr_energy = np.sqrt(curr_momentum**2 + proton_mass**2) - proton_mass
        dedx = prev_energy - curr_energy
        traj_energy += dedx
        print(dedx)
        prev_energy = curr_energy

    proton_init_4vec = proton_track.GetInitialMomentum()
    proton_init_ke = proton_init_4vec.E() - proton_init_4vec.M()
    proton_init_angle = proton_init_4vec.Vect().Angle(beam_angle) * 180.0 / np.pi

    h_proton_ke.Fill(proton_init_ke)
    h_proton_tcos.Fill(proton_init_angle, proton_init_ke)
    h_pr_smearing.Fill(edep_energy, proton_init_ke)

    print("Proton intial energy: {:.4f}".format(proton_init_ke))
    print("Deposited energy : {:.4f}".format(edep_energy))
    print("Trajectory energy: {:.4f}".format(traj_energy))

et = time.time()
print("Total NC1p: ", num_nc1p)
print("Total cont: ", num_cont)
print("Total time: ", time.strftime("%H:%M:%S", time.gmtime(et-st)))

output_file = RT.TFile("nc_elastic_output.root", "RECREATE")
h_proton_ke.Write()
h_proton_tcos.Write()
h_pr_smearing.Write()
