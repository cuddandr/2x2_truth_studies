import sys
import math
import ROOT as RT

GENIE_STATUS_DEF = {
    -1 : "kIStUndefined",
     0 : "kIStInitialState",              # generator-level initial state
     1 : "kIStStableFinalState",          # generator-level final state: particles to be tracked by detector-level MC
     2 : "kIStIntermediateState",
     3 : "kIStDecayedState",
    10 : "kIStCorrelatedNucleon",
    11 : "kIStNucleonTarget",
    12 : "kIStDISPreFragmHadronicState",
    13 : "kIStPreDecayResonantState",
    14 : "kIStHadronInTheNucleus",        # hadrons inside the nucleus marked for hadron transport modules to act on
    15 : "kIStFinalStateNuclearRemnant",  # low energy nuclear fragments entering the record collectively as a 'hadronic blob' pseudo-particle
    16 : "kIStNucleonClusterTarget",      # for composite nucleons before phase space decay
 }

class Inspector:
    """Class for a collection of methods to inspect and print information from an event in the
    TTree. Designed to be used in an interactive Python session.

    Reads edep-sim output and accesses the event TTree and GENIE pass-through informaion. Loads
    the first event by default.
    """
    def __init__(self, file_list):
        self.load_files(file_list)

    def load_files(self, file_list):
        """Load a list of edep-sim files for inspection. Adds each file to
        a TChain, and loads the GENIE pass-through information.

        Parameters
        ----------
        file_list : List of strings
            List of pathnames for each edep-sim file

        Returns
        -------
        None
        """
        print("Loading files...")
        self.edep_tree = RT.TChain("EDepSimEvents")
        self.grtk_tree = RT.TChain("DetSimPassThru/gRooTracker")
        # self.grtk_tree = RT.TChain("gRooTracker")

        for file in file_list:
            self.edep_tree.Add(file)
            self.grtk_tree.Add(file)

        print("Loading event 0 by default.")
        self.load_event(0)

    def load_event(self, num, verbose=False):
        """Load a specific event from the TTree.

        Parameters
        ----------
        num : int
            Event number to load.
        verbose : bool, optional
            Flag to call TTree::Show() on event after loading.

        Returns
        -------
        None
        """
        self.edep_tree.GetEntry(num)
        self.grtk_tree.GetEntry(num)
        if verbose:
            self.edep_tree.Show(num)

        self.Vtx  = self.edep_tree.Event.Primaries[0]
        self.Traj = self.edep_tree.Event.Trajectories
        self.get_neutrino()
        return

    def get_neutrino(self):
        """Extracts neutrino kinematics and PDG code from the GENIE tree. Information
        is stored as class members for use in other functions.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        genie_evt = self.grtk_tree
        for p in range(genie_evt.StdHepN):
            if genie_evt.StdHepStatus[p] != 0:
                continue

            if math.fabs(genie_evt.StdHepPdg[p]) not in [12, 14, 16]:
                continue

            self.nu_vec = RT.TLorentzVector(genie_evt.StdHepP4[p*4 + 0]*1000,
                                            genie_evt.StdHepP4[p*4 + 1]*1000,
                                            genie_evt.StdHepP4[p*4 + 2]*1000,
                                            genie_evt.StdHepP4[p*4 + 3]*1000)
            self.nu_pdg = genie_evt.StdHepPdg[p]

    def list_genie_stack(self, status_string=False):
        """Prints event and particle information from the GENIE record.

        Parameters
        ----------
        status_string : bool, optional
            Flag to print the particle status as a string instead of the enum value

        Returns
        -------
        None
        """
        genie_evt = self.grtk_tree
        print("EvtCode: ", genie_evt.EvtCode)

        for p in range(genie_evt.StdHepN):
            if status_string:
                p_status = GENIE_STATUS_DEF[genie_evt.StdHepStatus[p]]
                print("Status: {} | PDG: {:3d}".format(p_status, genie_evt.StdHepPdg[p]))
            else:
                print("Status: {:3d} | PDG: {:3d}".format(genie_evt.StdHepStatus[p], genie_evt.StdHepPdg[p]))


    def list_event_kinematics(self):
        """List event kinematics and information.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lep_vec = RT.TLorentzVector()
        lep_pdg = 0

        for p in self.Vtx.Particles:
            if math.fabs(p.GetPDGCode()) in [11, 12, 13, 14, 15, 16]:
                lep_vec = p.GetMomentum()
                lep_pdg = p.GetPDGCode()

        rec = self.Vtx.GetReaction()
        pos = self.Vtx.GetPosition().Vect()

        q0 = (self.nu_vec.E() - lep_vec.E()) / 1000.0
        q3 = (self.nu_vec.Vect() - lep_vec.Vect()).Mag() / 1000.0
        Q2 = (q3**2 - q0**2)

        self.list_neutrino()
        print("PDG: {:3d} | E: {:.3f}, P: {:.3f}".format(lep_pdg, lep_vec.E(), lep_vec.P()))
        print("Q^2: {:.3f}, q0: {:.3f}, q3: {:.3f} GeV".format(Q2, q0, q3))
        print("VTX: ({:.2f}, {:.2f}, {:.2f})".format(pos.X(), pos.Y(), pos.Z()))
        print("Reaction: {}".format(rec))

    def list_neutrino(self):
        """List neutrino information:
        - PDG of neutrino
        - Energy of neutrino

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("PDG: {:3d} | E: {:.3f}".format(self.nu_pdg, self.nu_vec.E()))

    def list_primaries(self):
        """List primary particles from the neutrino interaction.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        for p in self.Vtx.Particles:
            T = p.GetMomentum().E() - p.GetMomentum().M()
            # print("PDG: {:5d} {:12s} | TrkID: {:2d} E: {:.3f}".format(p.GetPDGCode(), p.GetName(), p.GetTrackId(), p.GetMomentum().E()))
            print("PDG: {:5d} {:12s} | TrkID: {:2d} T: {:.3f}".format(p.GetPDGCode(), p.GetName(), p.GetTrackId(), T))

    def list_parent(self, trk_id):
        """List the parent particle of a given track (or lists it as a primary particle).

        Parameters
        ----------
        track_id : int
            ID of the track to inspect

        Returns
        -------
        None
        """
        parent_id = self.Traj[trk_id].GetParentId()

        if parent_id == -1:
            print("Primary particle. Parent ID is -1")
            self.list_neutrino()
        else:
            trk = self.Traj[parent_id]
            self.trk_print(trk)

    def list_children(self, trk_id, filter=[]):
        """List all children track(s) of the given track.

        Parameters
        ----------
        track_id : int
            ID of the track to inspect
        filter : List of ints, optional
            PDG codes of particles/tracks to not print

        Returns
        -------
        None
        """
        for trk in self.Traj:
            parent_id = trk.GetParentId()
            if parent_id == trk_id and trk.GetPDGCode() not in filter:
                self.trk_print(trk)

    def list_ancestors(self, trk_id):
        """List all ancestor particles/tracks of the given track.
        Does not list neutrino.

        Parameters
        ----------
        track_id : int
            ID of the track to inspect

        Returns
        -------
        None
        """
        curr_track = self.Traj[trk_id]
        parent_id  = curr_track.GetParentId()
        print("\u2605", end = " ")
        self.trk_print(curr_track)

        while(parent_id != -1):
            curr_track = self.Traj[parent_id]
            parent_id  = curr_track.GetParentId()
            print("\u2BA1", end = " ")
            self.trk_print(curr_track)

    def find_particle(self, pdg_code):
        """List all tracks with the given PDG code.

        Parameters
        ----------
        pdg_code : int
            PDG code to search for

        Returns
        -------
        None
        """
        for trk in self.Traj:
            if trk.GetPDGCode() == pdg_code:
                self.trk_print(trk)

    def energy_deposit_trk(self, trk_id):
        """Get the total energy deposited along the given track as stored
        in SegmentDetectors. The energy deposit is included if the track
        is the primary contributor.

        Parameters
        ----------
        track_id : int
            ID of the track to get the energy deposit

        Returns
        -------
        reco_energy : the sum of energy deposited in each segment
        """
        reco_energy = 0
        segment_det = self.edep_tree.Event.SegmentDetectors
        for k,v in segment_det:
            for edep in v:
                prim_id = edep.GetPrimaryId()
                contrib = edep.GetContributors()
                # if prim_id != trk_id:
                    # continue
                # reco_energy += edep.GetEnergyDeposit()
                # if trk_id not in contrib:
                    # continue
                # if prim_id == trk_id or trk_id in contrib:
                # if trk_id in contrib:
                if trk_id == contrib[0]:
                    reco_energy += edep.GetEnergyDeposit()

        return reco_energy

    def trk_info(self, trk_id):
        """Print the following information for a given track:
        - PDG and Track ID for track and parent.
        - Energy, momentum, start and end positions.
        - Energy deposited along the track
        - If track was contained in 2x2 volume

        Parameters
        ----------
        track_id : int
            ID of the track to inspect

        Returns
        -------
        None
        """

        trk = self.Traj[trk_id]
        trk_4vec = trk.GetInitialMomentum()

        E = trk_4vec.E()
        P = trk_4vec.P()
        T = trk_4vec.E() - trk_4vec.M()

        parent_id = trk.GetParentId()
        if parent_id != -1:
            parent_name = self.Traj[parent_id].GetName()
        else:
            parent_name = "Primary"

        start = trk.Points[0].GetPosition()
        end   = trk.Points[-1].GetPosition()

        trk_edep = self.energy_deposit_trk(trk_id)

        print("PDG   : {:8s} {:5d} | TrkID: {:4d}".format(trk.GetName(), trk.GetPDGCode(), trk.GetTrackId()))
        print("Parent: {:14s} | TrkID: {:4d}".format(parent_name, parent_id))
        print("Energy = {:.3f}, P = {:.3f}, T = {:.3f}".format(E, P, T))
        print("Energy deposited: {:.4f}".format(trk_edep))
        print("Start : ({:.2f}, {:.2f}, {:.2f})".format(start.X(), start.Y(), start.Z()))
        print("End   : ({:.2f}, {:.2f}, {:.2f})".format(end.X(), end.Y(), end.Z()))
        print("Contained: {}".format(self.is_track_contained(trk)))

    def trk_print(self, trk):
        T = trk.GetInitialMomentum().E() - trk.GetInitialMomentum().M()
        # print("PDG: {:5d} {:8s} | TrkID: {:4d} E: {:.3f}".format(trk.GetPDGCode(), trk.GetName(), trk.GetTrackId(), trk.GetInitialMomentum().E()))
        print("PDG: {:5d} {:8s} | TrkID: {:4d} T: {:.3f}".format(trk.GetPDGCode(), trk.GetName(), trk.GetTrackId(), T))

    def is_point_contained(self, pos):
        """Checks if the given (X,Y,Z) point is within the 2x2 volume.

        Parameters
        ----------
        pos : List or vector of position in (X,Y,Z)

        Returns
        -------
        bool : True if contained, False if not contained
        """
        if abs(pos[0]) > 670: return False
        if abs(pos[1] - 430) > 670: return False
        if abs(pos[2]) > 670: return False
        return True

    def is_track_contained(self, trk):
        """Checks if the given track is contained in the 2x2 volume.

        Parameters
        ----------
        trk : Trajectory object

        Returns
        -------
        bool : True if contained, False if not contained
        """
        end_point = trk.Points[-1].GetPosition().Vect()
        return self.is_point_contained(end_point)

    def help(self):
        method_list = []

        for attribute in dir(Inspector):
            attribute_value = getattr(Inspector, attribute)
            if callable(attribute_value):
                if attribute.startswith('__') == False:
                    method_list.append(attribute)

        print("Available methods: ")
        print(method_list)
        print("Use help(inspector.method_name) to view more information for a given method.")
        print("Use help(inspector) to view all the class documentation at once.")

if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.exit("Requires one or more edep-sim output files as arguments!")
        # print("Requires one or more edep-sim output files as arguments!")
        # exit()

    filelist = [sys.argv[x] for x in range(1, len(sys.argv))]

    print("Instantiating the inspector...")
    print("Call inspector.help() for more information.")
    inspector = Inspector(filelist)
