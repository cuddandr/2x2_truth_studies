import sys
import ROOT
import numpy as np
from optparse import OptionParser

def skim_file(input_file_name, output_file_name):

    ## Open the input file
    edep_chain = ROOT.TChain("EDepSimEvents")
    gtrk_chain = ROOT.TChain("DetSimPassThru/gRooTracker")

    # for file in input_file_list:
    edep_chain.Add(input_file_name)
    gtrk_chain.Add(input_file_name)

    edep_chain.LoadTree(0)
    gtrk_chain.LoadTree(0)

    ## Make the skim file and tree
    skim_file = ROOT.TFile(output_file_name, "RECREATE")
    skim_edep = edep_chain.GetTree().CloneTree(0)
    skim_gtrk = gtrk_chain.GetTree().CloneTree(0)

    ## Count the number saved
    nsaved = 0

    ## Loop over events, decide if they're in the active region
    nevt = edep_chain.GetEntries()
    print("Skimming", nevt, "events from", input_file_name)

    for evt in range(nevt):

        if evt % (int(nevt/10)) == 0:
            print("Processed event: ", evt)

        edep_chain.GetEntry(evt)
        gtrk_chain.GetEntry(evt)

        vtx = edep_chain.Event.Primaries[0]
        num_vtx = len(edep_chain.Event.Primaries)
        primary_pdg = [x.GetPDGCode() for x in vtx.Particles]

        if not np.any(np.isin(np.abs(primary_pdg), [130, 310, 311, 321])):
            continue

        nsaved += 1
        skim_edep.Fill()
        skim_gtrk.Fill()

    ## Save output
    skim_edep.Write("EDepSimEvents")
    #skim_gtrk.Write("DetSimPassThru/gRooTracker")
    skim_gtrk.Write("gRooTracker")
    skim_file.Close()
    print("Saved", nsaved, "events to", output_file_name, "(%.3f)"%(nsaved/float(nevt)))
    return

if __name__ == '__main__':

    ## Get arguments
    parser = OptionParser()
    parser.add_option("-i", "--inFile",  action="store", type="string", dest="inFile")
    parser.add_option("-o", "--outFile", action="store", type="string", dest="outFile")
    (options, sys.argv[1:]) = parser.parse_args()

    # filelist = [sys.argv[x] for x in range(1, len(sys.argv))]
    ## Skim!
    skim_file(options.inFile, options.outFile)
