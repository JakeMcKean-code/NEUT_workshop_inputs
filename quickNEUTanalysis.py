import ROOT
import numpy as np
import sys 
from pathlib import Path

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)

MN = 938.918695
MN2 = MN * MN
PDG_PROTON = 2212
PDG_NEUTRON = 2112
PDG_MUON = 13
PDG_NUMU = 14

def four_vector(part):
    return [
        part.fP.E(),
        part.fP.Px(),
        part.fP.Py(),
        part.fP.Pz(),
    ]

def momentum2(part):
    px = part.fP.Px()
    py = part.fP.Py()
    pz = part.fP.Pz()
    return px*px + py*py + pz*pz

def find_largest_momentum_final_state_proton(neutev):
    best_proton = None
    best_p2 = -1.0

    Num_particles = neutev.Npart()
    for i in range(Num_particles):
        part = neutev.PartInfo(i)

        if part.fPID != PDG_PROTON:
            continue

        # skip initial-state or non-final-state particles
        if part.fStatus == 7:
            continue

        p2 = momentum2(part)

        if p2 > best_p2:
            best_p2 = p2
            best_proton = part

    return best_proton

def get_particle_info_by_pdg(neutev):
    neutrino        = None
    lepton          = None
    initial_nucleon = None
    proton          = None

    Num_particles = neutev.Npart()
    for i in range(Num_particles):
        part = neutev.PartInfo(i)

        # Incoming neutrino
        if part.fPID == PDG_NUMU:
            neutrino = part

        # Final-state lepton
        if part.fPID == PDG_MUON:
            lepton = part

        # Initial-state neutron
        if part.fPID == PDG_NEUTRON:
            initial_nucleon = part

        # Final-state proton 
        if Num_particles == 4:
            if part.fPID == PDG_PROTON:
                proton = part

    # This is for SRC events with the SF CCQE model
    if Num_particles > 4:
        proton = find_largest_momentum_final_state_proton(neutev)

    if neutrino is None:
        return None

    if lepton is None:
        return None

    if initial_nucleon is None:
        return None

    if proton is None:
        return None

    # Make four vectors with neutpart objects
    Kl_inc = four_vector(neutrino)
    Kl = four_vector(lepton)
    PA = four_vector(initial_nucleon)
    PN = four_vector(proton)

    return Kl_inc, Kl, PA, PN


def draw_histograms(histograms, pdf_name):
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    canvas.Print(pdf_name + "[")

    for hist in histograms:
        canvas.Clear()
        hist.Draw("EHIST")
        canvas.Print(pdf_name)
        canvas.Write(hist.GetName())

    canvas.Print(pdf_name + "]")


def quick_neut_analysis(input_file_name):
    fin = ROOT.TFile.Open(input_file_name, "READ")

    if not fin or fin.IsZombie():
        raise RuntimeError(f"Failed to read ROOT file: {input_file_name}")

    neuttree = fin.Get("neuttree")
    if not neuttree:
        raise RuntimeError(f'Failed to read "neuttree" from: {input_file_name}')

    nevs = neuttree.GetEntries()
    print(f"[INFO]: Reading tree with {nevs} NEUT events.")

    # Get scale factor from flux and event rate histograms
    flux_hist = fin.Get("flux_numu")
    event_hist = fin.Get("evtrt_numu")
    scale_factor = (
        event_hist.Integral("width") * 1.0e-38/(nevs * flux_hist.Integral("width"))
    )

    # Make histograms we want
    El_hist  = ROOT.TH1D("El_histo", r"d#sigma / dE_{l}", 75, 0, 3000)
    kl_hist  = ROOT.TH1D("kl_histo", r"d#sigma / dk_{l}", 75, 0, 3000)
    q2_hist  = ROOT.TH1D("Q2_histo", r"d#sigma / dQ^{2}", 30, 0, 3)
    thl_hist = ROOT.TH1D("thl_histo", r"d#sigma / dcos#theta_{l}", 21, -1, 1.1)
    hist_vec = [El_hist, kl_hist, q2_hist, thl_hist]

    # Set errors so we can plot statistical error bars
    for hist in hist_vec:
        hist.Sumw2()

    for ent_it in range(nevs):
        neuttree.GetEntry(ent_it)

        neutev = neuttree.vectorbranch
        npart = neutev.Npart()
        mode = neutev.Mode

        # CCQE only
        if mode != 1:
            continue

        # Event four vectors for particles
        result = get_particle_info_by_pdg(neutev)
        Kl_inc, Kl, PA, PN = result

        # Calculate kinematic variables 
        kl_mag = np.sqrt(Kl[1]**2 + Kl[2]**2 + Kl[3]**2)
        if kl_mag == 0:
            continue

        costhetal = Kl[3]/kl_mag
        q2 = (Kl_inc[0]**2 + kl_mag**2 - 2.0*Kl_inc[0]*kl_mag*costhetal)
        Q2 = (q2 - (Kl_inc[0]-Kl[0])**2)/1.0e6

        # Fill histograms with kinematics
        El_hist.Fill(Kl[0])
        kl_hist.Fill(kl_mag)
        q2_hist.Fill(Q2)
        thl_hist.Fill(costhetal)

    # Convert to differential cross section
    for hist in hist_vec:
        hist.Scale(scale_factor, "width")

    ROOT.gStyle.SetOptStat(0)

    input_path = Path(input_file_name)
    output_root_name = input_path.with_name(input_path.stem + "_analysed_py.root")
    pdf_name = input_path.with_name(input_path.stem + "_analysed_py.pdf")

    print(f"[INFO]: Writing output file: {output_root_name}")

    fout = ROOT.TFile.Open(output_root_name, "RECREATE")
    draw_histograms(hist_vec, pdf_name)

    for hist in hist_vec:
        hist.Write()

    fout.Close()
    fin.Close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"[ERROR]: Expected at least one input file.")
        print(f"[USAGE]: {sys.argv[0]} file1.root file2.root ...")
        sys.exit(1)

    for filename in sys.argv[1:]:
        print(f"\n[INFO]: Processing file: {filename}")
        quick_neut_analysis(filename)
