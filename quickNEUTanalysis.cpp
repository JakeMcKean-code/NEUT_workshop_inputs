#include "neutvect.h"
#include "neutpart.h"
#include "neutvtx.h"

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

void set_four_vector_part(double four_vector[], NeutPart* part)
{
  four_vector[0] = part->fP.E();
  four_vector[1] = part->fP.Px();
  four_vector[2] = part->fP.Py();
  four_vector[3] = part->fP.Pz();
}

NeutPart* FindLargestFinalStateProton(NeutVect* neutev) 
{
  NeutPart* best = nullptr;
  double best_p2 = -1.0;

  for (int i = 0; i < neutev->Npart(); ++i) {
    NeutPart* part = neutev->PartInfo(i);

    if (part->fPID != 2212) continue;   // proton only
    if (part->fStatus == 7) continue;   // skip non-final-state / initial-state particles

    double px = part->fP.Px();
    double py = part->fP.Py();
    double pz = part->fP.Pz();

    double p2 = px*px + py*py + pz*pz;

    if (p2 > best_p2) {
      best_p2 = p2;
      best = part;
    }
  }
  return best;
}

void get_particle_info_by_pdg(NeutVect* neutev,
    double Kl_inc[], double Kl[],
    double PA[], double PN[]) 
{
  NeutPart* neutrino = nullptr;
  NeutPart* lepton = nullptr;
  NeutPart* initial_nucleon = nullptr;

  for (int i = 0; i < neutev->Npart(); ++i) {
    NeutPart* part = neutev->PartInfo(i);

    // muon neutrino
    if (part->fPID == 14) {
      neutrino = part;
    }

    // Muon
    if (part->fPID == 13){ 
      lepton = part;
    }

    // Initial-state neutron
    if (part->fPID == 2112){ 
      initial_nucleon = part;
    }
  }

  // final-state proton
  NeutPart* proton = FindLargestFinalStateProton(neutev);

  if (!neutrino || !lepton || !initial_nucleon || !proton) {
    throw std::runtime_error("Could not identify required particles in event.");
  }

  set_four_vector_part(Kl_inc, neutrino);
  set_four_vector_part(Kl, lepton);
  set_four_vector_part(PA, initial_nucleon);
  set_four_vector_part(PN, proton);
}


void DrawHistograms(std::vector<TH1D*>& histograms,
                    const std::string& pdfName)
{
    TCanvas canvas("canvas","canvas",800,600);
    // Open multipage PDF
    canvas.Print((pdfName + "[").c_str());
    for (TH1D* hist : histograms) {
        canvas.Clear();
        hist->Draw("EHIST");
        // Add page to PDF
        canvas.Print(pdfName.c_str());
        // Write canvas into ROOT file
        canvas.Write(hist->GetName());
    }
    // Close PDF
    canvas.Print((pdfName + "]").c_str());
}


void quickNEUTanalysis(const char* inputFileName){

  // Error handling in case the root file cannot be read
  TFile fin(inputFileName,"READ");
  if(fin.IsZombie()){
    std::cout << "[ERROR]: Failed to read root file: " << inputFileName << std::endl;
    abort();
  }

  // Get pointer for neuttree object which is in the root file
  // Also do some error handling in case it is not there
  TTree *neuttree = 0;
  fin.GetObject("neuttree", neuttree);
  if(!neuttree){
    std::cout << "[ERROR]: Failed to read \"neuttree\" from: " << inputFileName << std::endl;
    abort();
  }

  // Get vector branch from neuttree (need to set the pointer to nullptr first)
  NeutVect *neutev = 0;
  neuttree->SetBranchAddress("vectorbranch", &neutev);

  // Get Num Entires for event loop
  Long64_t nevs = neuttree->GetEntries();
  std::cout << "[INFO]: Reading tree with " << nevs << " NEUT events." << std::endl;

  // Get the event rate and flux histos
  TH1D *fluxHist = (TH1D*)fin.Get("flux_numu");
  TH1D *eventHist = (TH1D*)fin.Get("evtrt_numu");

  // events -> diff xsec scale factor
  double scaleFactor = eventHist->Integral("width")*1E-38/(nevs*fluxHist->Integral("width")); 

  TH1D *El_hist = new TH1D("El_histo", "d#sigma / dE_{l}", 75, 0, 3000);
  TH1D *kl_hist = new TH1D("kl_histo", "d#sigma / dk_{l}", 75, 0, 3000);
  TH1D *q2_hist = new TH1D("Q2_histo", "d#sigma / dQ^{2}", 30, 0, 3);
  TH1D *thl_hist = new TH1D("thl_histo", "d#sigma / dcos#theta_{l}", 21, -1, 1.1);

  std::vector<TH1D*> hist_vec{El_hist, kl_hist, q2_hist, thl_hist};
  for (TH1D* hist : hist_vec){
    hist->Sumw2();
  }

  for (Long64_t ent_it = 0; ent_it < nevs; ++ent_it) {
    // four vectors for the event
    double Kl_inc[4];
    double Kl[4];
    double PA[4];
    double PN[4];

    // Get the event and event information (general)
    neuttree->GetEntry(ent_it);
    int npart = neutev->Npart();
    int mode = neutev->Mode;

    // Look at CCQE only for now
    if(mode != 1){continue;}

    try {
      get_particle_info_by_pdg(neutev, Kl_inc, Kl, PA, PN);
    } catch (const std::runtime_error& e) {
      std::cerr << "[WARNING]: Skipping event " << ent_it
                << ": " << e.what() << std::endl;
      continue;
    }

    // Get kinematic variables from event four vectors
    double kl_mag =  std::sqrt(Kl[1]*Kl[1] + Kl[2]*Kl[2] + Kl[3]*Kl[3]);
    double costhetal = Kl[3]/kl_mag;
    double q2 = Kl_inc[0]*Kl_inc[0] + kl_mag*kl_mag - 2.*Kl_inc[0]*kl_mag*costhetal; 
    double Q2 = (q2 - pow(Kl_inc[0]-Kl[0],2) )/1.E6; //GeV^2

    // Fill histograms 
    El_hist->Fill(Kl[0]);
    kl_hist->Fill(kl_mag);
    q2_hist->Fill(Q2);
    thl_hist->Fill(costhetal);
  }

  // Loop over the vector to access histogram and convert to diff xsec
  for (TH1D* histogram : hist_vec) {
      if (!histogram) {
          std::cerr << "Error: Histogram not found. Continuing." << std::endl;
          continue;
      }
      histogram->Scale(scaleFactor, "width");
  }
  gStyle->SetOptStat(0);

  // Get output filename from stem of input filename (remove .root from input)
  std::string inputFileNameStr = inputFileName;
  std::string stem = inputFileNameStr.substr(0, inputFileNameStr.find_last_of('.'));
  std::string outputRootName = stem + "_analysed_cpp.root";
  std::string pdfName        = stem + "_analysed_cpp.pdf";

  std::cout << "[INFO]: Writing output file: "
            << outputRootName << std::endl;

  TFile* outputFileRoot = TFile::Open(outputRootName.c_str(), "RECREATE");

  // Call the function to draw the histograms
  DrawHistograms(hist_vec, pdfName);
  outputFileRoot->Close();

  for (auto histogram : hist_vec) {
    delete histogram;
  }

}
