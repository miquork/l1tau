#define L1Tau_cxx
#include "L1Tau.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TLorentzVector.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <iostream>

bool useAbsEta = true;

double getEff(TProfile2D *p2, double pt, double eta) {

  double x = eta;
  // Adapt automatically from eta to |eta| binning
  if (p2->GetXaxis()->GetBinLowEdge(1)==0) x = fabs(eta);
  int i = p2->GetXaxis()->FindBin(x);
  int j = p2->GetYaxis()->FindBin(pt);
  double eff = p2->GetBinContent(i,j);

  return eff;
}

double getPre(TProfile2D *p2, double pt, double eta) {

  double x = eta;
  if (p2->GetXaxis()->GetBinLowEdge(1)==0) x = fabs(eta);
  int i = p2->GetXaxis()->FindBin(x);
  int j = p2->GetYaxis()->FindBin(pt);
  double pre = p2->GetBinContent(i,j);

  return pre;
}

TH2D *_h2BXvsRUN(0);
bool getFront(int run, int bx) {

  if (!_h2BXvsRUN) {
    TDirectory *curdir = gDirectory;
    //TFile *f = new TFile("../data/l1tau/maps/BXMaps_JetMET1_Run2023C-PromptNanoAODv11p9_v1-v1_NANOAOD.root","READ");
    TFile *f = new TFile("../data/l1tau/maps/Maps_2023.root","READ");
    assert(f && !f->IsZombie());
    curdir->cd();
    _h2BXvsRUN = (TH2D*)f->Get("simpleMap"); assert(_h2BXvsRUN);
  }
  
  // code from jecsys/OffsetTree.C
  TH2D *hbx = _h2BXvsRUN;
  
  // Fill BX at the front and end of bunch train (no early OOT PU)
  int irun = hbx->GetXaxis()->FindBin(run);
  int ibx = hbx->GetYaxis()->FindBin(bx);
  bool bxfront  = (hbx->GetBinContent(irun, ibx)!=0 &&
		   hbx->GetBinContent(irun, ibx-1)==0 &&
		   hbx->GetBinContent(irun, ibx-2)==0 &&
		   hbx->GetBinContent(irun, ibx-3)==0 &&
		   hbx->GetBinContent(irun, ibx-4)==0);
  /*
  bool bxback  =  (hbx->GetBinContent(irun, ibx)!=0 &&
		   hbx->GetBinContent(irun, ibx+1)==0 &&
		   hbx->GetBinContent(irun, ibx+2)==0 &&
		   hbx->GetBinContent(irun, ibx+3)==0 &&
		   hbx->GetBinContent(irun, ibx+4)==0);
  bool bxcenter = (hbx->GetBinContent(irun, ibx)!=0 &&
		   hbx->GetBinContent(irun, ibx-1)!=0 &&
		   hbx->GetBinContent(irun, ibx-2)!=0 &&
		   hbx->GetBinContent(irun, ibx+1)!=0 &&
		   hbx->GetBinContent(irun, ibx+2)!=0);
  */
  return bxfront;
} // getFront

void L1Tau::Loop()
{
//   In a ROOT session, you can do:
//      root> .L L1Tau.C
//      root> L1Tau t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();
   cout << endl;
   cout << "Looping over " << nentries << " entries for " << name << endl;
   
   TDirectory *curdir = gDirectory;
   TFile *fout = new TFile(Form("rootfiles/output-%s.root",name.c_str()),
			   "RECREATE");
   TFile *feff(0);
   //if (name=="23C") feff = new TFile("rootfiles/output-23C-refV2-3.root");
   //if (name=="23D") feff = new TFile("rootfiles/output-23D-refV2-3.root");
   if (name=="23C") feff = new TFile("rootfiles/output-23C-refV3.root");
   if (name=="23D") feff = new TFile("rootfiles/output-23D-refV3.root");
   assert(feff);

   TProfile2D *p2pre34ref = (TProfile2D*)feff->Get("p2pre34");
   assert(p2pre34ref);
   p2pre34ref->SetName(Form("p2pre34ref2_%s",name.c_str()));

   const char *cn = name.c_str();
   TProfile2D *p2pre32ref = (TProfile2D*)feff->Get("p2pre32");
   assert(p2pre32ref);
   p2pre32ref->SetName(Form("p2pre32ref_%s",name.c_str()));
   //
   TProfile2D *p2pre32to34ref = (TProfile2D*)feff->Get("p2pre32to34");
   assert(p2pre32to34ref);
   p2pre32to34ref->SetName(Form("p2pre32to34ref_%s",name.c_str()));
   
   TProfile2D *p2pre26ref = (TProfile2D*)feff->Get("p2pre26");
   assert(p2pre26ref);
   p2pre26ref->SetName(Form("p2pre26ref_%s",name.c_str()));
   
   TProfile2D *p2pre55ref = (TProfile2D*)feff->Get("p2pre55");
   assert(p2pre55ref);
   p2pre55ref->SetName(Form("p2pre55ref_%s",name.c_str()));

   //TProfile2D *p2pre90ref = (TProfile2D*)feff->Get("p2pre90"); // V2-3
   TProfile2D *p2pre90ref = (TProfile2D*)feff->Get("p2pre77p5t"); // V3
   assert(p2pre90ref);
   p2pre90ref->SetName(Form("p2pre90ref_%s",name.c_str()));
   //  
   TProfile2D *p2prejetref = (TProfile2D*)feff->Get("p2prejet");
   assert(p2prejetref);
   p2prejetref->SetName(Form("p2prejetref_%s",name.c_str()));
   
   fout->cd();
   
   TLorentzVector p4l1, p4l1tautag, p4l1jettag, p4l1tauprobe, p4l1jetprobe;
   TLorentzVector p4jet, p4tag, p4probe, p4jj, p4tt;
   TLorentzVector p4l1bx0, p4l1jetbx0, p4l1met, p4l1metnotag;
   TLorentzVector p4ditau;
   
   // Mjj binning from dijet mass search
   double vmjj[] =
     {1, 6, 16, 31, 50, 74, 103, 137, 176, 220, 270, 325, 386, 453, 526, 606, 693, 788, 890, 1000, 1118, 1246, 1383, 1530, 1687, 1856, 2037, 2231, 2438, 2659, 2895, 3147, 3416, 3704, 4010, 4337, 4686, 5058, 5455, 5877, 6328, 6808, 7320, 7866, 8447, 9067, 9726, 10430, 11179, 11977, 12827, 13732};
   double nmjj = sizeof(vmjj)/sizeof(vmjj[0])-1;

   // Inclusive jets pT binning
   double vpti[] = 
     {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
      2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
      4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
   double npti = sizeof(vpti)/sizeof(vpti[0])-1;

   // Widened inclusive jets pT binning above 2 TeV
   double vptw[] = 
     {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
      1327, 1410, 1497, 1588, 1684, 1784, 2000,
      2238, 2500, 2941, 4037, 6076, 7000};
   double nptw = sizeof(vptw)/sizeof(vptw[0])-1;

   // L1Jet binning (inclusive jets 10 GeV to 1.1 TeV, add 30 GeV point)
   double vptbx0[] = 
     {10, 12, 15, 18, 21, 24, 28, 30, 32, 37, 43, 49, 56, 64, 74, 84,
      97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
      507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101};
   double nptbx0 = sizeof(vptbx0)/sizeof(vptbx0[0])-1;

   // JEC L2Relative eta binning
   double veta[] =
     {-5.191,
      -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
      -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
      -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
      -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, 
      -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
      0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191};
   const int neta = sizeof(veta)/sizeof(veta[0])-1;

   // Widened JEC L2Relative eta binning
   double vetaw[] =
     {-5.191, -4.013, -3.314, -3.139, -2.5, -2.043,
      -1.74, -1.566, -1.305, -1.044, -0.783, -0.522,
      -0.261, 0, 0.261,
      0.522, 0.783, 1.044, 1.305, 1.566, 1.74, 2.043, 2.5,
      3.139, 4.013, 5.191};
   const int netaw = sizeof(vetaw)/sizeof(vetaw[0])-1;

   // Widened JEC L2Relative absolute |eta| binning
   double vetaw2[] =
     {0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.566, 1.74, 2.043, 2.5,
      3.139, 4.013, 5.191};
   const int netaw2 = sizeof(vetaw2)/sizeof(vetaw2[0])-1;

   // Controls for sample bias
   TH2D *hptl1metvsjet = new TH2D("hptl1metvsjet",";p_{T,L1Jet};p_{T,L1EtSum} (GeV);",72,0,180,36,0,90);//440,0,1100,440,1100);
   TH2D *hptl1metvsjet_notagjet = new TH2D("hptl1metvsjet_notagjet",";p_{T,L1Jet};p_{T,L1EtSum} (GeV);",72,0,180,36,0,90);//440,0,1100,440,1100);
   TH2D *hptl1metvstau = new TH2D("hptl1metvstau",";p_{T,L1Tau};p_{T,L1EtSum} (GeV);",72,0,180,36,0,90);//440,0,1100,440,1100);
   TH2D *hptl1tauvsjet = new TH2D("hptl1tauvsjet",";p_{T,L1Jet};p_{T,L1Tau} (GeV);",72,0,180,48,0,120);
   TH2D *hptl1jetvsoff = new TH2D("hptl1jetvsoff",";p_{T,jet};p_{T,L1Jet} (GeV);",nptw,vptw,72,0,180);
   TH2D *hptl1tauvsoff = new TH2D("hptl1tauvsoff",";p_{T,jet};p_{T,L1Tau} (GeV);",nptw,vptw,48,0,120);
   
   TH1D *hptl1met = new TH1D("hptl1met",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_wtagjet = new TH1D("hptl1met_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_notagjet = new TH1D("hptl1met_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_unp = new TH1D("hptl1met_unp",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_unp_wtagjet = new TH1D("hptl1met_unp_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_unp_notagjet = new TH1D("hptl1met_unp_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1met_bx1 = new TH1D("hptl1met_bx1",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_bx1_wtagjet = new TH1D("hptl1met_bx1_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_bx1_notagjet = new TH1D("hptl1met_bx1_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1met_pre = new TH1D("hptl1met_pre",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_pre_wtagjet = new TH1D("hptl1met_pre_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_pre_notagjet = new TH1D("hptl1met_pre_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1met_tag = new TH1D("hptl1met_tag",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_tag_wtagjet = new TH1D("hptl1met_tag_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1met_tag_notagjet = new TH1D("hptl1met_tag_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);

   
   TH1D *hptl1metnotag = new TH1D("hptl1metnotag",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_wtagjet = new TH1D("hptl1metnotag_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_notagjet = new TH1D("hptl1metnotag_notagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_unp = new TH1D("hptl1metnotag_unp",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_unp_wtagjet = new TH1D("hptl1metnotag_unp_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_unp_notagjet = new TH1D("hptl1metnotag_unp_notaget",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1metnotag_bx1 = new TH1D("hptl1metnotag_bx1",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_bx1_wtagjet = new TH1D("hptl1metnotag_bx1_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_bx1_notagjet = new TH1D("hptl1metnotag_bx1_notaget",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1metnotag_pre = new TH1D("hptl1metnotag_pre",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_pre_wtagjet = new TH1D("hptl1metnotag_pre_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_pre_notagjet = new TH1D("hptl1metnotag_pre_notaget",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1metnotag_tag = new TH1D("hptl1metnotag_tag",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_tag_wtagjet = new TH1D("hptl1metnotag_tag_wtagjet",";p_{T,L1EtSum} (GeV);",440,0,1100);
   TH1D *hptl1metnotag_tag_notagjet = new TH1D("hptl1metnotag_tag_notaget",";p_{T,L1EtSum} (GeV);",440,0,1100);

   
   //TH1D *hptl1methybrid = new TH1D("hptl1methybrid",";p_{T,L1EtSum} (GeV);",440,0,1100);
   //TH1D *hptl1methybrid_unp = new TH1D("hptl1methybrid_unp",";p_{T,L1EtSum} (GeV);",440,0,1100);

   TH1D *hptl1tau_all = new TH1D("hptl1tau_all",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_all_unp = new TH1D("hptl1tau_all_unp",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_all_bx1 = new TH1D("hptl1tau_all_bx1",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_all_pre = new TH1D("hptl1tau_all_pre",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_all_tag = new TH1D("hptl1tau_all_tag",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau = new TH1D("hptl1tau",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_unp = new TH1D("hptl1tau_unp",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_bx1 = new TH1D("hptl1tau_bx1",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_pre = new TH1D("hptl1tau_pre",";p_{T,L1tau} (GeV);",440,0,1100);
   TH1D *hptl1tau_tag = new TH1D("hptl1tau_tag",";p_{T,L1tau} (GeV);",440,0,1100);
   
   TH1D *hptl1jet = new TH1D("hptl1jet",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_unp = new TH1D("hptl1jet_unp",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_bx1 = new TH1D("hptl1jet_bx1",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_pre = new TH1D("hptl1jet_pre",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_tag = new TH1D("hptl1jet_tag",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_wtagjet = new TH1D("hptl1jet_wtagjet",";p_{T,L1jet} (GeV);",440,0,1100);
   TH1D *hptl1jet_notagjet = new TH1D("hptl1jet_notagjet",";p_{T,L1jet} (GeV);",440,0,1100);
   
   // L1_DoubleIsoTau34er2p1
   // Initial testing to wind suitable bin widths
   TProfile *p1eff = new TProfile("p1eff",";p_{T,probe} (GeV);L1Tau eff.;",
				  npti, vpti);
   TProfile *p1effw = new TProfile("p1effw",";p_{T,probe} (GeV);L1Tau eff.;",
				   nptw, vptw);
   TProfile *p1pre = new TProfile("p1pre",";p_{T,probe} (GeV);L1Tau prefiring;",
				  npti, vpti);
   TProfile *p1prew = new TProfile("p1prew",";p_{T,probe} (GeV);L1Tau prefiring;",
				   nptw, vptw);
   TProfile2D *p2eff = new TProfile2D("p2eff",";#eta_{jet,probe};"
				      "p_{T,probe} (GeV);L1Tau eff.;",
				      neta, veta, npti, vpti);
   TProfile2D *p2effw = new TProfile2D("p2effw",";#eta_{jet,probe};"
				       "p_{T,probe} (GeV);L1Tau eff.;",
				       netaw, vetaw, nptw, vptw);
   TProfile2D *p2effw2 = new TProfile2D("p2effw2",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Tau eff.;",
					netaw2, vetaw2, nptw, vptw);
   TProfile2D *p2pre = new TProfile2D("p2pre",";#eta_{jet,probe};"
				      "p_{T,probe} (GeV);L1Tau prefiring;",
				      neta, veta, npti, vpti);
   TProfile2D *p2prew = new TProfile2D("p2prew",";#eta_{jet,probe};"
				       "p_{T,probe} (GeV);L1Tau prefiring;",
				       netaw, vetaw, nptw, vptw);
   TProfile2D *p2prew2 = new TProfile2D("p2prew2",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Tau prefiring;",
					netaw2, vetaw2, nptw, vptw);

   TProfile *p1prebx0 = new TProfile("p1prebx0",";p_{T,L1jet,probe,BX0} (GeV);"
				     "L1Tau26 prefiring;",nptbx0, vptbx0);
   TProfile2D *p2prebx0 = new TProfile2D("p2prebx0",";|#eta_{L1jet,probe,BX0}|;"
					 "p_{T,L1jet,probe,BX0} (GeV);"
					 "L1Tau26 prefiring;",
					 netaw2, vetaw2, nptbx0, vptbx0);
   TH2D *h2jetbx0bx1 = new TH2D("h2jetbx0bx2",";p_{T,L1jet,probe,BX0} (GeV);p_{T,L1jet,probe,BX-1} (GeV);", nptbx0, vptbx0, nptw, vptw);
				
   
   // Proxy for L1MET90
   TProfile *p1pre75t = new TProfile("p1pre75t",";p_{T,probe} (GeV);"
				     "L1Tau75 prefiring;",nptw, vptw);
   TProfile2D *p2pre75t = new TProfile2D("p2pre75t",";|#eta_{jet,probe}|;"
					 "p_{T,probe} (GeV);L1Tau75 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre77p5t = new TProfile("p1pre77p5t",";p_{T,probe} (GeV);"
				     "L1Tau77.5 prefiring;",nptw, vptw);
   TProfile2D *p2pre77p5t = new TProfile2D("p2pre77p5t",";|#eta_{jet,probe}|;p_{T,probe} (GeV);L1Tau77.5 prefiring;",netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre80t = new TProfile("p1pre80t",";p_{T,probe} (GeV);"
				     "L1Tau80 prefiring;",nptw, vptw);
   TProfile2D *p2pre80t = new TProfile2D("p2pre80t",";|#eta_{jet,probe}|;"
					 "p_{T,probe} (GeV);L1Tau80 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   
   // Final efficiency maps (1D and 2D) for various triggers
   // L1_DoubleIsoTau34er2p1
   TProfile *p1pre34 = new TProfile("p1pre34",";p_{T,probe} (GeV);"
				    "L1Tau34 prefiring;",nptw, vptw);
   TProfile2D *p2pre34 = new TProfile2D("p2pre34",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Tau32 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   
   // L1_DoubleIsoTau32er2p1_Mass_Max80
   TProfile *p1pre32 = new TProfile("p1pre32",";p_{T,probe} (GeV);"
				    "L1Tau32 prefiring;",nptw, vptw);
   TProfile2D *p2pre32 = new TProfile2D("p2pre32",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Tau32 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre32to34 = new TProfile("p1pre32to34",";p_{T,probe} (GeV);"
					"L1Tau32-34 prefiring;",nptw, vptw);
   TProfile2D *p2pre32to34 = new TProfile2D("p2pre32to34",";|#eta_{jet,probe}|;p_{T,probe} (GeV);L1Tau32-34 prefiring;",netaw2, vetaw2, nptw, vptw);
   
   // L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dr0p5 (ditau leg, jet leg)
   TProfile *p1pre26 = new TProfile("p1pre26",";p_{T,probe} (GeV);"
				    "L1Tau26 prefiring;",nptw, vptw);
   TProfile2D *p2pre26 = new TProfile2D("p2pre26",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Tau26 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre55 = new TProfile("p1pre55",";p_{T,probe} (GeV);"
				    "L1Jet55 prefiring;",nptw, vptw);
   TProfile2D *p2pre55 = new TProfile2D("p2pre55",";|#eta_{jet,probe}|;"
					"p_{T,probe} (GeV);L1Jet55 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   
   // L1_ETMHF90 (leading jet leg, imbalancing opposite jet leg)
   TProfile *p1pre90 = new TProfile("p1pre90",";p_{T,probe} (GeV);"
				    "L1Jet90 prefiring;",nptw, vptw);
   TProfile2D *p2pre90 = new TProfile2D("p2pre90",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1Jet90 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   //
   TProfile *p1pre90h = new TProfile("p1pre90h",";p_{T,probe} (GeV);"
				     "L1MET90 prefiring;",nptw, vptw);
   TProfile2D *p2pre90h = new TProfile2D("p2pre90h",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1MET90 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   //
   TProfile *p1pre90_wtagjet = new TProfile("p1pre90_wtagjet",";p_{T,probe} (GeV);L1Jet90 prefiring;",nptw, vptw);
   TProfile2D *p2pre90_wtagjet = new TProfile2D("p2pre90_wtagjet",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1Jet90 prefiring;",netaw2, vetaw2, nptw, vptw);
   //
   TProfile *p1pre90_notagjet = new TProfile("p1pre90_notagjet",";p_{T,probe} (GeV);L1Jet90 prefiring;",nptw, vptw);
   TProfile2D *p2pre90_notagjet = new TProfile2D("p2pre90_notagjet",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1Jet90 prefiring;",netaw2, vetaw2, nptw, vptw);
   //
   TProfile *p1prejet = new TProfile("p1prejet",";p_{T,probe} (GeV);"
				     "L1Jet prefiring;",nptw, vptw);
   TProfile2D *p2prejet = new TProfile2D("p2prejet",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Jet prefiring;",
					 netaw2, vetaw2, nptw, vptw);

   // L1_SingleJet180
   TProfile *p1pre180jet = new TProfile("p1pre180jet",";p_{T,probe} (GeV);"
					"L1Jet180 prefiring;",nptw, vptw);
   TProfile2D *p2pre180jet = new TProfile2D("p2pre180jet",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1Jet180 prefiring;",netaw2, vetaw2, nptw, vptw);

   // L1_SingleTau120
   TProfile *p1pre120tau = new TProfile("p1pre120tau",";p_{T,probe} (GeV);"
				     "L1Tau120 prefiring;",nptw, vptw);
   TProfile2D *p2pre120tau = new TProfile2D("p2pre120tau",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1Tau120 prefiring;",netaw2, vetaw2, nptw, vptw);

   // Ditau
   TProfile *p1pre34ditau = new TProfile("p1pre34ditau",";p_{T,probe} (GeV);"
				       "DiTau34 prefiring;",nptw, vptw);
   TProfile2D *p2pre34ditau = new TProfile2D("p2pre34ditau",";|#eta_{jet,probe}|;p_{T,jet} (GeV);DiTau34 prefiring;",netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre32ditau = new TProfile("p1pre32ditau",";p_{T,probe} (GeV);"
					 "DiTau32 prefiring;",nptw, vptw);
   TProfile2D *p2pre32ditau = new TProfile2D("p2pre32ditau",";|#eta_{jet,probe}|;p_{T,jet} (GeV);DiTau32 prefiring;",netaw2, vetaw2, nptw, vptw);
   TProfile *p1prebx1 = new TProfile("p1prebx1",";p_{T,probe} (GeV);"
				      "BX1 prefiring;",nptw, vptw);
   TProfile2D *p2prebx1 = new TProfile2D("p2prebx1",";|#eta_{jet,probe}|;p_{T,jet} (GeV);BX1 prefiring;",netaw2, vetaw2, nptw, vptw);
   
   // L1_ETMHF90
   TProfile *p1pre90met = new TProfile("p1pre90met",";p_{T,probe} (GeV);"
				       "L1MET90 prefiring;",nptw, vptw);
   TProfile2D *p2pre90met = new TProfile2D("p2pre90met",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1MET90 prefiring;",netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre90metno180jet = new TProfile("p1pre90metno180jet",";p_{T,probe} (GeV);L1MET90 (not Jet180) prefiring;",nptw, vptw);
   TProfile2D *p2pre90metno180jet = new TProfile2D("p2pre90metno18jet",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1MET90 (not Jet180) prefiring;",netaw2, vetaw2, nptw, vptw);

   
   // L1_SingleTau120 + L1_SingleJet180
   TProfile *p1pre120and180 = new TProfile("p1pre120and180",";p_{T,probe} (GeV);L1Tau120+L1Jet180 prefiring;",nptw, vptw);
   TProfile2D *p2pre120and180 = new TProfile2D("p2pre120and180",";|#eta_{jet,probe}|;p_{T,jet} (GeV);L1Tau120+L1Jet180 prefiring;",netaw2, vetaw2, nptw, vptw);
   
   // Proxies for L1_SingleJet180, L1_ETMHF90, combinatorics
   TProfile *p1pre60m = new TProfile("p1pre60m",";p_{T,probe} (GeV);"
				     "L1MET60 prefiring;",nptw, vptw);
   TProfile2D *p2pre60m = new TProfile2D("p2pre60m",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1MET60 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre70m = new TProfile("p1pre70m",";p_{T,probe} (GeV);"
				    "L1MET70 prefiring;",nptw, vptw);
   TProfile2D *p2pre70m = new TProfile2D("p2pre70m",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1MET70 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre70 = new TProfile("p1pre70",";p_{T,probe} (GeV);"
				    "L1Jet70 prefiring;",nptw, vptw);
   TProfile2D *p2pre70 = new TProfile2D("p2pre70",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1Jet70 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre80 = new TProfile("p1pre80",";p_{T,probe} (GeV);"
				    "L1Jet80 prefiring;",nptw, vptw);
   TProfile2D *p2pre80 = new TProfile2D("p2pre80",";|#eta_{jet,probe}|;"
					"p_{T,jet} (GeV);L1Jet80 prefiring;",
					netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre80m = new TProfile("p1pre80m",";p_{T,probe} (GeV);"
				    "L1MET80 prefiring;",nptw, vptw);
   TProfile2D *p2pre80m = new TProfile2D("p2pre80m",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1MET80 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre110 = new TProfile("p1pre110",";p_{T,probe} (GeV);"
				     "L1Jet110 prefiring;",nptw, vptw);
   TProfile2D *p2pre110 = new TProfile2D("p2pre110",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Jet110 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre120 = new TProfile("p1pre120",";p_{T,probe} (GeV);"
				     "L1Jet120 prefiring;",nptw, vptw);
   TProfile2D *p2pre120 = new TProfile2D("p2pre120",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Jet110 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre130 = new TProfile("p1pre130",";p_{T,probe} (GeV);"
				     "L1Jet130 prefiring;",nptw, vptw);
   TProfile2D *p2pre130 = new TProfile2D("p2pre130",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Jet130 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   TProfile *p1pre140 = new TProfile("p1pre140",";p_{T,probe} (GeV);"
				     "L1Jet140 prefiring;",nptw, vptw);
   TProfile2D *p2pre140 = new TProfile2D("p2pre140",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Jet140 prefiring;",
					 netaw2, vetaw2, nptw, vptw);
   
   // Proxies for combinatorics
   TProfile *p1pretau = new TProfile("p1pretau",";p_{T,probe} (GeV);"
				     "L1Tau prefiring;",nptw, vptw);
   TProfile2D *p2pretau = new TProfile2D("p2pretau",";|#eta_{jet,probe}|;"
					 "p_{T,jet} (GeV);L1Tau prefiring;",
					 netaw2, vetaw2, nptw, vptw);


   // Estimates of prefiring efficiency
   TProfile *p34mjj = new TProfile("p34mjj",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *p32mjj = new TProfile("p32mjj",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *p32to34mjj = new TProfile("p32to34mjj",";M_{jj};Prefiring",
				       nmjj,vmjj);
   TProfile *p26mjj = new TProfile("p26mjj",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *p90mjj = new TProfile("p90mjj",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *ptotmjj = new TProfile("ptotmjj",";M_{jj};Prefiring",nmjj,vmjj);

   TProfile *p32mjjbx1 = new TProfile("p32mjjbx1",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *p34mjjbx1 = new TProfile("p34mjjbx1",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *p90mjjbx1 = new TProfile("p90mjjbx1",";M_{jj};Prefiring",nmjj,vmjj);
   TProfile *ptotmjjbx1 = new TProfile("ptotmjjbx1",";M_{jj};Prefiring",nmjj,vmjj);

   
   // other controls and plots
   TProfile *ptag = new TProfile("ptag",";p_{T,tag} (GeV);Fraction",nptw,vptw);
   TProfile *pprobe = new TProfile("pprobe",";p_{T,probe} (GeV);Fraction",
				   nptw,vptw);

   TH1D *hpttag0 = new TH1D("hpttag0",";p_{T,tag} (GeV);",nptw,vptw);
   TH1D *hptprobe0 = new TH1D("hptprobe0",";p_{T,probe} (GeV);",nptw,vptw);
   TH1D *hpttag1 = new TH1D("hpttag1",";p_{T,tag} (GeV);",nptw,vptw);
   TH1D *hptprobe1 = new TH1D("hptprobe1",";p_{T,probe} (GeV);",nptw,vptw);
   TH1D *hpttag2 = new TH1D("hpttag2",";p_{T,tag} (GeV);",nptw,vptw);
   TH1D *hptprobe2 = new TH1D("hptprobe2",";p_{T,probe} (GeV);",nptw,vptw);
   TH1D *hpttag = new TH1D("hpttag",";p_{T,tag} (GeV);N_{jet}",npti,vpti);
   TH1D *hptprobe = new TH1D("hptprobe",";p_{T,probe} (GeV);N_{jet}",npti,vpti);
   
   TH1D *hetatag = new TH1D("hetatag",";#eta_{tag};N_{jet}",neta,veta);
   TH1D *hetaprobe = new TH1D("hetaprobe",";#eta_{probe};N_{jet}",neta,veta);

   TH1D *hmjj = new TH1D("hmjj",";M_{jj} (GeV);N_{event}",nmjj,vmjj);
   TH1D *hmjjc = new TH1D("hmjjc",";M_{jj} (GeV);N_{event,corr}",nmjj,vmjj);
   TH1D *hmjj13 = new TH1D("hmjj13",";M_{jj} (GeV);N_{event}",nmjj,vmjj);
   TH1D *hmjj13c = new TH1D("hmjj13c",";M_{jj} (GeV);N_{event,corr}",nmjj,vmjj);

   TH2D *h2mtt = new TH2D("h2mtt",";p_{T,probe} (GeV);M_{#tau#tau} (GeV)",
			  nptw,vptw, 200,0,200);
   TH2D *h2mtt32 = new TH2D("h2mtt32",";p_{T,probe} (GeV);M_{#tau#tau} (GeV);",
			    nptw,vptw, 200,0,200);
   TH2D *h2mttr = new TH2D("h2mttr",";p_{T,probe} (GeV);M_{tt}/M_{80};",
			   nptw,vptw, 200,0,200);
   TH2D *h2mtt32r = new TH2D("h2mtt32r",";p_{T,probe} (GeV);M_{tt}/M_{80};",
			     nptw,vptw, 200,0,200);
   
   TH1D *hpt = new TH1D("hpt",";p_{T} (GeV);N_{jet}",npti,vpti);
   TH1D *hptc = new TH1D("hptc",";p_{T} (GeV);N_{jet,corr}",npti,vpti);
   TH1D *hpt13 = new TH1D("hpt13",";p_{T} (GeV);N_{jet}",npti,vpti);
   TH1D *hpt13c = new TH1D("hpt13c",";p_{T} (GeV);N_{jet,corr}",npti,vpti);

   
   fChain->SetBranchStatus("*",0);  // disable all branches

   fChain->SetBranchStatus("L1_UnprefireableEvent",1);
   fChain->SetBranchStatus("Flag_Run3",1);

   fChain->SetBranchStatus("run",1);
   fChain->SetBranchStatus("bunchCrossing",1);
   
   fChain->SetBranchStatus("nJet",1);
   fChain->SetBranchStatus("Jet_pt",1);
   fChain->SetBranchStatus("Jet_eta",1);
   fChain->SetBranchStatus("Jet_phi",1);
   fChain->SetBranchStatus("Jet_mass",1);
   fChain->SetBranchStatus("Jet_jetId",1);

   fChain->SetBranchStatus("nL1EtSum",1);
   fChain->SetBranchStatus("L1EtSum_bx",1);
   fChain->SetBranchStatus("L1EtSum_etSumType",1);
   fChain->SetBranchStatus("L1EtSum_pt",1);
   fChain->SetBranchStatus("L1EtSum_phi",1);

   fChain->SetBranchStatus("nL1Tau",1);
   fChain->SetBranchStatus("L1Tau_hwIso",1);
   fChain->SetBranchStatus("L1Tau_bx",1);
   fChain->SetBranchStatus("L1Tau_pt",1);
   fChain->SetBranchStatus("L1Tau_eta",1);
   fChain->SetBranchStatus("L1Tau_phi",1);

   fChain->SetBranchStatus("nL1Jet",1);
   fChain->SetBranchStatus("L1Jet_bx",1);
   fChain->SetBranchStatus("L1Jet_pt",1);
   fChain->SetBranchStatus("L1Jet_eta",1);
   fChain->SetBranchStatus("L1Jet_phi",1);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      // Patch missing variable
      L1_UnprefireableBX1 = getFront(run,bunchCrossing);
      L1_Prefireable = !(L1_UnprefireableEvent || L1_UnprefireableBX1);
      
      if (jentry%100000==0) cout << "." << flush;

      bool pass_njet = (nJet>=2); assert((int)Jet_pt->size()==nJet);
      bool pass_jetid = (nJet>=2 ? Jet_jetId[0]>=4 && Jet_jetId[1]>=4 : 0);

      if (nJet>2) p4jet.SetPtEtaPhiM((*Jet_pt)[2],Jet_eta[2],Jet_phi[2],
				     Jet_mass[2]);
      else p4jet.SetPtEtaPhiM(0,0,0,0);
      
      // Tag-and-probe selection
      for (int itag = 0; itag != 2 && pass_njet; ++itag) {

	int iprobe = (itag==0 ? 1 : 0);

	p4tag.SetPtEtaPhiM((*Jet_pt)[itag],Jet_eta[itag],Jet_phi[itag],
			   Jet_mass[itag]);
	p4probe.SetPtEtaPhiM((*Jet_pt)[iprobe],Jet_eta[iprobe],Jet_phi[iprobe],
			     Jet_mass[iprobe]);
	bool pass_trigmatch = (p4tag.Pt()>600.);
	bool pass_tageta = true;//(fabs(p4tag.Eta())<1.3);
	double eta = (useAbsEta ? fabs(p4probe.Eta()) : p4probe.Eta());

	// Find L1Tau matching tag jet, if any
	bool has_l1taumatch = false;
	p4l1tautag.SetPtEtaPhiM(0,Jet_eta[itag],Jet_phi[itag],0);
	for (int i = 0; i != nL1Tau && !has_l1taumatch; ++i) {
	  if (L1Tau_bx[i]==-1 && L1Tau_hwIso[i]>=1 && fabs(L1Tau_eta[i])<2.1) {
	    p4l1.SetPtEtaPhiM(L1Tau_pt[i],L1Tau_eta[i],L1Tau_phi[i],0.);
	    if (p4l1.DeltaR(p4tag)<0.4) {
	      p4l1tautag = p4l1;
	      has_l1taumatch = true;
	    }
	  }
	} // for i (tag)
	// Threshold could be 34, 32, 26, or 26 + Jet55 veto
	// Looks like mostly 32 needed, but 26 here for extra safety
	bool pass_l1tautag = (!has_l1taumatch || p4l1tautag.Pt()<26.);
	bool fail_l1tautag = (has_l1taumatch && p4l1tautag.Pt()>34. &&
			      L1_Prefireable);
	
	// Find L1Jet matching tag jet, if any
	bool has_l1jetmatch = false;
	p4l1jettag.SetPtEtaPhiM(0,Jet_eta[itag],Jet_phi[itag],0);
	for (int i = 0; i != nL1Jet && !has_l1jetmatch; ++i) {
	  if (L1Jet_bx[i]==-1) {
	    p4l1.SetPtEtaPhiM(L1Jet_pt[i],L1Jet_eta[i],L1Jet_phi[i],0.);
	    //L1Jet_mass[i]);
	    if (p4l1.DeltaR(p4tag)<0.4) {
	      p4l1jettag = p4l1;
	      has_l1jetmatch = true;
	    }
	  }
	} // for i (tag)
	// Jet veto for DiTau26+Jet55. Vetoing SingleJet180 useless
	bool pass_l1jettag = (!has_l1jetmatch || p4l1jettag.Pt()<55.);

	// Find L1Tau matching probe jet, if any
	bool has_l1tauprobe = false;
	p4l1tauprobe.SetPtEtaPhiM(0,Jet_eta[iprobe],Jet_phi[iprobe],0);
	for (int i = 0; i != nL1Tau && !has_l1tauprobe; ++i) {
	  if (L1Tau_bx[i]==-1 && L1Tau_hwIso[i]>=1 && fabs(L1Tau_eta[i])<2.1) {
	    p4l1.SetPtEtaPhiM(L1Tau_pt[i],L1Tau_eta[i],L1Tau_phi[i],0.);
	    if (p4l1.DeltaR(p4probe)<0.4) {
	      p4l1tauprobe = p4l1;
	      has_l1tauprobe = true;
	    }
	  }
	} // for i (probe)
	bool pass_l1pt120tau = (!has_l1tauprobe || p4l1tauprobe.Pt()<120.);
	bool pass_l1pt80t = (!has_l1tauprobe || p4l1tauprobe.Pt()<80.);
	bool pass_l1pt77p5t = (!has_l1tauprobe || p4l1tauprobe.Pt()<77.5);
	bool pass_l1pt75t = (!has_l1tauprobe || p4l1tauprobe.Pt()<75.);
	bool pass_l1pt34 = (!has_l1tauprobe || p4l1tauprobe.Pt()<34.);
	bool pass_l1pt32 = (!has_l1tauprobe || p4l1tauprobe.Pt()<32.);
	bool pass_l1pt32to34 = !(p4l1tauprobe.Pt()>32. &&
				 p4l1tauprobe.Pt()<34.);
	bool pass_l1pt26 = (!has_l1tauprobe || p4l1tauprobe.Pt()<26.);
	bool pass_l1tau = (!has_l1tauprobe);

	bool pass_ditau34 = !(p4l1tautag.Pt()>34. && p4l1tauprobe.Pt()>34.);
	p4ditau = (p4l1tautag + p4l1tauprobe);
	bool pass_ditau32 = !(p4l1tautag.Pt()>32. && p4l1tauprobe.Pt()>32. &&
			      p4ditau.M()<80.);
	// still to add ditau26+jet55

	// Find L1Jet matching probe jet, if any
	bool has_l1jetprobe = false;
	p4l1jetprobe.SetPtEtaPhiM(0,Jet_eta[iprobe],Jet_phi[iprobe],0);
	for (int i = 0; i != nL1Jet && !has_l1jetprobe; ++i) {
	  if (L1Jet_bx[i]==-1) {
	    p4l1.SetPtEtaPhiM(L1Jet_pt[i],L1Jet_eta[i],L1Jet_phi[i],0.);
	      //L1Jet_mass[i);
	    if (p4l1.DeltaR(p4probe)<0.4) {
	      p4l1jetprobe = p4l1;
	      has_l1jetprobe = true;
	    }
	  }
	} // for i (probe)
	bool pass_l1pt180 = (!has_l1jetprobe || p4l1jetprobe.Pt()<180.);
	bool pass_l1pt140 = (!has_l1jetprobe || p4l1jetprobe.Pt()<140.);
	bool pass_l1pt130 = (!has_l1jetprobe || p4l1jetprobe.Pt()<130.);
	bool pass_l1pt120 = (!has_l1jetprobe || p4l1jetprobe.Pt()<120.);
	bool pass_l1pt110 = (!has_l1jetprobe || p4l1jetprobe.Pt()<110.);
	bool pass_l1pt90 = (!has_l1jetprobe || p4l1jetprobe.Pt()<90.);
	bool pass_l1pt80 = (!has_l1jetprobe || p4l1jetprobe.Pt()<80.);
	bool pass_l1pt70 = (!has_l1jetprobe || p4l1jetprobe.Pt()<70.);
	bool pass_l1pt55 = (!has_l1jetprobe || p4l1jetprobe.Pt()<55.);
	bool pass_l1jet = (!has_l1jetprobe);

	// Find L1MET matching L1_ETMHF90 trigger
	bool has_l1met = false;
	p4l1met.SetPtEtaPhiM(0,0,0,0);
	for (int i = 0; i != nL1EtSum && !has_l1met; ++i) {
	  if (L1EtSum_bx[i]==-1 && L1EtSum_etSumType[i]==8) {
	    p4l1met.SetPtEtaPhiM(L1EtSum_pt[i],0.,L1EtSum_phi[i],0.);
	    has_l1met = true;
	  }
	}
	
	// Special shifted MET without tag object(s)
	if (has_l1jetmatch) p4l1metnotag = p4l1met + p4l1jettag;
	//else if (has_l1taumatch) p4l1metnotag = p4l1met + p4l1tautag;
	else p4l1metnotag = p4l1met;
	bool hasl1tag = has_l1jetmatch;// || has_l1taumatch); 

	//bool pass_l1met90h = (p4l1metnotag.Pt()<90.);
	bool pass_l1metnotag90 = (p4l1metnotag.Pt()<90.);
	bool pass_l1met90 = (p4l1met.Pt()<90.);
	bool pass_l1met80 = (p4l1met.Pt()<80. || p4l1met.Pt()>90.);
	bool pass_l1met70 = (p4l1met.Pt()<70. || p4l1met.Pt()>80.);
	bool pass_l1met60 = (p4l1met.Pt()<60. || p4l1met.Pt()>70.);
	
	// All prefiring triggers (still to implement ditau26_jet55)
	bool pass_bx1 = (pass_ditau32 && pass_ditau34 && pass_l1met90 &&
			 pass_l1pt180 && pass_l1pt120tau);

	
	// Tag-and-probe selection
	bool pass_tnp = (pass_trigmatch && pass_l1tautag && pass_l1jettag &&
			 pass_tageta && pass_jetid && Flag_Run3);
	bool fail_tnp = (pass_trigmatch && fail_l1tautag && 
			 pass_tageta && pass_jetid && Flag_Run3);

	// Single-leg selection
	bool pass_single = (pass_trigmatch && pass_tageta && pass_jetid &&
			    Flag_Run3 && L1_UnprefireableBX1);
	
	// Controls of tag-and-probe selection efficiencies
	if (Flag_Run3) {
	  ptag->Fill(p4tag.Pt(), pass_tnp ? 1 : 0);
	  pprobe->Fill(p4probe.Pt(), pass_tnp ? 1 : 0);

	  hpttag0->Fill(p4tag.Pt());
	  hptprobe0->Fill(p4probe.Pt());
	}
	if (pass_single) {
	  hpttag1->Fill(p4tag.Pt());
	  hptprobe1->Fill(p4probe.Pt());
	}
	if (pass_tnp) {
	  hpttag2->Fill(p4tag.Pt());
	  hptprobe2->Fill(p4probe.Pt());
	}

	if (pass_trigmatch && pass_tageta && pass_jetid && Flag_Run3) {
	  hptl1tau_all->Fill(p4l1tauprobe.Pt());
	  if (L1_UnprefireableEvent) hptl1tau_all_unp->Fill(p4l1tauprobe.Pt());
	  if (L1_UnprefireableBX1)   hptl1tau_all_bx1->Fill(p4l1tauprobe.Pt());
	  if (L1_Prefireable)        hptl1tau_all_pre->Fill(p4l1tauprobe.Pt());
	  if (fail_l1tautag)         hptl1tau_all_tag->Fill(p4l1tauprobe.Pt());
	}
      
	// Single leg efficiencies from unprefireable events
	if (pass_single) {
	  p2pre32ditau->Fill(eta, p4probe.Pt(), pass_ditau32 ? 0 : 1);
	  p2pre34ditau->Fill(eta, p4probe.Pt(), pass_ditau34 ? 0 : 1);
	  p2pre90met->Fill(eta, p4probe.Pt(), pass_l1met90 ? 0 : 1);
	  p2pre90metno180jet->Fill(eta, p4probe.Pt(), pass_l1met90 ? 0 : (!pass_l1pt180 ? 0 : 1));
	  p2pre120tau->Fill(eta,p4probe.Pt(),pass_l1pt120tau ? 0 : 1);
	  p2pre180jet->Fill(eta, p4probe.Pt(), pass_l1pt180 ? 0 : 1);
	  p2pre120and180->Fill(eta, p4probe.Pt(),
			       !pass_l1pt120tau && !pass_l1pt180 ? 1 : 0);
	  p2prebx1->Fill(eta, p4probe.Pt(), pass_bx1 ? 0 : 1);
	  
	  if (fabs(p4probe.Eta())<1.3) {
	    p1pre32ditau->Fill(p4probe.Pt(), pass_ditau32 ? 0 : 1);
	    p1pre34ditau->Fill(p4probe.Pt(), pass_ditau34 ? 0 : 1);
	    p1pre90met->Fill(p4probe.Pt(), pass_l1met90 ? 0 : 1);
	    p1pre90metno180jet->Fill(p4probe.Pt(), pass_l1met90 ? 0 : (!pass_l1pt180 ? 0 : 1));
	    p1pre120tau->Fill(p4probe.Pt(), pass_l1pt120tau ? 0 : 1);
	    p1pre180jet->Fill(p4probe.Pt(), pass_l1pt180 ? 0 : 1);
	    p1pre120and180->Fill(p4probe.Pt(),
				 !pass_l1pt120tau && !pass_l1pt180 ? 1 : 0);
	    p1prebx1->Fill(p4probe.Pt(), pass_bx1 ? 0 : 1);
	  }
	}

	// Checks of sample constitution
	if (fail_tnp) {

	  hptl1tau_tag->Fill(p4l1tauprobe.Pt());
	  hptl1jet_tag->Fill(p4l1jetprobe.Pt());

	  hptl1met_tag->Fill(p4l1met.Pt());
	  if (has_l1jetmatch) hptl1met_tag_wtagjet->Fill(p4l1met.Pt());
	  else                hptl1met_tag_notagjet->Fill(p4l1met.Pt());

	  hptl1metnotag_tag->Fill(p4l1metnotag.Pt());
	  if (hasl1tag) hptl1metnotag_tag_wtagjet->Fill(p4l1metnotag.Pt());
	  else          hptl1metnotag_tag_notagjet->Fill(p4l1metnotag.Pt());
	}
	
	// Two-leg (ditau, met) inefficiencies from 
	if (pass_tnp) {

	  // Controls of tag and probe pT, eta distributions
	  if (fabs(p4tag.Eta())<1.3)   hpttag->Fill(p4tag.Pt());
	  if (fabs(p4probe.Eta())<1.3) hptprobe->Fill(p4probe.Pt());
	  if (p4tag.Pt()>1200.)        hetatag->Fill(p4tag.Eta());
	  if (p4probe.Pt()>1200.)      hetaprobe->Fill(p4probe.Eta());

	  hptl1tau->Fill(p4l1tauprobe.Pt());
	  if (L1_UnprefireableEvent) hptl1tau_unp->Fill(p4l1tauprobe.Pt());
	  if (L1_UnprefireableBX1)   hptl1tau_bx1->Fill(p4l1tauprobe.Pt());
	  if (L1_Prefireable)        hptl1tau_pre->Fill(p4l1tauprobe.Pt());
	  
	  hptl1jet->Fill(p4l1jetprobe.Pt());
	  if (L1_UnprefireableEvent) hptl1jet_unp->Fill(p4l1jetprobe.Pt());
	  if (L1_UnprefireableBX1)   hptl1jet_bx1->Fill(p4l1jetprobe.Pt());
	  if (L1_Prefireable)        hptl1jet_pre->Fill(p4l1jetprobe.Pt());
	  
	  if (has_l1jetmatch)
	    hptl1jet_wtagjet->Fill(p4l1jetprobe.Pt());
	  else
	    hptl1jet_notagjet->Fill(p4l1jetprobe.Pt());

	  // Find L1Jet matching probe jet in BX0 also, if any
	  bool has_l1jetbx0 = false;
	  p4l1jetbx0.SetPtEtaPhiM(0,Jet_eta[iprobe],Jet_phi[iprobe],0);
	  for (int i = 0; i != nL1Jet && !has_l1jetbx0; ++i) {
	    if (L1Jet_bx[i]==0) {
	      p4l1bx0.SetPtEtaPhiM(L1Jet_pt[i],L1Jet_eta[i],L1Jet_phi[i],0.);
	      //L1Jet_mass[i);
	      if (p4l1bx0.DeltaR(p4probe)<0.4) {
		p4l1jetbx0 = p4l1bx0;
		has_l1jetbx0 = true;
	      }
	    }
	  } // for i (probe)
	  p2prebx0->Fill(p4l1jetbx0.Eta(),p4l1jetbx0.Pt(), pass_l1pt26 ? 0 : 1);
	  if (fabs(p4l1jetbx0.Eta())<1.3) {
	    p1prebx0->Fill(p4l1jetbx0.Pt(), pass_l1pt26 ? 0 : 1);
	    h2jetbx0bx1->Fill(p4l1jetbx0.Pt(), p4l1jetprobe.Pt());
	  }

	  hptl1metvsjet->Fill(p4l1jetprobe.Pt(),p4l1met.Pt());
	  hptl1metvstau->Fill(p4l1tauprobe.Pt(),p4l1met.Pt());
	  hptl1tauvsjet->Fill(p4l1jetprobe.Pt(),p4l1tauprobe.Pt());
	  hptl1jetvsoff->Fill(p4probe.Pt(),p4l1jetprobe.Pt());
	  hptl1tauvsoff->Fill(p4probe.Pt(),p4l1tauprobe.Pt());
	  
	  if (!has_l1jetmatch)
	    hptl1metvsjet_notagjet->Fill(p4l1jetprobe.Pt(),p4l1met.Pt());
	  
	  hptl1met->Fill(p4l1met.Pt());
	  if (has_l1jetmatch) hptl1met_wtagjet->Fill(p4l1met.Pt());
	  else                hptl1met_notagjet->Fill(p4l1met.Pt());
	  if (L1_UnprefireableEvent) {
	    hptl1met_unp->Fill(p4l1met.Pt());
	    if (has_l1jetmatch) hptl1met_unp_wtagjet->Fill(p4l1met.Pt());
	    else                hptl1met_unp_notagjet->Fill(p4l1met.Pt());
	  }
	  if (L1_UnprefireableBX1) {
	    hptl1met_bx1->Fill(p4l1met.Pt());
	    if (has_l1jetmatch) hptl1met_bx1_wtagjet->Fill(p4l1met.Pt());
	    else                hptl1met_bx1_notagjet->Fill(p4l1met.Pt());
	  }
	  if (L1_Prefireable) {
	    hptl1met_pre->Fill(p4l1met.Pt());
	    if (has_l1jetmatch) hptl1met_pre_wtagjet->Fill(p4l1met.Pt());
	    else                hptl1met_pre_notagjet->Fill(p4l1met.Pt());
	  }
	  
	  hptl1metnotag->Fill(p4l1metnotag.Pt());
	  if (hasl1tag) hptl1metnotag_wtagjet->Fill(p4l1metnotag.Pt());
	  else          hptl1metnotag_notagjet->Fill(p4l1metnotag.Pt());
	  if (L1_UnprefireableEvent) {
	    hptl1metnotag_unp->Fill(p4l1metnotag.Pt());
	    if (hasl1tag) hptl1metnotag_unp_wtagjet->Fill(p4l1metnotag.Pt());
	    else          hptl1metnotag_unp_notagjet->Fill(p4l1metnotag.Pt());
	  }
	  if (L1_UnprefireableBX1) {
	    hptl1metnotag_bx1->Fill(p4l1metnotag.Pt());
	    if (hasl1tag) hptl1metnotag_bx1_wtagjet->Fill(p4l1metnotag.Pt());
	    else          hptl1metnotag_bx1_notagjet->Fill(p4l1metnotag.Pt());
	  }
	  if (L1_Prefireable) {
	    hptl1metnotag_pre->Fill(p4l1metnotag.Pt());
	    if (hasl1tag) hptl1metnotag_pre_wtagjet->Fill(p4l1metnotag.Pt());
	    else          hptl1metnotag_pre_notagjet->Fill(p4l1metnotag.Pt());
	  }
	  
	  // Hybrid scheme to recover MET>90 with shift and double counting
	  //hptl1methybrid->Fill(p4l1met.Pt());
	  //if (p4l1metnotag.Pt()>90.) hptl1methybrid->Fill(p4l1met.Pt());
	  //if (L1_UnprefireableEvent) {
	  //hptl1methybrid_unp->Fill(p4l1met.Pt());
	  //if (p4l1metnotag.Pt()>90.) hptl1methybrid_unp->Fill(p4l1met.Pt());
	  //}
	  
	  
	  // Maps for initial testing
	  p2eff->Fill(p4probe.Eta(), p4probe.Pt(), pass_l1pt34 ? 1 : 0);
	  p2effw->Fill(p4probe.Eta(), p4probe.Pt(), pass_l1pt34 ? 1 : 0);
	  p2effw2->Fill(fabs(p4probe.Eta()), p4probe.Pt(), pass_l1pt34 ? 1 : 0);
	  p2pre->Fill(p4probe.Eta(), p4probe.Pt(), pass_l1pt34 ? 0 : 1);
	  p2prew->Fill(p4probe.Eta(), p4probe.Pt(), pass_l1pt34 ? 0 : 1);
	  p2prew2->Fill(fabs(p4probe.Eta()), p4probe.Pt(), pass_l1pt34 ? 0 : 1);

	  // Final prefiring maps
	  p2pre75t->Fill(eta,p4probe.Pt(),pass_l1pt75t ? 0 : 1);
	  p2pre77p5t->Fill(eta,p4probe.Pt(),
			   pass_l1pt77p5t ? 0 : 1);
	  p2pre80t->Fill(eta,p4probe.Pt(),pass_l1pt80t ? 0 : 1);
	  p2pre34->Fill(eta, p4probe.Pt(), pass_l1pt34 ? 0 : 1);
	  p2pre32->Fill(eta, p4probe.Pt(), pass_l1pt32 ? 0 : 1);
	  p2pre32to34->Fill(eta, p4probe.Pt(),
			    pass_l1pt32to34 ? 0 : 1);
	  p2pre26->Fill(eta, p4probe.Pt(), pass_l1pt26 ? 0 : 1);
	  p2pretau->Fill(eta, p4probe.Pt(), pass_l1tau ? 0 : 1);

	  //p2pre180->Fill(eta,p4probe.Pt(),pass_l1pt180 ? 0 : 1);
	  p2pre140->Fill(eta,p4probe.Pt(),pass_l1pt140 ? 0 : 1);
	  p2pre130->Fill(eta,p4probe.Pt(),pass_l1pt130 ? 0 : 1);
	  p2pre120->Fill(eta,p4probe.Pt(),pass_l1pt120 ? 0 : 1);
	  p2pre110->Fill(eta,p4probe.Pt(),pass_l1pt110 ? 0 : 1);
	  p2pre90->Fill(eta,p4probe.Pt(),pass_l1pt90 ? 0 : 1);
	  //p2pre90h->Fill(eta,p4probe.Pt(),pass_l1met90 ? 0 : 1);
	  //if (!pass_l1met90h) p2pre90h->Fill(eta, p4probe.Pt(), 0);
	  if (has_l1jetmatch) p2pre90h->Fill(eta,p4probe.Pt(),
					     pass_l1metnotag90 ? 0 : 1);
	  if (has_l1jetmatch)
	    p2pre90_wtagjet->Fill(eta,p4probe.Pt(),pass_l1pt90 ? 0 : 1);
	  else
	    p2pre90_notagjet->Fill(eta,p4probe.Pt(),pass_l1pt90 ? 0 : 1);
	  p2pre80m->Fill(eta,p4probe.Pt(),pass_l1met80 ? 0 : 1);
	  p2pre80->Fill(eta,p4probe.Pt(),pass_l1pt80 ? 0 : 1);
	  p2pre70m->Fill(eta,p4probe.Pt(),pass_l1met70 ? 0 : 1);
	  p2pre70->Fill(eta,p4probe.Pt(),pass_l1pt70 ? 0 : 1);
	  p2pre60m->Fill(eta,p4probe.Pt(),pass_l1met60 ? 0 : 1);
	  p2pre55->Fill(eta,p4probe.Pt(),pass_l1pt55 ? 0 : 1);
	  p2prejet->Fill(eta,p4probe.Pt(),pass_l1jet ? 0 : 1);
	  
	  if (fabs(p4probe.Eta())<1.3) {
	    // Maps for initial testing
	    p1eff->Fill(p4probe.Pt(), pass_l1pt34 ? 1 : 0);
	    p1effw->Fill(p4probe.Pt(), pass_l1pt34 ? 1 : 0);
	    p1pre->Fill(p4probe.Pt(), pass_l1pt34 ? 0 : 1);
	    p1prew->Fill(p4probe.Pt(), pass_l1pt34 ? 0 : 1);

	    // Final prefiring maps
	    p1pre75t->Fill(p4probe.Pt(), pass_l1pt75t ? 0 : 1);
	    p1pre77p5t->Fill(p4probe.Pt(), pass_l1pt77p5t ? 0 : 1);
	    p1pre80t->Fill(p4probe.Pt(), pass_l1pt80t ? 0 : 1);
	    p1pre34->Fill(p4probe.Pt(), pass_l1pt34 ? 0 : 1);
	    p1pre32->Fill(p4probe.Pt(), pass_l1pt32 ? 0 : 1);
	    p1pre32to34->Fill(p4probe.Pt(), pass_l1pt32to34 ? 0 : 1);
	    p1pre26->Fill(p4probe.Pt(), pass_l1pt26 ? 0 : 1);
	    p1pretau->Fill(p4probe.Pt(), pass_l1tau ? 0 : 1);

	    //p1pre180->Fill(p4probe.Pt(), pass_l1pt180 ? 0 : 1);
	    p1pre140->Fill(p4probe.Pt(), pass_l1pt140 ? 0 : 1);
	    p1pre130->Fill(p4probe.Pt(), pass_l1pt130 ? 0 : 1);
	    p1pre120->Fill(p4probe.Pt(), pass_l1pt120 ? 0 : 1);
	    p1pre110->Fill(p4probe.Pt(), pass_l1pt110 ? 0 : 1);
	    p1pre90->Fill(p4probe.Pt(), pass_l1pt90 ? 0 : 1);
	    //p1pre90h->Fill(p4probe.Pt(), pass_l1met90 ? 0 : 1);
	    //if (!pass_l1met90h) p1pre90h->Fill(p4probe.Pt(), 0);
	    if (has_l1jetmatch) p1pre90h->Fill(p4probe.Pt(),
					       pass_l1metnotag90 ? 0 : 1);
	    if (has_l1jetmatch)
	      p1pre90_wtagjet->Fill(p4probe.Pt(), pass_l1pt90 ? 0 : 1);
	    else
	      p1pre90_notagjet->Fill(p4probe.Pt(), pass_l1pt90 ? 0 : 1);
	    p1pre80m->Fill(p4probe.Pt(), pass_l1met80 ? 0 : 1);
	    p1pre80->Fill(p4probe.Pt(), pass_l1pt80 ? 0 : 1);
	    p1pre70m->Fill(p4probe.Pt(), pass_l1met70 ? 0 : 1);
	    p1pre70->Fill(p4probe.Pt(), pass_l1pt70 ? 0 : 1);
	    p1pre60m->Fill(p4probe.Pt(), pass_l1met60 ? 0 : 1);
	    p1pre55->Fill(p4probe.Pt(), pass_l1pt55 ? 0 : 1);
	    p1prejet->Fill(p4probe.Pt(), pass_l1jet ? 0 : 1);
	  }
	} // pass tag-and-probe selection

	// Use only one permutation for Mjj and inclusive 2-jet
	if (itag==0 && pass_jetid) {

	  // Calculate dijet mass
	  p4jj = p4tag + p4probe;
	  double mjj = p4jj.M();
	  double deta = fabs(p4tag.Eta()-p4probe.Eta());

	  // Calculate ditau mass
	  p4tt = p4l1tautag + p4l1tauprobe;
	  double mtt = p4tt.M();
	  
	  // Estimate trigger efficiency (1-pre) for pairs
	  // Two-leg trigger needs both prefiring to fail:
	  //   pre_both = pre_tag * pre_probe
	  //   eff = max(0.05, 1 - pre_both)
	  // where minimum of 0.05 is set by unprefireable rate

	  // Estimate prefiring probability for 34 GeV
	  double pret = getPre(p2pre34ref,p4tag.Pt(),p4tag.Eta());
	  double prep = getPre(p2pre34ref,p4probe.Pt(),eta);
	  double eff34 = max(0.05, 1 - pret * prep);

	  // Estimate prefiring probability for 32 GeV + mass below 80 GeV 
	  double pre32t = getPre(p2pre32ref,p4tag.Pt(),p4tag.Eta());
	  double pre32p = getPre(p2pre32ref,p4probe.Pt(),eta);
	  // M^2 = 2*pT1*pT2*(cosh(eta1-eta2) - cos(phi1-phi2))
	  // Estimate upper end prefire using pT2=pT1=pTmin=32 GeV
	  double mass80 = sqrt(2*32.*32.*(cosh(p4tag.Eta()-p4probe.Eta()) -
					  cos(p4tag.Phi()-p4probe.Phi())));
	  double pre80m = (mass80<80. ? 1 : 0);
	  // Two-leg needs both prefiring, and small enough mass
	  double eff32no80m = max(0.05, 1 - pre32t * pre32p);
	  double eff32 = max(0.05, 1 - pre32t * pre32p * pre80m);

	  // Estimate extra rate of 32 GeV above 34 GeV, two permutations
	  double pre32to34t_1 = getPre(p2pre32to34ref,p4tag.Pt(),p4tag.Eta());
	  double pre32p_1 = getPre(p2pre32ref,p4probe.Pt(),eta);
	  double pre32t_2 = getPre(p2pre32ref,p4tag.Pt(),p4tag.Eta());
	  double pre32to34p_2 = getPre(p2pre32to34ref,p4probe.Pt(),eta);
	  double eff32to34 = max(0.05, 1
				 - pre32to34t_1 * pre32p_1 * pre80m
				 - pre32t_2 * pre32to34p_2 * pre80m
				 + pre32to34t_1 * pre32to34p_2 * pre80m);
	  
	  // Estimate prefiring probability for 26 GeV + jet of 55 GeV
	  // Need three permutations
	  double pre26t_1 = getPre(p2pre26ref,p4tag.Pt(),p4tag.Eta());
	  double pre26p_1 = getPre(p2pre26ref,p4probe.Pt(),eta);
	  double pre55j_1 = getPre(p2pre55ref,p4jet.Pt(),p4jet.Eta());
	  double pre26t_2 = getPre(p2pre26ref,p4tag.Pt(),p4tag.Eta());
	  double pre55p_2 = getPre(p2pre55ref,p4probe.Pt(),eta);
	  double pre26j_2 = getPre(p2pre26ref,p4jet.Pt(),p4jet.Eta());
	  double pre55t_3 = getPre(p2pre55ref,p4tag.Pt(),p4tag.Eta());
	  double pre26p_3 = getPre(p2pre26ref,p4probe.Pt(),eta);
	  double pre26j_3 = getPre(p2pre26ref,p4jet.Pt(),p4jet.Eta());
	  // Probability of firing any triplet is small, so ignore firing two
	  double eff26 = max(0.05, 1 
			     - pre26t_1 * pre26p_1 * pre55j_1
			     - pre26t_2 * pre55p_2 * pre26j_2
			     - pre55t_3 * pre26p_3 * pre26j_3
			     );

	  // Estimate prefiring for L1_ETMHF90
	  // One jet above 90, another zero; two permutations
	  // Correct pre90 for missing events without L1Jet on tag side
	  double pre90t_1 = getPre(p2pre90ref,p4tag.Pt(),p4tag.Eta());
	  double prejtp_1 = getPre(p2prejetref,p4probe.Pt(),eta);
	  double prejtt_2 = getPre(p2prejetref,p4tag.Pt(),p4tag.Eta());
	  double pre90p_2 = getPre(p2pre90ref,p4probe.Pt(),eta);

	  // Not sure if right, but seems to work (~correct shape and size)
	  //double eff90 = max(0.05, 1
	  //		     - pre90t_1/pre0jt_1*pre0jp_1
	  //		     - pre0jt_2*pre90p_2/pre0jp_2);
	  // Not working like this? V2-3
	  //double eff90 = max(0.05, 1
	  //		     - pre90t_1*(1-prejtp_1)
	  //		     - (1-prejtt_2)*pre90p_2);
	  // How about this, V3?
	  //double eff90 = max(0.05,1 - pre90t_1 - pre90p_2);
	  // MET90 estimate was matched to BX1 rate so just max of tag/probe?
	  //double eff90 = max(0.05,1 - max(pre90t_1,pre90p_2));
	  // MET90 as geometric average of tag and probe?
	  //double eff90 = max(0.05, 1 - sqrt(pre90t_1*pre90p_2));
	  // MET90 as arithmetic average of tag and probe?
	  double eff90 = max(0.05, 1 - 0.5*(pre90t_1+pre90p_2));
	
	  // Total efficiency: (jet180+tau120)+ met90 + (ditau34+ditau32to34)
	  // Ignore cross terms, because met90 and ditau determined from
	  // tag-and-probe sample that has cross-terms already removed. Their
	  // parameterizations are therefore for purely additive rate
	  double pretot = (1-eff90) + ((1-eff34)+(1-eff32to34));
	  //  - (1-eff90) * ((1-eff34)+(1-eff32to34));
	  //- (1-eff90) * ((1-eff32)+(1-eff34));
	  double efftot = max(0.05, 1 - pretot);
	  
	  
	  // Monitor ditau mass vs pT(probe)
	  // NB: could do profile of fraction with Mjj>80 GeV
	  h2mtt->Fill(p4probe.Pt(), mass80);
	  h2mttr->Fill(p4probe.Pt(), mass80/80.);
	  h2mtt32->Fill(p4probe.Pt(), mass80, eff32no80m);
	  h2mtt32r->Fill(p4probe.Pt(), mass80/80., eff32no80m);
	  
	  // Fill histograms before and after prefiring corrections
	  hmjj->Fill(mjj);
	  if (eff34>0) hmjjc->Fill(mjj, 1./eff34);
	  if (deta<1.3) {
	    hmjj13->Fill(mjj);
	    if (eff34>0) hmjj13c->Fill(mjj, 1./eff34);
	  }

	  hpt->Fill(p4tag.Pt());
	  hpt->Fill(p4probe.Pt());
	  if (eff34>0) {
	    hptc->Fill(p4tag.Pt(), 1./eff34);
	    hptc->Fill(p4probe.Pt(), 1./eff34);
	  }
	  
	  if (fabs(p4tag.Eta())  <1.3) hpt13->Fill(p4tag.Pt());
	  if (fabs(p4probe.Eta())<1.3) hpt13->Fill(p4probe.Pt());
	  if (eff34>0) {
	    if (fabs(p4tag.Eta())  <1.3) hpt13c->Fill(p4tag.Pt(), 1./eff34);
	    if (fabs(p4probe.Eta())<1.3) hpt13c->Fill(p4probe.Pt(), 1./eff34);
	  }

	  // Efficiency estimate for dijet mass
	  if (deta<1.3) {
	    // correct for biased eta-pT distribution in data with weight
	    double w = (efftot>0.05 ? 1./efftot : 1./0.05); 
	    p34mjj->Fill(mjj, 1-eff34, w);
	    p32mjj->Fill(mjj, 1-eff32, w);
	    p32to34mjj->Fill(mjj, 1-eff32to34, w);
	    p26mjj->Fill(mjj, 1-eff26, w);

	    p90mjj->Fill(mjj, 1-eff90, w);

	    ptotmjj->Fill(mjj, 1-efftot, w);

	    if (pass_single) {
	      p34mjjbx1->Fill(mjj, pass_ditau34 ? 0 : 1);
	      p32mjjbx1->Fill(mjj, pass_ditau32 ? 0 : 1);
	      p90mjjbx1->Fill(mjj, pass_l1met90 ? 0 : 1);
	      ptotmjjbx1->Fill(mjj, pass_bx1 ? 0 : 1);
	    }
	  } // deta<1.3
	} // tag=0
	
      } // itag
   } // jentry
   
   fout->Write();
   fout->Close();
   curdir->cd();
} // Loop
