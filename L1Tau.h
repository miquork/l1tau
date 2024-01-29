//////////////////////////////////////////////////////////
// This class has been automatically generated on (and then modded)
// Mon Dec 11 15:11:02 2023 by ROOT version 6.26/06
// from TTree Events/Events
// found on file: rootfiles/taus_Run2023C_PromptNanoAODv12_v3.root
//////////////////////////////////////////////////////////

#ifndef L1Tau_h
#define L1Tau_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "ROOT/RVec.hxx"
//#include "ROOT/RVec.hxx"
#include "c++/v1/vector"


class L1Tau {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   string name;
  
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
  static const int nMaxL1Tau = 100;
   Int_t           nL1Tau;
   Short_t         L1Tau_hwIso[nMaxL1Tau];   //[nL1Tau]
   Short_t         L1Tau_bx[nMaxL1Tau];   //[nL1Tau]
   Float_t         L1Tau_eta[nMaxL1Tau];   //[nL1Tau]
   Float_t         L1Tau_phi[nMaxL1Tau];   //[nL1Tau]
   Float_t         L1Tau_pt[nMaxL1Tau];   //[nL1Tau]

  static const int nMaxL1Jet = 100;
   Int_t           nL1Jet;
   Short_t         L1Jet_bx[nMaxL1Jet];   //[nL1Jet]
   Float_t         L1Jet_eta[nMaxL1Jet];   //[nL1Jet]
   Float_t         L1Jet_phi[nMaxL1Jet];   //[nL1Jet]
   Float_t         L1Jet_pt[nMaxL1Jet];   //[nL1Jet]

  static const int nMaxL1EtSum = 100;
   Int_t           nL1EtSum;
   Short_t         L1EtSum_bx[nMaxL1EtSum];   //[nL1EtSum]
   Int_t           L1EtSum_etSumType[nMaxL1EtSum];   //[nL1EtSum]
   Float_t         L1EtSum_phi[nMaxL1EtSum];   //[nL1EtSum]
   Float_t         L1EtSum_pt[nMaxL1EtSum];   //[nL1EtSum]
  
  static const int nMaxJet = 200;
   Int_t           nJet;
   Float_t         Jet_eta[nMaxJet];   //[nJet]
   Float_t         Jet_phi[nMaxJet];   //[nJet]
  //Float_t         Jet_pt[nMaxJet];   //[nJet]
  // Branch with newly re-calculated JEC.
  vector<float>   *Jet_pt;
   Float_t         Jet_mass[nMaxJet];   //[nJet]
   UChar_t         Jet_jetId[nMaxJet];   //[nJet]
   UChar_t         Jet_nConstituents[nMaxJet];   //[nJet]
  //Float_t         Jet_rawFactor[nMaxJet];   //[nJet]
   vector<float>   *Jet_rawFactor;
  
   Bool_t          HLT_PFJet500;
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing;

   Bool_t          Flag_Run3;
   Bool_t          L1_UnprefireableEvent;
   Bool_t          L1_UnprefireableBX1;
   Bool_t          L1_Prefireable;
  
   // List of branches
   TBranch        *b_nL1Tau;   //!
   TBranch        *b_L1Tau_hwIso;   //!
   TBranch        *b_L1Tau_bx;   //!
   TBranch        *b_L1Tau_eta;   //!
   TBranch        *b_L1Tau_phi;   //!
   TBranch        *b_L1Tau_pt;   //!

   TBranch        *b_nL1Jet;   //!
   TBranch        *b_L1Jet_bx;   //!
   TBranch        *b_L1Jet_eta;   //!
   TBranch        *b_L1Jet_phi;   //!
   TBranch        *b_L1Jet_pt;   //!

   TBranch        *b_nL1EtSum;   //!
   TBranch        *b_L1EtSum_bx;   //!
   TBranch        *b_L1EtSum_etSumType;   //!
   TBranch        *b_L1EtSum_phi;   //!
   TBranch        *b_L1EtSum_pt;   //!
  
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_rawFactor;   //!

   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!

   TBranch        *b_Flag_Run3;   //!
   TBranch        *b_L1_UnprefireableEvent;   //!

  L1Tau(TTree *tree=0, string name="X");
   virtual ~L1Tau();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef L1Tau_cxx
L1Tau::L1Tau(TTree *tree, string _name) : fChain(0), name(_name)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfiles/taus_Run2023C_PromptNanoAODv12_v3.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootfiles/taus_Run2023C_PromptNanoAODv12_v3.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

L1Tau::~L1Tau()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t L1Tau::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t L1Tau::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void L1Tau::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jet_pt = 0;
  
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nL1Tau", &nL1Tau, &b_nL1Tau);
   fChain->SetBranchAddress("L1Tau_hwIso", L1Tau_hwIso, &b_L1Tau_hwIso);
   fChain->SetBranchAddress("L1Tau_bx", L1Tau_bx, &b_L1Tau_bx);
   fChain->SetBranchAddress("L1Tau_eta", L1Tau_eta, &b_L1Tau_eta);
   fChain->SetBranchAddress("L1Tau_phi", L1Tau_phi, &b_L1Tau_phi);
   fChain->SetBranchAddress("L1Tau_pt", L1Tau_pt, &b_L1Tau_pt);

   fChain->SetBranchAddress("nL1Jet", &nL1Jet, &b_nL1Jet);
   fChain->SetBranchAddress("L1Jet_bx", L1Jet_bx, &b_L1Jet_bx);
   fChain->SetBranchAddress("L1Jet_eta", L1Jet_eta, &b_L1Jet_eta);
   fChain->SetBranchAddress("L1Jet_phi", L1Jet_phi, &b_L1Jet_phi);
   fChain->SetBranchAddress("L1Jet_pt", L1Jet_pt, &b_L1Jet_pt);

   fChain->SetBranchAddress("nL1EtSum", &nL1EtSum, &b_nL1EtSum);
   fChain->SetBranchAddress("L1EtSum_bx", L1EtSum_bx, &b_L1EtSum_bx);
   fChain->SetBranchAddress("L1EtSum_etSumType", L1EtSum_etSumType, &b_L1EtSum_etSumType);
   fChain->SetBranchAddress("L1EtSum_phi", L1EtSum_phi, &b_L1EtSum_phi);
   fChain->SetBranchAddress("L1EtSum_pt", L1EtSum_pt, &b_L1EtSum_pt);
   
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   //fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt); // old JEC
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt); // new JEC
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   //fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_rawFactor", &Jet_rawFactor, &b_Jet_rawFactor);
   
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);

   fChain->SetBranchAddress("Flag_Run3", &Flag_Run3, &b_Flag_Run3);
   fChain->SetBranchAddress("L1_UnprefireableEvent", &L1_UnprefireableEvent, &b_L1_UnprefireableEvent);

   
   Notify();
}

Bool_t L1Tau::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void L1Tau::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t L1Tau::Cut(Long64_t entry)
{
  if (entry) {}; // suppress warning
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef L1Tau_cxx
