// Compare 19Dec2023 vs 22Sep2023
// Only load a few branches (event ID, jet pt, eta, phi)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "ROOT/RVec.hxx"
//#include "ROOT/RVec.hxx"
#include "c++/v1/vector"


#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStopwatch.h"

//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

//#include "compareLite/settings.h"
#include "tools.h"
#include "../tdrstyle_mod22.C"

#include <iostream>

using namespace std;
using namespace tools;

class evtid {
private:
  UInt_t run_, ls_;
  ULong64_t evt_;
public :
  evtid() : run_(0), ls_(0), evt_(0) {}
  evtid(UInt_t run, UInt_t ls, ULong64_t evt) : run_(run), ls_(ls), evt_(evt) {}
  bool operator()(evtid const& a, evtid const& b) const {
    if (a.run_ < b.run_) return true;
    if (a.run_ > b.run_) return false;
    if (a.ls_  < b.ls_)  return true;
    if (a.ls_  > b.ls_)  return false;
    return (a.evt_ < b.evt_);
  }
  UInt_t run() const { return run_; }
  UInt_t lbn() const { return ls_; }
  ULong64_t evt() const { return evt_; }
};

int nparl1_ = 9;
Double_t funcL1(Double_t *x, Double_t *p) {

  double val = 0;
  double y = fabs(x[0]);
  if (y<0.6)           val = p[0] + (p[1]-p[0])*(y-0.0)/(0.6-0.0);
  if (y>=0.6 && y<1.3) val = p[1] + (p[2]-p[1])*(y-0.6)/(1.3-0.6);
  if (y>=1.3 && y<2.4) val = p[2] + (p[3]-p[2])*(y-1.3)/(2.5-1.3);
  if (y>=2.4 && y<2.6) val = p[3] + (p[4]-p[3])*(y-2.4)/(2.6-2.4);
  if (y>=2.6 && y<2.8) val = p[4] + (p[5]-p[4])*(y-2.6)/(2.8-2.6);
  if (y>=2.8 && y<3.0) val = p[5] + (p[6]-p[5])*(y-2.8)/(3.0-2.8);
  if (y>=3.0 && y<3.2) val = p[6] + (p[7]-p[6])*(y-3.0)/(3.2-3.0);
  if (y>=3.2 && y<3.4) val = p[7] + (p[8]-p[7])*(y-3.2)/(3.4-3.2);
  if (y>=3.4)          val = p[8];

  return val;
}

void compareLite(string run="2023D") {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  // Book tree tA (19Dec2023)
  TChain *c_tA = new TChain("Events");
  cout << "A is 19Dec2023" << endl;
  if (run=="2023Cv123")
    c_tA->AddFile("../data/l1tau/19Dec2023/Run2023C_taus.root");
  if (run=="2023D")
    c_tA->AddFile("../data/l1tau/19Dec2023/Run2023D_taus.root");
  
  // Set branches to sort events
  TBranch *b_run_tA, *b_lbn_tA, *b_evt_tA;
  UInt_t run_tA, lbn_tA;
  ULong64_t evt_tA;
  c_tA->SetBranchAddress("run",&run_tA,&b_run_tA);
  c_tA->SetBranchAddress("luminosityBlock",&lbn_tA,&b_lbn_tA);
  c_tA->SetBranchAddress("event",&evt_tA,&b_evt_tA);

  cout << "Sort TA entries" << endl << flush;
  map<evtid, pair<Long64_t, Long64_t>, evtid> mtA;
  Long64_t nentries = c_tA->GetEntries();//Fast();
  cout << "..Processing " << nentries << " entries" << endl << flush;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = c_tA->LoadTree(jentry);
    if (ientry < 0) break;
    //nb = c_tA->GetEntry(jentry);   nbytes += nb;
    b_run_tA->GetEntry(ientry);
    b_lbn_tA->GetEntry(ientry);
    b_evt_tA->GetEntry(ientry);
    assert(run_tA);
    assert(lbn_tA);
    assert(evt_tA);

    mtA[evtid(run_tA, lbn_tA, evt_tA)]
      = pair<Long64_t, Long64_t>(jentry, ientry);

    if (jentry%1000000==0) cout << "." << flush;
  }
  cout << endl;
  cout << "Found " << mtA.size() << " unique entries" << endl;
  cout << endl;


  // Book TB tree (22Sep2023)
  TChain *c_tB = new TChain("Events");
  //cout << "B is 22Sep2023" << endl;
  cout << "B is Prompt23" << endl;
  if (run=="2023Cv123") {
    c_tB->Add("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv11p9_v1-v1_NANOAOD_taus.root");
    c_tB->Add("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v2-v2_NANOAOD_taus.root");
    c_tB->Add("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v2-v4_NANOAOD_taus.root");
    c_tB->AddFile("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v3-v1_NANOAOD_taus.root");
  }
  if (run=="2023D") {
    c_tB->AddFile("../data/l1tau/Summer22Prompt23/Run2023D-PromptReco-v1_NANOAOD_taus.root");
    c_tB->AddFile("../data/l1tau/Summer22Prompt23/Run2023D-PromptReco-v2_NANOAOD_taus.root");
  }

  // Set branches to sort events
  TBranch *b_run_tB, *b_lbn_tB, *b_evt_tB;
  UInt_t run_tB, lbn_tB;
  ULong64_t evt_tB;
  c_tB->SetBranchAddress("run",&run_tB,&b_run_tB);
  c_tB->SetBranchAddress("luminosityBlock",&lbn_tB,&b_lbn_tB);
  c_tB->SetBranchAddress("event",&evt_tB,&b_evt_tB);

  cout << "Sort TB entries" << endl << flush;
  map<evtid, pair<Long64_t, Long64_t>, evtid> mtB;
  nentries = c_tB->GetEntries();//Fast();
  cout << "..Processing " << nentries << " entries" << endl << flush;

  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = c_tB->LoadTree(jentry);
    if (ientry < 0) break;
    //nb = c_tB->GetEntry(jentry);   nbytes += nb;
    b_run_tB->GetEntry(ientry);
    b_lbn_tB->GetEntry(ientry);
    b_evt_tB->GetEntry(ientry);
    assert(run_tB);
    assert(lbn_tB);
    assert(evt_tB);

    mtB[evtid(run_tB, lbn_tB, evt_tB)]
      = pair<Long64_t, Long64_t>(jentry, ientry);

    if (jentry%1000000==0) cout << "." << flush;
  }
  cout << endl;
  cout << "Found " << mtB.size() << " unique entries" << endl;
  cout << endl;

  // Use sorted events to find matching ones
  cout << "Match TA to TB events" << endl;
  map<Long64_t, Long64_t> mAtoB;
  Long64_t ntotA = mtA.size();
  int failsAtoB(0);
  typedef map<evtid, pair<Long64_t, Long64_t> > IT;
  for (IT::const_iterator it = mtA.begin(); it != mtA.end(); ++it) {

    if (mtB.find(it->first)!=mtB.end()) {
      assert(mAtoB.find(it->second.first)==mAtoB.end());
      mAtoB[it->second.first] = mtB[it->first].first;
    }
    else if (++failsAtoB<10) {
      
      cout << "For tA("<<it->first.run()<<","<<it->first.lbn()<<","
	   <<it->first.evt()<<"), did not find a matching entry in tB" << endl;
    }
  }
  cout << "Found " << mAtoB.size() << " matching entries" << endl;

  map<Long64_t, Long64_t> mBtoA;
  Long64_t ntotB = mtB.size();
  int failsBtoA(0);
  for (IT::const_iterator it = mtB.begin(); it != mtB.end(); ++it) {

    if (mtA.find(it->first)!=mtA.end()) {
      assert(mBtoA.find(it->second.first)==mBtoA.end());
      mBtoA[it->second.first] = mtA[it->first].first;
    }
    else if (++failsBtoA<10) {
      
      cout << "For tB("<<it->first.run()<<","<<it->first.lbn()<<","
	   <<it->first.evt()<<"), did not find a matching entry in tA" << endl;
    }
  }
  cout << "Found " << mBtoA.size() << " matching entries" << endl;

  // Set branches needed from trees
  TBranch *b_flag_tA, *b_flag_tB;
  TBranch *b_jttrg_tA, *b_jttrg_tB;
  TBranch *b_njt_tA, *b_njt_tB;
  TBranch *b_jtpt_tA, *b_jtpt_tB;
  TBranch *b_jteta_tA, *b_jteta_tB;
  TBranch *b_jtphi_tA, *b_jtphi_tB;
  TBranch *b_jtidtight_tA, *b_jtidtight_tB;
  TBranch *b_jtjes_tA, *b_jtjes_tB;

  Bool_t flag_tA, flag_tB;
  Bool_t jttrg_tA, jttrg_tB;
  Int_t njt_tA, njt_tB;
  //Int_t npv_tA(0), npv_tB(0);
  //
  const int njt = 100;
  //Float_t jtpt_tA[njt], jtpt_tB[njt];
  vector<float> *jtpt_tA, *jtpt_tB;
  Float_t jteta_tA[njt], jteta_tB[njt];
  Float_t jtphi_tA[njt], jtphi_tB[njt];
  UChar_t jtidtight_tA[njt], jtidtight_tB[njt];
  //Float_t jtjes_tA[njt], jtjes_tB[njt];
  vector<float> *jtjes_tA, *jtjes_tB;


  c_tA->SetBranchAddress("Flag_Run3",&flag_tA,&b_flag_tA);
  c_tB->SetBranchAddress("Flag_Run3",&flag_tB,&b_flag_tB);
  c_tA->SetBranchAddress("HLT_PFJet500",&jttrg_tA,&b_jttrg_tA);
  c_tB->SetBranchAddress("HLT_PFJet500",&jttrg_tB,&b_jttrg_tB);
  //
  c_tA->SetBranchAddress("nJet",&njt_tA,&b_njt_tA);
  c_tB->SetBranchAddress("nJet",&njt_tB,&b_njt_tB);
  jtpt_tA = 0, jtpt_tB = 0; // Avoid crash in GetEntry
  c_tA->SetBranchAddress("Jet_pt",&jtpt_tA,&b_jtpt_tA);
  c_tB->SetBranchAddress("Jet_pt",&jtpt_tB,&b_jtpt_tB);
  c_tA->SetBranchAddress("Jet_eta",jteta_tA,&b_jteta_tA);
  c_tB->SetBranchAddress("Jet_eta",jteta_tB,&b_jteta_tB);
  c_tA->SetBranchAddress("Jet_phi",jtphi_tA,&b_jtphi_tA);
  c_tB->SetBranchAddress("Jet_phi",jtphi_tB,&b_jtphi_tB);
  c_tA->SetBranchAddress("Jet_jetId",jtidtight_tA,&b_jtidtight_tA);
  c_tB->SetBranchAddress("Jet_jetId",jtidtight_tB,&b_jtidtight_tB);
  jtjes_tA = 0, jtjes_tB = 0; // Avoid crash in GetEntry
  c_tA->SetBranchAddress("Jet_rawFactor",&jtjes_tA,&b_jtjes_tA);
  c_tB->SetBranchAddress("Jet_rawFactor",&jtjes_tB,&b_jtjes_tB);

  const int nsample = 1;//100; // 0.5h
  const float frac = 1;//0.1;//0.01;
  cout << "Pairing TA and TB" << endl << flush;
  cout << "Sampling 1/"<<nsample<<" of events" << endl << flush;
  cout << "Keeping first "<<frac*100<<"% of events" << endl << flush;

  c_tA->SetBranchStatus("*",0);  // disable all branches
  c_tA->SetBranchStatus("Flag_Run3",1);
  c_tA->SetBranchStatus("run",1);
  c_tA->SetBranchStatus("luminosityBlock",1);
  c_tA->SetBranchStatus("event",1);
  c_tA->SetBranchStatus("HLT_PFJet500",1);
  //
  c_tA->SetBranchStatus("nJet",1);
  c_tA->SetBranchStatus("Jet_pt",1);
  c_tA->SetBranchStatus("Jet_eta",1);
  c_tA->SetBranchStatus("Jet_phi",1);
  c_tA->SetBranchStatus("Jet_jetId",1);
  c_tA->SetBranchStatus("Jet_rawFactor",1);

  c_tB->SetBranchStatus("*",0);  // disable all branches
  c_tB->SetBranchStatus("Flag_Run3",1);
  c_tB->SetBranchStatus("run",1);
  c_tB->SetBranchStatus("luminosityBlock",1);
  c_tB->SetBranchStatus("event",1);
  c_tB->SetBranchStatus("HLT_PFJet500",1);
  //
  c_tB->SetBranchStatus("nJet",1);
  c_tB->SetBranchStatus("Jet_pt",1);
  c_tB->SetBranchStatus("Jet_eta",1);
  c_tB->SetBranchStatus("Jet_phi",1);
  c_tB->SetBranchStatus("Jet_jetId",1);
  c_tB->SetBranchStatus("Jet_rawFactor",1);

  // Results of interest
  // pT binning from JEC?
  /*
  double vx[] =
    {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967,
     1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
     2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450};
  int nx = sizeof(vx)/sizeof(vx[0])-1;
  double vxw[] =
      {10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468,
     548, 638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3450};
  int nxw = sizeof(vxw)/sizeof(vxw[0])-1;
  */
  double vx[] =
    {10, 15, 21, 28, 37, 49, 64, 84, 114, 153, 196, 245, 300, 395, 468,
     548, 638, 790, 967, 1172, 1410, 1684, 2000, 2500, 3000, 3500,
     4000, 4500, 5000, 6000};
  int nx = sizeof(vx)/sizeof(vx[0])-1;
  
  // Open file for outputting results
  TFile *f = new TFile(Form("compareLite/compareLite_%s.root",run.c_str()),
		       "RECREATE");

  // Tag-and-probe method
  TProfile *pta_tp = new TProfile("pta_tp",";p_{T,tag};p_{T,A}",nx,vx);
  TProfile *pa_tp = new TProfile("pa_tp",";p_{T,tag};p_{T,A}/p_{T,tag}",nx,vx);
  TProfile *pb_tp = new TProfile("pb_tp",";p_{T,tag};p_{T,B}/p_{T,tag}",nx,vx);
  TProfile *pd_tp = new TProfile("pd_tp",";p_{T,tag};"
				 "(p_{T,B}-p_{T,A})/p_{T,tag}",nx,vx);
  TH2D *h2a_tp =new TH2D("h2a_tp",";p_{T,tag};p_{T,A}/p_{T,tag}",nx,vx,400,-1,3);
  TH2D *h2b_tp =new TH2D("h2b_tp",";p_{T,tag};p_{T,B}/p_{T,tag}",nx,vx,400,-1,3);
  TH2D *h2d_tp =new TH2D("h2d_tp",";p_{T,tag};"
			 "0.5*(p_{T,B}-p_{T,A})/p_{T,tag}",nx,vx,600,-3,3);
  
  // Direct match method in two variants
  TProfile *pta_dm = new TProfile("pta_dm",";p_{T,A};p_{T,A}",nx,vx);
  TProfile *pa_dm = new TProfile("pa_dm",";p_{T,A};p_{T,B}/p_{T,A}",nx,vx);
  TProfile *ptb_dm = new TProfile("ptb_dm",";p_{T,B};p_{T,A}",nx,vx);
  TProfile *pb_dm = new TProfile("pb_dm",";p_{T,B};p_{T,A}/p_{T,B}",nx,vx);
  TProfile *ptd_dm = new TProfile("ptd_dm",";(p_{T,B}+p_{T,A})/2;"
				  "p_{T,A}",nx,vx);
  TProfile *pd_dm = new TProfile("pd_dm",";(p_{T,B}+p_{T,A})/2;"
				 "(p_{T,B}-p_{T,A})/(p_{T,B}+p_{T,A})",nx,vx);
  TH2D *h2a_dm = new TH2D("h2a_dm",";p_{T,A};p_{T,B}/p_{T,A}",nx,vx,400,-1,3);
  TH2D *h2b_dm = new TH2D("h2b_dm",";p_{T,B};p_{T,A}/p_{T,B}",nx,vx,400,-1,3);
  TH2D *h2d_dm = new TH2D("h2d_dm",";(p_{T,B}+p_{T,A})/2;"
			  "(p_{T,B}-p_{T,A})/(p_{T,B}+p_{T,A})",nx,vx,600,-3,3);

  curdir->cd();
  
  TStopwatch t;
  t.Start();

  int nev = 0;
  int nmatch = 0;
  int nj = 0;
  for (map<Long64_t,Long64_t>::const_iterator it = mAtoB.begin();
       it != mAtoB.end(); ++it) {

    if (++nev%10000==0) cout << "." << flush;
    if (nev%nsample!=0) continue;
    if (nev>frac*ntotA) continue;

    if (nev==1 || nev==10000 || nev==100000 || nev==1000000 || nev==5000000){
      cout << endl
	   << Form("Processed %ld events (%1.1f%%) in %1.0f sec. ETA:",
		   (long int)nev, 100.*nev/ntotA,
		   t.RealTime()) << endl;
      TDatime now; now.Set(now.Convert()+t.RealTime()*ntotA/nev);
      now.Print();
      t.Continue();
    }

    ++nmatch;
    
    // Load matching entries
    Long64_t jentrytA = it->first;
    Long64_t jentrytB = it->second;
    
    if (jentrytA<0 || jentrytA>=ntotA) continue;
    Long64_t ientrytA = c_tA->LoadTree(jentrytA);
    if (ientrytA < 0) break;
    c_tA->GetEntry(jentrytA);
    if (jentrytB<0 || jentrytB>=ntotB) continue;
    Long64_t ientrytB = c_tB->LoadTree(jentrytB);
    if (ientrytA < 0) break;
    c_tB->GetEntry(jentrytB);
    
    
    // Loop over two leading jets to find probe pairs
    for (int i = 0; i != min(2,int(njt_tA)); ++i) {
      
      double pt = (*jtpt_tA)[i];
      double eta = jteta_tA[i];
      double phi = jtphi_tA[i];
      double jes = (1-(*jtjes_tA)[i]);
      
      bool hasmatch = false;
      for (int j = 0; j != min(2,int(njt_tB)) && !hasmatch; ++j) {
	
	double ptB = (*jtpt_tB)[j];
	double etaB = jteta_tB[j];
	double phiB = jtphi_tB[j];
	double jesB = (1-(*jtjes_tB)[j]);
	
	// Match probe jets with deltaR<R/cone2
	double dr = tools::oplus(delta_eta(eta,etaB),delta_phi(phi,phiB));
	if (dr < 0.20) {
	  
	  hasmatch = true;
	  double ptave = 0.5 * (ptB + pt);

	  // Tag selection
	  bool istp(false);
	  double pttag(0);
	  if (i<2 && j<2 && njt_tA>1 && njt_tB>1) {

	    // Tag is the other one of the two leading jets
	    int k = (i==0 ? 1 : 0);
	    int l = (j==0 ? 1 : 0);
	    
	    double pttagA = (*jtpt_tA)[k];
	    double pttagB = (*jtpt_tB)[l];
	    pttag = 0.5 * (pttagA+pttagB);
	    
	    double etaTA = jtphi_tA[k];
	    double etaTB = jtphi_tB[l];
	    double etatag = 0.5*(etaTA + etaTB);

	    double phiTA = jtphi_tA[k];
	    double phiTB = jtphi_tB[l];
	    double dphiTA = delta_phi(phiTA,phi);
	    double dphiTB = delta_phi(phiTB,phiB);
	    double dphitag = 0.5*(dphiTA + dphiTB);
	    
	    double alphaTA = (njt_tA>2 ? (*jtpt_tA)[2]/pttagA : 0);
	    double alphaTB = (njt_tB>2 ? (*jtpt_tB)[2]/pttagB : 0);
	    double alphatag = (alphaTA>0 && alphaTB>0 ?
			       0.5*(alphaTA+alphaTB) : max(alphaTA,alphaTB));
	    
	    istp = (dphitag>2.7 && alphatag<0.3 && fabs(etatag)<1.3 &&
		    pt/pttag>0.45 && ptB/pttag>0.45 &&
		    pttag>600.);
	  } // tag selection

	  // Look at good probe jets in barrel in good events
	  if (fabs(eta)<1.3 && fabs(etaB)<1.3 &&
	      jtidtight_tA[i]>=4 && jtidtight_tB[j]>=4 &&
	      flag_tA && flag_tB &&
	      pt > 0.5*ptB && ptB > 0.5*pt
	    //dr < 0.10) {
	    ) {

	    // Tag-and-probe method
	    if (istp) {
	      pta_tp->Fill(pttag, pt);
	      pa_tp->Fill(pttag, pt / pttag);
	      pb_tp->Fill(pttag, ptB / pttag);
	      pd_tp->Fill(pttag, 0.5*(ptB-pt) / pttag);
	      
	      h2a_tp->Fill(pttag, pt / pttag);
	      h2b_tp->Fill(pttag, ptB / pttag);
	      h2d_tp->Fill(pttag, 0.5*(ptB-pt) / pttag);
	    } //istp

	    // Direct matching method
	    pta_dm->Fill(pt, pt);
	    pa_dm->Fill(pt, ptB / pt);
	    ptb_dm->Fill(ptB, pt);
	    pb_dm->Fill(ptB, pt / ptB);
	    ptd_dm->Fill(ptave, pt);
	    pd_dm->Fill(ptave, 0.5*(ptB-pt) / ptave);
	    
	    h2a_dm->Fill(pt, ptB / pt);
	    h2b_dm->Fill(ptB, pt / ptB);
	    h2d_dm->Fill(ptave, 0.5*(ptB-pt) / ptave);
	  } // good barrel probe
	} // dr match
      } // for j
    } // for i
  } // for events
  cout << endl << "Found " << nmatch << " matching events"// << endl;
       << " of which " << nj << " had same number of jets" << endl;
    
  cout << "Output stored to " << f->GetName() << endl;
  f->Write();
  f->Close();

  t.Stop();
  cout << "Processing used " << t.CpuTime() << "s CPU time ("
       << t.CpuTime()/3600. << "h)" << endl;
  cout << "Processing used " << t.RealTime() << "s real time ("
       << t.RealTime()/3600. << "h)" << endl;
  cout << endl << endl;

}
