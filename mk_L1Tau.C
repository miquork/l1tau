#include "L1Tau.h"
//#include "drawL1Tau.C"
//#include "drawL1TauSelection.C"

R__LOAD_LIBRARY(L1Tau_C)
R__LOAD_LIBRARY(drawL1Tau_C)
R__LOAD_LIBRARY(drawL1TauSelection_C)


void mk_L1Tau() {

  gROOT->ProcessLine(".L L1Tau.C+g");
  gROOT->ProcessLine(".L drawL1Tau.C+g");

  TChain *c23c = new TChain("Events");
  //c23c->AddFile("../data/l1tau/Summer22Prompt23/Run2023B-PromptNanoAODv11p9_v1-v1_NANOAOD_taus.root");
  c23c->AddFile("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv11p9_v1-v1_NANOAOD_taus.root");
  c23c->AddFile("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v2-v2_NANOAOD_taus.root");
  c23c->AddFile("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v2-v4_NANOAOD_taus.root");
  c23c->AddFile("../data/l1tau/Summer22Prompt23/Run2023C-PromptNanoAODv12_v3-v1_NANOAOD_taus.root");
  L1Tau tauc(c23c,"23C");
  tauc.Loop();

  TChain *c23d = new TChain("Events");
  c23d->AddFile("../data/l1tau/Summer22Prompt23/Run2023D-PromptReco-v1_NANOAOD_taus.root");
  c23d->AddFile("../data/l1tau/Summer22Prompt23/Run2023D-PromptReco-v2_NANOAOD_taus.root");
  L1Tau taud(c23d,"23D");
  taud.Loop();


  // Draw efficiency parameterizations in 1D and 2D
  drawL1Tau("IsoTau34");
  drawL1Tau("IsoTau32");
  drawL1Tau("IsoTau26");

  drawL1Tau("Jet55");
  drawL1Tau("Jet90");

  // These cannot be reliably estimated from tag-and-probe sample
  // L1_ETMHF90 starts firing too often when L1Jet>90
  //drawL1Tau("Jet110");
  //drawL1Tau("Jet140");
  //drawL1Tau("Jet180");

  // Draw sample selection basics and Mjj prefiring estimates
  drawL1TauSelection("23C");
  drawL1TauSelection("23D");
}
