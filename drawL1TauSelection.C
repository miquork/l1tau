// Purpose: draw L1Tau event selection
// 1) hpttag, hptprobe for event counts
// 2) ptag,pprobe for phase space and selection efficiency
//    - pT>600 GeV
// 3) hptl1tau, hptl1jet so show safe phase space for efficiencies
//    - L1Tau<70 GeV, L1Jet<80(?); lines for MHT90, Tau120, Jet180 trigger
// 4) efficiencies vs Mjj
#include "TProfile.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2D.h"

#include "tdrstyle_mod22.C"

void drawL1TauSelection(string name = "23C") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f0 = new TFile("rootfiles/unprefireable.root","READ");
  assert(f0 && !f0->IsZombie());
  TH1D *hr = (TH1D*)f0->Get("hr"); assert(hr);
  TH1D *hr23d = (TH1D*)f0->Get("hr23d"); assert(hr23d);

  const char *cn = name.c_str();
  TFile *f = new TFile(Form("rootfiles/output-%s.root",cn),"READ");
  assert(f && !f->IsZombie());

  TH1D *hpttag0 = (TH1D*)f->Get("hpttag0"); assert(hpttag0);
  TH1D *hptprobe0 = (TH1D*)f->Get("hptprobe0"); assert(hptprobe0);
  TH1D *hpttag1 = (TH1D*)f->Get("hpttag1"); assert(hpttag1);
  TH1D *hptprobe1 = (TH1D*)f->Get("hptprobe1"); assert(hptprobe1);
  TH1D *hpttag2 = (TH1D*)f->Get("hpttag2"); assert(hpttag2);
  TH1D *hptprobe2 = (TH1D*)f->Get("hptprobe2"); assert(hptprobe2);

  TProfile *ptag = (TProfile*)f->Get("ptag"); assert(ptag);
  TProfile *pprobe = (TProfile*)f->Get("pprobe"); assert(pprobe);

  TH1D *hptl1jet = (TH1D*)f->Get("hptl1jet_pre"); assert(hptl1jet);
  TH1D *hptl1tau = (TH1D*)f->Get("hptl1tau_pre"); assert(hptl1tau);
  TH1D *hptl1met = (TH1D*)f->Get("hptl1met_pre"); assert(hptl1met);

  TH1D *hptl1jet_tag = (TH1D*)f->Get("hptl1jet_tag"); assert(hptl1jet_tag);
  TH1D *hptl1tau_tag = (TH1D*)f->Get("hptl1tau_tag"); assert(hptl1tau_tag);
  TH1D *hptl1met_tag = (TH1D*)f->Get("hptl1met_tag"); assert(hptl1met_tag);
  
  TH1D *hptl1jet_unp = (TH1D*)f->Get("hptl1jet_unp"); assert(hptl1jet_unp);
  TH1D *hptl1tau_unp = (TH1D*)f->Get("hptl1tau_unp"); assert(hptl1tau_unp);
  TH1D *hptl1met_unp = (TH1D*)f->Get("hptl1met_unp"); assert(hptl1met_unp);

  TH1D *hptl1jet_bx1 = (TH1D*)f->Get("hptl1jet_bx1"); assert(hptl1jet_bx1);
  TH1D *hptl1tau_bx1 = (TH1D*)f->Get("hptl1tau_bx1"); assert(hptl1tau_bx1);
  TH1D *hptl1met_bx1 = (TH1D*)f->Get("hptl1met_bx1"); assert(hptl1met_bx1);

  TH1D *hptl1metnotag = (TH1D*)f->Get("hptl1metnotag"); assert(hptl1metnotag);
  TH1D *hptl1metnotag_wtagjet = (TH1D*)f->Get("hptl1metnotag_wtagjet");
  assert(hptl1metnotag_wtagjet);
  TH1D *hptl1metnotag_notagjet = (TH1D*)f->Get("hptl1metnotag_notagjet");
  assert(hptl1metnotag_notagjet);
  //
  TH1D *hptl1met_wtagjet = (TH1D*)f->Get("hptl1met_wtagjet");
  assert(hptl1met_wtagjet);
  TH1D *hptl1met_notagjet = (TH1D*)f->Get("hptl1met_notagjet");
  assert(hptl1met_notagjet);
  
  TProfile *p34 = (TProfile*)f->Get("p1pre34"); assert(p34);
  TProfile *p32 = (TProfile*)f->Get("p1pre32"); assert(p32);
  TProfile *p32to34 = (TProfile*)f->Get("p1pre32to34"); assert(p32to34);
  TProfile *p26 = (TProfile*)f->Get("p1pre26"); assert(p26);
  TProfile *p55 = (TProfile*)f->Get("p1pre55"); assert(p55);
  TProfile *p90 = (TProfile*)f->Get("p1pre90"); assert(p90);
  TProfile *pjet = (TProfile*)f->Get("p1prejet"); assert(pjet);
  TProfile *ptau = (TProfile*)f->Get("p1pretau"); assert(ptau);

  TProfile *p34mjj = (TProfile*)f->Get("p34mjj"); assert(p34mjj);
  TProfile *p32mjj = (TProfile*)f->Get("p32mjj"); assert(p32mjj);
  TProfile *p32to34mjj = (TProfile*)f->Get("p32to34mjj"); assert(p32to34mjj);
  TProfile *p26mjj = (TProfile*)f->Get("p26mjj"); assert(p26mjj);
  TProfile *p90mjj = (TProfile*)f->Get("p90mjj"); assert(p90mjj);
  TProfile *ptotmjj = (TProfile*)f->Get("ptotmjj"); assert(ptotmjj);
  TProfile *ptotmjjbx1 = (TProfile*)f->Get("ptotmjjbx1"); assert(ptotmjjbx1);

  TH2D *h2metvstau = (TH2D*)f->Get("hptl1metvstau"); assert(h2metvstau);
  TH2D *h2metvsjet = (TH2D*)f->Get("hptl1metvsjet"); assert(h2metvsjet);
  TH2D *h2tauvsjet = (TH2D*)f->Get("hptl1tauvsjet"); assert(h2tauvsjet);
  
  TProfile *p90met = (TProfile*)f->Get("p1pre90met"); assert(p90met);
  TProfile *p90metno180jet = (TProfile*)f->Get("p1pre90metno180jet");
  assert(p90metno180jet);
  //TProfile *p75tau = (TProfile*)f->Get("p1pre75t"); assert(p75tau);
  //TProfile *p80tau = (TProfile*)f->Get("p1pre80t"); assert(p80tau);
  TProfile *p80tau = (TProfile*)f->Get("p1pre77p5t"); assert(p80tau);
  TProfile *p130jet = (TProfile*)f->Get("p1pre130"); assert(p130jet);

  TProfile *p180jet = (TProfile*)f->Get("p1pre180jet"); assert(p180jet);
  TProfile *p120tau = (TProfile*)f->Get("p1pre120tau"); assert(p120tau);
  TProfile *p120tau180jet = (TProfile*)f->Get("p1pre120and180");
  assert(p120tau180jet);
  
  curdir->cd();
  
  TH1D *h1 = tdrHist("h1","Events",0.5,5e6,"p_{T} (GeV)",100,7000);
  TH1D *h1d = tdrHist("h1d","Ratio",0,1,"p_{T} (GeV)",100,7000);
  lumi_136TeV = Form("Run20%s",cn);
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h1,h1d,8,11);

  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(hpttag0,"HIST",kNone,kBlue+1,kSolid,-1,kNone);
  tdrDraw(hptprobe0,"HIST",kNone,kRed+1,kSolid,-1,kNone);
  tdrDraw(hpttag2,"HIST",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hptprobe2,"HIST",kNone,kRed,kSolid,-1,kNone);

  TLegend *leg1 = tdrLeg(0.65,0.88-0.05*4,0.80,0.88);
  leg1->AddEntry(hpttag0," Tag (all)","L");
  leg1->AddEntry(hptprobe0," Probe (all)","L");
  leg1->AddEntry(hpttag2," Tag (pass)","L");
  leg1->AddEntry(hptprobe2," Probe (pass)","L");

  
  c1->cd(2);
  gPad->SetLogx();

  tdrDraw(ptag,"HIST",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(pprobe,"HIST",kNone,kRed,kSolid,-1,kNone);

  c1->SaveAs(Form("pdf/drawL1TauSelection_hpt_%s.pdf",cn));
  

  TH1D *h2 = tdrHist("h2","Events",0.5,5e6,"Probe L1 p_{T} in BX-1 (GeV)",
		     10,1100);
  TH1D *h2d = tdrHist("h2d","Ratio",0,2,"Probe L1 p_{T} in BX-1 (GeV)",10,1100);
  TCanvas *c2 = tdrDiCanvas("c2",h2,h2d,8,11);

  c2->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D *hptl1met_unp_w = (TH1D*)hptl1met_unp->Clone("hptl1met_unp_w");
  double ku = hptl1met->Integral()/hptl1met_unp->Integral();
  hptl1jet_unp->Scale(ku);
  hptl1tau_unp->Scale(ku);
  hptl1met_unp->Scale(ku);

  TH1D *hptl1met_bx1_w = (TH1D*)hptl1met_bx1->Clone("hptl1met_bx1_w");
  double kbx1 = hptl1met->Integral()/hptl1met_bx1->Integral();
  hptl1jet_bx1->Scale(kbx1);
  hptl1tau_bx1->Scale(kbx1);
  hptl1met_bx1->Scale(kbx1);
  
  TLine *l = new TLine();
  l->SetLineColor(kGray+2);
  l->SetLineStyle(kDashed);
  l->DrawLine(10,kbx1,1100,kbx1);
  l->SetLineStyle(kDotted);
  l->DrawLine(10,ku,1100,ku);

  /*
  // Try returning these once 1st BX unprefireabla available
  tdrDraw(hptl1jet_unp,"E",kNone,kBlue-9,kSolid,-1,kNone);
  tdrDraw(hptl1met_unp,"E",kNone,kGreen-9,kSolid,-1,kNone);
  tdrDraw(hptl1tau_unp,"E",kNone,kRed-9,kSolid,-1,kNone);
  */
  // Try returning these once 1st BX unprefireabla available
  tdrDraw(hptl1jet_bx1,"E",kNone,kBlue-9,kSolid,-1,kNone);
  //tdrDraw(hptl1met_bx1,"E",kNone,kGreen-9,kSolid,-1,kNone);
  tdrDraw(hptl1tau_bx1,"E",kNone,kRed-9,kSolid,-1,kNone);

  double kf = 27.;
  hptl1jet_tag->Scale(kf);
  hptl1met_tag->Scale(kf);
  hptl1tau_tag->Scale(kf);
  tdrDraw(hptl1jet_tag,"HIST",kNone,kBlue+2,kSolid,-1,kNone);
  //tdrDraw(hptl1met_tag,"HIST",kNone,kGreen+4,kSolid,-1,kNone);
  tdrDraw(hptl1tau_tag,"HIST",kNone,kRed+2,kSolid,-1,kNone);
  
  tdrDraw(hptl1jet,"HIST",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hptl1met,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hptl1tau,"HIST",kNone,kRed,kSolid,-1,kNone);

  TF1 *f1j = new TF1("f1j","[0]*pow(x,[1])",40,80);
  f1j->SetParameters(1e11,-5);
  hptl1jet->Fit(f1j,"RN");
  f1j->SetLineColor(kBlue);
  f1j->DrawClone("SAME");
  f1j->SetLineStyle(kDotted);
  f1j->SetRange(27.5,1100);
  f1j->Draw("SAME");

  TF1 *f1m = new TF1("f1t","[0]*pow(x,[1])",50.,90.);
  f1m->SetParameters(1e11,-5);
  hptl1met->Fit(f1m,"RN");
  f1m->SetLineColor(kGreen+2);
  f1m->DrawClone("SAME");
  f1m->SetLineStyle(kDotted);
  f1m->SetRange(25.,200);
  f1m->Draw("SAME");

  TF1 *f1t = new TF1("f1t","[0]*pow(x,[1])",27.5,75);
  f1t->SetParameters(1e11,-5);
  hptl1tau->Fit(f1t,"RN");
  f1t->SetLineColor(kRed);
  f1t->DrawClone("SAME");
  f1t->SetLineStyle(kDotted);
  f1t->SetRange(22.5,300);
  f1t->Draw("SAME");

  l->SetLineColor(kBlue);
  l->DrawLine(180,0.5,180,1e3);
  //l->DrawLine(90,0.5,90,2e4);
  l->SetLineColor(kRed);
  l->DrawLine(120,0.5,120,4e3);//1e2);
  //l->DrawLine(90,0.5,90,5e2);
  l->DrawLine(32,0.5,32,5e5);
  l->SetLineColor(kGreen+2);
  l->DrawLine(90,0.5,90,2e4);

  TLegend *leg2 = tdrLeg(0.38,0.88-0.05*4,0.53,0.88);
  leg2->SetHeader("Prefireable events, tag tau veto");
  leg2->AddEntry(hptl1jet," L1Jet","L");
  leg2->AddEntry(hptl1met," L1EtSum","L");
  leg2->AddEntry(hptl1tau," L1Tau","L");

  /*
  TLegend *leg2b = tdrLeg(0.58,0.73-0.05*4,0.73,0.73);
  leg2b->SetHeader(Form("Unprefireable #times %1.0f",k));
  leg2b->AddEntry(hptl1jet_unp," L1Jet","LE");
  leg2b->AddEntry(hptl1met_unp," L1EtSum","LE");
  leg2b->AddEntry(hptl1tau_unp," L1Tau","LE");
  */
  TLegend *leg2b = tdrLeg(0.58,0.53,0.73,0.53+0.05*3);
  leg2b->SetHeader(Form("BX1 #times %1.1f",kbx1));
  leg2b->AddEntry(hptl1jet_bx1," L1Jet","LE");
  //leg2b->AddEntry(hptl1met_bx1," L1EtSum","LE");
  leg2b->AddEntry(hptl1tau_bx1," L1Tau","LE");

  TLegend *leg2c = tdrLeg(0.68,0.38,0.83,0.38+0.05*3);
  leg2c->SetHeader(Form("Tag fail #times %1.1f",kf));
  leg2c->AddEntry(hptl1jet_tag," L1Jet","LE");
  leg2c->AddEntry(hptl1tau_tag," L1Tau","LE");
  
  gPad->RedrawAxis();
  
  c2->cd(2);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(10,1,1100,1);
  l->SetLineColor(kBlue);
  l->DrawLine(180,0.0,180,12.0);
  //l->DrawLine(90,0.5,90,2.0);
  l->SetLineColor(kRed);
  l->DrawLine(120,0.,120,2.0);
  //l->DrawLine(90,0.,90,0.5);
  l->DrawLine(32,0.0,32,2.0);
  l->SetLineColor(kGreen+2);
  l->DrawLine(90,0.,90,2.0);

  /*
  // Try returning these once 1st BX unprefireable available
  TH1D *hrjet_unp = (TH1D*)hptl1jet->Clone("hrjet_unp");
  hrjet_unp->Divide(f1j);
  //tdrDraw(hrjet_unp,"HIST",kNone,kBlue-9,kSolid,-1,kNone);

  TH1D *hrtau_unp = (TH1D*)hptl1tau_unp->Clone("hrtau_unp");
  hrtau_unp->Divide(f1t);
  //tdrDraw(hrtau_unp,"HIST",kNone,kRed-9,kSolid,-1,kNone);

  TH1D *hrmet_unp = (TH1D*)hptl1met_unp->Clone("hrmet_unp");
  hrmet_unp->Divide(f1m);
  //tdrDraw(hrmet_unp,"HIST",kNone,kGreen-9,kSolid,-1,kNone);
  */

  TH1D *hrjet_bx1 = (TH1D*)hptl1jet_bx1->Clone("hrjet_bx1");
  hrjet_bx1->Divide(f1j);
  tdrDraw(hrjet_bx1,"HIST",kNone,kBlue-9,kSolid,-1,kNone);

  TH1D *hrtau_bx1 = (TH1D*)hptl1tau_bx1->Clone("hrtau_bx1");
  hrtau_bx1->Divide(f1t);
  tdrDraw(hrtau_bx1,"HIST",kNone,kRed-9,kSolid,-1,kNone);

  TH1D *hrmet_bx1 = (TH1D*)hptl1met_bx1->Clone("hrmet_bx1");
  hrmet_bx1->Divide(f1m);
  //tdrDraw(hrmet_bx1,"HIST",kNone,kGreen-9,kSolid,-1,kNone);

  TH1D *hrjet_tag = (TH1D*)hptl1jet_tag->Clone("hrjet_tag");
  hrjet_tag->Divide(f1j);
  tdrDraw(hrjet_tag,"HIST",kNone,kBlue+2,kSolid,-1,kNone);

  TH1D *hrtau_tag = (TH1D*)hptl1tau_tag->Clone("hrtau_tag");
  hrtau_tag->Divide(f1t);
  tdrDraw(hrtau_tag,"HIST",kNone,kRed+2,kSolid,-1,kNone);

  TH1D *hrmet_tag = (TH1D*)hptl1met_tag->Clone("hrmet_tag");
  hrmet_tag->Divide(f1m);
  //tdrDraw(hrmet_tag,"HIST",kNone,kGreen+4,kSolid,-1,kNone);
  
  TH1D *hrjet = (TH1D*)hptl1jet->Clone("hrjet");
  hrjet->Divide(f1j);
  tdrDraw(hrjet,"HIST",kNone,kBlue,kSolid,-1,kNone);

  TH1D *hrtau = (TH1D*)hptl1tau->Clone("hrtau");
  hrtau->Divide(f1t);
  tdrDraw(hrtau,"HIST",kNone,kRed,kSolid,-1,kNone);

  TH1D *hrmet = (TH1D*)hptl1met->Clone("hrmet");
  hrmet->Divide(f1m);
  tdrDraw(hrmet,"HIST",kNone,kGreen+2,kSolid,-1,kNone);

  gPad->RedrawAxis();
  c2->SaveAs(Form("pdf/drawL1TauSelection_hptl1_%s.pdf",cn));


  // Look deeper into L1EtSum to patch prefiring
  TH1D *h6 = tdrHist("h6","Events",0.5,5e6,"L1EtSum p_{T} in BX-1 (GeV)",
		     10,1100);
  TH1D *h6d = tdrHist("h6d","Ratio",0,2,"L1EtSum p_{T} in BX-1 (GeV)",10,1100);
  TCanvas *c6 = tdrDiCanvas("c6",h6,h6d,8,11);

  c6->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  
  //double k = hptl1met->Integral()/hptl1met_unp->Integral();
  hptl1met_unp_w->Rebin(4);
  hptl1met_unp_w->Scale(ku/4.);

  //double k = hptl1met->Integral()/hptl1met_bx1->Integral();
  hptl1met_bx1_w->Rebin(4);
  hptl1met_bx1_w->Scale(kbx1/4.);

  l->SetLineColor(kGray+2);
  l->SetLineStyle(kDashed);
  l->DrawLine(10,kbx1,1100,kbx1);
  l->SetLineStyle(kDotted);
  l->DrawLine(10,ku,1100,ku);
  
  tdrDraw(hptl1met_unp_w,"E",kNone,kOrange+2,kSolid,-1,kNone);
  tdrDraw(hptl1met_bx1_w,"E",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hptl1met,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
  tdrDraw(hptl1met_wtagjet,"HIST",kNone,kBlue,kSolid,-1,kNone);

  //tdrDraw(hptl1metnotag,"HIST",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hptl1metnotag_wtagjet,"HIST",kNone,kMagenta+2,kSolid,-1,kNone);
  //tdrDraw(hptl1metnotag_notagjet,"HIST",kNone,kOrange+2,kSolid,-1,kNone);
  
  l->SetLineColor(kGreen+2);
  l->DrawLine(90,0.5,90,2e4);

  TLegend *leg6 = tdrLeg(0.38,0.88-0.05*3,0.53,0.88);
  leg6->AddEntry(hptl1met," L1EtSum (all events)","L");
  leg6->AddEntry(hptl1met_wtagjet," L1EtSum (tag has L1Jet)","L");
  //leg6->AddEntry(hptl1metnotag," L1EtSum+L1Jet (all)","L");
  leg6->AddEntry(hptl1metnotag_wtagjet," L1EtSum + tag L1Jet","L");
  //leg6->AddEntry(hptl1metnotag_notagjet," L1EtSum+L1Jet (no tag)","L");


  TLegend *leg6b = tdrLeg(0.58,0.53,0.73,0.53+0.05*2);
  leg6b->SetHeader(Form("Unprefireable #times %1.0f",ku));
  leg6b->AddEntry(hptl1met_unp_w," L1EtSum","LE");

  TLegend *leg6c = tdrLeg(0.58,0.43,0.73,0.43+0.05*2);
  leg6c->SetHeader(Form("BX1 #times %1.1f",kbx1));
  leg6c->AddEntry(hptl1met_bx1_w," L1EtSum","LE");

  gPad->RedrawAxis();
  
  c6->cd(2);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(10,1,1100,1);
  l->SetLineColor(kGreen+2);
  l->DrawLine(90,0.,90,2.0);

  gPad->RedrawAxis();
  c6->SaveAs(Form("pdf/drawL1TauSelection_hptl1met_%s.pdf",cn));
  

  TH1D *h3 = tdrHist("h3","Fraction in BX-1",0,1,"p_{T,probe} (GeV)",600,4000);
  TCanvas *c3 = tdrCanvas("c3",h3,8,0,kSquare);

  // Correct for missing some of L1_ETMHF90
  // Only get these when jet on tag side
  // But we throw away events with (tau26 or) jet55 on tag side
  // Latter more frequent and likely contains tau26 case
  // Update: consider this as additive rate over single jets and ditaus
  TH1D *h90c = p90->ProjectionX("h90c");
  for (int i = 1; i != h90c->GetNbinsX()+1; ++i) {
    double fjet = pjet->GetBinContent(i);// - p55->GetBinContent(i);
    if (h90c->GetBinContent(i)!=0 && fjet>0) {
      h90c->SetBinContent(i, h90c->GetBinContent(i)/fjet);
      h90c->SetBinError(i, h90c->GetBinError(i)/fjet);
    }
  }

  // Updated evaluations
  TProfile *p90combo = (TProfile*)p90met->Clone("p90combo");
  p90combo->Add(p130jet);
  //p90combo->Add(p75tau);
  p90combo->Add(p80tau);

  /*
  tdrDraw(pjet,"Pz",kFullSquare,kGray+1);
  tdrDraw(p90,"Pz",kFullSquare,kRed);
  tdrDraw(p55,"Pz",kOpenSquare,kGreen+1);
  tdrDraw(p26,"Pz",kOpenSquare,kGreen+3);
  tdrDraw(p32,"Pz",kFullCircle,kMagenta+1);
  tdrDraw(p32to34,"Pz",kFullDiamond,kMagenta+1);
  tdrDraw(p34,"Pz",kFullCircle,kOrange+1);

  tdrDraw(h90c,"Pz",kOpenSquare,kRed);
  tdrDraw(p90combo,"Pz",kFullDiamond,kBlack);

  TLegend *leg3 = tdrLeg(0.18,0.90-8*0.03,0.43,0.90);
  leg3->SetTextSize(0.030);
  leg3->AddEntry(pjet,"L1Jet","PLE");
  leg3->AddEntry(p26,"L1Tau26","PLE");
  leg3->AddEntry(p55,"L1Jet55","PLE");
  leg3->AddEntry(p32,"L1Tau32","PLE");
  leg3->AddEntry(p34,"L1Tau34","PLE");
  //leg3->AddEntry(h90c,"L1MHT90 (corr.)","PLE");
  //leg3->AddEntry(p90,"L1MHT90 (raw)","PLE");
  leg3->AddEntry(h90c,"L1Jet90/L1Jet","PLE");
  leg3->AddEntry(p90,"L1JEt90","PLE");
  leg3->AddEntry(p32to34,"L1Tau32-34","PLE");
  leg3->AddEntry(p90combo,"L1MET90 (est.)","PLE");
  */
  
  tdrDraw(pjet,"Pz",kFullSquare,kGray+1);
  tdrDraw(ptau,"Pz",kOpenSquare,kGray+2);
  tdrDraw(p55,"Pz",kOpenSquare,kGreen+1);
  tdrDraw(p26,"Pz",kOpenSquare,kGreen+3);
  tdrDraw(p32,"Pz",kFullCircle,kMagenta+1);
  tdrDraw(p34,"Pz",kOpenCircle,kOrange+1);
  //tdrDraw(p130jet,"Pz",kFullDiamond,kRed);
  tdrDraw(p80tau,"Pz",kFullDiamond,kRed);

  TLegend *leg3 = tdrLeg(0.18,0.90-8*0.03,0.43,0.90);
  leg3->SetTextSize(0.030);
  leg3->AddEntry(pjet,"L1Jet","PLE");
  leg3->AddEntry(ptau,"L1Tau","PLE");
  leg3->AddEntry(p26,"L1Tau26","PLE");
  leg3->AddEntry(p55,"L1Jet55","PLE");
  leg3->AddEntry(p32,"L1Tau32","PLE");
  leg3->AddEntry(p34,"L1Tau34","PLE");
  //leg3->AddEntry(p130jet,"L1Jet130 (MET90)","PLE");
  leg3->AddEntry(p80tau,"L1Tau77.5","PLE");
  leg3->AddEntry(p80tau,"(MET90)","");

  // Add fits
  double ptmin(600), ptmax(4000);//2941);
  TF1 *f2c = new TF1("f2c","[0]*0.5*(1+erf((0.5*x*[1]-32)/([2]*32)))+"
		     "[3]*0.5*(1+erf((0.5*x*[4]-32)/([5]*32)))",
		     ptmin,ptmax);
  f2c->SetParameters(0.1613,0.03737,0.2247, 0.1593,0.06843,0.4608); // IsoTau32
  f2c->SetParLimits(0,0,1);
  f2c->SetParLimits(3,0,1);
  //
  pjet->Fit(f2c,"QRN");
  f2c->SetLineColor(kGray+1);
  f2c->DrawClone("SAME");
  //
  ptau->Fit(f2c,"QRN");
  f2c->SetLineColor(kGray+2);
  f2c->DrawClone("SAME");
  //
  p26->Fit(f2c,"QRN");
  f2c->SetLineColor(kGreen+3);
  f2c->DrawClone("SAME");
  //
  p55->Fit(f2c,"QRN");
  f2c->SetLineColor(kGreen+1);
  f2c->DrawClone("SAME");
  //
  p32->Fit(f2c,"QRN");
  f2c->SetLineColor(kMagenta+1);
  f2c->DrawClone("SAME");
  //
  p34->Fit(f2c,"QRN");
  f2c->SetLineColor(kOrange+1);
  f2c->DrawClone("SAME");
  //
  //f2c->FixParameter(0,f2c->GetParameter(0));
  f2c->FixParameter(3,f2c->GetParameter(3));
  if (name=="23D") f2c->FixParameter(3,0);//f2c->GetParameter(3));
  p80tau->Fit(f2c,"QRN");
  f2c->SetLineColor(kRed);
  f2c->DrawClone("SAME");
  
  gPad->RedrawAxis();
  c3->SaveAs(Form("pdf/drawL1TauSelection_pre_%s.pdf",cn));

  h3->GetYaxis()->SetRangeUser(1e-3,2e1);
  gPad->SetLogy();
  c3->SaveAs(Form("pdf/drawL1TauSelection_pre_%s_log.pdf",cn));
  
  TH1D *h5 = tdrHist("h5","Prefiring",0,1,"M_{jj} (GeV)",1200,8000);
  TCanvas *c5 = tdrCanvas("c5",h5,8,0,kSquare);

  tdrDraw(ptotmjjbx1,"Pz",kNone,kGray+1);
  tdrDraw(ptotmjj,"Pz",kFullSquare,kBlack);
  tdrDraw(p90mjj,"Pz",kOpenSquare,kRed);
  //tdrDraw(p26mjj,"Pz",kOpenSquare,kGreen+2);
  tdrDraw(p32mjj,"Pz",kOpenCircle,kMagenta+1);
  //tdrDraw(p32to34mjj,"Pz",kOpenDiamond,kMagenta+1);
  //tdrDraw(p34mjj,"Pz",kOpenCircle,kOrange+1);

  // chi2 / NDF = 12.8 / 10
  TF1 *f1_22FG_23BCv3_23Cv4 = new TF1("f1_22FG_23BCv3_23Cv4","[p0]*0.5*(1+erf((0.5*x*[p1]-34)/([p2]*34)))",30,8000);
  double pars_f1_22FG_23BCv3_23Cv4[3] = {0.5, 0.012, 0.4588};
  f1_22FG_23BCv3_23Cv4->SetParameters(pars_f1_22FG_23BCv3_23Cv4);
  // par1 = 0.0120 +/- 0.0002

  f1_22FG_23BCv3_23Cv4->SetLineColor(kGray+2);//kBlack);
  //f1_22FG_23BCv3_23Cv4->SetLineWidth(2);
  if (name=="23C") f1_22FG_23BCv3_23Cv4->Draw("SAME");
  if (name=="23C") tdrDraw(hr,"Pz",kNone,kGray+2); //hr->SetLineWidth(1);

  // chi2 / NDF = 10.0 / 8
  TF1 *f1_23D = new TF1("f1_23D","[p0]*0.5*(1+erf((0.5*x*[p1]-34)/([p2]*34)))",30,8000);
  double pars_f1_23D[3] = {0.5, 0.006347, 0.4588};
  f1_23D->SetParameters(pars_f1_23D);
  f1_23D->SetLineColor(kGreen+2);//kBlack);
  f1_23D->SetLineWidth(2);
  if (name=="23D") f1_23D->Draw("SAME");
  if (name=="23D") tdrDraw(hr23d,"Pz",kNone,kGreen+2); hr23d->SetLineWidth(2);


  //TLegend *leg5 = tdrLeg(0.18,0.90-8*0.030,0.43,0.90);
  TLegend *leg5 = tdrLeg(0.18,0.90-5*0.030,0.43,0.90);
  leg5->SetTextSize(0.030);
  leg5->AddEntry(hr,"Unprefireable","PLE");
  leg5->AddEntry(f1_22FG_23BCv3_23Cv4,"Unpref. fit","L");
  leg5->AddEntry(ptotmjj,"L1_* (TnPSim)","PLE");
  leg5->AddEntry(p90mjj,"L1_ETMHF90 (TnPSim)","PLE");
  //leg5->AddEntry(p34mjj,"L1_DoubleIsoTau34er2p1 (TnPSim)","PLE");
  //leg5->AddEntry(p32mjj,"(L1_DoubleIsoTau32er2p1_Mass_Max80)","PLE");
  leg5->AddEntry(p32mjj,"L1_DoubleIsoTau32er2p1_Mass_Max80 (TnPSim)","PLE");
  //leg5->AddEntry(p26mjj,"(L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5)","PLE");
  //leg5->AddEntry(p32to34mjj,"32 vs 34 extra rate (TnPSim)","PLE");
  
  gPad->RedrawAxis();
  c5->SaveAs(Form("pdf/drawL1TauSelection_mjj_%s.pdf",cn));

  h5->GetYaxis()->SetRangeUser(1e-3,1e1);
  gPad->SetLogy();
  c5->SaveAs(Form("pdf/drawL1TauSelection_mjj_%s_log.pdf",cn));


  TH1D *h7 = tdrHist("h7","L1MET (GeV)",0,90,"L1Jet,L1Tau (GeV)",24,180);
  TCanvas *c7 = tdrCanvas("c7",h7,8,0,kSquare);

  //tdrDraw(h2metvstau,"COLZ",kNone,kRed);
  //tdrDraw(h2metvsjet,"SAME BOX",kNone,kBlue,kSolid,-1,kNone);
  //h2metvstau->Draw("SAME COLZ");
  h2metvsjet->Draw("SAME COLZ");
  h2metvstau->Draw("SAME BOX");
  h2metvstau->Scale(50.);
  h2metvstau->SetLineColor(kRed-9);
  gPad->SetLogz();

  TProfile *pmetvsjet = h2metvsjet->ProfileX("pmetvsjet");
  //tdrDraw(pmetvsjet,"Pz",kFullCircle,kBlue+1);
  int ix1 = h2metvsjet->GetXaxis()->FindBin(30.);
  int ix2 = h2metvsjet->GetXaxis()->GetNbins();
  TProfile *pjetvsmet = h2metvsjet->ProfileY("pjetvsmet",ix1,ix2);
  TGraphErrors *gjetvsmet = new TGraphErrors();
  for (int i = 1; i != pjetvsmet->GetNbinsX()+1; ++i) {
    if (pjetvsmet->GetBinContent(i)) {
      gjetvsmet->SetPoint(i-1, pjetvsmet->GetBinContent(i),
			  pjetvsmet->GetBinCenter(i));
    }
  }
  tdrDraw(gjetvsmet,"Pz",kFullCircle,kBlue);

  l->SetLineStyle(kSolid);
  l->SetLineColor(kBlack);//kBlue+1);
  l->DrawLine(90,0,90,90);
  l->SetLineColor(kBlue);
  l->DrawLine(130,0,130,90);

  TProfile *pmetvstau = h2metvstau->ProfileX("pmetvstau");
  //tdrDraw(pmetvstau,"Pz",kFullCircle,kRed+1);
  int jx1 = h2metvstau->GetXaxis()->FindBin(25.);
  int jx2 = h2metvstau->GetXaxis()->GetNbins();
  TProfile *ptauvsmet = h2metvstau->ProfileY("pjetvstau",jx1,jx2);
  TGraphErrors *gtauvsmet = new TGraphErrors();
  for (int i = 1; i != ptauvsmet->GetNbinsX()+1; ++i) {
    if (ptauvsmet->GetBinContent(i)) {
      gtauvsmet->SetPoint(i-1, ptauvsmet->GetBinContent(i),
			  ptauvsmet->GetBinCenter(i));
    }
  }
  tdrDraw(gtauvsmet,"Pz",kFullCircle,kRed);

  l->SetLineColor(kRed);
  //l->DrawLine(75,0,75,90);
  l->DrawLine(77.5,0,77.5,90);

  h2metvstau->SetFillColor(kNone);
  h2metvstau->SetFillStyle(kNone);
  h2metvsjet->SetFillColor(kBlue);
  //TLegend *leg7 = tdrLeg(0.55,0.15,0.75,0.15+0.05*4);
  TLegend *leg7 = tdrLeg(0.63,0.15,0.78,0.15+0.035*4);
  leg7->SetTextSize(0.035);
  leg7->AddEntry(h2metvstau,"L1MET vs L1Tau","F");
  leg7->AddEntry(h2metvsjet,"L1MET vs L1Jet","F");
  leg7->AddEntry(gtauvsmet,"#LTL1Tau#GT vs L1MET","P");
  leg7->AddEntry(gjetvsmet,"#LTL1Jet#GT vs L1MET","P");
  
  gPad->RedrawAxis();
  c7->SaveAs(Form("pdf/drawL1TauSelection_metvstau_%s.pdf",cn));

  TH1D *h8 = tdrHist("h8","Fraction in BX-1",0,1,"p_{T,probe}",600,4000);
  //TH1D *h8d = tdrHist("h8d","Ratio",0,2,"p_{T,probe}",600,4000);
  //TCanvas *c8 = tdrDiCanvas("c8",h8,h8d,8,0);
  TCanvas *c8 = tdrCanvas("c8",h8,8,0,kSquare);

  c8->cd(1);

  tdrDraw(p90met,"Pz",kFullSquare,kGreen+2); //p90met->SetMarkerSize(0.8);
  //tdrDraw(p90metno180jet,"Pz",kOpenSquare,kGreen+2); // small change
  tdrDraw(p130jet,"Pz",kFullCircle,kBlue);
  tdrDraw(p80tau,"Pz",kFullDiamond,kRed);
  //tdrDraw(p75tau,"Pz",kFullDiamond,kRed);
  //tdrDraw(p90combo,"Pz",kFullDiamond,kBlack);

  TF1 *f90 = new TF1("f90","[p0]*0.5*(1+erf((0.5*x*[p1]-[p3])/([p2]*[p3])))+[p4]*x",600,8000);
  //double pars_90[5] = {0.5, 0.012, 0.4588, 34.,0};
  //double pars_90[5] = {0.5, 0.012, 0.4588, 130.,0};
  //double pars_90[5] = {0.5, 0.012, 0.4588, 90.,0};
  double pars_90[5] = {0.2, 0.012, 0.4588, 90.,0};
  f90->SetParameters(pars_90);
  //f90->FixParameter(4,0);
  //p130jet->Fit(f90,"QRN");
  //p90combo->Fit(f90,"QRN");
  p90met->Fit(f90,"QRN");
  f90->SetLineColor(kGreen+2);
  f90->DrawClone("SAME");

  //f90->SetParameters(pars_90);
  //p75tau->Fit(f90,"QRN");
  p80tau->Fit(f90,"QRN");
  f90->SetLineColor(kRed);
  f90->DrawClone("SAME");
  
  //f90->SetParameters(pars_90);
  p130jet->Fit(f90,"QRN");
  f90->SetLineColor(kBlue);
  f90->SetLineWidth(2); // favorite model for 23C
  f90->DrawClone("SAME");
  f90->SetLineWidth(1);
  
  p90met->Fit(f90,"QRN");
  f90->SetLineColor(kGreen+2);
  f90->DrawClone("SAME");

  //TLegend *leg8 = tdrLeg(0.20,0.88-0.05*4,0.35,0.88);
  TLegend *leg8 = tdrLeg(0.20,0.88-0.05*3,0.35,0.88);
  leg8->AddEntry(p90met,"MET90 (BX1 unprefireable)","PLE");
  leg8->AddEntry(p130jet,"Jet130 (TnP)","PLE");
  //leg8->AddEntry(p80tau,"Tau80 (TnP)","PLE");
  leg8->AddEntry(p80tau,"Tau77.5 (TnP)","PLE");
  //leg8->AddEntry(p75tau,"Tau75 (TnP)","PLE");
  //leg8->AddEntry(p90combo,"Jet130+Tau75 (TnP)","PLE");

  // Add upper limit on p90met
  TProfile *p90metplus1 = (TProfile*)p90met->Clone("p90metplus1");
  for (int i = 1; i != p90met->GetNbinsX()+1; ++i) {
    if (p90met->GetBinContent(i)==0) {
      p90metplus1->Fill(p90met->GetBinCenter(i), 1);
    }
  }
  tdrDraw(p90metplus1,"Pz",kOpenSquare,kGreen+2);
  p90metplus1->Fit(f90,"QRN");
  f90->SetLineColor(kGreen+2-9);
  f90->DrawClone("SAME");
  
  // Add upper limit on p130jet
  TProfile *p130jetplus1 = (TProfile*)p130jet->Clone("p130jetplus1");
  for (int i = 1; i != p130jet->GetNbinsX()+1; ++i) {
    if (p130jet->GetBinContent(i)==0) {
      p130jetplus1->Fill(p130jet->GetBinCenter(i), 1);
    }
  }
  tdrDraw(p130jetplus1,"Pz",kOpenCircle,kBlue);
  p130jetplus1->Fit(f90,"QRN");
  f90->SetLineColor(kBlue-9);
  f90->DrawClone("SAME");
  
  // Add upper limit on p80tau
  TProfile *p80tauplus1 = (TProfile*)p80tau->Clone("p80tauplus1");
  for (int i = 1; i != p80tau->GetNbinsX()+1; ++i) {
    if (p80tau->GetBinContent(i)==0) {
      p80tauplus1->Fill(p80tau->GetBinCenter(i), 1);
    }
  }
  tdrDraw(p80tauplus1,"Pz",kOpenDiamond,kRed);
  p80tauplus1->Fit(f90,"QRN");
  f90->SetLineColor(kRed-9);
  f90->DrawClone("SAME");
  
  c8->SaveAs(Form("pdf/drawL1TauSelection_met90_%s.pdf",cn));

  c8->cd(1);
  h8->GetYaxis()->SetRangeUser(1e-5,1e1);
  gPad->SetLogy();
  c8->SaveAs(Form("pdf/drawL1TauSelection_met90_%s_log.pdf",cn));


  TH1D *h9 = tdrHist("h9","L1Tau (GeV)",0,120,"L1Jet (GeV)",0,180);
  TCanvas *c9 = tdrCanvas("c9",h9,8,11,kSquare);

  h2tauvsjet->Draw("SAME COLZ");
  gPad->SetLogz();

  int ky1 = h2tauvsjet->GetYaxis()->FindBin(20.);
  int ky2 = h2tauvsjet->GetXaxis()->GetNbins();
  TProfile *ptauvsjet = h2tauvsjet->ProfileX("ptauvsjet",ky1,ky2);
  tdrDraw(ptauvsjet,"Pz",kFullCircle,kRed);
  
  int kx1 = h2tauvsjet->GetXaxis()->FindBin(30.);
  int kx2 = h2tauvsjet->GetXaxis()->GetNbins();
  TProfile *pjetvstau = h2tauvsjet->ProfileY("pjetvstau",kx1,kx2);
  TGraphErrors *gjetvstau = new TGraphErrors();
  for (int i = 1; i != pjetvstau->GetNbinsX()+1; ++i) {
    if (pjetvstau->GetBinContent(i)) {
      gjetvstau->SetPoint(i-1, pjetvstau->GetBinContent(i),
			  pjetvstau->GetBinCenter(i));
    }
  }
  tdrDraw(gjetvstau,"Pz",kFullCircle,kBlue);

  l->SetLineStyle(kSolid);
  //l->SetLineColor(kBlack);
  //l->DrawLine(0,90,180,90);
  //l->DrawLine(90,0,90,180);
  l->SetLineColor(kRed);
  //l->DrawLine(0,75,180,75);
  l->DrawLine(0,77.5,180,77.5);
  l->SetLineColor(kBlue);
  l->DrawLine(130,0,130,120);

  //TLegend *leg9 = tdrLeg(0.63,0.15,0.78,0.15+0.035*4);
  TLegend *leg9 = tdrLeg(0.20,0.78-0.035*4,0.35,0.78);
  leg9->SetTextSize(0.035);
  leg9->AddEntry(h2tauvsjet,"L1Tau vs L1Jet","F");
  leg9->AddEntry(ptauvsjet,"#LTL1Tau#GT vs L1Jet","P");
  leg9->AddEntry(gjetvstau,"#LTL1Jet#GT vs L1Tau","P");
  
  gPad->RedrawAxis();
  c9->SaveAs(Form("pdf/drawL1TauSelection_tauvsjet_%s.pdf",cn));

  
  TH1D *h10 = tdrHist("h10","Fraction in BX-1",0,1,"p_{T,probe}",600,4000);
  //TH1D *h10d = tdrHist("h10d","Ratio",0,2,"p_{T,probe}",600,4000);
  //TCanvas *c10 = tdrDiCanvas("c10",h10,h10d,8,0);
  TCanvas *c10 = tdrCanvas("c10",h10,8,0,kSquare);

  c10->cd(1);

  TH1D *h30jets = pjet->ProjectionX("h30jets");
  h30jets->Scale(pow(180./30.,-5+2));
  TH1D *h55jets = p55->ProjectionX("h55jets");
  h55jets->Scale(pow(180./55.,-5+2));
  TH1D *h130jets = p130jet->ProjectionX("h130jets");
  h130jets->Scale(pow(180./130.,-5+2));
  TH1D *h80taus = p80tau->ProjectionX("h80taus");
  h80taus->Scale(pow(120./77.5,-5+3));//2));
  //h80taus->Scale(pow(120./80.,-5+3));//2));
  //TH1D *h75taus = p75tau->ProjectionX("h75taus");
  //h75taus->Scale(pow(120./75.,-5+2));

  // Add upper limit on p90met
  TProfile *p180jetplus1 = (TProfile*)p180jet->Clone("p180jetplus1");
  for (int i = 1; i != p180jet->GetNbinsX()+1; ++i) {
    if (p180jet->GetBinContent(i)==0) {
      p180jetplus1->Fill(p180jet->GetBinCenter(i), 1);
    }
  }
  
  TH1D *h130jetplus1s = p130jetplus1->ProjectionX("h130jetplus1s");
  h130jetplus1s->Scale(pow(180./130.,-5+2));
  TH1D *h80tauplus1s = p80tauplus1->ProjectionX("h80tauplus1s");
  h80tauplus1s->Scale(pow(120./77.5,-5+3));//2));
  
  tdrDraw(p180jet,"Pz",kFullSquare,kGreen+2);
  tdrDraw(h130jets,"Pz",kFullCircle,kBlue);
  //tdrDraw(h55jets,"Pz",kOpenCircle,kGreen+2);
  //tdrDraw(h30jets,"Pz",kOpenDiamond,kGray+2);
  tdrDraw(h80taus,"Pz",kFullDiamond,kRed);
  //tdrDraw(h75taus,"Pz",kOpenDiamond,kRed);

  tdrDraw(p180jetplus1,"Pz",kOpenSquare,kGreen+2);
  tdrDraw(h130jetplus1s,"Pz",kOpenCircle,kBlue);
  tdrDraw(h80tauplus1s,"Pz",kOpenDiamond,kRed);
  
  TF1 *f180 = new TF1("f180","[p0]*0.5*(1+erf((0.5*x*[p1]-[p3])/([p2]*[p3])))+[p4]*x",600,8000);
  //double pars_180[5] = {0.5, 0.012, 0.4588, 34.,0};
  //double pars_180[5] = {0.5, 0.012, 0.4588, 180.,0};
  double pars_180[5] = {0.1, 0.012, 0.4588, 180.,0};
  f180->SetParameters(pars_180);
  //f180->FixParameter(4,0);
  p180jet->Fit(f180,"QRN");
  f180->SetLineColor(kGreen+2);
  //f180->DrawClone("SAME");

  h130jets->Fit(f180,"QRN");
  f180->SetLineColor(kBlue);
  f180->SetLineWidth(2); // favorite for 23C
  f180->DrawClone("SAME");
  f180->SetLineWidth(1);

  //h75taus->Fit(f180,"QRN");
  h80taus->Fit(f180,"QRN");
  f180->SetLineColor(kRed);
  f180->DrawClone("SAME");

  p180jet->Fit(f180,"QRN");
  f180->SetLineColor(kGreen+2);
  f180->DrawClone("SAME");
  
  TLegend *leg10 = tdrLeg(0.20,0.88-0.05*3,0.35,0.88);
  leg10->AddEntry(p180jet,"Jet180 (BX1 unprefireable)","PLE");
  leg10->AddEntry(p130jet,"Jet130 (TnP scaled)","PLE");
  leg10->AddEntry(p80tau,"Tau77.5 (TnP scaled)","PLE");
  
  c10->SaveAs(Form("pdf/drawL1TauSelection_jet180_%s.pdf",cn));

  c10->cd(1);
  h10->GetYaxis()->SetRangeUser(1e-5,1e1);
  gPad->SetLogy();
  c10->SaveAs(Form("pdf/drawL1TauSelection_jet180_%s_log.pdf",cn));
}
