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

#include "tdrstyle_mod22.C"

void drawL1TauSelection(string name = "23D") {

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

  TProfile *ptag = (TProfile*)f->Get("ptag"); assert(ptag);
  TProfile *pprobe = (TProfile*)f->Get("pprobe"); assert(pprobe);

  TH1D *hptl1jet = (TH1D*)f->Get("hptl1jet"); assert(hptl1jet);
  TH1D *hptl1tau = (TH1D*)f->Get("hptl1tau"); assert(hptl1tau);

  TProfile *p34 = (TProfile*)f->Get("p1pre34"); assert(p34);
  TProfile *p32 = (TProfile*)f->Get("p1pre32"); assert(p32);
  TProfile *p32to34 = (TProfile*)f->Get("p1pre32to34"); assert(p32to34);
  TProfile *p26 = (TProfile*)f->Get("p1pre26"); assert(p26);
  TProfile *p55 = (TProfile*)f->Get("p1pre55"); assert(p55);
  TProfile *p90 = (TProfile*)f->Get("p1pre90"); assert(p90);
  TProfile *pjet = (TProfile*)f->Get("p1prejet"); assert(pjet);

  TProfile *p34mjj = (TProfile*)f->Get("p34mjj"); assert(p34mjj);
  TProfile *p32mjj = (TProfile*)f->Get("p32mjj"); assert(p32mjj);
  TProfile *p32to34mjj = (TProfile*)f->Get("p32to34mjj"); assert(p32to34mjj);
  TProfile *p26mjj = (TProfile*)f->Get("p26mjj"); assert(p26mjj);
  TProfile *p90mjj = (TProfile*)f->Get("p90mjj"); assert(p90mjj);
  TProfile *ptotmjj = (TProfile*)f->Get("ptotmjj"); assert(ptotmjj);
  
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
  tdrDraw(hpttag1,"HIST",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hptprobe1,"HIST",kNone,kRed,kSolid,-1,kNone);

  TLegend *leg1 = tdrLeg(0.65,0.88-0.05*4,0.80,0.88);
  leg1->AddEntry(hpttag0," Tag (all)","L");
  leg1->AddEntry(hptprobe0," Probe (all)","L");
  leg1->AddEntry(hpttag1," Tag (pass)","L");
  leg1->AddEntry(hptprobe1," Probe (pass)","L");

  
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
  tdrDraw(hptl1jet,"HIST",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hptl1tau,"HIST",kNone,kRed,kSolid,-1,kNone);

  TF1 *f1j = new TF1("f1j","[0]*pow(x,[1])",40,80);
  f1j->SetParameters(1e11,-5);
  hptl1jet->Fit(f1j,"RN");
  f1j->SetLineColor(kBlue);
  f1j->DrawClone("SAME");
  f1j->SetLineStyle(kDotted);
  f1j->SetRange(27.5,1100);
  f1j->Draw("SAME");

  TF1 *f1t = new TF1("f1t","[0]*pow(x,[1])",27.5,75);
  f1t->SetParameters(1e11,-5);
  hptl1tau->Fit(f1t,"RN");
  f1t->SetLineColor(kRed);
  f1t->DrawClone("SAME");
  f1t->SetLineStyle(kDotted);
  f1t->SetRange(22.5,300);
  f1t->Draw("SAME");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kBlue);
  l->DrawLine(180,0.5,180,1e3);
  l->DrawLine(90,0.5,90,2e4);
  l->SetLineColor(kRed);
  l->DrawLine(120,0.5,120,1e2);
  l->DrawLine(90,0.5,90,5e2);
  l->DrawLine(32,0.5,32,5e5);

  TLegend *leg2 = tdrLeg(0.753,0.88-0.05*2,0.88,0.88);
  leg2->AddEntry(hptl1jet," L1Jet","L");
  leg2->AddEntry(hptl1tau," L1Tau","L");
  
  c2->cd(2);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(10,1,1100,1);
  l->SetLineColor(kBlue);
  l->DrawLine(180,0.0,180,12.0);
  l->DrawLine(90,0.5,90,2.0);
  l->SetLineColor(kRed);
  l->DrawLine(120,0.,120,2.0);
  l->DrawLine(90,0.,90,0.5);
  l->DrawLine(32,0.0,32,2.0);

  
  TH1D *hrjet = (TH1D*)hptl1jet->Clone("hrjet");
  hrjet->Divide(f1j);
  tdrDraw(hrjet,"HIST",kNone,kBlue,kSolid,-1,kNone);

  TH1D *hrtau = (TH1D*)hptl1tau->Clone("hrtau");
  hrtau->Divide(f1t);
  tdrDraw(hrtau,"HIST",kNone,kRed,kSolid,-1,kNone);

  gPad->RedrawAxis();
  c2->SaveAs(Form("pdf/drawL1TauSelection_hptl1_%s.pdf",cn));


  TH1D *h3 = tdrHist("h3","Fraction in BX-1",0,1,"p_{T,probe} (GeV)",600,4000);
  TCanvas *c3 = tdrCanvas("c3",h3,8,0,kSquare);

  // Correct for missing some of L1_ETMHF90
  // Only get these when jet on tag side
  // But we throw away events with (tau26 or) jet55 on tag side
  // Latter more frequent and likely contains tau26 case
  TH1D *h90c = p90->ProjectionX("h90c");
  for (int i = 1; i != h90c->GetNbinsX()+1; ++i) {
    double fjet = pjet->GetBinContent(i);// - p55->GetBinContent(i);
    if (h90c->GetBinContent(i)!=0 && fjet>0) {
      h90c->SetBinContent(i, h90c->GetBinContent(i)/fjet);
      h90c->SetBinError(i, h90c->GetBinError(i)/fjet);
    }
  }
  
  tdrDraw(pjet,"Pz",kFullSquare,kGray+1);
  tdrDraw(p90,"Pz",kFullSquare,kRed);
  tdrDraw(p55,"Pz",kOpenSquare,kGreen+1);
  tdrDraw(p26,"Pz",kOpenSquare,kGreen+3);
  tdrDraw(p32,"Pz",kFullCircle,kMagenta+1);
  tdrDraw(p32to34,"Pz",kFullDiamond,kMagenta+1);
  tdrDraw(p34,"Pz",kFullCircle,kOrange+1);

  tdrDraw(h90c,"Pz",kOpenSquare,kRed);

  TLegend *leg3 = tdrLeg(0.18,0.90-8*0.03,0.43,0.90);
  leg3->SetTextSize(0.030);
  leg3->AddEntry(pjet,"L1Jet","PLE");
  leg3->AddEntry(p26,"L1Tau26","PLE");
  leg3->AddEntry(p55,"L1Jet55","PLE");
  leg3->AddEntry(p32,"L1Tau32","PLE");
  leg3->AddEntry(p34,"L1Tau34","PLE");
  leg3->AddEntry(h90c,"L1MHT90 (corr.)","PLE");
  leg3->AddEntry(p90,"L1MHT90 (raw)","PLE");
  leg3->AddEntry(p32to34,"L1Tau32-34","PLE");
  
  gPad->RedrawAxis();
  c3->SaveAs(Form("pdf/drawL1TauSelection_pre_%s.pdf",cn));

  h3->GetYaxis()->SetRangeUser(1e-3,2e1);
  gPad->SetLogy();
  c3->SaveAs(Form("pdf/drawL1TauSelection_pre_%s_log.pdf",cn));
  
  TH1D *h5 = tdrHist("h5","Prefiring",0,1,"M_{jj} (GeV)",1200,8000);
  TCanvas *c5 = tdrCanvas("c5",h5,8,0,kSquare);

  tdrDraw(ptotmjj,"Pz",kFullSquare,kBlack);
  tdrDraw(p90mjj,"Pz",kOpenSquare,kRed);
  tdrDraw(p26mjj,"Pz",kOpenSquare,kGreen+2);
  tdrDraw(p32mjj,"Pz",kOpenCircle,kMagenta+1);
  tdrDraw(p32to34mjj,"Pz",kOpenDiamond,kMagenta+1);
  tdrDraw(p34mjj,"Pz",kOpenCircle,kOrange+1);

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


  TLegend *leg5 = tdrLeg(0.18,0.90-8*0.030,0.43,0.90);
  leg5->SetTextSize(0.030);
  leg5->AddEntry(hr,"Unprefireable","PLE");
  leg5->AddEntry(f1_22FG_23BCv3_23Cv4,"Unpref. fit","L");
  leg5->AddEntry(ptotmjj,"L1_* (TnPSim)","PLE");
  leg5->AddEntry(p90mjj,"L1_ETMHF90 (TnPSim)","PLE");
  leg5->AddEntry(p34mjj,"L1_DoubleIsoTau34er2p1 (TnPSim)","PLE");
  leg5->AddEntry(p32mjj,"(L1_DoubleIsoTau32er2p1_Mass_Max80)","PLE");
  leg5->AddEntry(p26mjj,"(L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5)","PLE");
  leg5->AddEntry(p32to34mjj,"32 vs 34 extra rate (TnPSim)","PLE");
  
  gPad->RedrawAxis();
  c5->SaveAs(Form("pdf/drawL1TauSelection_mjj_%s.pdf",cn));

  h5->GetYaxis()->SetRangeUser(1e-3,1e1);
  gPad->SetLogy();
  c5->SaveAs(Form("pdf/drawL1TauSelection_mjj_%s_log.pdf",cn));
}
