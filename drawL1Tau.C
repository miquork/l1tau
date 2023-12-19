// Purpose: draw L1Tau trigger inefficiency
#include "TFile.h"
#include "TProfile.h"
#include "TLine.h"
#include "TF1.h"

#include "tdrstyle_mod22.C"

void printPars(TF1 *f1, string name) {
  cout << Form("// chi2 / NDF = %1.1f / %d\n",
	       f1->GetChisquare(),f1->GetNDF());
  cout << "TF1 *"<<name<<" = new TF1(\""<<name<<"\",\""
       <<f1->GetExpFormula() << "\",30,8000);" << endl;
  cout << "double pars_"<<name<<"[" <<f1->GetNpar() << "] = {";
  for (int i = 0; i != f1->GetNpar(); ++i) {
    cout << Form("%s%1.4g", (i==0 ? "" : ", "), f1->GetParameter(i));
  }
  cout << "};" << endl;
  cout << name << "->SetParameters(pars_"<<name<<");" << endl;
} // printPars

void drawL1Tau(string path = "IsoTau34") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *fc = new TFile("rootfiles/output-23C.root","READ");
  assert(fc && !fc->IsZombie());
  
  TFile *fd = new TFile("rootfiles/output-23D.root","READ");
  assert(fd && !fd->IsZombie());

  curdir->cd();

  TProfile *p1c(0), *p2c(0), *p1d(0), *p2d(0);
  const char *cp = path.c_str();
  if (path=="IsoTau34") {
    p1c = (TProfile*)fc->Get("p1prew"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2prew2"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1prew"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2prew2"); assert(p2d);
  
  }
  else if (path=="IsoTau32") {
    p1c = (TProfile*)fc->Get("p1pre32"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre32"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre32"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre32"); assert(p2d);
  }
  else if (path=="IsoTau26") {
    p1c = (TProfile*)fc->Get("p1pre26"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre26"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre26"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre26"); assert(p2d);
  }
  else if (path=="Jet55") {
    p1c = (TProfile*)fc->Get("p1pre55"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre55"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre55"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre55"); assert(p2d);
  }
  else if (path=="Jet110") {
    p1c = (TProfile*)fc->Get("p1pre110"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre110"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre110"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre110"); assert(p2d);
  }
  else if (path=="Jet90") {
    p1c = (TProfile*)fc->Get("p1pre90"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre90"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre90"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre90"); assert(p2d);
  }
  else if (path=="Jet140") {
    p1c = (TProfile*)fc->Get("p1pre140"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre140"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre140"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre140"); assert(p2d);
  }
  else if (path=="Jet180") {
    p1c = (TProfile*)fc->Get("p1pre180"); assert(p1c);
    p2c = (TProfile*)fc->Get("p2pre180"); assert(p2c);
    p1d = (TProfile*)fd->Get("p1pre180"); assert(p1d);
    p2d = (TProfile*)fd->Get("p2pre180"); assert(p2d);
  }
  else assert(false);

  TH1D *h1c = p1c->ProjectionX("h1c");
  TH1D *h1d = p1d->ProjectionX("h1d");

  double ptmin = 600;
  double ptmax = 2941;
  TF1 *f1c = new TF1("f1c","[0]*0.5*(1+erf((0.5*x*[1]-34)/([2]*34)))",
		     ptmin,ptmax);
  f1c->SetParameters(0.209, 0.06136, 0.4795); // IsoTau32
  f1c->SetParLimits(0,0,1);
  h1c->Fit(f1c,"QRN");

  TF1 *f2c = new TF1("f2c","[0]*0.5*(1+erf((0.5*x*[1]-34)/([2]*34)))+"
		     "[3]*0.5*(1+erf((0.5*x*[4]-34)/([5]*34)))",
		     ptmin,ptmax);
  f2c->SetParameters(0.1613,0.03737,0.2247, 0.1593,0.06843,0.4608); // IsoTau32
  f2c->SetParLimits(0,0,1);
  f2c->SetParLimits(3,0,1);
  h1c->Fit(f2c,"QRN");

  TF1 *f1d = new TF1("f1d","[0]*0.5*(1+erf((0.5*x*[1]-34)/([2]*34)))",
		     ptmin,ptmax);
  f1d->SetParameters(0.1202, 0.04038, 0.3529); // IsoTau32
  f1d->SetParLimits(0,0,1);
  h1d->Fit(f1d,"QRN");

  TF1 *f2d = new TF1("f2d","[0]*0.5*(1+erf((0.5*x*[1]-34)/([2]*34)))+"
		     "[3]*0.5*(1+erf((0.5*x*[4]-34)/([5]*34)))",
		     ptmin,ptmax);
  f2d->SetParameters(0.1946,0.03221,0.355, 0.01585,0.05468,0.3112); // IsoTau32
  f2d->SetParLimits(0,0,1);
  f2d->SetParLimits(3,0,1);
  h1d->Fit(f2d,"QRN");
  
  lumi_136TeV = "Run3, 2023Cv123 + 2023D";
  extraText = "Private";
  TH1D *h = tdrHist("h",Form("%s BX-1 prefiring probability",cp),1e-3,1.2,
		    "p_{T,probe} (GeV)",300.,6076);
  h->GetXaxis()->SetRangeUser(300,4000.);
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(300,1,4000,1);
  l->SetLineColor(kRed+2);
  l->SetLineStyle(kDotted);
  l->DrawLine(600,1e-3,600,0.08);
  l->DrawLine(f1c->GetXmax(),1e-3,f1c->GetXmax(),0.8);
  
  tdrDraw(h1d,"Pz",kFullCircle,kBlack); h1d->SetLineWidth(3);
  tdrDraw(h1c,"Pz",kOpenCircle,kGray+2); h1c->SetLineWidth(2);


  f2c->SetLineColor(kRed);
  f2c->SetRange(300,4000.);
  f2c->Draw("SAME");

  f1c->SetLineColor(kBlue);
  f1c->SetRange(300,4000.);
  f1c->Draw("SAME");

  f2d->SetLineWidth(2);
  f2d->SetLineColor(kRed);
  f2d->SetRange(300,4000.);
  f2d->Draw("SAME");

  f1d->SetLineColor(kBlue);
  f1d->SetRange(300,4000.);
  f1d->Draw("SAME");

  
  TLegend *leg = tdrLeg(0.36,0.89-0.045*4,0.61,0.89);
  leg->AddEntry(h1d,"2023D","PLE");
  leg->AddEntry(h1c,"2023Cv123","PLE");
  leg->AddEntry(f2c,"2-comp. erf","L");
  leg->AddEntry(f1c,"1-comp. erf","L");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  tex->DrawLatex(0.18,0.65,"|#eta_{probe}| < 1.3");
  tex->DrawLatex(0.18,0.58,"L2L3Res_V3");
  
  gPad->RedrawAxis();
  
  
  TH1D *h2c = tdrHist("h2c","p_{T,probe} (GeV)",15,4500,
		      "|#eta_{probe}|",0,5.2);
  TCanvas *c2c = new TCanvas("c2c","c2c",600,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);

  h2c->Draw();
  p2c->Draw("SAME COLZ");
  p2c->GetZaxis()->SetRangeUser(0,1);
  p2c->GetZaxis()->SetTitleOffset(1.1);
  p2c->GetZaxis()->SetTitleFont(42);
  p2c->GetZaxis()->SetTitleSize(0.045);
  p2c->GetZaxis()->SetRangeUser(0.,1);
  p2c->SetZTitle(Form("%s BX-1 prefiring probability",cp));

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0,600.,5.2,600.);
  l->SetLineStyle(kDotted);
  l->DrawLine(-1.305,15,-1.305,4500.);
  l->DrawLine(+1.305,15,+1.305,4500.);

  tex->DrawLatex(0.50,0.89,"2023Cv123");
  tex->DrawLatex(0.50,0.84,cp);
  
  gPad->RedrawAxis();

  TH1D *h2d = tdrHist("h2d","p_{T,probe} (GeV)",15,4500,
		      "|#eta_{probe}|",0.,5.2);
  TCanvas *c2d = new TCanvas("c2d","c2d",600,600);
  gPad->SetLogy();
  gPad->SetRightMargin(0.15);

  h2d->Draw();
  p2d->Draw("SAME COLZ");
  p2d->GetZaxis()->SetRangeUser(0,1);
  p2d->GetZaxis()->SetTitleOffset(1.1);
  p2d->GetZaxis()->SetTitleFont(42);
  p2d->GetZaxis()->SetTitleSize(0.045);
  p2d->GetZaxis()->SetRangeUser(0.,1);
  p2d->SetZTitle(Form("%s BX-1 prefiring probability",cp));
  
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0,600.,5.2,600.);
  l->SetLineStyle(kDotted);
  l->DrawLine(-1.305,15,-1.305,4500.);
  l->DrawLine(+1.305,15,+1.305,4500.);

  tex->DrawLatex(0.50,0.89,"2023D");
  tex->DrawLatex(0.50,0.84,cp);

  gPad->RedrawAxis();

  c1->SaveAs(Form("pdf/drawL1Tau_1D_%s.pdf",cp));

  // Log scale
  c1->SetLogx();
  c1->SetLogy();
  h1c->GetXaxis()->SetRangeUser(300,4000);
  h1d->GetXaxis()->SetRangeUser(300,4000);
  h->GetXaxis()->SetRangeUser(300,4000);
  h->GetYaxis()->SetRangeUser(5e-4,1e1);
  c1->Update();
  c1->SaveAs(Form("pdf/drawL1Tau_1D_log_%s.pdf",cp));
  
  c2c->SaveAs(Form("pdf/drawL1Tau_2D_%s_23C.pdf",cp));
  c2d->SaveAs(Form("pdf/drawL1Tau_2D_%s_23D.pdf",cp));


  printPars(f1c,"23Cv123_1c");
  printPars(f2c,"23Cv123_2c");
  printPars(f1d,"23D_1c");
  printPars(f2d,"23D_2c");
}
