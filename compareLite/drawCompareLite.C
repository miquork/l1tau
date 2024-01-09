// Purpose: Draw tag-and-probe and direct matching results for reco-reco match
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMultiGraph.h"

#include "../tdrstyle_mod22.C"

TGraphErrors *tagandprobe(TProfile *pt, TProfile *pa, TProfile *pb,
			  TProfile *pd) {

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != pa->GetNbinsX()+1; ++i) {
    //double x = pa->GetBinContent(i) * pa->GetBinCenter(i); // <A/tag>*<tag>~<A>
    double x = pt->GetBinContent(i); // <A>
    double ex = pt->GetBinError(i); // d<A>
    if (pa->GetBinContent(i)!=0 && pb->GetBinContent(i)!=0) {
      double y = pb->GetBinContent(i) / pa->GetBinContent(i);
      double ey = 2*pd->GetBinError(i); // x2?
      //double ex = pa->GetBinError(i) * pa->GetBinCenter(i);
      int n = g->GetN();
      g->SetPoint(n, x, y);
      g->SetPointError(n, ex, ey);
    }
  } // for i
  return g;
} // tagandprobe

TGraphErrors *directmatch(TProfile *pt, TProfile *p, bool invert = false) {

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != p->GetNbinsX()+1; ++i) {
    //double x = p->GetBinCenter(i); // <A>, or <B>
    //if (invert) x = p->GetBinContent(i) * x; // <A/B>*<B>~<A>
    double x = pt->GetBinCenter(i); // <A>
    double ex = pt->GetBinError(i); // d<A>
    if (p->GetBinContent(i)!=0) {
      double y = p->GetBinContent(i); // <A/B> or <B/A>
      if (invert) y = 1./y; // 1/<A/B>~<B>/<A>
      double ey = p->GetBinError(i);
      //double ex = p->GetBinError(i) * x;
      int n = g->GetN();
      g->SetPoint(n, x, y);
      g->SetPointError(n, ex, ey);
    }
  } // for i
  return g;
} // directmatch

TGraphErrors *directaverage(TProfile *pt, TProfile *p) {

  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != p->GetNbinsX()+1; ++i) {
    if (p->GetBinContent(i)!=0) {

      //double x = p->GetBinCenter(i); // 0.5*(B+A)
      double y = p->GetBinContent(i); // (B-A)/(B+A)
      // 2x  = B+A
      // 2xy = B-A
      // => 2B = 2x+2xy, 2A = 2x-2xy
      // B/A = (x+xy)/(x-xy) = (1+y)/(1-y)

      //x = x-x*y; // <A>
      double x = pt->GetBinContent(i); // <A>
      double ex = p->GetBinContent(i); // d<A>
      y = (1+y)/(1-y);  // <B>/<A>
      double ey = 2*p->GetBinError(i); // x2?
      //double ex = p->GetBinError(i) * x;
      int n = g->GetN();
      g->SetPoint(n, x, y);
      g->SetPointError(n, ex, ey);
    }
  } // for i
  return g;
} // directmatch

void drawCompareLite(string run = "2023D") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *crun = run.c_str();
  TFile *f = new TFile(Form("compareLite/compareLite_%s.root",crun),"READ");
  assert(f && !f->IsZombie());

  string sA = "19Dec2023";
  string sB = "Prompt23";//"22Sep";
  const char *cA = sA.c_str();
  const char *cB = sB.c_str();
  
  double xmin = 600;
  double xmax = 4000;//3500;

  // Background canvas
  TH1D *h = tdrHist("h","#LTp_{T}(B)#GT / #LTp_{T}(A)#GT", 0.98,1.13,
		    "#LTp_{T}(A)#GT (GeV)",xmin,xmax);
  if (run=="2023Cv123") {}
  if (run=="2023D") { h->GetYaxis()->SetRangeUser(0.92,1.07); }
  
  lumi_136TeV = run.c_str();//"2023Cv123";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  
  // Tag-and-probe
  TProfile *pta_tp = (TProfile*)f->Get("pta_tp"); assert(pta_tp);
  TProfile *pa_tp = (TProfile*)f->Get("pa_tp"); assert(pa_tp);
  TProfile *pb_tp = (TProfile*)f->Get("pb_tp"); assert(pb_tp);
  TProfile *pd_tp = (TProfile*)f->Get("pd_tp"); assert(pd_tp);
  TGraphErrors *g = tagandprobe(pta_tp,pa_tp,pb_tp,pd_tp);

  // Direct match
  TProfile *pta_dm = (TProfile*)f->Get("pta_dm"); assert(pta_dm);
  TProfile *pa_dm = (TProfile*)f->Get("pa_dm"); assert(pa_dm);
  TGraphErrors *ga = directmatch(pta_dm,pa_dm);
  TProfile *ptb_dm = (TProfile*)f->Get("ptb_dm"); assert(ptb_dm);
  TProfile *pb_dm = (TProfile*)f->Get("pb_dm"); assert(pb_dm);
  TGraphErrors *gb = directmatch(ptb_dm,pb_dm,true);
  TProfile *ptd_dm = (TProfile*)f->Get("ptd_dm"); assert(ptd_dm);
  TProfile *pd_dm = (TProfile*)f->Get("pd_dm"); assert(pd_dm);
  TGraphErrors *gd = directaverage(ptd_dm,pd_dm);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(xmin,1,xmax,1);
  //l->DrawLine(xmin,0.99,xmax,0.99);
  //l->DrawLine(xmin,1.01,xmax,1.01);

  TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.045);
  t->DrawLatex(0.19,0.75,"|#eta| < 1.3");
  t->DrawLatex(0.19,0.25,Form("B: %s",cB));
  t->DrawLatex(0.19,0.20,Form("A: %s",cA));

  tdrDraw(ga,"Pz",kOpenTriangleDown,kBlue);
  tdrDraw(gd,"Pz",kOpenDiamond,kGreen+2);
  tdrDraw(gb,"Pz",kOpenTriangleUp,kRed);
  tdrDraw(g,"Pz",kFullCircle,kBlack);

  // SPRH -3% variation
  //if (!fhh) fhh = new TF1("fhh","[p0]+[p1]*pow(x/[p2],[p3])/(1+pow(x/[p2],[p3]))*(1-pow(x/[p2],-[p3]))",15,4500);
  //fhh->SetParameters(-0.7938, -0.5798, 396.1, 1.412);
  
  //TF1 *f1 = new TF1("f1","[0]+[1]*x/3000.",600,3000.);
  TF1 *f1 = new TF1("f1","[0]+[1]*0.01*(-0.798-0.5798*pow(x/396.1,1.412)/(1+pow(x/396.1,1.412))*(1-pow(x/396.1,-1.412)))+[2]*-0.1*x/3000.",600,3300);
  f1->SetParameters(1,1,1);
  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(g);
  mg->Add(gd);
  f1->SetLineColor(kGreen+2);
  f1->SetLineWidth(2);
  f1->SetParameters(1.08,-0.05);
  mg->Fit(f1,"RN");
  f1->Draw("SAME");

  TLegend *leg = new TLegend(0.40,0.90-0.05*4,0.65,0.90, "", "brNDC");
  leg->SetBorderSize(0); leg->SetFillStyle(kNone); leg->SetTextSize(0.045);
  leg->AddEntry(g, "Tag-and-probe", "PLE");
  leg->AddEntry(ga, "Direct match vs A", "PLE");
  leg->AddEntry(gd, "Direct match vs (A+B)/2", "PLe");
  leg->AddEntry(gb, "Direct match vs B", "PLE");
  leg->Draw();

  c1->SaveAs(Form("pdf/drawCompareLite_%s_vs_%s_TnP_%s.pdf",cB,cA,crun));
} // drawTnP
  

