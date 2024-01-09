{

  //gROOT->ProcessLine(".exception");

  gROOT->ProcessLine(".L compareLite/tools.C+g");
  // Link JEC libraries that are modified to work stand-alone
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  //gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L compareLite/compareLite.C+g");
  gROOT->ProcessLine(".L compareLite/drawCompareLite.C+g");

  compareLite("2023Cv123");
  //compareLite("2023D");
  drawCompareLite("2023Cv123");
  //drawCompareLite("2023D");
}
