//////////////////////////////////////////////////////////////////////////
//
// First exercise for the CMS DAS
// 
// Creating a dataset, fitting saving into a Workspace
//
// pdf = gauss(x,m,s) + exp(x,tau)
//
//
// 2014 - Mario Pelliccioni, Luca Lista
// 
/////////////////////////////////////////////////////////////////////////
//gSystem->Load("libRooFit");


using namespace RooFit;

void fitNoSignalArbHist(TString fitHisto)
{
// how to:?
// make a function for the spectrum (constBG + cuka1 + cuka2 + nika1 + nika2) and fit it to the no current data
// then fix the parameter for the constBG and add a signal gaussian with fixed width and mean and variable gain
// add the uncertainty of the constBG as nuisance parameter
// also fit the cuka1, cuka2, nika1, nika2 parameters new!

// or go full simultaneous fit and fit both histograms simultaneously while leaving the slope?? in common

  // get the histograms

  TFile *fIN = new TFile("energyHistograms.root");
  
  //TH1F *withH = (TH1F*)fIN->Get("withCurrentSum");
  //TH1F *noH   = (TH1F*)fIN->Get("noCurrentSmallSum"); 

  TH1F *fitH   = (TH1F*)fIN->Get(fitHisto);

  //noH->Rebin(25);
  //withH->Rebin(25);
  fitH->Rebin(25);

  //withH->Draw();


  //First, define the observable for the analysis
  RooRealVar energy("energy","energy",7000.,8600.);


  //Construct the P.D.F.

  //CuKa1
  
  RooRealVar meanCuKa1("meanCuKa1","mean of Cu Ka1 gaussian",8047.78,8040.,8080.);
  RooRealVar sigmaCuKa("sigmaCuKa","width of Cu Ka1 gaussian",75.,70.,90.);
  RooGaussian gaussCuKa1("gaussCuKa1","Cu Ka1 PDF",energy,meanCuKa1,sigmaCuKa);
 
  RooRealVar cuKa1N("cuKa1N","cu Ka1 Events",15000,0,100000);

  //Cuka2

  RooRealVar CuKa2Diff("CuKa2Diff","diff Ka1 - Ka2",19.95,19.,20.);
  RooRealVar CuKa2Ratio("CuKa2Ratio","ratio Ka1 / Ka2",0.51,0.,1.);

  RooGenericPdf meanCuKa2("meanCuKa2","diff Cu Ka1 - Ka2  PDF","meanCuKa1 - CuKa2Diff",RooArgSet(meanCuKa1,CuKa2Diff));
  RooGaussian gaussCuKa2("gaussCuKa2","Cu Ka2 PDF",energy,meanCuKa2,sigmaCuKa); 

  RooGenericPdf cuKa2N("CuPdfRatio","ratio Cu Ka1 / Ka2  PDF","cuKa1N*CuKa2Ratio",RooArgSet(cuKa1N,CuKa2Ratio));

  //NiKa1

  RooRealVar meanNiKa1("meanNiKa1","mean of Ni Ka1 gaussian",7478.15,7470.,7500.);
  RooRealVar sigmaNiKa("sigmaNiKa","width of Ni Ka1 gaussian",70.,50.,90.);
  RooGaussian gaussNiKa1("gaussNiKa1","Ni Ka1 PDF",energy,meanNiKa1,sigmaNiKa); 
 
  RooRealVar niKa1N("niKa1N","Nickel Ka1 Events",200,0,1000);

  //Nika2

  RooRealVar NiKa2Diff("NiKa2Diff","diff Ka1 - Ka2",17.26,17.,18.);
  RooRealVar NiKa2Ratio("NiKa2Ratio","ratio Ka1 / Ka2",0.51,0.,1.);

  RooGenericPdf meanNiKa2("meanNiKa2","diff Ni Ka1 - Ka2  PDF","meanNiKa1 - NiKa2Diff",RooArgSet(meanNiKa1,NiKa2Diff));
  RooGaussian gaussNiKa2("gaussNiKa2","Cu Ka2 PDF",energy,meanNiKa2,sigmaNiKa);

  RooGenericPdf niKa2N("NiPdfRatio","ratio Ni Ka1 / Ka2  PDF","niKa1N*NiKa2Ratio",RooArgSet(niKa1N,NiKa2Ratio)); 

  //Background
  

  RooRealVar backC("backC","background constant",500,0,100000);
  RooRealVar backSl("backSl","background Slope",0,-1,1);

  
  //RooUniform backgF("backgF","background PDF",energy);
  RooChebychev backgF("backgF","Background",energy,RooArgSet(backSl));

  //Now define the background P.D.F, a simple exponential
 // RooRealVar tau("tau","exponential function parameter",-0.05,-10.,-0.001);
 // RooExponential exponential("exponential","Background PDF",mass,tau);

  //Now construct the total PDF. We need to define the number of signal and background events in the model
  //
  //for UL calculation use Nsig = 5, Nbkg = 100
  //for mH calculation use Nsig = 50, Nbkg = 450
  //for systematics inclusion use Nsig = 20, Nbkg = 100, also, with the width set to 5 GeV!!
  //RooRealVar NoCurrentN("NoCurrentN","Number of events without current",1000.,0.,2000.);
  //RooRealVar Nbkg("Nbkg","Number of background events",450.,0.,1000.);

  RooAddPdf PDFtot("PDFtot","PDFtot",RooArgList(gaussCuKa1,gaussCuKa2,gaussNiKa1,gaussNiKa2,backgF),RooArgList(cuKa1N,cuKa2N,niKa1N,niKa2N,backC));

  CuKa2Diff.setConstant(kTRUE);
  CuKa2Ratio.setConstant(kTRUE);
  NiKa2Diff.setConstant(kTRUE);
  NiKa2Ratio.setConstant(kTRUE);
  
  // set the forbidden stuff to constant first

  //Nsig.setConstant(kTRUE); // set to 0 for no current data -> but apparently this must not be set to zero for the model that is later used for the current data!!
  //meanForbidden.setConstant(kTRUE);

  //Now generate a sample with the total PDF
  //RooDataSet *data = PDFtot.generate(RooArgSet(energy),10000,Extended(1));

  // or take the real data

  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  //RooDataHist noDH("noDH","noDH",energy,Import(*noH));
  //RooDataHist withDH("withDH","withDH",energy,Import(*withH));
  RooDataHist fitDH("fitDH","fitDH",energy,Import(*fitH));
  

  //Now fit the PDF to the data
  //For low statistics, fix the mean and the width
  //mean.setConstant(kTRUE);
  //sigma.setConstant(kTRUE);
  PDFtot.fitTo(fitDH);

  
  //Now plot the data and the fitted PDF
  RooPlot *energyFrame = energy.frame(50);
  //noDH.plotOn(energyFrame);
  fitDH.plotOn(energyFrame);
  PDFtot.plotOn(energyFrame);

  //One can also plot the single components of the total PDF, like the background component
  PDFtot.plotOn(energyFrame,Components(backgF),LineStyle(kDashed),LineColor(kRed));
  PDFtot.plotOn(energyFrame,Components(gaussCuKa1),LineStyle(kDashed),LineColor(kBlack));
  PDFtot.plotOn(energyFrame,Components(gaussCuKa2),LineStyle(kDashed),LineColor(kYellow));
  PDFtot.plotOn(energyFrame,Components(gaussNiKa1),LineStyle(kDashed),LineColor(kGreen));

  //Actually plot the result
  TCanvas c1;
  c1.cd();
  energyFrame->Draw();
  c1.SaveAs("fitArbHisto.gif");

  // Print values of mean and sigma (that now reflect fitted values and errors, unless you fixed them)
  meanCuKa1.Print();
  meanCuKa2.Print();
  sigmaCuKa.Print();


  //Now save the data and the PDF into a Workspace, for later use for statistical analysis
  //and save the workspace into a ROOT file
/*
  RooWorkspace* w = new RooWorkspace("w") ;
  w->import(noDH);
  w->import(withDH);
  w->import(PDFtot);
 

  TFile fOut("fitNoCurrent.root","RECREATE");
  fOut.cd();
  w->Write();
  fOut.Close();
  fIN->Close();
*/
}

