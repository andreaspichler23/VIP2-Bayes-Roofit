
// fit the with current data with the function with signal, but with the restriction of the background parameters (back constant, back slope) from the data without current

using namespace RooFit;
using namespace RooStats;

void fitCurrentwithSignal(){
      
  //withW->Print();
  TFile fno("fitNoSignalNoCurrent.root");
  fno.cd();
  RooWorkspace *noW = (RooWorkspace*)fno.Get("w");
  
  noW->Print();
  
  RooRealVar *noBackC = noW->var("backC");
  RooRealVar *noBackSl = noW->var("backSl");

  
  Double_t noBackCDouble = noBackC->getValV();
  Double_t noBackSlDouble = noBackSl->getValV();
  
  //First, define the observable for the analysis
  RooRealVar energy("energy","energy",7000.,8600.);
  //RooRealVar energy("energy","energy",7600.,7850.);


  //--------------------------------------------------Construct the P.D.F. with a forbidden PEP signal component ------------------------------------

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
  
  RooChebychev backgF("backgF","Background",energy,RooArgSet(backSl));
  
  // PEP violating tranistion

  RooRealVar meanForbidden("meanForbidden","mean of the forbidden tranistion", 7729, 7728.,7730.);
  RooGaussian gaussForbidden("gaussForbidden","Forbidden pdf",energy,meanForbidden,sigmaCuKa);

  RooRealVar Nsig("Nsig","signal Events",0,-10000.,10000);
  
  
  
  RooAddPdf PDFwith("PDFwith","PDFwith",RooArgList(gaussCuKa1,gaussCuKa2,gaussNiKa1,gaussNiKa2,backgF,gaussForbidden),RooArgList(cuKa1N,cuKa2N,niKa1N,niKa2N,backC,Nsig));
  
  CuKa2Diff.setConstant(kTRUE);
  CuKa2Ratio.setConstant(kTRUE);
  NiKa2Diff.setConstant(kTRUE);
  NiKa2Ratio.setConstant(kTRUE);
  
  // -------------------------------- set the background parameter to the values of the data without current -----------
  
  
  backC.setVal(noBackCDouble);
  backC.setConstant(kTRUE); // for no current data:   3.18384e+04
  backSl.setVal(noBackSlDouble);
  backSl.setConstant(kTRUE); // for no current data:   -5.17610e-02
  
  backC.Print();

  meanForbidden.setConstant(kTRUE);
  
  
  // -------------------------------   end of pdf definition -------------------------------------------------------------
  
  
  
  noBackC->Print();
  noBackSl->Print();
 
 
  
  RooAbsData *withAD = noW->data("withDH");
  RooAbsData *noAD = noW->data("noDH");// data with current is also stored in the workspace with the fit without current
  
  PDFwith.fitTo(*withAD);
  
  backC.Print();
  
  RooPlot *energyFrame = noW->var("energy")->frame();
  
//  TAxis *xaxis = energyFrame->GetXaxis();
//  xaxis->SetRangeUser(7600,7800);
//  
//  TAxis *yaxis = energyFrame->GetYaxis();
//  yaxis->SetRangeUser(580,700);
  
  withAD->plotOn(energyFrame);
  PDFwith.plotOn(energyFrame);
  
  PDFwith.plotOn(energyFrame,Components(gaussForbidden),LineStyle(kDashed),LineColor(kOrange));
  PDFwith.plotOn(energyFrame,Components(backgF),LineStyle(kDashed),LineColor(kRed));
  PDFwith.plotOn(energyFrame,Components(gaussCuKa1),LineStyle(kDashed),LineColor(kBlack));
  PDFwith.plotOn(energyFrame,Components(gaussCuKa2),LineStyle(kDashed),LineColor(kYellow));
  PDFwith.plotOn(energyFrame,Components(gaussNiKa1),LineStyle(kDashed),LineColor(kGreen));
  
  
  //Actually plot the result
  TCanvas c1;
  c1.cd();
  gPad->SetLogy();
  energyFrame->Draw();
  c1.SaveAs("fitCurrentwithSignal.gif");
  
  //Now save the data and the PDF into a Workspace, for later use for statistical analysis
  //and save the workspace into a ROOT file
  RooWorkspace* w = new RooWorkspace("w") ;
  w->import(*withAD);
  w->import(*noAD);
  w->import(PDFwith);
  //w->import(*energyFrame);
  
  TFile fOut("fitWithCurrentWithSignal-ConstBG.root","RECREATE");
  fOut.cd();
  energyFrame->Write();
  w->Write();
  fOut.Close();
  fno.Close();
  
  cout << endl << "Number of signal counts: " << endl;
  Nsig.Print();
  cout << endl;
    
}
//TFile *f = new TFile("fitWithCurrentWithSignal-ConstBG.root")
// f->ls()
//RooPlot *pl = (RooPlot*)f->Get("frame_energy_3268c50") change the name by the outcome of ls()