
// fit the with current data with the function with signal, but with the restriction of the background parameters (back constant, back slope) from the data without current

using namespace RooFit;
using namespace RooStats;

void bayesianAnalyser_NuisancePar(){
      
  //withW->Print();
  //TFile *fHisto = new TFile("energyHistograms.root");
  
    //First, define the observable for the analysis
  RooRealVar energy("energy","energy",7000.,8600.);
  
//  TH1F *withH = (TH1F*)fHisto->Get("withCurrentSum");
//  TH1F *noH   = (TH1F*)fHisto->Get("noCurrentSmallSum"); 
//
//  noH->Rebin(reBin);
//  withH->Rebin(reBin);
// 
//  // Create a binned dataset that imports contents of TH1 and associates its contents to observable energy
//  RooDataHist noDH("noDH","noDH",energy,Import(*noH)) ;
//  RooDataHist withDH("withDH","withDH",energy,Import(*withH)) ;
  
  
  
  TFile fFit("fitWithCurrentWithSignal-ConstBG.root"); // in this file the background is set to the values without current
  fFit.cd();
  RooWorkspace *wConst = (RooWorkspace*)fFit.Get("w");
  
  RooAbsData *withAD = wConst->data("withDH");
  RooAbsData *noAD = wConst->data("noDH");
  
  
  wConst->Print();
  
  // backC no current 3.18384e+04   2.52693e+02 -> so the relative error is around 0,008 = 0.8 %
  // backSl      -5.17610e-02   1.09214e-02
  
  RooRealVar *backC = wConst->var("backC"); // backC and backSl are the parameters from the fit without current
  RooRealVar *backSl = wConst->var("backSl");


  RooWorkspace* nuisW = new RooWorkspace("nuisW"); 

 // ------------------------------- Constructing the pdf with a forbidden signal component ------------------------------
  
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
  

  //the variables for the background function are defined earlier
  
  RooChebychev backgF("backgF","Background",energy,RooArgSet(*backSl));
  
  // PEP violating tranistion

  RooRealVar meanForbidden("meanForbidden","mean of the forbidden tranistion", 7729, 7728.,7730.);
  RooGaussian gaussForbidden("gaussForbidden","Forbidden pdf",energy,meanForbidden,sigmaCuKa);

  RooRealVar Nsig("Nsig","signal Events",0,-10000.,10000);
  
  
  // -------------------------------- now defining background as a nuisance parameter with an uncertainty  and the PDFs ------------------------------------
  
  
  RooRealVar backC_alpha("backC_alpha","Dimension of systematic variation of number of BG events",1.,0.01,10.);
  RooFormulaVar backC_nuis("backC_nuis","@0*@1",RooArgList(*backC,backC_alpha)); 

  RooRealVar one("one","one",1.);
  RooRealVar backC_syst("backC_syst","The systematic uncertainty on Nbkg",0.008); // defines the uncertainty 
  RooGaussian constr_backC("constr_backC","Background uncertainty constraint",backC_alpha,one,backC_syst); // gaussian distribution of backC_alpha with mean = 1 and sigma = backC_syst
  
  RooAddPdf PDFtot_nuis_unconstr("PDFtot","PDFtot",RooArgList(gaussCuKa1,gaussCuKa2,gaussNiKa1,gaussNiKa2,backgF,gaussForbidden),RooArgList(cuKa1N,cuKa2N,niKa1N,niKa2N,backC_nuis,Nsig)); // here backC_nuis
  RooProdPdf PDFtot_nuis("PDFtot_nuis","PDFtot_nuis",RooArgList(PDFtot_nuis_unconstr,constr_backC)); // this is to make backC_nuis gaussian distributed with mean backC and sigma = backC/100
  
  
  //RooAddPdf PDFtot_nuis("PDFtot_nuis","PDFtot_nuis",RooArgList(gaussCuKa1,gaussCuKa2,gaussNiKa1,gaussNiKa2,backgF,gaussForbidden),RooArgList(cuKa1N,cuKa2N,niKa1N,niKa2N,*backC,Nsig)); for testing with 
  //the constant BG value
  
  
  
  // ----------------------------set some things constant --------------------------------------------------------------
  
  CuKa2Diff.setConstant(kTRUE);
  CuKa2Ratio.setConstant(kTRUE);
  NiKa2Diff.setConstant(kTRUE);
  NiKa2Ratio.setConstant(kTRUE);
  
  
//backC        3.45985e+04   2.62129e+02 for no current histogram fit without signal gaussian
  //backC->setVal();
  //backSl->setVal(backSlDouble);
  //backSl->setConstant(kTRUE);
  
  //backC->Print();
  backC->setConstant(kTRUE);
  backSl->setConstant(kTRUE);
  meanForbidden.setConstant(kTRUE);
  
  meanCuKa1.setConstant(kTRUE);
  meanNiKa1.setConstant(kTRUE);
  
  
  
  // -------------------------------   end of pdf definition - start of definition os nuisance parameters -------------------------------------------------------------
 /* 
   //Assume an uncertainty on the number of background events
  //Construct the uncertainty with a lognormal assumption
  RooRealVar Nbkg_alpha("Nbkg_alpha","Dimension of systematic variation of Nbkg",1.,0.01,10.);
  RooFormulaVar Nbkg_nuis("Nbkg_nuis","@0*@1",RooArgList(*Nbkg,Nbkg_alpha));

  //Now prepare a gaussian for the nuisance parameter, to be multiplied to the total PDF
  RooRealVar one("one","one",1.);
  RooRealVar Nbkg_syst("Nbkg_syst","The systematic uncertainty on Nbkg",0.3);    //10% uncertainty
  RooGaussian constr_Nbkg("constr_Nbkg","Background uncertainty constraint",one,Nbkg_alpha,Nbkg_syst); // i think "one" should be mean, the variable should be Nbkg_alpha, and sigma = 0.3
// this  

  //Now construct the total PDF
  RooAddPdf PDFtot_nuis_unconstr("PDFtot_nuis_unconstr","PDFtot_nuis_unconstr",RooArgList(*gauss,*exponential),RooArgList(*Nsig,Nbkg_nuis));// the exponential bg function is multiplied by
// Nbgk_nuis = Nbgk * Nbgk_alpha , but Nbgk_alpha == 1 the starting value
 
  //Now add the gaussian constraint to the total PDF
  RooProdPdf PDFtot_nuis("PDFtot_nuis","PDFtot_nuis",RooArgList(PDFtot_nuis_unconstr,constr_Nbkg)); // the pdf of signal + bg * Nbgk_alpha is multiplied by: a gaussian for !Nbgk_alpha! -> 
  * so only this variable is affected by the gaussian; with mean 1 and sigma 0.3 -> this is the number which is multiplied to Nbgk to get Nbgk_nuis -> so finally Nbgk_nuis is gaussian distributed with mean Nbgk

  PDFtot_nuis.fitTo(*w->data("PDFtotData"), RooFit::Constrain(RooArgSet(Nbkg_alpha)),Extended(1));
 */


  
  //RooAbsData *withAD = noW->data("withDH");
  //RooAbsData *noAD = noW->data("noDH");// data with current is also stored in the workspace with the fit without current
  
  PDFtot_nuis.fitTo(*withAD);
  //PDFtot_nuis.fitTo(withDH, RooFit::Constrain(RooArgSet(backC_alpha)),Extended(1)); // Constrain(const RooArgSet&pars) – Apply constraints to listed parameters in likelihood using internal constrains in p.d.f
  
  //backC->Print();
  
  nuisW->import(PDFtot_nuis);
  
   //Configure the model, we need both the S+B and the B only models
  ModelConfig sbModel;
  sbModel.SetWorkspace(*nuisW);
  sbModel.SetPdf("PDFtot_nuis"); // in this model, only the background parameters are constant and the mean of the forbidden transition
  sbModel.SetName("S+B Model");
  RooRealVar* poi = nuisW->var("Nsig");
  poi->setRange(0.,700.);  //i think this also sets the lower and upper limit for the prior!! 
  sbModel.SetParametersOfInterest(*poi);

  ModelConfig *bModel = (ModelConfig*) sbModel.Clone();
  bModel->SetPdf("PDFtot_nuis");
  bModel->SetName("B model: poi=0");      
  poi->setVal(0);// for the background only model, nsig = 0
  bModel->SetSnapshot(*poi);
  
  nuisW->factory("Uniform::prior(Nsig)");
  sbModel.SetPriorPdf(*nuisW->pdf("prior"));

  //Construct the bayesian calculator
  BayesianCalculator bc(*(wConst->data("withDH")), sbModel); // rooabsdata but still as data histogram somehow
  bc.SetConfidenceLevel(0.997);
  bc.SetLeftSideTailFraction(0.); // set the fraction of probability content on the left tail Central limits use 0.5 (default case) for upper limits it is 0 and 1 for lower limit
  SimpleInterval* bcInt = bc.GetInterval();
  
  cout << "################" << endl;

  //Now let's see what the bayesian limit is
  cout << "Bayesian upper limit on Nsig = " << bcInt->UpperLimit() << endl;
  
  Nsig.Print();
  backC->Print();
  backSl->Print();
  cuKa1N.Print();
  meanCuKa1.Print();
  sigmaCuKa.Print();

  // plot now the result of the scan 
  //First the CLs
  //Then the Bayesian posterior
  RooPlot *bcPlot = bc.GetPosteriorPlot();

  // plot in a new canvas with style
  TCanvas dataCanvas("dataCanvas");
  bcPlot->Draw();
  dataCanvas.SaveAs("bayesianAnalysis_Nuisance.gif");
    
}
