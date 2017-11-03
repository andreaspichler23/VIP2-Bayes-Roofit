//////////////////////////////////////////////////////////////////////////
//
// Third exercise for the CMS DAS
// 
// ULs with the CLs and BayesianCalculator methods
//
// pdf = gauss(x,m,s) + exp(x,tau)
//
//
// 2014 - Mario Pelliccioni, Luca Lista
// 
/////////////////////////////////////////////////////////////////////////

//gSystem->Load("libRooFit");
//gSystem->Load("libRooStats");

using namespace RooFit;
using namespace RooStats;

void bayesianAnalyser()
{

  //Open the rootfile and get the workspace from the exercise_0
  TFile fIn("fitWithCurrentWithSignal-ConstBG.root"); // in this file the background is set to the values without current
  fIn.cd();
  RooWorkspace *w = (RooWorkspace*)fIn.Get("w");
  
  w->Print();

  //You can set constant parameters that are known
  //If you leave them floating, the fit procedure will determine their uncertainty
  RooRealVar *nobackC = w->var("backC");
  nobackC->setConstant(kTRUE);
  //nobackC->Print();
  
  RooRealVar *nobackSl = w->var("backSl");
  nobackSl->setConstant(kTRUE);
  //nobackSl->Print();
  
  RooRealVar *cuKa1N = w->var("cuKa1N");
  RooRealVar *meanCuKa1 = w->var("meanCuKa1");
  RooRealVar *sigmaCuKa = w->var("sigmaCuKa");
  RooRealVar *Nsig = w->var("Nsig");
  
  

  //Configure the model, we need both the S+B and the B only models
  ModelConfig sbModel;
  sbModel.SetWorkspace(*w);
  sbModel.SetPdf("PDFwith"); // in this model, only the background parameters are constant and the mean of the forbidden transition
  sbModel.SetName("S+B Model");
  RooRealVar* poi = w->var("Nsig");
  poi->setRange(0.,700.);  //i think this also sets the lower and upper limit for the prior!! 
  sbModel.SetParametersOfInterest(*poi);

  ModelConfig *bModel = (ModelConfig*) sbModel.Clone();
  bModel->SetPdf("PDFwith");
  bModel->SetName("B model: poi=0");      
  poi->setVal(0);// for the background only model, nsig = 0
  bModel->SetSnapshot(*poi);

  
 // --------------------------------  start the frequentist calculation  -----------------------------------

//  FrequentistCalculator fc(*(w->data("withDH")), *bModel, sbModel);
//  fc.SetToys(100,100); // toys for null, toys for alternative hypothesis
//
//  //Create hypotest inverter passing the desired calculator 
//  HypoTestInverter calc(fc);
//
//  // set confidence level (e.g. 95% upper limits)
//  calc.SetConfidenceLevel(0.95);
//
//  //use CLs
//  calc.UseCLs(true);
//
//  //reduce the noise
//  calc.SetVerbose(false);
//
//  //configure ToyMC Samler
//  ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
//
//  //profile likelihood test statistics 
//  ProfileLikelihoodTestStat profll(*(sbModel.GetPdf()));
//  //for CLs (bounded intervals) use one-sided profile likelihood
//  profll.SetOneSided(true);
//
//  //set the test statistic to use 
//  toymcs->SetTestStatistic(&profll);
//
//  int npoints = 3;  // number of points to scan
//  // min and max (better to choose smaller intervals)
//  double poimin = poi->getMin();
//  double poimax = poi->getMax();
//
//  cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << endl;
//  calc.SetFixedScan(npoints,poimin,poimax);
//  
//  HypoTestInverterResult * r = calc.GetInterval();
//  double upperLimit = r->UpperLimit();

// ------------------------- start the bayesian calculator -----------------------------------------------------  

  //Example using the BayesianCalculator
  //Now we also need to specify a prior in the ModelConfig
  //To be quicker, we'll use the PDF factory facility of RooWorkspace
  //NB!! For simplicity, we are using a flat prior, but this doesn't mean it's the best choice!
  w->factory("Uniform::prior(Nsig)");
  sbModel.SetPriorPdf(*w->pdf("prior"));

  //Construct the bayesian calculator
  BayesianCalculator bc(*(w->data("withDH")), sbModel); // rooabsdata but still as data histogram somehow
  bc.SetConfidenceLevel(0.997);
  bc.SetLeftSideTailFraction(0.); // set the fraction of probability content on the left tail Central limits use 0.5 (default case) for upper limits it is 0 and 1 for lower limit
  SimpleInterval* bcInt = bc.GetInterval();
/*
  //Now let's print the result of the two methods
  //First the CLs
  cout << "################" << endl;
  cout << "The computed CLs upper limit is: " << upperLimit << endl;

  //Compute expected limit
  cout << "Expected upper limits, using the B (alternate) model : " << endl;
  cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << endl;
  cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
  cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
*/
  cout << "################" << endl;

  //Now let's see what the bayesian limit is
  cout << "Bayesian upper limit on Nsig = " << bcInt->UpperLimit() << endl;
  
  Nsig->Print();
  nobackC->Print();
  nobackSl->Print();
  cuKa1N->Print();
  meanCuKa1->Print();
  sigmaCuKa->Print();

  // plot now the result of the scan 
  //First the CLs
  //HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",r);
  //Then the Bayesian posterior
  RooPlot *bcPlot = bc.GetPosteriorPlot();

  // plot in a new canvas with style
  TCanvas dataCanvas("dataCanvas");
  //dataCanvas.Divide(2,1);
  //dataCanvas.SetLogy(false);
  dataCanvas.cd(1);
  //plot->Draw("2CL");
  //dataCanvas.cd(2);
  bcPlot->Draw();

  dataCanvas.SaveAs("bayesianAnalysis.gif");

  TFile *outF = new TFile("bayesFile.root","recreate");
  outF->cd();
  
  bcPlot->Write();
  
  outF->Close();
  
}

//TFile *f = new TFile("bayesFile.root")
//f->ls()
//RooPlot *pl = (RooPlot*)f->Get("frame_Nsig_4d1d430")