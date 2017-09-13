//////////////////////////////////////////////////////////////////////////
//
// Second exercise for the CMS DAS
// 
// Hypothesis testing with the Profile Likelihood method
//
// pdf = gauss(x,m,s) + exp(x,tau)
//
//
// 2014 - Mario Pelliccioni, Luca Lista
// 
/////////////////////////////////////////////////////////////////////////

//gSystem->Load("libRooFit");
//gSystem->Load("libRooStats");

// How likely is it that we measure this spectrum when the null hypothesis (Nsig = 0, no signal) is true?

using namespace RooFit;
using namespace RooStats;

void getPValue()
{
  // get the histogram with current

  //TFile *fIN = new TFile("energyHistograms.root");
  
  //TH1F *withH = (TH1F*)fIN->Get("withCurrentSum");
  //withH->Rebin(25);

  //Open the rootfile and get the workspace from the previous exercise (fitNoCurrent.C)
  TFile fIn("fitnoCurrent.root");
  fIn.cd();
  RooWorkspace *w = (RooWorkspace*)fIn.Get("w");

  //You can set constant parameters that are known
  //If you leave them floating, the fit procedure will determine their uncertainty
  //w->var("mean")->setConstant(kTRUE);
  //w->var("sigma")->setConstant(kTRUE);
  //w->var("tau")->setConstant(kTRUE);
  //w->var("Nbkg")->setConstant(kTRUE);
  w->var("backC")->setConstant(kTRUE);
  w->var("backSl")->setConstant(kTRUE);

  w->var("cuKa1N")->Print();

  //Set the RooModelConfig and let it know what the content of the workspace is about
  ModelConfig model;
  model.SetWorkspace(*w);
  model.SetPdf("PDFtot");

  // here we explicitly set the value of the parameters for the null.  
  // We want no signal contribution, eg. Nsig = 0
  // How likely is it that we measure this spectrum when the null hypothesis (Nsig = 0, no signal) is true?
  RooRealVar* Nsig = w->var("Nsig");
  RooArgSet poi(*Nsig); // poi = parameters of interest: paramter which i want the p-value for; nuisance parameters will be all other parameters of the model (all but the ones to be set constant?)
  RooArgSet *nullParams = (RooArgSet*) poi.snapshot(); // makes a clone of poi 
  nullParams->setRealValue("Nsig",0.); // set the parameter for the null hypothesis

  //Build the profile likelihood calculator
  ProfileLikelihoodCalculator plc; 

  plc.SetData(*(w->data("withDH"))); // virtual void RooStats::CombinedCalculator::SetData 	( 	RooAbsData &  	data	) 	
  // RooAbsData * RooWorkspace::data 	( 	const char *  	name	): Retrieve dataset (binned or unbinned) with given name. A null pointer is returned if not found. 
  plc.SetModel(model);
  plc.SetParameters(poi); // specify the parameters of interest in the interval 
  plc.SetNullParameters(*nullParams); // set parameter values for the null if using a common PDF 
  cout << " I AM HERE " << endl << endl << endl;

  // We get a HypoTestResult out of the calculator, and we can query it.
  HypoTestResult* htr = plc.GetHypoTest();
  cout << " I AM HERE 2" << endl << endl << endl;
  cout << "-------------------------------------------------" << endl;
  cout << "The p-value for the null is " << htr->NullPValue() << endl;
  cout << "Corresponding to a signifcance of " << htr->Significance() << endl;
  cout << "-------------------------------------------------\n\n" << endl;

}

