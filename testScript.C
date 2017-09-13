

using namespace RooFit;
using namespace RooStats;

void testScript(){

//Open the rootfile and get the workspace from the exercise_0
  TFile fIn("fitNoCurrent.root");
  fIn.cd();
  RooWorkspace *w = (RooWorkspace*)fIn.Get("w");
  w->Print();

  RooRealVar *cuKa1N = w->var("cuKa1N");
  RooRealVar *Nsig = w->var("Nsig");
  RooRealVar *energy = w->var("energy");

  RooAbsData *noAD = w->data("noDH");
  RooDataHist noDH = RooDataHist("noDH","noDH",*energy,*noAD);

  RooAbsData *withAD = w->data("withDH");
  RooDataHist withDH = RooDataHist("withDH","withDH",*energy,*withAD);

  Nsig->Print();
  cuKa1N->Print();
   
  RooPlot *frame = w->var("energy")->frame() ; // roorealvar->frame(): Create a new RooPlot on the heap with a drawing frame initialized for this object, but no plot contents. 

  
  noDH.plotOn(frame); // has to be plotted first so that the pdf gets automatically normalized to amount of counts
  withDH.plotOn(frame,MarkerColor(kRed));
  w->pdf("PDFtot")->plotOn(frame) ; 
  

 //Actually plot the result
  TCanvas c1;
  c1.cd();
  frame->Draw();
  c1.SaveAs("testScript.gif");

  //RooAbsData *withC = w->data("withDH");

  //withC->Draw();

}
