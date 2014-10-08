#include "alphaTree.C"
#include <iostream>
void makeQcdRate2D(){

  std::vector<TString> qcdFiles;
  qcdFiles.push_back("/home/matthew/trees/QCD3050_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD5080_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD80120_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD120170_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD170300_2014-10-07.root");

  std::map<TString,TH2D*> rateTotal;

  TFile * testOut = (TFile*) new TFile("qcdRate2D.root","recreate");

  for (int i = 0; i < qcdFiles.size(); i++)
  {

    std::cout << "File: " << qcdFiles.at(i) << std::endl; 

    TFile * fIn = (TFile *) new TFile(qcdFiles.at(i));
    TTree * tree = (TTree *) fIn->Get("MakeL1Trees/Ntuple");
    testOut->cd();
    alphaTree t(tree);
    if (i ==0)
    {
      rateTotal  = t.plotRate2dQcd("");
    }
    else 
    {
      std::map<TString,TH2D*> rates  = t.plotRate2dQcd(qcdFiles.at(i));
      for (std::map<TString,TH2D*>::iterator rate = rates.begin(); rate != rates.end();rate++)
      {
	TH2D * temp = (TH2D*) (rate->second)->Clone();
	rateTotal[rate->first]->Add(temp);
	delete temp;
      }
    }
  }
  for (std::map<TString,TH2D*>::iterator rate = rateTotal.begin(); rate != rateTotal.end();rate++)
  {
    std::cout << rate->first << std::endl;
    rate->second->Write();
  }
}

void makeQcdRate1D(){

  std::vector<TString> qcdFiles;
  qcdFiles.push_back("/home/matthew/trees/QCD3050_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD5080_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD80120_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD120170_2014-10-07.root");
  qcdFiles.push_back("/home/matthew/trees/QCD170300_2014-10-07.root");

  std::map<TString,TH1D*> rateTotal;

  TFile * testOut = (TFile*) new TFile("qcdRate1D.root","recreate");
  for (int i = 0; i < qcdFiles.size(); i++)
  {

    std::cout << "File: " << qcdFiles.at(i) << std::endl; 

    TFile * fIn = (TFile *) new TFile(qcdFiles.at(i));
    TTree * tree = (TTree *) fIn->Get("MakeL1Trees/Ntuple");
    testOut->cd();
    alphaTree t(tree);
    if (i == 0)
    {
      rateTotal  = t.plotRate1dQcd("");
    }
    else
    {
      std::map<TString,TH1D*> rates  = t.plotRate1dQcd(qcdFiles.at(i));
      for (std::map<TString,TH1D*>::iterator rate = rates.begin(); rate != rates.end();rate++)
      {
	TH1D * temp = (TH1D*) (rate->second)->Clone();
	rateTotal[rate->first]->Add(temp);
	delete temp;
      }
    }
  }
  for (std::map<TString,TH1D*>::iterator rate = rateTotal.begin(); rate != rateTotal.end();rate++)
  {
    std::cout << rate->first << std::endl;
    rate->second->Write();
  }

}
