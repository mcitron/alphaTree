#define alphaTree_cxx
#include "alphaTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

int N_50NS_BUNCHES=1368;
int N_25NS_BUNCHES=2508;
int LHC_FREQUENCY=11246;
TH2D * makeCumu(TH2D * input, double norm);
TH1D * makeCumu1D(TH1D * input, double norm);
TH1 * makeCumuInv1D(TH1 * input, double norm);
TH2 * makeCumuInv2D(TH2 * input, double norm);
int ZB_XSECTION=LHC_FREQUENCY*N_25NS_BUNCHES;
void alphaTree::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L alphaTree.C
  //      Root > alphaTree t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //

  //     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch


}
void alphaTree::findRate()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::map<TString,std::pair<int,bool> > cuts;
  cuts["Original"] = std::make_pair(0,false);
  cuts["Update"] =std::make_pair(0,false);
  cuts["dPhiNew"] = std::make_pair(0,false);
  cuts["met70"] = std::make_pair(0,false);
  cuts["mht60Ht100Full"] = std::make_pair(0,false);
  cuts["met80"] = std::make_pair(0,false);
  cuts["single"] = std::make_pair(0,false);
  cuts["mhtDivHtBin1"] = std::make_pair(0,false);
  cuts["mhtDivHtBin2"] =std::make_pair(0,false);
  cuts["mhtDivHtBin3"] =std::make_pair(0,false);
  cuts["mhtDivHtBin4"] =std::make_pair(0,false);
  cuts["mhtDivHtBin5"] =std::make_pair(0,false);
  cuts["mhtDivHtBin6"] =std::make_pair(0,false);
  cuts["mhtDivHtMet"] =std::make_pair(0,false);
  cuts["dPhiFull"] =std::make_pair(0,false);
  cuts["dPhiNewFull"] =std::make_pair(0,false);
  cuts["met70Full"] =std::make_pair(0,false);
  cuts["met60CrossDJ60MetFull"] =std::make_pair(0,false);
  cuts["met80Full"] =std::make_pair(0,false);
  cuts["mhtDivHtBin1Full"] = std::make_pair(0,false);
  cuts["mhtDivHtBin2Full"] = std::make_pair(0,false);
  cuts["mhtDivHtBin3Full"] = std::make_pair(0,false);
  cuts["mhtDivHtBin4Full"] =std::make_pair(0,false);
  cuts["mhtDivHtBin5Full"] =std::make_pair(0,false);
  cuts["mhtDivHtBin6Full"] =std::make_pair(0,false);
  cuts["mhtDivHtMetFull"] =std::make_pair(0,false);


  cuts["metMht60Ht100Full"]= std::make_pair(0,false);
  cuts["metMhtDivHtFull"]= std::make_pair(0,false);




  nentries = 10000000;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullmenuOrig = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    cuts["Original"].second = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || (jetPtsUct->size()>1 && jetPtsUct->at(1) >= 60 && metPtUct >= 60) || metPtUct > 70;
    cuts["Update"].second = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
    cuts["dPhi"].second = (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 30 && htUct  >= 125. && abs(jetDPhi12Uct/0.3409) <= 7);
    cuts["dPhiNew"].second =(jetPtsUct->size() > 1 && jetPtsUct->at(0) >120 && jetPtsUct->at(1) >= 30 && htUct  >= 125. && abs(jetDPhi12Uct/0.3409) <= 7);
    cuts["met70"].second = (metPtUct >= 70.);
    cuts["met80"].second = (metPtUct >= 80.);
    cuts["met60CrossDJ60MetFull"].second = fullmenuOrig || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 60 && metPtUct >= 60) || metPtUct >= 70;
    cuts["single"].second = (jetPtsUct->size() > 0 && jetPtsUct->at(0) >= 200.);
    cuts["mhtDivHtBin1"].second = (mhtDivHtUct >= 0.30 && htUct >= 125.);
    cuts["mhtDivHtBin2"].second = (mhtDivHtUct >= 0.3 && htUct >= 110.);
    cuts["mhtDivHtBin3"].second = (mhtDivHtUct >= 0.3 && htUct >= 125.) || (mhtDivHtUct >= 0.7 && htUct >= 70.);
    cuts["mhtDivHtBin4"].second = (mhtDivHtUct >= 0.3 && htUct >= 85.);
    cuts["mhtDivHtBin4"].second = (mhtDivHtUct >= 0.2 && htUct >= 120. ) || (mhtDivHtUct >= 0.3 && htUct >= 125.);
    cuts["mhtDivHtBin5"].second = (mhtDivHtUct >= 0.1 && htUct >= 150. ) || (mhtDivHtUct >= 0.3 && htUct >= 125.);
    cuts["mhtDivHtBin6"].second = (mhtDivHtUct >= 0.18 && htUct >= 125. );
    cuts["mhtDivHtMet"].second = (metPtUct >=70.) || (mhtDivHtUct >= 0.3 && htUct >= 125.);

    cuts["dPhiNewFull"].second = (jetPtsUct->size() > 1 && jetPtsUct->at(0) > 120 && jetPtsUct->at(1) >= 30 && htUct  >= 125. && abs(jetDPhi12Uct/0.3409) <= 7 ) || fullmenu;
    cuts["met70Full"].second = (metPtUct >= 70.) || fullmenuOrig;
    cuts["mht60Ht100Full"].second = (mhtPtUct >= 60. && htUct >= 100) || fullmenu;
    cuts["met80Full"].second = (metPtUct >= 80.) || fullmenu;
    cuts["mhtDivHtBin1Full"].second = (mhtDivHtUct >= 0.3 && htUct >= 125.) || fullmenu;
    cuts["mhtDivHtBin2Full"].second = (mhtDivHtUct >= 0.3 && htUct >= 110.) || fullmenu;
    cuts["mhtDivHtBin3Full"].second = (mhtDivHtUct >= 0.3 && htUct >= 125.) || (mhtDivHtUct >= 0.7 && htUct >= 70.) || fullmenu;
    cuts["mhtDivHtBin4Full"].second = (mhtDivHtUct >= 0.2 && htUct >= 120. ) || (mhtDivHtUct >= 0.3 && htUct >= 125.) || fullmenu;
    cuts["mhtDivHtBin5Full"].second = (metPtUct >=70.)|| (mhtDivHtUct >= 0.1 && htUct >= 150. ) || (mhtDivHtUct >= 0.3 && htUct >= 125.) || fullmenu;
    cuts["mhtDivHtMetFull"].second = (metPtUct >=70.) || (mhtDivHtUct >= 0.3 && htUct >= 125.) || fullmenu;
    cuts["mhtDivHtBin6Full"].second = (mhtDivHtUct >= 0.18 && htUct >= 125. ) || fullmenu;

    cuts["metMht60Ht100Full"].second = (mhtPtUct >= 60. && htUct >= 100) || fullmenu || metPtUct >= 70;;
    cuts["metMhtDivHtFull"].second = (mhtDivHtUct >= 0.30 && htUct >= 125.) || fullmenu || metPtUct >= 70;
    cuts["metOriginal"].second = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct;

    for (std::map<TString,std::pair<int,bool> >::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
    {
      //  std::cout << iCut->second.first << std::endl;
      if (((iCut)->second).second) (iCut->second).first++;
    }


    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,std::pair<int,bool> >::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
  {
    std::cout << iCut->first << ": " << ZB_XSECTION*(double)(iCut->second).first/(nentries)<< std::endl;
  }
}

void alphaTree::findEff()
{
  if (fChain == 0) return;

  int jetPtMin = 50;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  std::map<TString,bool> mult;
  mult["2Jet"] =  false;
  mult["3Jet"] = false;
  mult["ge4Jet"] = false;
  std::map<TString,std::pair<int,bool> > analysisBins;

  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
  {
    analysisBins["Bin1"+iMult->first] =std::make_pair(0,false);
    analysisBins["Bin2"+iMult->first] =std::make_pair(0,false);
    analysisBins["Bin3"+iMult->first] =std::make_pair(0,false);
  }
  analysisBins["NoSelection"] =std::make_pair(0,false);

  std::map<TString,std::pair<int,bool> > cuts;
  for (std::map<TString,std::pair<int,bool> >::iterator iBin = analysisBins.begin(); iBin != analysisBins.end();iBin++)
  {
    cuts["original"+iBin->first] =std::make_pair(0,false);
    cuts["ht"+iBin->first] =std::make_pair(0,false);
    cuts["mhtDivHt"+iBin->first] = std::make_pair(0,false);
    cuts["mht"+iBin->first] = std::make_pair(0,false);
    cuts["met60Full"+iBin->first] =std::make_pair(0,false);
    cuts["met70Full"+iBin->first] =std::make_pair(0,false);
    cuts["metMhtDivHt"+iBin->first] = std::make_pair(0,false);
    cuts["metMht"+iBin->first] = std::make_pair(0,false);
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    mult["2Jet"] = multiplicity50Gen5 ==2;
    mult["3Jet"] = multiplicity50Gen5 ==3;
    mult["ge4Jet"] = multiplicity50Gen5 > 3;

    for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
    {
      analysisBins["Bin1"+iMult->first].second = (htGen5 >= 200. && htGen5 < 300. && alphaTGen5 > 0.65) && iMult->second; 
      analysisBins["Bin2"+iMult->first].second = (htGen5 >= 300. && htGen5 < 400. && alphaTGen5 > 0.60) && iMult->second; 
      analysisBins["Bin3"+iMult->first].second = (htGen5 >= 400. && htGen5 < 500. && alphaTGen5 > 0.55) && iMult->second;
    }
    mult["1Jet"] = multiplicity50Gen5==1;
    analysisBins["NoSelection"].second = true;

    analysisBins["Bin11Jet"].second = (htGen5 >= 200. && htGen5 < 300.) && multiplicity50Gen5 ==1;
    analysisBins["Bin21Jet"].second = (htGen5 >= 300. && htGen5 < 400.) && multiplicity50Gen5 ==1;
    analysisBins["Bin31Jet"].second = (htGen5 >= 400. && htGen5 < 500.) && multiplicity50Gen5 ==1; 

    bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullmenuMET70 = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct > 70; 
    bool fullmenuOrig = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 

    for (std::map<TString,std::pair<int,bool> >::iterator iBin = analysisBins.begin(); iBin != analysisBins.end();iBin++)
    {

      cuts["mhtDivHt"+iBin->first].second = ((mhtDivHtUct >= 0.4 && htUct >= 112.) || fullmenu) && (iBin->second).second;
      cuts["mht"+iBin->first].second =((mhtPtUct >= 56 && htUct >= 112) || fullmenu) && (iBin->second).second;
      cuts["original"+iBin->first].second =(fullmenuOrig) && (iBin->second).second;
      cuts["ht"+iBin->first].second =(htUct >= 175) && (iBin->second).second;
      cuts["met60Full"+iBin->first].second = (fullmenu || metPtUct >= 60) && (iBin->second).second;
      cuts["met70Full"+iBin->first].second = (fullmenuOrig || metPtUct >= 70) && (iBin->second).second;
      cuts["metMhtDivHt"+iBin->first].second = ((mhtDivHtUct >= 0.3 && htUct >= 125.) || fullmenuMET70) && (iBin->second).second;
      cuts["metMht"+iBin->first].second = ((mhtPtUct >= 44 && htUct >= 125.) || fullmenuMET70) && (iBin->second).second;

      if((iBin->second).second) {iBin->second.first++;}
    }

    for (std::map<TString,std::pair<int,bool> >::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
    {
      if (((iCut)->second).second){ (iCut->second).first++;}
    }

    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }

  for (std::map<TString,std::pair<int,bool> >::iterator iBin = analysisBins.begin(); iBin != analysisBins.end();iBin++)
  {
    std::cout << "original"+iBin->first+": " << (double) cuts["original"+iBin->first].first/iBin->second.first << "\t"<<cuts["original"+iBin->first].first << "\t" << iBin->second.first << std::endl;
    std::cout << "ht"+iBin->first+": " << (double) cuts["ht"+iBin->first].first/iBin->second.first << "\t"<<cuts["ht"+iBin->first].first << "\t" << iBin->second.first << std::endl;
    std::cout << "met70Full"+iBin->first+": " << (double) cuts["met70Full"+iBin->first].first/iBin->second.first <<"\t"<<cuts["met70Full"+iBin->first].first << "\t" << iBin->second.first <<  std::endl;
    std::cout << "met60Full"+iBin->first+": " << (double) cuts["met60Full"+iBin->first].first/iBin->second.first <<"\t"<<cuts["met60Full"+iBin->first].first << "\t" << iBin->second.first <<  std::endl;
    std::cout << "mhtDivHt"+iBin->first+": " << (double) cuts["mhtDivHt"+iBin->first].first/iBin->second.first << "\t"<<cuts["mhtDivHt"+iBin->first].first << "\t" << iBin->second.first << std::endl;
    std::cout << "mht"+iBin->first+": " << (double) cuts["mht"+iBin->first].first/iBin->second.first << "\t"<<cuts["mht"+iBin->first].first << "\t" << iBin->second.first << std::endl;
    std::cout << "metMhtDivHt"+iBin->first+": " << (double) cuts["metMhtDivHt"+iBin->first].first/iBin->second.first << "\t"<<cuts["metMhtDivHt"+iBin->first].first << "\t" << iBin->second.first << std::endl;
    std::cout << "metMht"+iBin->first+": " << (double) cuts["metMht"+iBin->first].first/iBin->second.first << "\t"<<cuts["metMht"+iBin->first].first << "\t" << iBin->second.first << std::endl;

    std::cout << std::endl;
  }
}
void alphaTree::plotOfflineEff()
{
  TFile* fOut = new TFile("offlineEff.root","recreate");
  std::map<TString,bool> cuts;


  cuts["ht160R"] = false;
  cuts["ht150R"] = false;
  cuts["ht140R"] = false;
  cuts["ht130R"] = false;
  cuts["ht120R"] = false;

  cuts["ht160RMet"] = false;
  cuts["ht150RMet"] = false;
  cuts["ht140RMet"] = false;
  cuts["ht130RMet"] = false;
  cuts["ht120RMet"] = false;

  std::map<TString,bool> mult;
  mult["1Jet"] = false;
  mult["2Jet"] =  false;
  mult["3Jet"] = false;
  mult["ge4Jet"] = false;

  int jetPtMin = 50;
  std::map<TString,bool> hardlead;
  hardlead["Hard"] = false;
  hardlead["noHard"] = false;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::map <TString,TEfficiency *> offlineEff;
  std::map <TString,TDirectory *> dir;

  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
  {
    //dir[iMult->first] = fOut->mkdir(iMult->first);
    for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
    {
      for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
      {
	offlineEff[iHard->first+iCut->first+iMult->first] = new TEfficiency(iCut->first+iMult->first+iHard->first+"Bin1",";HT/GeV;alphaT/GeV",10,0,500,10,0,1);  
	//l1Turn[iCut->first+iMult->first]->SetDirectory(dir);
      }
    } 
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    bool veto = false;
    for (unsigned int i = 0; i < jetPtsAllUct->size(); i++)
    {
      if (abs(jetEtasAllUct->at(i)) > 3) {veto = true; break;}
    }
    if (ientry < 0) break;

    if (jetPtsGen5->size() > 1)
    {
      hardlead["noHard"] = jetPtsGen5->at(1) > jetPtMin;
      hardlead["Hard"] = jetPtsGen5->at(1) >= 100;
    }
    else
    {
      hardlead["noHard"] = false;
      hardlead["Hard"] = false;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    bool offlineCut = true;
    /*  
	for (int j = 0; j < jetEtasAllGen5->size(); j++)
	{
	if (jetPtsAllGen5->at(j) >= 50 && fabs(jetEtasAllGen5->at(j)) > 3.){jetVeto = false; break;}
	}
	*/
    ///////FULLMENUS WITHOUT JET TRIGGERS
    bool fullmenu = htUct >=200.;// || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullmenuMet = htUct >=200. || metPtUct >= 70; //|| (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 

    cuts["ht160R"] = fullmenu || htUct >= 160;
    cuts["ht150R"] = fullmenu || (htUct >= 150 && mhtDivHtUct >= 0.17); 
    cuts["ht140R"] = fullmenu || (htUct >= 140 && mhtDivHtUct >= 0.25);
    cuts["ht130R"] = fullmenu || (htUct >= 130 && mhtDivHtUct >= 0.34);
    cuts["ht120R"] = fullmenu || (htUct >= 120 && mhtDivHtUct >= 0.40);

    cuts["ht160RMet"] = fullmenuMet || (htUct >= 160);
    cuts["ht150RMet"] = fullmenuMet || (htUct >= 150 && mhtDivHtUct >= 0.17 );
    cuts["ht140RMet"] = fullmenuMet || (htUct >= 140 && mhtDivHtUct >= 0.25 );
    cuts["ht130RMet"] = fullmenuMet || (htUct >= 130 && mhtDivHtUct >= 0.32 );
    cuts["ht120RMet"] = fullmenuMet || (htUct >= 120 && mhtDivHtUct >= 0.36 );

    for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
    {
      mult["1Jet"] = multiplicity50Gen5 ==1;
      mult["2Jet"] = multiplicity50Gen5 ==2 && offlineCut&& iHard->second;
      mult["3Jet"] = multiplicity50Gen5 ==3 &&offlineCut&& iHard->second;
      mult["ge4Jet"] = multiplicity50Gen5 >= 4 && offlineCut&& iHard->second;

      for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
      {
	if(iMult->second)
	{
	  for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
	  {
	    offlineEff[iHard->first+iCut->first+iMult->first]->Fill(iCut->second,htGen5,alphaTGen5);
	  } 
	}
      }
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
  {
    for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
    {
      dir[iHard->first+iMult->first] = fOut->mkdir(iHard->first+iMult->first);
      dir[iHard->first+iMult->first]->cd();
      for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
      {
	TH2 * passed =  (TH2*) offlineEff[iHard->first+iCut->first+iMult->first]->GetCopyPassedHisto();
	TH2 * total = (TH2*) offlineEff[iHard->first+iCut->first+iMult->first]->GetCopyTotalHisto();

	passed = makeCumuInv2D(passed,1.0);
	total = makeCumuInv2D(total,1.0);

	TH2 * eff = (TH2*)passed->Clone();
	eff->Divide(total);
	eff->SetName(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	eff->SetTitle(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	eff->Write();
	delete passed;
	delete total;
      }
    } 
  }
}


void alphaTree::plotTurnOn()
{
  TFile* fOut = new TFile("turnonPlots.root","recreate");
  std::map<TString,bool> cuts;

  cuts["Original"] = false;
  cuts["Update"] = false;
  cuts["dPhi"] = false;
  cuts["dPhiNew"] = false;

  cuts["met70"] = false;
  cuts["met62"] = false;
  cuts["met60"] = false;
  cuts["ht175"] = false;

  cuts["ht200"] = false;
  cuts["met70Full"] = false;
  cuts["met62Full"] = false;
  cuts["met60Full"] = false;

  cuts["mhtDivHt70"] = false;
  cuts["mhtDivHt70Full"] =false; 
  cuts["metMhtDivHt70"] = false;
  cuts["metMhtDivHt70Full"] =false;

  cuts["metMhtDivHt130"] = false;
  cuts["mhtDivHt130"] = false;
  cuts["mhtDivHt130OR70"] =false;
  cuts["mhtDivHt130Full"] =false;
  cuts["mht60Full"] =false;
  cuts["metMht60Full"] =false;
  cuts["mhtDivHt130OR70Full"] = false;

  cuts["metMhtDivHt130Full"] =false;
  cuts["metMhtDivHt130O70RFull"] = false;

  cuts["dPhiFull"] = false;
  cuts["dPhiNewFull"] = false;

  std::map<TString,bool> mult;
  mult["1Jet"] = false;
  mult["2Jet"] =  false;
  mult["3Jet"] = false;
  mult["ge4Jet"] = false;

  int jetPtMin = 50;
  std::map<TString,bool> hardlead;
  hardlead["Hard"] = false;
  hardlead["noHard"] = false;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  std::map <TString,TEfficiency *> l1TurnBin1;
  std::map <TString,TEfficiency *> l1TurnBin2;
  std::map <TString,TDirectory *> dir;

  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
  {
    //dir[iMult->first] = fOut->mkdir(iMult->first);
    for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
    {
      for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
      {
	l1TurnBin1[iHard->first+iCut->first+iMult->first] = new TEfficiency(iCut->first+iMult->first+iHard->first+"Bin1",";MHT/GeV;",40,0,1000);  
	l1TurnBin2[iHard->first+iCut->first+iMult->first] = new TEfficiency(iCut->first+iMult->first+iHard->first+"Bin2",";MHT/GeV;",40,0,1000);  
	//l1Turn[iCut->first+iMult->first]->SetDirectory(dir);
      }
    } 
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    bool veto = false;
    for (unsigned int i = 0; i < jetPtsAllUct->size(); i++)
    {
      if (abs(jetEtasAllUct->at(i)) > 3) {veto = true; break;}
    }
    if (ientry < 0) break;

    if (jetPtsGen5->size() > 1)
    {
      hardlead["noHard"] = jetPtsGen5->at(1) > jetPtMin;
      hardlead["Hard"] = jetPtsGen5->at(1) >= 100;
    }
    else
    {
      hardlead["noHard"] = false;
      hardlead["Hard"] = false;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    bool offlineCut = true;
    /*  
	for (int j = 0; j < jetEtasAllGen5->size(); j++)
	{
	if (jetPtsAllGen5->at(j) >= 50 && fabs(jetEtasAllGen5->at(j)) > 3.){jetVeto = false; break;}
	}
	*/
    ///////FULLMENUS WITHOUT JET TRIGGERS
    bool fullmenu = htUct >=200.;// || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullmenuOrig = htUct >=175.; //|| (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    cuts["Original"] = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
    cuts["Update"] = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
    cuts["dPhi"] = (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 30 && htUct  >= 125. && abs(jetDPhi12Uct/0.3409) <= 7);
    cuts["dPhiNew"] = (jetPtsUct->size() > 1 && jetPtsUct->at(0) > 120 && jetPtsUct->at(1)>30 && htUct > 125);

    cuts["met70"] = (metPtUct >= 70.);
    cuts["met62"] = (metPtUct >= 62.);
    cuts["met60"] = (metPtUct >= 60.);

    cuts["ht175"] = (htUct >= 175.);
    cuts["ht200"] = (htUct >= 200.);
    cuts["met70Full"] = (metPtUct > 70) || fullmenuOrig;
    cuts["met62Full"] = (metPtUct > 62) || fullmenu;
    cuts["met60Full"] = (metPtUct > 60) || fullmenu;
    cuts["mht60Full"] = (mhtPtUct > 60 && htUct > 100) || fullmenu;

    cuts["mhtDivHt70"] = (mhtDivHtUct >= 0.70 && htUct >= 70.);
    cuts["mhtDivHt70Full"] = (mhtDivHtUct >= 0.70 && htUct >= 70.) || fullmenu;
    cuts["metMhtDivHt70"] = (metPtUct >= 70.) || (mhtDivHtUct >= 0.70 && htUct >= 70.);
    cuts["metMhtDivHt70Full"] = (metPtUct >= 70.) || (mhtDivHtUct >= 0.70 && htUct >= 70.) || fullmenu;

    cuts["mhtDivHt130"] = (mhtDivHtUct >= 0.30 && htUct >= 125.);
    cuts["metMhtDivHt130"] = (mhtDivHtUct >= 0.30 && htUct >= 125.) || (metPtUct >= 70);
    cuts["mhtDivHt130OR70"] = (mhtDivHtUct >= 0.30 && htUct >= 125.) || (mhtDivHtUct >= 0.70 && htUct >= 70.);
    cuts["mhtDivHt130Full"] = (mhtDivHtUct >= 0.30 && htUct >= 125.) || fullmenu;
    cuts["mhtDivHt130OR70Full"] = (mhtDivHtUct >= 0.30 && htUct >= 125.) || (mhtDivHtUct >= 0.70 && htUct >= 70.) || fullmenu; 

    cuts["metMhtDivHt130Full"] = ((metPtUct >= 70.) || (mhtDivHtUct >= 0.30 && htUct >= 125.)) || fullmenu;
    cuts["metMhtDivHt130O70RFull"] = ((metPtUct >= 70.) || (mhtDivHtUct >= 0.30 && htUct >= 125.) || (mhtDivHtUct >= 0.70 && htUct >= 70.)) || fullmenu;
    cuts["metMht60Full"] = (mhtPtUct > 60 && htUct > 100) || fullmenu || metPtUct > 70;

    cuts["dPhiFull"] = (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 30 && htUct  >= 125. && abs(jetDPhi12Uct/0.3409) <= 7 ) || fullmenu;
    cuts["dPhiNewFull"] = (jetPtsUct->size() > 1 && jetPtsUct->at(0) > 120 && jetPtsUct->at(1) > 30 && htUct > 125) || fullmenu;

    for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
    {
      mult["1Jet"] = multiplicity50Gen5 ==1;
      mult["2Jet"] = multiplicity50Gen5 ==2 && offlineCut&& iHard->second;
      mult["3Jet"] = multiplicity50Gen5 ==3 &&offlineCut&& iHard->second;
      mult["ge4Jet"] = multiplicity50Gen5 >= 4 && offlineCut&& iHard->second;

      for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
      {
	if(iMult->second)
	{
	  for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
	  {
	    if(htGen5 >= 200 && htGen5 <= 300)
	    {
	      l1TurnBin1[iHard->first+iCut->first+iMult->first]->Fill(iCut->second,mhtPtGen5);
	    }
	    if(htGen5 > 300 && htGen5 <= 400)
	    {
	      l1TurnBin2[iHard->first+iCut->first+iMult->first]->Fill(iCut->second,mhtPtGen5);
	    }
	  } 
	}
      }
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
  {
    for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
    {
      dir[iHard->first+iMult->first] = fOut->mkdir(iHard->first+iMult->first+"Bin1");
      dir[iHard->first+iMult->first]->cd();
      for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
      {
	//l1Turn[iHard->first+iCut->first+iMult->first]->Write();
	TGraphAsymmErrors * temp = l1TurnBin1[iHard->first+iCut->first+iMult->first]->CreateGraph(); 
	temp->SetName(iHard->first+iCut->first+iMult->first+ "_Eff");
	temp->Write();

	TH1 * passed = l1TurnBin1[iHard->first+iCut->first+iMult->first]->GetCopyPassedHisto();
	TH1 * total = l1TurnBin1[iHard->first+iCut->first+iMult->first]->GetCopyTotalHisto();
	passed = makeCumuInv1D(passed,1.0);
	total = makeCumuInv1D(total,1.0);
	TGraphAsymmErrors * effCumulativeGraph = new TGraphAsymmErrors(passed,total); 
	effCumulativeGraph->SetName(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	effCumulativeGraph->SetTitle(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	effCumulativeGraph->GetXaxis()->SetTitle("MHT/GeV");
	effCumulativeGraph->Write();
	delete passed;
	delete total;
      }

      dir[iHard->first+iMult->first] = fOut->mkdir(iHard->first+iMult->first+"Bin2");
      dir[iHard->first+iMult->first]->cd();
      for (std::map<TString,bool>::iterator iCut = cuts.begin(); iCut != cuts.end();iCut++)
      {
	TGraphAsymmErrors * temp2 = l1TurnBin2[iHard->first+iCut->first+iMult->first]->CreateGraph(); 
	temp2->SetName(iHard->first+iCut->first+iMult->first+ "_Eff");
	temp2->Write();

	TH1 * passed2 = l1TurnBin2[iHard->first+iCut->first+iMult->first]->GetCopyPassedHisto();
	TH1 * total2 = l1TurnBin2[iHard->first+iCut->first+iMult->first]->GetCopyTotalHisto();
	passed2 = makeCumuInv1D(passed2,1.0);
	total2 = makeCumuInv1D(total2,1.0);
	TGraphAsymmErrors * effCumulativeGraph2 = new TGraphAsymmErrors(passed2,total2); 
	effCumulativeGraph2->SetName(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	effCumulativeGraph2->SetTitle(iHard->first+iCut->first+iMult->first+ "_CumuEff");
	effCumulativeGraph2->GetXaxis()->SetTitle("MHT/GeV");
	effCumulativeGraph2->Write();
	delete passed2;
	delete total2;
      }

    } 
  }
}

void alphaTree::regionStudy()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  TH2D * test = new TH2D("valid",";my H_{T};true H_{T} Cut (GeV)",400,0.,1000.,400,0.,1000.);
  std::cout << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;
    double sum = 0.;


    for (int i =0; i< regionEt->size();i++ )
    {
      sum += int(regionEt->at(i));
    }
    test->Fill(int(sum),htUct);
  }
  test->Draw("colz");
}

void alphaTree::dhtStudy()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  TH2D * test2D = new TH2D("valid",";dH_{T}/H_{T};mH_{T}/H_{T} Cut",400,0.,1.,400,0.,1.);
  //TH1D * test = new TH1D("valid",";my H_{T};",400,0.,1000.);
  //TH1D * test2 = new TH1D("gen",";my H_{T};",400,0.,1000.);
  std::cout << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    double ht=0.;
    double mhtx=0.;
    double mhty=0.;
    double mht=0.;
    double dht=0.;
    double dpt=0.;
    double minpt=0.;
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;
    if (ientry < 0) break;

    if (jetPtsUct->size() > 0 && jetPtsUct->at(0) > 0) 
    {
      dht = 2*jetPtsUct->at(0);
      dpt = jetPtsUct->at(0);
      for (int i =0; i< jetPtsUct->size();i++ )
      {
	if(jetPtsUct->at(i) > 30)
	{
	  ht+=jetPtsUct->at(i);
	  mhtx+=jetPtsUct->at(i)*cos(jetPhisUct->at(i));
	  mhty+=jetPtsUct->at(i)*sin(jetPhisUct->at(i));
	}
	if(i < 3) dht -= jetPtsUct->at(i);
	if (i != 0 && jetPtsUct->at(i) > 1.) {minpt=jetPtsUct->at(i);}
      }

      mht = sqrt(mhtx*mhtx+mhty*mhty);
      dpt -= minpt;
      //if(jetPtsUct->size() <= 2 || (jetPtsUct->size() > 2 && jetPtsUct->at(2) < 30)) 
      test2D->Fill(fabs(dpt/jetPtsUct->at(0)),mht/ht);
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  test2D->Draw("colz");
  //test->SetLineColor(2);
  //  //test2->Draw("same");
  //  
}
void alphaTree::findCalo()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int counter = 0;
  nentries = 1E6;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;
    if (ientry < 0) break;

    double ht40 = 0.;
    for (int i = 0; i < jetPtsCalo->size();i++) {if (jetPtsCalo->at(i) >= 40) {ht40 += jetPtsCalo->at(i);};} 
    if (htUct > 200 && ht40 > 200) counter++;
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  std::cout << std::endl;
  std::cout << (double) ZB_XSECTION * counter / nentries << std::endl;
}
void alphaTree::onOffHt()
{
  TFile* fOut = new TFile("eff2D.root","recreate");
  fOut->cd();
  //TH2D * effNum = new TH2D("effNum",";Offline;Online",40,0.,1000.,40,0.,1000.);
  //TH2D * effDenom = new TH2D("effDenom",";Offline;Online",40,0.,1000.,40,0.,1000.);

  std::map<TString,bool> mult;
  mult["1Jet"] = false;
  mult["2Jet"] =  false;
  mult["3Jet"] = false;
  mult["ge4Jet"] = false;
  std::map<TString,TH2D*> num;
  std::map<TString,TH2D*> numCumu;
  std::map<TString,TH2D*> denom;
  std::map<TString,TH2D*> denomCumu;
  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end(); iMult++)
  {
    num[iMult->first] = new TH2D("effNum"+iMult->first,";Offline H_{T};Online H_{T} Cut (GeV)",100,0.,1000.,100,0.,1000.);
    denom[iMult->first] = new TH2D("effDenom"+iMult->first,";Offline H_{T};Online H_{T} Cut  (GeV)",100,0.,1000.,100,0.,1000.);
    numCumu[iMult->first] = new TH2D("effNumCumu"+iMult->first,";Offline H_{T};Online H_{T} Cut (GeV)",100,0.,1000.,100,0.,1000.);
    denomCumu[iMult->first] = new TH2D("effDenomCumu"+iMult->first,";Offline H_{T}; Online H_{T} Cut (GeV)",100,0.,1000.,100,0.,1000.);
  }

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  /*std::map<int,TEfficiency *> result;
    for (int i =0; i <= 1000;i+=100)
    {
    result[i] = new TEfficiency("temp",";pT;",10,0.,1000.);
    } */ 
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    bool offlineCut = alphaTGen5 > 0.65;
    ///Mult cuts
    mult["1Jet"] = jetPtsGen5->at(1) < 50; //&& mhtPtGen5 > 150;
    mult["2Jet"] = (jetPtsGen5->at(1) > 50 && jetPtsGen5->at(2) < 50);//&&offlineCut;
    mult["3Jet"] = (jetPtsGen5->at(2) >= 50 && jetPtsGen5->at(3) < 50) ;//&&offlineCut;
    mult["ge4Jet"] = (jetPtsGen5->at(3) >= 50);// && offlineCut;

    for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
    {
      int nBinsY = num[iMult->first]->GetNbinsY();
      int nBinsX = num[iMult->first]->GetNbinsX();
      if(iMult->second)
      {
	for (int binY = 0; binY <= nBinsY+1; binY++)
	{
	  denom[iMult->first]->Fill(htGen5,denom[iMult->first]->GetYaxis()->GetBinCenter(binY));
	  if(htUct >= (num[iMult->first]->GetYaxis()->GetBinLowEdge(binY)))
	  {
	    //std::cout << effDenom->GetXaxis()->GetBinCenter(binX) << " " << htGen5 << std::endl;
	    num[iMult->first]->Fill(htGen5,num[iMult->first]->GetYaxis()->GetBinCenter(binY));
	  } 
	  for (int binX = 0; binX <= nBinsX+1; binX++)
	  {
	    if(htGen5 >= (num[iMult->first]->GetXaxis()->GetBinLowEdge(binX)))
	    {
	      denomCumu[iMult->first]->Fill(denomCumu[iMult->first]->GetXaxis()->GetBinCenter(binX),denomCumu[iMult->first]->GetYaxis()->GetBinCenter(binY));
	      if(htUct >= (numCumu[iMult->first]->GetYaxis()->GetBinLowEdge(binY)))
	      {
		numCumu[iMult->first]->Fill(numCumu[iMult->first]->GetXaxis()->GetBinCenter(binX),numCumu[iMult->first]->GetYaxis()->GetBinCenter(binY));
	      } 

	    }
	  }
	}
      }
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end(); iMult++)
  {
    num[iMult->first]->Sumw2();
    denom[iMult->first]->Sumw2(); 
    num[iMult->first]->Divide(denom[iMult->first]);
    num[iMult->first]->GetXaxis()->SetRangeUser(0,500);
    num[iMult->first]->GetYaxis()->SetRangeUser(0,500);
    num[iMult->first]->Write();
    numCumu[iMult->first]->Sumw2();
    denomCumu[iMult->first]->Sumw2(); 
    numCumu[iMult->first]->Divide(denomCumu[iMult->first]);
    numCumu[iMult->first]->Write();
  }

}
void alphaTree::plotRateNint()
{
  TFile* fOut = new TFile("ratePlotsNint.root","recreate");
  fOut->cd();
  std::map<TString,bool > cuts;
  cuts["original"] =false;
  cuts["ht"] =false;
  cuts["mhtDivHt"] = false;
  cuts["mht"] = false;
  cuts["met60Full"] =false;
  cuts["met70Full"] =false;
  cuts["metMhtDivHt"] = false;
  cuts["metMht"] = false;
  cuts["alphaTCalo"] = false;
  cuts["MHTCalo"] = false;

  std::map<TString,TH1D *> num;
  TH1D * denom = new TH1D("total",";nInt;",20,20,80);
  for (std::map<TString,bool>::iterator iCut = cuts.begin();iCut != cuts.end(); iCut++)
  {
    num[iCut->first] = new TH1D(iCut->first,";nInt;",20,20,80);
  }

  Long64_t nentries = fChain->GetEntriesFast();
  nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;
  int test = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;
    //for (std::map<TString,TH1D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
    bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullmenuMET70 = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct > 70; 
    bool fullmenuOrig = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 

    cuts["mhtDivHt"] = ((mhtDivHtUct >= 0.4 && htUct >= 112.) ) ;
    cuts["mht"] =((mhtPtUct >= 56 && htUct >= 112) ) ;
    cuts["original"] =(fullmenuOrig) ;
    cuts["ht"] =(htUct >= 175) ;
    cuts["met60Full"] = (metPtUct >= 60) ;
    cuts["met70Full"] = (metPtUct >= 70) ;
    cuts["metMhtDivHt"] = ((mhtDivHtUct >= 0.3 && htUct >= 125.)) ;
    cuts["metMht"] = ((mhtPtUct >= 44 && htUct >= 125)) ;

    cuts["alphaTCalo"] = alphaTCalo >= 0.57 && htCalo >= 200;
    cuts["MHTCalo"] = mhtPtCalo >= 130 && htCalo >= 200;

    denom->Fill(nInt);
    for (std::map<TString,bool>::iterator iCut = cuts.begin();iCut != cuts.end(); iCut++)
    {
      if (iCut->second) num[iCut->first]->Fill(nInt,ZB_XSECTION);
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,bool>::iterator iCut = cuts.begin();iCut != cuts.end(); iCut++)
  {
    num[iCut->first]->Divide(denom);
    num[iCut->first]->Write();
  }
}
void alphaTree::plotRate1D()
{
  TFile* fOut = new TFile("ratePlots1D.root","recreate");
  fOut->cd();
  std::map<TString,TH1D*> rateVar;
  rateVar["HTUct"] = new TH1D("HtUctRate",";HT;",500,0,500);
  rateVar["metUct"] = new TH1D("metUctRate",";HT;",500,0,500);
  rateVar["metFullMenuUct"] = new TH1D("metUctFullMenuRate",";HT;",500,0,500);
  rateVar["metNewFullMenuUct"] = new TH1D("metUctNewFullMenuRate",";HT;",500,0,500);
  rateVar["MHTUct"] = new TH1D("MHtUctRate",";MHT;",500,0,500);
  rateVar["jet1Uct"] = new TH1D("Jet1UctRate",";pt;",500,0,500);
  rateVar["jet2Uct"] = new TH1D("Jet2UctRate",";pt;",500,0,500);
  rateVar["jet3Uct"] = new TH1D("Jet3UctRate",";pt;",500,0,500);
  rateVar["jet4Uct"] = new TH1D("Jet4UctRate",";pt;",500,0,500);

  rateVar["HTS2Global"] = new TH1D("HtS2GlobalRate",";HT;",500,0,500);
  rateVar["MHTS2Global"] = new TH1D("MHtS2GlobalRate",";MHT;",500,0,500);
  rateVar["jet1S2Global"] = new TH1D("Jet1S2GlobalRate",";pt;",500,0,500);
  rateVar["jet2S2Global"] = new TH1D("Jet2S2GlobalRate",";pt;",500,0,500);
  rateVar["jet3S2Global"] = new TH1D("Jet3S2GlobalRate",";pt;",500,0,500);
  rateVar["jet4S2Global"] = new TH1D("Jet4S2GlobalRate",";pt;",500,0,500);
  rateVar["alphaTCalo"] = new TH1D("alphaTCalo",";AT;",20,0,1);
  rateVar["htCalo"] = new TH1D("htCalo",";Ht;",100,0,1000);

  Long64_t nentries = fChain->GetEntriesFast();
  nentries = 10000000;
  Long64_t nbytes = 0, nb = 0;
  int test = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;
    //for (std::map<TString,TH1D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
    double ht40 = 0.;
    for (int i = 0; i < jetPtsCalo->size();i++) {if (jetPtsCalo->at(i) >= 40) {ht40 += jetPtsCalo->at(i);};} 
    /*
       if ((jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct >= 70 || (jetPtsUct->size() > 1 && jetPtsUct->at(1)>=60 && metPtUct>60))
       {
       rateVar["HTUct"]->Fill(1000);
       }
       else
       {
       } 
       */
    rateVar["HTUct"]->Fill(htUct);
    rateVar["htCalo"]->Fill(ht40);
    rateVar["MHTUct"]->Fill(mhtPtUct);  
    rateVar["metUct"]->Fill(metPtUct);  
    bool fullMenu = htUct >=175. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 100.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    bool fullMenuNew = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    if(jetPtsUct->size() > 0) rateVar["jet1Uct"]->Fill(jetPtsUct->at(0)); 
    /*
       if(htUct > 155. ||(jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct >= 70 || (jetPtsUct->size() > 1 && jetPtsUct->at(1)>=60 && metPtUct>60))
       {
       rateVar["jet2Uct"]->Fill(1000); 
       }
       else
       {
       }
       */
    if(jetPtsUct->size() > 1) rateVar["jet2Uct"]->Fill(jetPtsUct->at(1)); 
    if(jetPtsUct->size() > 2) rateVar["jet3Uct"]->Fill(jetPtsUct->at(2)); 
    if(jetPtsUct->size() > 3) rateVar["jet4Uct"]->Fill(jetPtsUct->at(3)); 

    rateVar["HTS2Global"]->Fill(htS2Global);
    rateVar["MHTS2Global"]->Fill(mhtPtS2Global);  

    if(jetPtsS2Global->size() > 0) rateVar["jet1S2Global"]->Fill(jetPtsS2Global->at(0)); 
    if(jetPtsS2Global->size() > 1) rateVar["jet2S2Global"]->Fill(jetPtsS2Global->at(1)); 
    if(jetPtsS2Global->size() > 2) rateVar["jet3S2Global"]->Fill(jetPtsS2Global->at(2)); 
    if(jetPtsS2Global->size() > 3) rateVar["jet4S2Global"]->Fill(jetPtsS2Global->at(3)); 
    if(fullMenu)
    {
      rateVar["metFullMenuUct"]->Fill(1000);
    }
    else
    {
      rateVar["metFullMenuUct"]->Fill(metPtUct);
    }
    if(fullMenuNew)
    {
      rateVar["metNewFullMenuUct"]->Fill(1000);
    }
    else
    {
      rateVar["metNewFullMenuUct"]->Fill(metPtUct);
    }
    if(htCalo > 400) {rateVar["alphaTCalo"]->Fill(alphaTCalo);}

    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,TH1D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
  {
    TH1D * temp = makeCumu1D(iVar->second,(double)nentries/ZB_XSECTION);
    temp->SetName(iVar->first+"_rate");
    temp->Write();
  }
}

void alphaTree::plotEff1D()
{
  TFile* fOut = new TFile("effPlots1D.root","recreate");
  fOut->cd();

  std::map<TString,TH1D*> effVarUct;
  effVarUct["HT"] = new TH1D("HtUctRate",";HT;",500,0,500);
  effVarUct["MHT"] = new TH1D("MHtUctRate",";MHT;",500,0,500);
  effVarUct["jet1"] = new TH1D("Jet1UctRate",";pt;",500,0,500);
  effVarUct["jet2"] = new TH1D("Jet2UctRate",";pt;",500,0,500);
  effVarUct["jet3"] = new TH1D("Jet3UctRate",";pt;",500,0,500);
  effVarUct["jet4"] = new TH1D("Jet4UctRate",";pt;",500,0,500);

  //std::map<TString,TH1D*> effVar;
  std::map<TString,TH1D*> effVarS2;
  effVarS2["HT"] = new TH1D("HtRateS2Global",";HT;",500,0,500);
  effVarS2["MHT"] = new TH1D("MHtRateS2Global",";MHT;",500,0,500);
  effVarS2["jet1"] = new TH1D("Jet1RateS2Global",";pt;",500,0,500);
  effVarS2["jet2"] = new TH1D("Jet2RateS2Global",";pt;",500,0,500);
  effVarS2["jet3"] = new TH1D("Jet3RateS2Global",";pt;",500,0,500);
  effVarS2["jet4"] = new TH1D("Jet4RateS2Global",";pt;",500,0,500);

  std::map<TString,int> offlineCut;
  offlineCut["HT"] = 0;
  offlineCut["MHT"] = 0;
  offlineCut["jet1"] = 0;
  offlineCut["jet2"] = 0;
  offlineCut["jet3"] = 0;
  offlineCut["jet4"] = 0;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  int test = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   

    if (ientry < 0) break;
    //for (std::map<TString,TH1D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
    if (htGen5 >= 400) {offlineCut["HT"]++; effVarUct["HT"]->Fill(htUct); effVarS2["HT"]->Fill(htS2Global);}

    if (mhtPtGen5 >= 130) {offlineCut["MHT"]++; effVarUct["MHT"]->Fill(mhtPtUct);effVarS2["MHT"]->Fill(mhtPtS2Global);}

    if (jetPtsGen5->size() > 0 && jetPtsGen5->at(0) >= 200) {
      offlineCut["jet1"]++; 
      if(jetPtsUct->size() > 0) {effVarUct["jet1"]->Fill(jetPtsUct->at(0));}
      if(jetPtsS2Global->size() > 0) {effVarS2["jet1"]->Fill(jetPtsS2Global->at(0));}
    }

    if (jetPtsGen5->size() > 1 && jetPtsGen5->at(1) >= 150) {
      offlineCut["jet2"]++; 
      if(jetPtsUct->size() > 1) {effVarUct["jet2"]->Fill(jetPtsUct->at(1));}
      if(jetPtsS2Global->size() > 1) {effVarS2["jet2"]->Fill(jetPtsS2Global->at(1));}
    }
    if (jetPtsGen5->size() > 2 && jetPtsGen5->at(2) >= 100) {
      offlineCut["jet3"]++; 
      if(jetPtsUct->size() > 2) {effVarUct["jet3"]->Fill(jetPtsUct->at(2));}
      if(jetPtsS2Global->size() > 2) {effVarS2["jet3"]->Fill(jetPtsS2Global->at(2));}
    }
    if (jetPtsGen5->size() > 3 && jetPtsGen5->at(3) >= 50) {
      offlineCut["jet4"]++; 
      if(jetPtsUct->size() > 3) {effVarUct["jet4"]->Fill(jetPtsUct->at(3));}
      if(jetPtsS2Global->size() > 3) {effVarS2["jet4"]->Fill(jetPtsS2Global->at(3));}
    }


    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,TH1D*>::iterator iVar = effVarUct.begin(); iVar != effVarUct.end(); iVar++)
  {
    TH1D * temp = makeCumu1D(iVar->second,offlineCut[iVar->first]);
    temp->SetName(iVar->first+"Uct_eff");
    temp->Write();
  }
  for (std::map<TString,TH1D*>::iterator iVar = effVarS2.begin(); iVar != effVarS2.end(); iVar++)
  {
    TH1D * temp = makeCumu1D(iVar->second,offlineCut[iVar->first]);
    temp->SetName(iVar->first+"S2Global_eff");
    temp->Write();
  }
}

std::map<TString,TH2D *> alphaTree::plotRate2dQcd(TString tempstring2)
{

  std::map<TString,TH2D*> rateVar;
  rateVar["htMhtUct"] = new TH2D(tempstring2+"mhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);

  Long64_t nentries =fChain->GetEntriesFast();
  nentries = 1000;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;

    rateVar["htMhtUct"]->Fill(htUct,mhtDivHtUct);
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,TH2D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
  {
    TH2D * temp = makeCumu(iVar->second, (double)nentries/(weight *1.4E-2) ); /// 1.4 E34 * 1E -36 to convert from pb to cm2
    iVar->second = (TH2D*) temp->Clone();
    delete temp;
  }
  return rateVar;
}
std::map<TString,TH1D*> alphaTree::plotRate1dQcd(TString temp)
{
  std::map<TString,TH1D*> rateVar;
  rateVar["HTUct"] = new TH1D(temp+"HtUctRate",";HT;",500,0,500);
  rateVar["jet1Uct"] = new TH1D(temp+"Jet1UctRate",";pt;",500,0,500);
  rateVar["jet2Uct"] = new TH1D(temp+"Jet2UctRate",";pt;",500,0,500);
  rateVar["jet3Uct"] = new TH1D(temp+"Jet3UctRate",";pt;",500,0,500);
  rateVar["jet4Uct"] = new TH1D(temp+"Jet4UctRate",";pt;",500,0,500);

  Long64_t nentries = fChain->GetEntriesFast();
  //nentries = 1000000;
  Long64_t nbytes = 0, nb = 0;
  int test = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;
    rateVar["HTUct"]->Fill(htUct);
    if(jetPtsUct->size() > 0) rateVar["jet1Uct"]->Fill(jetPtsUct->at(0)); 
    if(jetPtsUct->size() > 1) rateVar["jet2Uct"]->Fill(jetPtsUct->at(1)); 
    if(jetPtsUct->size() > 2) rateVar["jet3Uct"]->Fill(jetPtsUct->at(2)); 
    if(jetPtsUct->size() > 3) rateVar["jet4Uct"]->Fill(jetPtsUct->at(3)); 

    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,TH1D*>::iterator iVar = rateVar.begin(); iVar != rateVar.end(); iVar++)
  {
    TH1D * temp = makeCumu1D(iVar->second, (double)nentries/(weight *1.4E-2));
    iVar->second = (TH1D*) temp->Clone();
    delete temp;
  }
  return rateVar;
}

void alphaTree::plotRate()
{
  TFile* fOut = new TFile("ratePlots.root","recreate");
  fOut->cd();
  TH2D* ratesMhtDivHt = new TH2D("mhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);
  TH2D * fullRateMhtDivHt = new TH2D("fullMhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);
  TH2D * fullRateMht = new TH2D("fullMht",";Ht;mhtDivHt",50,0,200.,50,0.,200.);
  TH2D * fullRateMETMht = new TH2D("fullMETMht",";Ht;mhtDivHt",50,0,200.,50,0.,200.);
  TH2D * fullRateMETMhtDivHt = new TH2D("fullMETMhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);
  TH2D * addRateMhtDivHt = new TH2D("addMhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);
  TH2D * addFullRateMhtDivHt = new TH2D("addFullMhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);
  TH2D * addFullRateMETMhtDivHt = new TH2D("addFullMETMhtDivHt",";Ht;mhtDivHt",50,0,200.,20,0.,1.);

  TH2D* ratesDphi = new TH2D("dphi",";jet2pT;dPhi12",50,0.,200.,10,-0.5,9.5);
  TH2D * fullRateDphi = new TH2D("fulldphi",";jet2pT;dPhi12",50,0.,200.,10,-0.5,9.5);

  Long64_t nentries =fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;

    ratesMhtDivHt->Fill(htUct,mhtDivHtUct);
    if(!(htUct > 125 && mhtDivHtUct > 0.3)) addRateMhtDivHt->Fill(htUct,mhtDivHtUct);
    if(htUct >= 125.)
    {
      if(jetPtsUct->size() > 1) {ratesDphi->Fill(jetPtsUct->at(1),9-abs((jetDPhi12Uct)/0.3409));}
    }

    bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.); 
    if (fullmenu)
    {
      fullRateMhtDivHt->Fill(1000.,1000.);
      fullRateMht->Fill(1000.,1000.);
      fullRateDphi->Fill(1000.,1000.);
      //if(!(htUct > 125 && mhtDivHtUct > 0.3)) addFullRateMhtDivHt->Fill(1000.,1000.);
    }
    else 
    {
      fullRateMhtDivHt->Fill(htUct,mhtDivHtUct);
      fullRateMht->Fill(htUct,mhtPtUct);
      if(htUct >= 125. && jetPtsUct->size() > 1) {fullRateDphi->Fill(jetPtsUct->at(1),9-abs((jetDPhi12Uct)/0.3409));}
      if(!(htUct > 125 && mhtDivHtUct > 0.3))  addFullRateMhtDivHt->Fill(htUct,mhtDivHtUct);
    }

    bool fullmenuMET = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct > 70.; 
    if (fullmenuMET)
    {
      fullRateMETMhtDivHt->Fill(1000.,1000.);
      fullRateMETMht->Fill(1000.,1000.);
      //if(!(htUct > 125 && mhtDivHtUct > 0.3)) addFullRateMhtDivHt->Fill(1000.,1000.);
    }
    else 
    {
      fullRateMETMhtDivHt->Fill(htUct,mhtDivHtUct);
      fullRateMETMht->Fill(htUct,mhtPtUct);
      if(!(htUct > 125 && mhtDivHtUct > 0.3))  addFullRateMETMhtDivHt->Fill(htUct,mhtDivHtUct);
    }
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  TH2D * rates2 = makeCumu(ratesMhtDivHt, (double) nentries/ZB_XSECTION);
  TH2D * rates3 = makeCumu(ratesDphi,(double)  nentries/ZB_XSECTION);
  TH2D * rates4 = makeCumu(addRateMhtDivHt, (double) nentries/ZB_XSECTION);
  TH2D * rates5 = makeCumu(fullRateMhtDivHt,(double)  nentries/ZB_XSECTION);
  TH2D * rates6 = makeCumu(addFullRateMhtDivHt, (double) nentries/ZB_XSECTION);
  TH2D * rates7 = makeCumu(fullRateDphi, (double) nentries/ZB_XSECTION);
  TH2D * rates8 = makeCumu(fullRateMETMhtDivHt, (double) nentries/ZB_XSECTION);
  TH2D * rates9 = makeCumu(addFullRateMETMhtDivHt,(double)  nentries/ZB_XSECTION);
  TH2D * rates10 = makeCumu(fullRateMht, (double) nentries/ZB_XSECTION);
  TH2D * rates11 = makeCumu(fullRateMETMht, (double) nentries/ZB_XSECTION);
  //temp.Write();
  rates2->Write();
  rates3->Write();
  rates4->Write();
  rates5->Write();
  rates6->Write();
  rates7->Write();
  rates8->Write();
  rates9->Write();
  rates10->Write();
  rates11->Write();
}
void alphaTree::plotPtForwardCut()
{

  TH1D * histo = new TH1D("jetPtVeto",";pT/GeV;",20,0,1000);  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    nb = fChain->GetEntry(jentry); nbytes+=nb;   
    if (ientry < 0) break;
    bool veto = false;
    for (unsigned int i = 0; i < jetPtsAllUct->size(); i++)
    {
      if (jetPtsAllUct->at(i) > 50 && abs(jetEtasAllUct->at(i)) > 3) {veto = true; break;}
    }
    for (unsigned int i = 0; i < jetPtsUct->size(); i++)
    {
      if (!veto) histo->Fill(jetPtsUct->at(i));
    }
  }
  histo->Draw();
}

void alphaTree::plotEff()
{
  TFile* fOut = new TFile("effPlots.root","recreate");
  fOut->cd();
  std::map<TString,bool> analysisBins;
  analysisBins["Bin1"] =false; 
  analysisBins["Bin2"] = false;
  analysisBins["Bin3"] = false;
  analysisBins["Bin4"] = false;
  std::map<TString,bool> monoBins;
  monoBins["Bin1"] =false; 
  monoBins["Bin2"] = false;
  monoBins["Bin3"] = false;

  std::map<TString,bool> mult;
  //mult["1Jet"] = false;
  mult["2Jet"] =  false;
  mult["3Jet"] = false;
  mult["ge4Jet"] = false;

  std::map<TString,bool> hardlead;
  hardlead["Hard"] = false;
  hardlead["noHard"] = false;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  std::map <TString,TH2D *> numMhtDivHt;
  std::map <TString,TH2D *> numMhtFull;
  std::map <TString,TH2D *> numMhtFullMet;
  std::map <TString,TH2D *> numMhtDivHtAdd;
  std::map <TString,TH2D *> numMhtDivHtFull;
  std::map <TString,TH2D *> numMhtDivHtFullAdd;
  std::map <TString,TH2D *> numMhtDivHtFullMet;
  std::map <TString,TH2D *> numMhtDivHtFullMetAdd;
  std::map <TString,TH2D *> numDphi;
  std::map <TString,double> denomEntries;

  std::map <TString,TDirectory *> dir;
  int jetPtMin = 50;
  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
  {
    //dir[iMult->first] = fOut->mkdir(iMult->first);
    for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
    {
      for (std::map<TString,bool>::iterator iCut = analysisBins.begin(); iCut != analysisBins.end();iCut++)
      {
	numMhtDivHt[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHt",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	numMhtFull[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtFull",";hTUct;mhtUct",50,0.,200.,50,0.,200.); 
	numMhtDivHtAdd[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtAdd",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	numMhtDivHtFull[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtFull",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	numMhtFullMet[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtFullMet",";hTUct;mhtUct",50,0.,200.,50,0.,200.); 
	numMhtDivHtFullAdd[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtFullAdd",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	numMhtDivHtFullMet[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtFullMet",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	numMhtDivHtFullMetAdd[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtFullMetAdd",";hTUct;mhtDivHtUct",50,0.,200.,20,0.,1.); 
	//	numMhtDivHtDhtDivHt[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"mhtDivHtDhtDivHt",";dhtDivHtUct;mhtDivHtUct",100,0,1,100,0,1.); 
	numDphi[iCut->first+iMult->first+iHard->first] = new TH2D(iHard->first+iCut->first+iMult->first+"dphi",";jet2pT;dPhi12",50,0.,200.,10,-0.5,9.5); 
	denomEntries[iCut->first+iMult->first+iHard->first] = 0.;
	//l1Turn[iCut->first+iMult->first]->SetDirectory(dir);
      }
    } 
  }
  for (std::map<TString,bool>::iterator iCut = monoBins.begin(); iCut != monoBins.end();iCut++)
  {
    numMhtDivHt["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHt",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numMhtDivHtAdd["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtAdd",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numMhtFull["1Jet"+iCut->first] = new TH2D(iCut->first+"1JetmhtFull",";hTUct;mhtUct",50,0.,200.,50,0.,200.); 
    numMhtFullMet["1Jet"+iCut->first] = new TH2D(iCut->first+"1JetmhtFullMet",";hTUct;mhtUct",50,0.,200.,50,0.,200.); 
    numMhtDivHtFull["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtFull",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numMhtDivHtFullAdd["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtFullAdd",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numMhtDivHtFullMet["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtMet",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numMhtDivHtFullMetAdd["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtFullMetAdd",";hTUct;mhtDivHtUct",50,0,200,20,0,1.);
    numDphi["1Jet"+iCut->first]= new TH2D(iCut->first+"1Jetdphi",";jet2pT;dPhi12",300,0,300,10,-0.5,9.5);
    //numMhtDivHtDhtDivHt["1Jet"+iCut->first]= new TH2D(iCut->first+"1JetmhtDivHtDhtDivHt",";dhtDivHtUct;mhtDivHtUct",100,0,1,100,0,1.);
    denomEntries["1Jet"+iCut->first]=0;
  }
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //mult["1Jet" && iHard->second;] = jetPtsGen5->at(1) < 50;
    if(jetPtsGen5->size() > 1)
    {
      hardlead["noHard"] = jetPtsGen5->at(1) > jetPtMin;
      hardlead["Hard"] = jetPtsGen5->at(1) >= 90;
    }
    else
    {
      hardlead["noHard"] = false;
      hardlead["Hard"] = false;

    }
    for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
    {
      mult["2Jet"] = (iHard->second &&  multiplicity50Gen5 == 2);
      mult["3Jet"] = (iHard->second&&   multiplicity50Gen5 ==3);
      mult["ge4Jet"] = (iHard->second&& multiplicity50Gen5 >= 4);

      //Simplified A bins
      analysisBins["Bin1"] = (htGen5 >= 200. && htGen5 < 300. && alphaTGen5 > 0.65); 
      analysisBins["Bin2"] = (htGen5 >= 300. && htGen5 < 400. && alphaTGen5 > 0.60); 
      analysisBins["Bin3"] = (htGen5 >= 400. && htGen5 < 500. && alphaTGen5 > 0.55);
      analysisBins["Bin4"] = (htGen5 >= 180. && alphaTGen5 > 0.6); 


      for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
      {
	if(iMult->second)
	{
	  for (std::map<TString,bool>::iterator iCut = analysisBins.begin(); iCut != analysisBins.end();iCut++)
	  {
	    if (iCut->second) 
	    {
	      //std::cout << jetPtsGen5->at(3) << iMult->first << std::endl;
	      numMhtDivHt[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
	      bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
	      if (fullmenu)
	      {
		numMhtDivHtFull[iCut->first+iMult->first+iHard->first]->Fill(1000.,1000.);
		numMhtFull[iCut->first+iMult->first+iHard->first]->Fill(1000.,1000.);
	      }
	      else
	      {
		numMhtDivHtFull[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
		numMhtFull[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtPtUct);
		if(!(htUct > 125 && mhtDivHtUct > 0.3))   numMhtDivHtFullAdd[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
	      }
	      bool fullmenuMet = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct > 70;
	      if (fullmenuMet)
	      {
		numMhtDivHtFullMet[iCut->first+iMult->first+iHard->first]->Fill(1000.,1000.);
		numMhtFullMet[iCut->first+iMult->first+iHard->first]->Fill(1000.,1000.);
	      }
	      else
	      {
		numMhtDivHtFullMet[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
		numMhtFullMet[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtPtUct);
		if(!(htUct > 125 && mhtDivHtUct > 0.3))   numMhtDivHtFullMetAdd[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
	      }


	      if (!(htUct > 125 && mhtDivHtUct > 0.3)) numMhtDivHtAdd[iCut->first+iMult->first+iHard->first]->Fill(htUct,mhtDivHtUct);
	      //numMhtDivHtDhtDivHt[iCut->first+iMult->first+iHard->first]->Fill(dhtUct/htUct,mhtDivHtUct);
	      if (jetPtsUct->size() > 1 && htUct>125)
	      {
		numDphi[iCut->first+iMult->first+iHard->first]->Fill(jetPtsUct->at(1),9-abs((jetDPhi12Uct)/0.3409));
	      }
	      denomEntries[iCut->first+iMult->first+iHard->first]+=1.;
	      /*   int nXbins = fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetNbinsX();
		   int nYbins = fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetNbinsY();
		   bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
		   for (int xbins = 0; xbins <= nXbins+1; xbins++)
		   {
		   for (int ybins = 0; ybins <= nYbins+1; ybins++)
		   {
		   if((htUct >=(fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetXaxis()->GetBinLowEdge(xbins)) && mhtDivHtUct >= (fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetYaxis()->GetBinLowEdge(ybins))) || fullmenu )
		   {
		   fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->Fill(fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetXaxis()->GetBinCenter(xbins),fullEffMhtDivHt[iCut->first+iMult->first+iHard->first]->GetYaxis()->GetBinCenter(ybins),1.);
		   }
		   }
		   }
		   */
	    }
	  } 
	}
      }
    }
    monoBins["Bin1"] = (htGen5 >= 200. && htGen5 < 300.);
    monoBins["Bin2"] = (htGen5 >= 300. && htGen5 < 400.);
    monoBins["Bin3"] = (htGen5 >= 400. && htGen5 < 500.);
    if(multiplicity50Gen5==1) 
    {
      //Simplified A bins
      for (std::map<TString,bool>::iterator iCut = monoBins.begin(); iCut != monoBins.end();iCut++)
      {

	if(iCut->second)
	{
	  numMhtDivHt["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct);
	  if (!(htUct > 125 && mhtDivHtUct > 0.3)) numMhtDivHtAdd["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct); 
	  //numMhtDivHtDhtDivHt["1Jet"+iCut->first]->Fill(dhtUct/htUct,mhtDivHtUct);
	  if(jetPtsUct->size()>1 && htUct > 125)
	  {
	    numDphi["1Jet"+iCut->first]->Fill(jetPtsUct->at(1),9-abs((jetDPhi12Uct)/0.3409));
	  }
	  bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
	  if (fullmenu)
	  {
	    numMhtDivHtFull["1Jet"+iCut->first]->Fill(1000.,1000.);
	    numMhtFull["1Jet"+iCut->first]->Fill(1000.,1000.);
	  }
	  else
	  {
	    numMhtDivHtFull["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct);
	    numMhtFull["1Jet"+iCut->first]->Fill(htUct,mhtPtUct);
	    if(!(htUct > 125 && mhtDivHtUct > 0.3))   numMhtDivHtFullAdd["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct);
	  }
	  bool fullmenuMet = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.) || metPtUct > 70;
	  if (fullmenuMet)
	  {
	    numMhtDivHtFullMet["1Jet"+iCut->first]->Fill(1000.,1000.);
	    numMhtFullMet["1Jet"+iCut->first]->Fill(1000.,1000.);
	  }
	  else
	  {
	    numMhtDivHtFullMet["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct);
	    numMhtFullMet["1Jet"+iCut->first]->Fill(htUct,mhtPtUct);
	    if(!(htUct > 125 && mhtDivHtUct > 0.3))   numMhtDivHtFullMetAdd["1Jet"+iCut->first]->Fill(htUct,mhtDivHtUct);
	  }
	  /*int nXbins = fullEffMhtDivHt["1Jet"+iCut->first]->GetNbinsX();
	    int nYbins = fullEffMhtDivHt["1Jet"+iCut->first]->GetNbinsY();
	    bool fullmenu = htUct >=200. || (jetPtsUct->size() > 1 && jetPtsUct->at(1) >= 120.) || (jetPtsUct->size() > 3 && jetPtsUct->at(3) >= 60.);
	    for (int xbins = 0; xbins <= nXbins+1; xbins++)
	    {
	    for (int ybins = 0; ybins <= nYbins+1; ybins++)
	    {
	    if((htUct >=(fullEffMhtDivHt["1Jet"+iCut->first]->GetXaxis()->GetBinLowEdge(xbins)) && mhtDivHtUct >= (fullEffMhtDivHt["1Jet"+iCut->first]->GetYaxis()->GetBinLowEdge(ybins))) || fullmenu )
	    {
	    fullEffMhtDivHt["1Jet"+iCut->first]->Fill(fullEffMhtDivHt["1Jet"+iCut->first]->GetXaxis()->GetBinCenter(xbins),fullEffMhtDivHt["1Jet"+iCut->first]->GetYaxis()->GetBinCenter(ybins),1.);
	    }
	    }
	    }*/
	  denomEntries["1Jet"+iCut->first]+=1;
	}
      }
    } 
    if (jentry%10000 == 0) std::cout << std::setprecision(4) << jentry*100./nentries << "%     " << "\r" <<std::flush;
  }
  for (std::map<TString,bool>::iterator iMult = mult.begin(); iMult != mult.end();iMult++)
  {
    for (std::map<TString,bool>::iterator iHard = hardlead.begin(); iHard != hardlead.end();iHard++)
    {
      dir[iMult->first+iHard->first] = fOut->mkdir(iMult->first+iHard->first);
      dir[iMult->first+iHard->first]->cd();
      for (std::map<TString,bool>::iterator iCut = analysisBins.begin(); iCut != analysisBins.end();iCut++)
      {
	//numMhtDivHt[iCut->first+iMult->first]->Sumw2();
	//numMhtDivHt[iCut->first+iMult->first]->Divide(denom[iMult->first]);
	//std::cout << iMult->first << numMhtDivHt[iCut->first+iMult->first]->GetEntries() << " " << denomEntries[iCut->first+iMult->first] << std::endl;
	TH2D * temp= makeCumu(numDphi[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp2= makeCumu(numMhtDivHt[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp3= makeCumu(numMhtDivHtAdd[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp4= makeCumu(numMhtDivHtFull[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp5= makeCumu(numMhtDivHtFullAdd[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp6= makeCumu(numMhtDivHtFullMet[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp7= makeCumu(numMhtDivHtFullMetAdd[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp8= makeCumu(numMhtFull[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	TH2D * temp9= makeCumu(numMhtFullMet[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	//TH2D * temp3= makeCumu(numMhtDivHtDhtDivHt[iCut->first+iMult->first+iHard->first],denomEntries[iCut->first+iMult->first+iHard->first]);
	temp->Write();
	temp2->Write();
	temp3->Write();
	temp4->Write();
	temp5->Write();
	temp6->Write();
	temp7->Write();
	temp8->Write();
	temp9->Write();
	//temp3->Write();
      } 
    }
  }
  //std::cout << iMult->first << num[iCut->first+iMult->first]->GetEntries() << " " << denomEntries[iCut->first+iMult->first] << std::endl;
  dir["1Jet"] = fOut->mkdir("1Jet");
  dir["1Jet"]->cd();
  for (std::map<TString,bool>::iterator iCut = monoBins.begin(); iCut != monoBins.end();iCut++)
  {
    TH2D * temp= makeCumu(numDphi["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp2= makeCumu(numMhtDivHt["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp3= makeCumu(numMhtDivHtAdd["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp4= makeCumu(numMhtDivHtFull["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp5= makeCumu(numMhtDivHtFullAdd["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp6= makeCumu(numMhtDivHtFullMet["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp7= makeCumu(numMhtDivHtFullMetAdd["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    TH2D * temp8= makeCumu(numMhtFullMet["1Jet"+iCut->first],denomEntries["1Jet"+iCut->first]);
    temp->Write();
    temp2->Write();
    temp3->Write();
    temp4->Write();
    temp5->Write();
    temp6->Write();
    temp7->Write();
    temp8->Write();
  }
}

TH2D * makeCumu(TH2D * input, double norm){

  TH2D * output = new TH2D(*input);
  //output=input;
  //int norm = input->GetEntries();
  //output->SetBinContent(0,1.);
  int nXbins = input->GetNbinsX();
  int nYbins = input->GetNbinsY();
  for (int xbins = 0; xbins <= nXbins+1; xbins++)
  {
    int dummy = 0;
    //int dummy = input->GetBinContent(xbins,nYbins+1);
    //output->SetBinContent(xbins,nYbins+1,dummy);

    for (int ybins = 0; ybins <= nYbins+1; ybins++)
    {
      dummy += input->GetBinContent(xbins,nYbins+1-ybins);
      output->SetBinContent(xbins,nYbins+1-ybins,dummy);
    }
  } 
  for (int ybins = 0; ybins <= nYbins+1; ybins++)
  {
    int dummy=0;
    //int dummy = input->GetBinContent(nXbins+1,ybins);
    //output->SetBinContent(nXbins+1,ybins,(double)dummy/norm);
    for (int xbins = 0; xbins <= nXbins+1; xbins++)
    {
      dummy += output->GetBinContent(nXbins+1-xbins,ybins);
      output->SetBinContent(nXbins+1-xbins,ybins,(double)dummy/norm);

      //output->SetBinContent(xbins,ybins,ZB_XSECTION*(1-(double)dummy/norm));
    }
  } 

  return output;

}

TH1D * makeCumu1D(TH1D * input, double norm){

  TH1D * output = (TH1D*)input->Clone();
  int nXbins = input->GetNbinsX();
  int dummy=0;
  for (int xbins = 0; xbins <= nXbins+1; xbins++)
  {
    dummy += output->GetBinContent(nXbins+1-xbins);
    output->SetBinContent(nXbins+1-xbins,(double)dummy/norm);

  }
  return output;

}
TH1 * makeCumuInv1D(TH1 * input, double norm){

  TH1 * output = (TH1*)input->Clone();
  int nXbins = input->GetNbinsX();
  int dummy=0;
  for (int xbins = 0; xbins <= nXbins+1; xbins++)
  {
    dummy += output->GetBinContent(nXbins+1-xbins);
    output->SetBinContent(nXbins+1-xbins,(double)dummy/norm);

  }
  return output;

}
TH2 * makeCumuInv2D(TH2 * input, double norm){

  TH2 * output = (TH2*)input->Clone();
  int nYbins = input->GetNbinsY();
  int nXbins = input->GetNbinsX();
  for (int xbins = 0; xbins <= nXbins+1; xbins++)
  {
    int dummy=0;
    for (int ybins = 0; ybins <= nYbins+1; ybins++)
    {
      dummy += output->GetBinContent(xbins,nYbins+1-ybins);
      output->SetBinContent(xbins,nYbins+1-ybins,(double)dummy/norm);

    }
  }
  return output;

}
