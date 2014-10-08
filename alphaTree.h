//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 21 09:57:05 2014 by ROOT version 5.34/07
// from TTree Ntuple/Ntuple
// found on file: ../output/TTbar_2014-07-16.root
//////////////////////////////////////////////////////////

#ifndef alphaTree_h
#define alphaTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class alphaTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t jetVeto50;
   Bool_t jetVeto30;
   Float_t         mhtPtUct;
   Float_t         mhtPhiUct;
   Float_t         metPtUct;
   Float_t         metPhiUct;
   Float_t         htUct;
   Float_t         etUct;
   Float_t         dhtUct;
   UInt_t          multiplicityUct;
   Float_t         mhtDivHtUct;
   Float_t         alphaTUct;
   Float_t         jetDPhi12Uct;
   vector<float>   *jetPtsUct;
   vector<float>   *jetPhisUct;
   vector<float>   *jetEtasUct;
   vector<float>   *jetPtsAllUct;
   vector<float>   *jetPhisAllUct;
   vector<float>   *jetEtasAllUct;
   Float_t         mhtPtS2Global;
   Float_t         mhtPhiS2Global;
   Float_t         metPtS2Global;
   Float_t         metPhiS2Global;
   Float_t         htS2Global;
   Float_t         etS2Global;
   Float_t         dhtS2Global;
   UInt_t          multiplicityS2Global;
   Float_t         mhtDivHtS2Global;
   Float_t         alphaTS2Global;
   Float_t         jetDPhi12S2Global;
   vector<float>   *jetPtsS2Global;
   vector<float>   *jetPhisS2Global;
   vector<float>   *jetEtasS2Global;
   vector<float>   *jetPtsAllS2Global;
   vector<float>   *jetPhisAllS2Global;
   vector<float>   *jetEtasAllS2Global;
   Float_t         mhtPtGct;
   Float_t         mhtPhiGct;
   Float_t         metPtGct;
   Float_t         metPhiGct;
   Float_t         htGct;
   Float_t         etGct;
   Float_t         dhtGct;
   UInt_t          multiplicityGct;
   Float_t         mhtDivHtGct;
   Float_t         alphaTGct;
   Float_t         jetDPhi12Gct;
   vector<float>   *jetPtsGct;
   vector<float>   *jetPhisGct;
   vector<float>   *jetEtasGct;
   vector<float>   *jetPtsAllGct;
   vector<float>   *jetPhisAllGct;
   vector<float>   *jetEtasAllGct;
   Float_t         mhtPtGen5;
   Float_t         mhtPhiGen5;
   Float_t         metPtGen5;
   Float_t         metPhiGen5;
   Float_t         htGen5;
   Float_t         etGen5;
   Float_t         dhtGen5;
   UInt_t          multiplicityGen5;
   UInt_t          multiplicity50Gen5;
   Float_t         mhtDivHtGen5;
   Float_t         alphaTGen5;
   Float_t         jetDPhi12Gen5;
   vector<float>   *jetPtsGen5;
   vector<float>   *jetPhisGen5;
   vector<float>   *jetEtasGen5;
   Float_t         mhtPtCalo;
   Float_t         mhtPhiCalo;
   Float_t         metPtCalo;
   Float_t         metPhiCalo;
   Float_t         htCalo;
   Float_t         etCalo;
   Float_t         dhtCalo;
   UInt_t          multiplicityCalo;
   Float_t         mhtDivHtCalo;
   Float_t         alphaTCalo;
   Float_t         jetDPhi12Calo;
   vector<float>   *jetPtsCalo;
   vector<float>   *jetPhisCalo;
   vector<float>   *jetEtasCalo;
   vector<float>   *jetPtsAllCalo;
   vector<float>   *jetPhisAllCalo;
   vector<float>   *jetEtasAllCalo;
   Float_t         mhtPtPf;
   Float_t         mhtPhiPf;
   Float_t         metPtPf;
   Float_t         metPhiPf;
   Float_t         htPf;
   Float_t         etPf;
   Float_t         dhtPf;
   UInt_t          multiplicity50Pf;
   Float_t         mhtDivHtPf;
   Float_t         alphaTPf;
   Float_t         jetDPhi12Pf;
   vector<float>   *jetPtsPf;
   vector<float>   *jetPhisPf;
   vector<float>   *jetEtasPf;
   vector<float>   *jetPtsAllPf;
   vector<float>   *jetPhisAllPf;
   vector<float>   *jetEtasAllPf;
   vector<float>   *jetMuonsPf;
   UInt_t          MuonNumber;
   UInt_t          run;
   UInt_t          nInt;
   Float_t         weight;
   UInt_t          lumi;
   ULong64_t       event;
   vector<double> *regionEt;
   // List of branches
   TBranch        *b_regionEt; //!
   TBranch        *b_jetVeto50; //!
   TBranch        *b_jetVeto30; //!
   TBranch        *b_mhtPtUct;   //!
   TBranch        *b_mhtPhiUct;   //!
   TBranch        *b_metPtUct;   //!
   TBranch        *b_metPhiUct;   //!
   TBranch        *b_htUct;   //!
   TBranch        *b_etUct;   //!
   TBranch        *b_dhtUct;   //!
   TBranch        *b_multiplicityUct;   //!
   TBranch        *b_mhtDivHtUct;   //!
   TBranch        *b_alphaTUct;   //!
   TBranch        *b_jetDPhi12Uct;   //!
   TBranch        *b_jetPtsUct;   //!
   TBranch        *b_jetPhisUct;   //!
   TBranch        *b_jetEtasUct;   //!
   TBranch        *b_jetPtsAllUct;   //!
   TBranch        *b_jetPhisAllUct;   //!
   TBranch        *b_jetEtasAllUct;   //!
   TBranch        *b_mhtPtS2Global;   //!
   TBranch        *b_mhtPhiS2Global;   //!
   TBranch        *b_metPtS2Global;   //!
   TBranch        *b_metPhiS2Global;   //!
   TBranch        *b_htS2Global;   //!
   TBranch        *b_etS2Global;   //!
   TBranch        *b_dhtS2Global;   //!
   TBranch        *b_multiplicityS2Global;   //!
   TBranch        *b_mhtDivHtS2Global;   //!
   TBranch        *b_alphaTS2Global;   //!
   TBranch        *b_jetDPhi12S2Global;   //!
   TBranch        *b_jetPtsS2Global;   //!
   TBranch        *b_jetPhisS2Global;   //!
   TBranch        *b_jetEtasS2Global;   //!
   TBranch        *b_jetPtsAllS2Global;   //!
   TBranch        *b_jetPhisAllS2Global;   //!
   TBranch        *b_jetEtasAllS2Global;   //!
   TBranch        *b_mhtPtGct;   //!
   TBranch        *b_mhtPhiGct;   //!
   TBranch        *b_metPtGct;   //!
   TBranch        *b_metPhiGct;   //!
   TBranch        *b_htGct;   //!
   TBranch        *b_etGct;   //!
   TBranch        *b_dhtGct;   //!
   TBranch        *b_multiplicityGct;   //!
   TBranch        *b_mhtDivHtGct;   //!
   TBranch        *b_alphaTGct;   //!
   TBranch        *b_jetDPhi12Gct;   //!
   TBranch        *b_jetPtsGct;   //!
   TBranch        *b_jetPhisGct;   //!
   TBranch        *b_jetEtasGct;   //!
   TBranch        *b_jetPtsAllGct;   //!
   TBranch        *b_jetPhisAllGct;   //!
   TBranch        *b_jetEtasAllGct;   //!
   TBranch        *b_mhtPtGen5;   //!
   TBranch        *b_mhtPhiGen5;   //!
   TBranch        *b_metPtGen5;   //!
   TBranch        *b_metPhiGen5;   //!
   TBranch        *b_htGen5;   //!
   TBranch        *b_etGen5;   //!
   TBranch        *b_dhtGen5;   //!
   TBranch        *b_multiplicityGen5;   //!
   TBranch        *b_multiplicity50Gen5;   //!
   TBranch        *b_mhtDivHtGen5;   //!
   TBranch        *b_alphaTGen5;   //!
   TBranch        *b_jetDPhi12Gen5;   //!
   TBranch        *b_jetPtsGen5;   //!
   TBranch        *b_jetPhisGen5;   //!
   TBranch        *b_jetEtasGen5;   //!
   TBranch        *b_mhtPtCalo;   //!
   TBranch        *b_mhtPhiCalo;   //!
   TBranch        *b_metPtCalo;   //!
   TBranch        *b_metPhiCalo;   //!
   TBranch        *b_htCalo;   //!
   TBranch        *b_etCalo;   //!
   TBranch        *b_dhtCalo;   //!
   TBranch        *b_multiplicityCalo;   //!
   TBranch        *b_mhtDivHtCalo;   //!
   TBranch        *b_alphaTCalo;   //!
   TBranch        *b_jetDPhi12Calo;   //!
   TBranch        *b_jetPtsCalo;   //!
   TBranch        *b_jetPhisCalo;   //!
   TBranch        *b_jetEtasCalo;   //!
   TBranch        *b_jetPtsAllCalo;   //!
   TBranch        *b_jetPhisAllCalo;   //!
   TBranch        *b_jetEtasAllCalo;   //!
   TBranch        *b_mhtPtPf;   //!
   TBranch        *b_mhtPhiPf;   //!
   TBranch        *b_metPtPf;   //!
   TBranch        *b_metPhiPf;   //!
   TBranch        *b_htPf;   //!
   TBranch        *b_etPf;   //!
   TBranch        *b_dhtPf;   //!
   TBranch        *b_multiplicity50Pf;   //!
   TBranch        *b_mhtDivHtPf;   //!
   TBranch        *b_alphaTPf;   //!
   TBranch        *b_jetDPhi12Pf;   //!
   TBranch        *b_jetPtsPf;   //!
   TBranch        *b_jetPhisPf;   //!
   TBranch        *b_jetEtasPf;   //!
   TBranch        *b_jetPtsAllPf;   //!
   TBranch        *b_jetPhisAllPf;   //!
   TBranch        *b_jetEtasAllPf;   //!
   TBranch        *b_jetMuonsPf;   //!
   TBranch        *b_run;   //!
   TBranch        *b_nint;   //!
   TBranch        *b_weight;   //!

   alphaTree(TTree *tree=0);
   virtual ~alphaTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     findRate();
   virtual void     findCalo();
   virtual void regionStudy();
   virtual void     plotPtForwardCut();
   virtual void     findEff();
   virtual void     dhtStudy();
   virtual void     onOffHt();
   virtual void     plotTurnOn();
   virtual std::map<TString,TH1D*> plotRate1dQcd(TString temp);
   virtual std::map<TString,TH2D *> plotRate2dQcd(TString tempstring2);
   virtual void     plotRate();
   virtual void     plotRate1D();
   virtual void     plotRateNint();
   virtual void     plotEff1D();
   virtual void     plotEff();
   virtual void     plotOfflineEff();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef alphaTree_cxx
alphaTree::alphaTree(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../output/TTbar_2014-07-16.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../output/TTbar_2014-07-16.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("../output/TTbar_2014-07-16.root:/MakeL1Trees");
    dir->GetObject("Ntuple",tree);

  }
  Init(tree);
}

alphaTree::~alphaTree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t alphaTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t alphaTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void alphaTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  jetPtsUct = 0;
  jetPhisUct = 0;
  jetEtasUct = 0;
  jetPtsAllUct = 0;
  jetPhisAllUct = 0;
  jetEtasAllUct = 0;
  jetPtsS2Global = 0;
  jetPhisS2Global = 0;
  jetEtasS2Global = 0;
  jetPtsAllS2Global = 0;
  jetPhisAllS2Global = 0;
  jetEtasAllS2Global = 0;
  jetPtsGct = 0;
  jetPhisGct = 0;
  jetEtasGct = 0;
  jetPtsAllGct = 0;
  jetPhisAllGct = 0;
  jetEtasAllGct = 0;
  jetPtsGen5 = 0;
  jetPhisGen5 = 0;
  jetEtasGen5 = 0;
  jetPtsCalo = 0;
  jetPhisCalo = 0;
  jetEtasCalo = 0;
  jetPtsAllCalo = 0;
  jetPhisAllCalo = 0;
  jetEtasAllCalo = 0;
  jetPtsPf = 0;
  jetPhisPf = 0;
  jetEtasPf = 0;
  jetPtsAllPf = 0;
  jetPhisAllPf = 0;
  jetEtasAllPf = 0;
  jetMuonsPf = 0;
  regionEt = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  //fChain->SetBranchAddress("regionEt", &regionEt, &b_regionEt);
  //fChain->SetBranchAddress("jetVeto50",&jetVeto50,&b_jetVeto50);
  //fChain->SetBranchAddress("jetVeto30",&jetVeto30,&b_jetVeto30);
  fChain->SetBranchAddress("mhtPtUct", &mhtPtUct, &b_mhtPtUct);
  fChain->SetBranchAddress("mhtPhiUct", &mhtPhiUct, &b_mhtPhiUct);
  fChain->SetBranchAddress("metPtUct", &metPtUct, &b_metPtUct);
  fChain->SetBranchAddress("metPhiUct", &metPhiUct, &b_metPhiUct);
  fChain->SetBranchAddress("htUct", &htUct, &b_htUct);
  fChain->SetBranchAddress("etUct", &etUct, &b_etUct);
  fChain->SetBranchAddress("dhtUct", &dhtUct, &b_dhtUct);
  fChain->SetBranchAddress("multiplicityUct", &multiplicityUct, &b_multiplicityUct);
  fChain->SetBranchAddress("mhtDivHtUct", &mhtDivHtUct, &b_mhtDivHtUct);
  fChain->SetBranchAddress("alphaTUct", &alphaTUct, &b_alphaTUct);
  fChain->SetBranchAddress("jetDPhi12Uct", &jetDPhi12Uct, &b_jetDPhi12Uct);
  fChain->SetBranchAddress("jetPtsUct", &jetPtsUct, &b_jetPtsUct);
  fChain->SetBranchAddress("jetPhisUct", &jetPhisUct, &b_jetPhisUct);
  fChain->SetBranchAddress("jetEtasUct", &jetEtasUct, &b_jetEtasUct);
  fChain->SetBranchAddress("jetPtsAllUct", &jetPtsAllUct, &b_jetPtsAllUct);
  fChain->SetBranchAddress("jetPhisAllUct", &jetPhisAllUct, &b_jetPhisAllUct);
  fChain->SetBranchAddress("jetEtasAllUct", &jetEtasAllUct, &b_jetEtasAllUct);

  fChain->SetBranchAddress("mhtPtS2Global", &mhtPtS2Global, &b_mhtPtS2Global);
  fChain->SetBranchAddress("mhtPhiS2Global", &mhtPhiS2Global, &b_mhtPhiS2Global);
  //fChain->SetBranchAddress("metPtS2Global", &metPtS2Global, &b_metPtS2Global);
  //fChain->SetBranchAddress("metPhiS2Global", &metPhiS2Global, &b_metPhiS2Global);
  fChain->SetBranchAddress("htS2Global", &htS2Global, &b_htS2Global);
  //fChain->SetBranchAddress("etS2Global", &etS2Global, &b_etS2Global);
  fChain->SetBranchAddress("dhtS2Global", &dhtS2Global, &b_dhtS2Global);
  fChain->SetBranchAddress("multiplicityS2Global", &multiplicityS2Global, &b_multiplicityS2Global);
  fChain->SetBranchAddress("mhtDivHtS2Global", &mhtDivHtS2Global, &b_mhtDivHtS2Global);
  fChain->SetBranchAddress("alphaTS2Global", &alphaTS2Global, &b_alphaTS2Global);
  fChain->SetBranchAddress("jetDPhi12S2Global", &jetDPhi12S2Global, &b_jetDPhi12S2Global);
  fChain->SetBranchAddress("jetPtsS2Global", &jetPtsS2Global, &b_jetPtsS2Global);
  fChain->SetBranchAddress("jetPhisS2Global", &jetPhisS2Global, &b_jetPhisS2Global);
  fChain->SetBranchAddress("jetEtasS2Global", &jetEtasS2Global, &b_jetEtasS2Global);
//  fChain->SetBranchAddress("jetPtsAllS2Global", &jetPtsAllS2Global, &b_jetPtsAllS2Global);
//  fChain->SetBranchAddress("jetPhisAllS2Global", &jetPhisAllS2Global, &b_jetPhisAllS2Global);
//  fChain->SetBranchAddress("jetEtasAllS2Global", &jetEtasAllS2Global, &b_jetEtasAllS2Global);

  fChain->SetBranchAddress("mhtPtGct", &mhtPtGct, &b_mhtPtGct);
  fChain->SetBranchAddress("mhtPhiGct", &mhtPhiGct, &b_mhtPhiGct);
  fChain->SetBranchAddress("metPtGct", &metPtGct, &b_metPtGct);
  fChain->SetBranchAddress("metPhiGct", &metPhiGct, &b_metPhiGct);
  fChain->SetBranchAddress("htGct", &htGct, &b_htGct);
  fChain->SetBranchAddress("etGct", &etGct, &b_etGct);
  fChain->SetBranchAddress("dhtGct", &dhtGct, &b_dhtGct);
  fChain->SetBranchAddress("multiplicityGct", &multiplicityGct, &b_multiplicityGct);
  fChain->SetBranchAddress("mhtDivHtGct", &mhtDivHtGct, &b_mhtDivHtGct);
  fChain->SetBranchAddress("alphaTGct", &alphaTGct, &b_alphaTGct);
  fChain->SetBranchAddress("jetDPhi12Gct", &jetDPhi12Gct, &b_jetDPhi12Gct);
  fChain->SetBranchAddress("jetPtsGct", &jetPtsGct, &b_jetPtsGct);
  fChain->SetBranchAddress("jetPhisGct", &jetPhisGct, &b_jetPhisGct);
  fChain->SetBranchAddress("jetEtasGct", &jetEtasGct, &b_jetEtasGct);
  fChain->SetBranchAddress("jetPtsAllGct", &jetPtsAllGct, &b_jetPtsAllGct);
  fChain->SetBranchAddress("jetPhisAllGct", &jetPhisAllGct, &b_jetPhisAllGct);
  fChain->SetBranchAddress("jetEtasAllGct", &jetEtasAllGct, &b_jetEtasAllGct);
  fChain->SetBranchAddress("mhtPtGen5", &mhtPtGen5, &b_mhtPtGen5);
  fChain->SetBranchAddress("mhtPhiGen5", &mhtPhiGen5, &b_mhtPhiGen5);
  fChain->SetBranchAddress("metPtGen5", &metPtGen5, &b_metPtGen5);
  fChain->SetBranchAddress("metPhiGen5", &metPhiGen5, &b_metPhiGen5);
  fChain->SetBranchAddress("htGen5", &htGen5, &b_htGen5);
  fChain->SetBranchAddress("etGen5", &etGen5, &b_etGen5);
  fChain->SetBranchAddress("dhtGen5", &dhtGen5, &b_dhtGen5);
  fChain->SetBranchAddress("multiplicityGen5", &multiplicityGen5, &b_multiplicityGen5);
  fChain->SetBranchAddress("multiplicity50Gen5", &multiplicity50Gen5, &b_multiplicity50Gen5);
  fChain->SetBranchAddress("mhtDivHtGen5", &mhtDivHtGen5, &b_mhtDivHtGen5);
  fChain->SetBranchAddress("alphaTGen5", &alphaTGen5, &b_alphaTGen5);
  fChain->SetBranchAddress("jetDPhi12Gen5", &jetDPhi12Gen5, &b_jetDPhi12Gen5);
  fChain->SetBranchAddress("jetPtsGen5", &jetPtsGen5, &b_jetPtsGen5);
  fChain->SetBranchAddress("jetPhisGen5", &jetPhisGen5, &b_jetPhisGen5);
  fChain->SetBranchAddress("jetEtasGen5", &jetEtasGen5, &b_jetEtasGen5);
  fChain->SetBranchAddress("mhtPtCalo", &mhtPtCalo, &b_mhtPtCalo);
  fChain->SetBranchAddress("mhtPhiCalo", &mhtPhiCalo, &b_mhtPhiCalo);
  fChain->SetBranchAddress("metPtCalo", &metPtCalo, &b_metPtCalo);
  fChain->SetBranchAddress("metPhiCalo", &metPhiCalo, &b_metPhiCalo);
  fChain->SetBranchAddress("htCalo", &htCalo, &b_htCalo);
  fChain->SetBranchAddress("etCalo", &etCalo, &b_etCalo);
  fChain->SetBranchAddress("dhtCalo", &dhtCalo, &b_dhtCalo);
  fChain->SetBranchAddress("multiplicityCalo", &multiplicityCalo, &b_multiplicityCalo);
  fChain->SetBranchAddress("mhtDivHtCalo", &mhtDivHtCalo, &b_mhtDivHtCalo);
  fChain->SetBranchAddress("alphaTCalo", &alphaTCalo, &b_alphaTCalo);
  fChain->SetBranchAddress("jetDPhi12Calo", &jetDPhi12Calo, &b_jetDPhi12Calo);
  fChain->SetBranchAddress("jetPtsCalo", &jetPtsCalo, &b_jetPtsCalo);
  fChain->SetBranchAddress("jetPhisCalo", &jetPhisCalo, &b_jetPhisCalo);
  fChain->SetBranchAddress("jetEtasCalo", &jetEtasCalo, &b_jetEtasCalo);
  //fChain->SetBranchAddress("jetPtsAllCalo", &jetPtsAllCalo, &b_jetPtsAllCalo);
  //fChain->SetBranchAddress("jetPhisAllCalo", &jetPhisAllCalo, &b_jetPhisAllCalo);
  //fChain->SetBranchAddress("jetEtasAllCalo", &jetEtasAllCalo, &b_jetEtasAllCalo);
  fChain->SetBranchAddress("mhtPtPf", &mhtPtPf, &b_mhtPtPf);
  fChain->SetBranchAddress("mhtPhiPf", &mhtPhiPf, &b_mhtPhiPf);
  fChain->SetBranchAddress("metPtPf", &metPtPf, &b_metPtPf);
  fChain->SetBranchAddress("metPhiPf", &metPhiPf, &b_metPhiPf);
  fChain->SetBranchAddress("htPf", &htPf, &b_htPf);
  fChain->SetBranchAddress("etPf", &etPf, &b_etPf);
  fChain->SetBranchAddress("dhtPf", &dhtPf, &b_dhtPf);
  fChain->SetBranchAddress("multiplicity50Pf", &multiplicity50Pf, &b_multiplicity50Pf);
  fChain->SetBranchAddress("mhtDivHtPf", &mhtDivHtPf, &b_mhtDivHtPf);
  fChain->SetBranchAddress("alphaTPf", &alphaTPf, &b_alphaTPf);
  fChain->SetBranchAddress("jetDPhi12Pf", &jetDPhi12Pf, &b_jetDPhi12Pf);
  fChain->SetBranchAddress("jetPtsPf", &jetPtsPf, &b_jetPtsPf);
  fChain->SetBranchAddress("jetPhisPf", &jetPhisPf, &b_jetPhisPf);
  fChain->SetBranchAddress("jetEtasPf", &jetEtasPf, &b_jetEtasPf);
  //fChain->SetBranchAddress("jetPtsAllPf", &jetPtsAllPf, &b_jetPtsAllPf);
  //fChain->SetBranchAddress("jetPhisAllPf", &jetPhisAllPf, &b_jetPhisAllPf);
  //fChain->SetBranchAddress("jetEtasAllPf", &jetEtasAllPf, &b_jetEtasAllPf);
  fChain->SetBranchAddress("jetMuonsPf", &jetMuonsPf, &b_jetMuonsPf);
  fChain->SetBranchAddress("MuonNumber", &MuonNumber, &b_run);
  fChain->SetBranchAddress("nInt", &nInt, &b_nint);
  fChain->SetBranchAddress("weight", &weight, &b_weight);
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("lumi", &lumi, &b_run);
  fChain->SetBranchAddress("event", &event, &b_run);
  Notify();
}

Bool_t alphaTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void alphaTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t alphaTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef alphaTree_cxx
