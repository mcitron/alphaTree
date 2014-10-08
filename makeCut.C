void makeCut()
{
  TFile * rateFile=new TFile("./ratePlots.root");
  TFile * effFile=new TFile("./effPlots.root");
  TFile * output = new TFile("cutOut.root","recreate");

  std::vector<TString> jetMult;
  //jetMult.push_back("1Jet");
  jetMult.push_back("2Jet");
  jetMult.push_back("3Jet");
  jetMult.push_back("ge4Jet");
  //TH2F * effPlot=effFile->Get("mhtDivHtUct_htUct/mhtDivHtUct_htUct_htPf > 200 && htPf < 275 && alphaTPf > 0.65 && jetPtsPf[1]>100_Eff");
  std::vector<TString> hardLead;
  hardLead.push_back("Hard");
  hardLead.push_back("noHard");
  std::vector<TString> cutType;
  cutType.push_back("mhtDivHtFull");
  cutType.push_back("mhtFull");
  cutType.push_back("mhtDivHtFullMet");
  cutType.push_back("mhtFullMet");
  // cutType.push_back("mhtDivHtFull");
  // cutType.push_back("mhtFull");
  //cutType.push_back("dphi");
  std::vector<TString> rateType;
  rateType.push_back("fullMhtDivHt");
  rateType.push_back("fullMht");
  rateType.push_back("fullMETMhtDivHt");
  rateType.push_back("fullMETMht");
  // rateType.push_back("fullMhtDivHt");
  //rateType.push_back("fullMht");
  //rateType.push_back("dphi");
  //
  std::vector<int> rateMin;
  rateMin.push_back(14000);
  rateMin.push_back(14000);
  rateMin.push_back(18000);
  rateMin.push_back(18000);
  std::vector<int> rateMax;
  rateMax.push_back(16000);
  rateMax.push_back(16000);
  rateMax.push_back(20000);
  rateMax.push_back(20000);

  output->cd();
  int i = 0;
  for(std::vector<TString>::iterator cut = cutType.begin(); cut != cutType.end(); cut++)
  {
    TDirectory * cutDir = output->mkdir(*cut);
    cutDir->cd();
    TH2F * ratePlot=rateFile->Get(rateType.at(i));
    for(std::vector<TString>::iterator mult = jetMult.begin(); mult != jetMult.end(); mult++)
    {
      TDirectory * numDir = cutDir->mkdir(*mult);
      numDir->cd();
      for(std::vector<TString>::iterator hard = hardLead.begin(); hard != hardLead.end(); hard++)
      {
	std::vector<TString> pfDivString;
	pfDivString.push_back(*mult+*hard+"/"+*hard+"Bin1"+*mult+*cut);
	pfDivString.push_back(*mult+*hard+"/"+*hard+"Bin2"+*mult+*cut);
	pfDivString.push_back(*mult+*hard+"/"+*hard+"Bin3"+*mult+*cut);
	for (std::vector<TString>::iterator eff = pfDivString.begin(); eff != pfDivString.end(); eff++)
	{

	  std::cout << *eff << std::endl;
	  TH2F * effPlot=effFile->Get(*eff);
	  effPlot->SetMaximum(1);
	  //effPlot->Write();
	  TH2F * effPlotCut = makeCut(effPlot,ratePlot,rateMin.at(i),rateMax.at(i));
	  //effPlotCut->SetTitle("cut");
	  TString min; min.Form("%d",rateMin.at(i));
	  TString max; max.Form("%d",rateMax.at(i));
          TString name = min + " to " + max + " Hz";
	  effPlotCut->SetMaximum(1.0);
	  effPlotCut->SetTitle(name);
	  effPlotCut->Write();

	}
      }

    }
    i++;
  }

  TDirectory * numDir = cutDir->mkdir("1Jet");
  numDir->cd();
  i = 0;
  /*
     for(std::vector<TString>::iterator cut = cutType.begin(); cut != cutType.end(); cut++)
     {
     TH2F * ratePlot=(TH2F*)rateFile->Get(rateType.at(i));
     i++;
     std::cout << "1Jet/Bin11Jet"+*cut << std::endl;
     TH2F * effPlot=(TH2F*)effFile->Get("1Jet/Bin11Jet"+*cut);
     effPlot->SetMaximum(1);
     TH2F * effPlotCut = makeCut(effPlot,ratePlot,13000,17000);
     effPlotCut->SetMaximum(1.0);
     effPlotCut->Write();
     }
     */
  output->Close();
}
TH2F * makeCut(TH2F * eff, TH2F * rate,double lowCut,double highCut){

  TH2F * output = new TH2F(*eff);

  int nXbins = eff->GetNbinsX();
  int nYbins = eff->GetNbinsY();
  for (int xbins = 0; xbins <= nXbins+1; xbins++)
  {
    for (int ybins = 0; ybins <= nYbins+1; ybins++)
    {
      output->SetBinContent(xbins,ybins,0.);
      if(rate->GetBinContent(xbins,ybins) > lowCut && rate->GetBinContent(xbins,ybins) < highCut){
	output->SetBinContent(xbins,ybins,eff->GetBinContent(xbins,ybins));
      }
      //output->SetBinContent(xbins,ybins,ZB_XSECTION*(1-(double)dummy/norm));
    }
  }

  return output;
}
