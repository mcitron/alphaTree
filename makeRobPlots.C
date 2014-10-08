makeRobPlots()
{
  //TFile * fIn = new TFile("140929/turnonPlotsT2TT.root");
  TFile * fIn = new TFile("turnonPlots.root");
  
  std::vector<TString> njets;
  // njets.push_back("1");
  njets.push_back("2");
  njets.push_back("3");
  njets.push_back("ge4");
  std::vector<TString> hard;
  //hard.push_back("Hard");
  hard.push_back("noHard");
  std::vector<TString> bins;
  bins.push_back("Bin1");
  bins.push_back("Bin2");

  for (std::vector<TString>::iterator iBin = bins.begin(); iBin != bins.end(); iBin++)
  {
    for (std::vector<TString>::iterator iHard = hard.begin(); iHard != hard.end(); iHard++)
    {
      for (std::vector<TString>::iterator iNum = njets.begin(); iNum !=njets.end(); iNum++)
      {
	//TGraphAsymmErrors * ours = fIn->Get(*iHard+""+*iNum+"Jet/"+*iHard+"mhtDivHt130Full"+*iNum+"Jet_Eff");
	TGraphAsymmErrors * oursMet = fIn->Get(*iHard+""+*iNum+"Jet"+*iBin+"/"+*iHard+"metMhtDivHt130Full"+*iNum+"Jet_CumuEff");
	TGraphAsymmErrors * orig = fIn->Get(*iHard+""+*iNum+"Jet"+*iBin+"/"+*iHard+"met70Full"+*iNum+"Jet_CumuEff");
	TGraphAsymmErrors * mht = fIn->Get(*iHard+""+*iNum+"Jet"+*iBin+"/"+*iHard+"metMht60Full"+*iNum+"Jet_CumuEff");
	TGraphAsymmErrors * met = fIn->Get(*iHard+""+*iNum+"Jet"+*iBin+"/"+*iHard+"met60Full"+*iNum+"Jet_CumuEff");

	TCanvas * c = new TCanvas("Turn on Curve (PF MHT) "+*iNum+*iHard+*iBin,"",600,600);
	c->cd();
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle(*iNum+" Jet"+" "+*iHard+" "+*iBin);
	//ours->SetLineColor(1);
	oursMet->SetLineColor(4);
	mht->SetLineColor(45);
	orig->SetLineColor(6);
	met->SetLineColor(2);

	//mg->Add(ours,"lp");
	mg->Add(oursMet,"lp");
	mg->Add(mht,"lp");
	mg->Add(orig,"lp");
	mg->Add(met,"lp");

	mg->Draw("a");
	mg->GetXaxis()->SetRangeUser(0,500);
	mg->GetXaxis()->SetTitle("MHT PF (GeV)");
	//TLegend* leg = new TLegend(0.5,0.1,0.6,0.4);
	TLegend* leg = new TLegend(0.5,0.1,0.6,0.4);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.02);
	//leg->AddEntry(ours,"mHTOvHT0.32 HT125 (full menu)","lep");
	//leg->AddEntry(dom,"dPhi8 Jet2Pt30 Ht125 (full menu)","lep");
	//leg->AddEntry(dom,"Jet1Pt120 Jet2Pt30 Ht125 (full menu)","lep");
	leg->AddEntry(orig,"MET70 || HT175","lep");
	leg->AddEntry(oursMet,"MET70 || mHTOvHT0.32 HT125 || HT 200","lep");
	leg->AddEntry(mht,"MET70 || mht 60 HT100 || HT 200","lep");
	leg->AddEntry(met,"MET60 || HT200","lep");
	leg->Draw("L");
	c->SaveAs("talkplots/TT25TurnOn"+*iNum+*iHard+*iBin+"JetCumu.png");
      }
    }
  }
}
