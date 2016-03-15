{
  TFile *f = new TFile("AnalysisResults.root");
  TTree *t = f->Get("T");

  Int_t n=t->GetEntries();

  TH1F *hinttime = new TH1F("hinttime",";interaction time (ps)",200,-500,500);
  TH1F *htofstarttime = new TH1F("htofstarttime",";interaction time (ps)",200,-500,500);
  TH1F *hpro_hr = new TH1F("hpro_hr","T0-TOF res < 40 ps;T0-TOF (ps)",200,-500,500);

  TProfile *hpro = new TProfile("hpro",";interaction time (ps);T0-TOF (ps)",500,-1000,1000);
  TProfile *hpro2 = new TProfile("hpro2",";interaction time (ps);eff T0-TOF",500,-1000,1000);

  htofstarttime->SetLineColor(2);
  hpro_hr->SetLineColor(4);

  Float_t inttimeprev = 1000000;
  Float_t inttime,tofstarttime,tofstarttimeres;

  for(Int_t i=0;i<n;i++){
    t->GetEvent(i);

    inttime = t->GetLeaf("InteractionTime")->GetValue();
    tofstarttime = t->GetLeaf("StartTime")->GetValue();
    tofstarttimeres = t->GetLeaf("StartTimeRes")->GetValue();

    if(inttime != inttimeprev){
      hpro2->Fill(inttime,tofstarttime!=0);
      hinttime->Fill(inttime);
      htofstarttime->Fill(tofstarttime);
      if(tofstarttime != 0) hpro->Fill(inttime,tofstarttime);
      if(tofstarttimeres < 40)  hpro_hr->Fill(tofstarttime);
      inttimeprev = inttime;
    }
  }




  hinttime->Draw();
  htofstarttime->Draw("SAME");
  hpro_hr->Draw("SAME");
  new TCanvas();
  hpro->Draw();
  new TCanvas();
  hpro2->Draw();
}
