#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLeaf.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

void analisi_tree_1hit(){ //faccio gli istogrammi dal Tree T creato nel file CheckESD.C
  gROOT->Reset();
  gStyle->SetOptStat(0012);
  gStyle->SetOptFit(0111);
  
  // definire istogrammi (ricordarsi di fare il write nel file successivamente)
  //TH1F *hdeltat = new TH1F("hdeltat","inside the pad (cl_{1}) - cluster along x;t_{1} - t_{2} (ps)",400,-500,500);
  //TH1F *hdeltax = new TH1F("hdeltax","inside the pad (cl_{1}) - cluster along x;#Deltax_{1} - #Deltax_{2} (cm)",100,-10,10);
  TH1F *hch = new TH1F("hch","inside the pad (cl_{1}) - cluster along x;ch_{1} - ch_{2}",500,0,500);
  //TH2F *pxt = new TH2F("pxt","inside the pad (cl_{1}) - cluster along x ;t_{1} - t_{2} (ps);cl_{1} #Deltax (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4); // 24.4 ps quantizzazione TDC
  //TH2F *pzt = new TH2F("pzt","inside the pad (cl_{1}) - cluster along z ;t_{1} - t_{2} (ps);cl_{1} #Deltaz (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4); // 24.4 ps quantizzazione TDC
  //CFC:Istogramma numero cluster(qua coincidono con le hit)
  //TH1F *hnc = new TH1F("hnc","Number of Hits",12,0.,12.);
  //Istogramma  tempi TOF
  TH1F *ht1 = new TH1F("ht1","TOF's Time ",41,0.,30000.);
  
  //Istogramma 1D per i residui
  TH1F *hresx1 = new TH1F("hresx1","Residui x1",100,-10.,10.);
  TH1F *hresz1 = new TH1F("hresz1","Residui z1",100,-10.,10.);
  TH1F *hresdist1 = new TH1F("hresdist1","Residui sqrt(x1^2+z1^2)",100,-10.,10.);
  
  //Istogramma 2D per i residui
  TH2F *h2resxz1 = new TH2F("h2resxz1" , "Residui dx1 e dz1", 100, -4. , 4. , 100 , -4. , 4.);
  
  //Istogramma tempi meno tempi attesi
  TH1F *ht1_texp = new TH1F("ht1_texp","",100,-2000.,2000.);
  TH1F *ht1_texptot = new TH1F("ht1_texptot","",100,-2000.,2000.);
  
  //TH1F *hexp_time_pi = new TH1F("hexp_time_pi","hexp_time_pi",1000,0.,30000.);
  
  TProfile *hprofx  = new TProfile("hprofx","Profile t1-t_exp_pi vs dx1",26, -2.,2.);
  TProfile *hprofxcorr  = new TProfile("hprofxcorr","Profile t1-t_exp_pi corr vs dx1",26, -2.,2.);
  
  TProfile *hprofz  = new TProfile("hprofz","Profile t1-t_exp_pi vs dz1",26, -3.,3.);
  
  TProfile *hprofd=new TProfile("hprofd","Profile t1-t_exp_pi vs d",26, 0.,4.);
  
  TH1F *htbest = new TH1F("htbest","htbest",100,-2000.,2000.);
  
  TH1F *htcorrtw= new TH1F("htcorrtw","htcorrtw",100,-2000.,2000.);
  
  //TH2F *h2dxdzt1_texp=new TH2F("h2dxdzt1_texp","h2dxdzt1_texp",30,-10.25.,10.25,50,-10.75.,10.75.);
  //TH2F *h2dxdzt1_texp_dummy=new TH2F("h2dxdzt1_texp_dummy","h2dxdzt1_texp_dummy",30,-10.25.,10.25,50,-10.75.,10.75.);
  TH2F *h2dxdzt1_texp=new TH2F("h2dxdzt1_texp","h2dxdzt1_texp",20,-1.25,1.25,10,-1.75,1.75);
  TH2F *h2dxdzt1_texp_dummy=new TH2F("h2dxdzt1_texp_dummy","h2dxdzt1_texp_dummy",20,-1.25,1.25,10,-1.75,1.75);
  
  //Istogramma 2D distanza eff vs delay time
  TH2F *h2t1_texp_deff= new TH2F("h2t1_texp_deff" , "Delay time 0 vs deff", 50, 0. , 10.,100,-400.,400.);
  TH2F *h2t1_texp_deff_tw= new TH2F("h2t1_texp_deff_tw" , "Delay time 0 corr tw vs deff", 50, 0. , 10.,100,-400.,400.);
  
  TProfile *hprofdeff=new TProfile("hprofdeff","Profile t1-t_exp_pi vs deff",50, -3.,10.);
  
  // Utilizzo TOT.
  TH2F *h2t1_texp_TOT= new TH2F("h2t1_texp_TOT" , "Delay time vs TOT", 50, 0. , 100.,100,-400.,400.);
  TProfile *hproft1_texp_TOT = new TProfile("hproft1_texp_TOT","Profile t1-t_exp_pi vs TOT",50, 0. , 100.);
  TH2F *h2t1_deff_TOT= new TH2F("h2t1_deff_TOT" , "deff vs TOT", 50, 0. , 100.,50, 0. , 10.);
  TH2F *h2t1_TOT_deff= new TH2F("h2t1_TOT_deff" , "TOT vs deff",50, 0. , 10., 50, 0. , 100.);
  
  
  TFile *f = new TFile("AnalysisResults.root");
  TTree *T = (TTree*)f->Get("T"); //in generale . (e non freccia) se Tfile è un oggetto e NON un puntatore(*)
  
  
  //Varibili tree "T"
  //Int_t nevento;
  //Int_t ntracks;
  Int_t ncluster;
  Float_t tempo[100];//con start time sottratto
  Float_t DeltaX[100];
  Float_t DeltaZ[100];
  Int_t ChannelTOF[100];
  Float_t impulso_trasv;
  Float_t exp_time_pi[100];
  Float_t L[100];
  Float_t TOT[100];
  Float_t res[3];
  Int_t charge;
  Float_t phi,eta;
  Float_t secAngle;
  Float_t cval[5];
  Float_t thetay;
  Float_t StartTime,StartTimeRes;
  Float_t z;
  
  //T->Branch("nevento",&nevento,"nevento/I");
  //T->SetBranchAddress("ntracks",&ntracks);
  T->SetBranchAddress("ncluster",&ncluster);
  T->SetBranchAddress("tempo",tempo);
  T->SetBranchAddress("DeltaX",DeltaX);
  T->SetBranchAddress("DeltaZ",DeltaZ);
  T->SetBranchAddress("ChannelTOF",ChannelTOF);
  T->SetBranchAddress("impulso_trasv",&impulso_trasv);
  T->SetBranchAddress("exp_time_pi",exp_time_pi);
  T->SetBranchAddress("L",L);
  T->SetBranchAddress("TOT",TOT);
  T->SetBranchAddress("res",res);
  T->SetBranchAddress("charge",&charge);
  T->SetBranchAddress("phi",&phi);
  T->SetBranchAddress("eta",&eta);
  T->SetBranchAddress("secPhi",&secAngle);
  T->SetBranchAddress("cval",cval);
  T->SetBranchAddress("thetay",&thetay);
  T->SetBranchAddress("StartTime",&StartTime);
  T->SetBranchAddress("StartTimeRes",&StartTimeRes);
  
  
  Int_t nentries = (Int_t)T->GetEntries();
  
  for(Int_t i=0;i<nentries;i++) {
    T->GetEntry(i);
    
    for(Int_t ip=0;ip<ncluster;ip++)
      tempo[ip] -= StartTime;
    
    if(ncluster == 1){
      if(impulso_trasv>0.8 && impulso_trasv<1.){ // serve per gli exp time
	//if( TMath::Abs(DeltaZ[0])<1.75){
	
	
	if(TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ ){
	  if(TMath::Abs(DeltaX[0])<1.5){
	    //if((ChannelTOF[0]/48)%2==0){z=1.75+DeltaZ[0];}
	    //if((ChannelTOF[0]/48)%2 == 1){z=1.75-DeltaZ[0];}
	    if((ChannelTOF[0]/48)%2==0){z=DeltaZ[0];}
	    if((ChannelTOF[0]/48)%2 == 1){z= -DeltaZ[0];}
	    Float_t w=tempo[0]- exp_time_pi[0];
	    h2dxdzt1_texp->Fill(DeltaX[0],z, w);
	    h2dxdzt1_texp_dummy->Fill(DeltaX[0],z);
	    
	    if( TMath::Abs(DeltaZ[0])<1.75 && TMath::Abs(DeltaX[0])<1.25){ //sto selezionando che il pad matchato sia dove passa la traccia, poichè altrimenti nel grafico con d e delay time vi sarebbero degli effetti di bordo
			  
	      Float_t t1=tempo[0];
	      ht1->Fill(t1);
	      ht1->GetXaxis()->SetTitle("t1(ps)");
	      
	      //cout<<(ChannelTOF[0]/48)%2<<endl;
	      
	      
	      Float_t d;
	      
	      if((ChannelTOF[0]/48)%2 == 0){
		d=sqrt(DeltaX[0]*DeltaX[0]+(1.75+DeltaZ[0])*(1.75+DeltaZ[0])); //da quello che ho capito il pad 0 è quello in alto quindi ha punto d raccolta in alto
		/*
		  if((ChannelTOF[0]/96)==((ChannelTOF[0]+48)/96)) //if fatto per capire quale sia quallo in alto e quale quello in basso poichè hanno una raccolta di carica differente
		  {
		  cout<<"ciao, io ,pad 0, sono quello sopra"<<endl;
		  }
		*/
	      }
	      if((ChannelTOF[0]/48)%2 == 1) {
		d=sqrt(DeltaX[0]*DeltaX[0]+(1.75-DeltaZ[0])*(1.75-DeltaZ[0])); //da quello che ho capito il pad 1 è quello in basso quindi ha punto d raccolta in basso
		/*//if fatto per capire quale sia quallo in alto e quale quello in basso poichè hanno una raccolta di carica differente
		  
		if((ChannelTOF[0]/96)==((ChannelTOF[0]+48)/96))
		{
		cout<<"ciao, io ,pad 1, sono quello sopra"<<endl;
		}
		*/
	      }
	      
	      h2resxz1->Fill(DeltaX[0],DeltaZ[0]);
	      h2resxz1->GetXaxis()->SetTitle("Dx1 (cm)");
	      h2resxz1->GetYaxis()->SetTitle("Dz1 (cm)");
	      
	      hresx1->Fill(DeltaX[0]);
	      hresx1->GetXaxis()->SetTitle("Dx1 (cm)");
	      
	      hresz1->Fill(DeltaZ[0]);
	      hresz1->GetXaxis()->SetTitle("Dz1 (cm)");
	      
	      hresdist1->Fill(sqrt(DeltaX[0]*DeltaX[0]+DeltaZ[0]*DeltaZ[0]));
	      hresdist1->GetXaxis()->SetTitle("sqrt(Dx^2+Dz^2) (cm)");
	      
	      
              Float_t res1 = DeltaX[0]*DeltaX[0] + DeltaZ[0]*DeltaZ[0];
	      
              //Float_t posx = (DeltaX[0])*(ChannelTOF[0]-ChannelTOF[1]);
              
	      // hdeltax->Fill(posx);
              //hdeltax->GetXaxis()->SetTitle("posx (cm)");
	      //
	      //              pxt->Fill(tempo[0], DeltaX[0]);
	      //              pxt->GetXaxis()->SetTitle("t (cm)");
	      //              pxt->GetYaxis()->SetTitle("dx (cm)");
	      
              
	      //              pzt->Fill(tempo[0], DeltaZ[0]);
	      //              pzt->GetXaxis()->SetTitle("t(cm)");
	      //              pzt->GetYaxis()->SetTitle("Dz1 (cm)");
	      
              hch->Fill(ChannelTOF[0]);
              hch->GetXaxis()->SetTitle("ch");
	      
              Float_t tw1=tempo[0]- exp_time_pi[0];
	      
              ht1_texp->Fill(tw1);
              ht1_texp->GetXaxis()->SetTitle("t_{1} - t_{exp #pi} (ps)");
	      
	      //cerco correlazione tra TOT e delay time
              h2t1_texp_TOT->Fill(TOT[0],tw1);
              h2t1_texp_TOT->GetXaxis()->SetTitle("TOT ");
              h2t1_texp_TOT->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
              hproft1_texp_TOT->Fill(TOT[0],tw1 ,1);
              hproft1_texp_TOT->GetXaxis()->SetTitle("TOT");
              hproft1_texp_TOT->GetYaxis()->SetTitle("t_{0}-t_{exp #pi} (ps)");
	      
	      /////////////////////PROFILE  1
	      
              //riempio istogrammma, poi lo fitto fuori dal loop
	      
	      //               hprofx->Fill(DeltaX[0],tw1 ,1);
	      //               hprofx->GetXaxis()->SetTitle("dx (cm)");
	      //               hprofx->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
	      //        
	      //               hprofz->Fill(DeltaZ[0],tw1 ,1);
	      //               hprofz->GetXaxis()->SetTitle("dz (cm)");
	      //               hprofz->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
	      
	      
	      hprofd->Fill(d,tw1 ,1);
	      hprofd->GetXaxis()->SetTitle("d (cm)");
	      hprofd->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
	      
	      
	      
	      //}
	    }
	  }
	}
      }
    }
  }
    
  
  h2dxdzt1_texp->Divide(h2dxdzt1_texp_dummy);
  h2dxdzt1_texp->SetOption("LEGO2");
  h2dxdzt1_texp->GetXaxis()->SetTitle("dx (cm)");
  h2dxdzt1_texp->GetYaxis()->SetTitle("dz (cm)");
  
  //Double_t par[0]=0;
  //Double_t par[1]=0;
  
  TF2 *myfunc=new TF2("myfunc"," [0] * sqrt(x*x+[1]*[1]*(1.75+y)*(1.75+y))+[2]", -1.25,1.25,-1.75,1.75);
  myfunc->SetParameter(1,0.5);
  h2dxdzt1_texp->Fit(myfunc,"R");
  Double_t vel1, alfa;
  vel1= pow(myfunc->GetParameter(0),-1);
  alfa= myfunc->GetParameter(1);
  cout<<"Ciao, io sono la velocità in 3d "<<vel1<<endl;
  
  for(Int_t i=0;i<nentries;i++){
    T->GetEntry(i);
    
    for(Int_t ip=0;ip<ncluster;ip++)
      tempo[ip] -= StartTime;

    if(ncluster == 1) {
      if(impulso_trasv>0.8 && impulso_trasv<1.){ // serve per gli exp time
	if(TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ ){
	  
	  if(TMath::Abs(DeltaX[0])<1.5){
	    if((ChannelTOF[0]/48)%2==0){z=DeltaZ[0];}
	    if((ChannelTOF[0]/48)%2==1){z= -DeltaZ[0];}
	    Float_t deff=sqrt(DeltaX[0]*DeltaX[0]+ alfa*alfa * (1.75+z)*(1.75+z)*(z>-1.75));
	    Float_t tw1tot=tempo[0]- exp_time_pi[0];
	    h2t1_texp_deff->Fill(deff,tw1tot);
	    h2t1_texp_deff->GetXaxis()->SetTitle("deff (cm)");
	    h2t1_texp_deff->GetYaxis()->SetTitle("delay time (ps)");
            
	    hprofdeff ->Fill(deff,tw1tot ,1);
	    hprofdeff->GetXaxis()->SetTitle("deff (cm)");
	    hprofdeff->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
            
	    h2t1_deff_TOT->Fill(TOT[0],deff);
	    h2t1_deff_TOT->GetXaxis()->SetTitle("TOT (ps)");
	    h2t1_deff_TOT->GetYaxis()->SetTitle("deff (cm)");
            
	    h2t1_TOT_deff->Fill(deff, TOT[0]);
	    h2t1_TOT_deff->GetXaxis()->SetTitle("deff (cm)");
	    h2t1_TOT_deff->GetYaxis()->SetTitle("TOT (ps)");
	           
	  }
	}
      }
    }
  }
  
  TF1 *fs = new TF1("fs","gaus",-400.,400.);
  h2t1_texp_deff->FitSlicesY(fs);
  TH1D *h2t1_texp_deff_1 = (TH1D *) gDirectory->FindObject("h2t1_texp_deff_1");
  h2t1_texp_deff_1->Draw("same");
  //h2t1_texp_deff_2->Draw("same");
  
  //   TF1 *f1 = new TF1("f1","pol1",0., 2.15);
  //    h2t1_texp_deff_1->Fit("f1","R");
  //    Double_t vel2;
  //    vel2 = pow(f1->GetParameter(1),-1);
  //    cout<<"Ciao, io sono la velocità in 2d "<<vel2<<endl;
  
  TF1 *f1 = new TF1("f1","[0]*TMath::Min(x,1.75)+[1]",0., 2.8); //Min(x,1.75) , l'1.75 è stato scelto guardando il fit, come meglio sembrava interploare. Per ora, per quanto ne so, è solo una coincidenza che coincida con pad z
  h2t1_texp_deff_1->Fit(f1,"R");
  Double_t p0,vel2, offset_tw;
  p0=f1->GetParameter(0);
  vel2= pow(p0,-1);
  offset_tw=f1->GetParameter(1);
  cout<<"Ciao, io sono la velocità in 2d "<<vel2<<endl;
  
  for(Int_t i=0;i<nentries;i++){
    T->GetEntry(i);
    
    for(Int_t ip=0;ip<ncluster;ip++)
      tempo[ip] -= StartTime;
    
    if(ncluster == 1){
      if(impulso_trasv>0.8 && impulso_trasv<1.){ // serve per gli exp time
	if(TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ ){                   
	  if(TMath::Abs(DeltaX[0])<1.5){
	    if((ChannelTOF[0]/48)%2==0){z=DeltaZ[0];}
	    if((ChannelTOF[0]/48)%2==1){z= -DeltaZ[0];}
            
	    Float_t deff=sqrt(DeltaX[0]*DeltaX[0]+ alfa*alfa * (1.75+z)*(1.75+z)*(z>-1.75));
            
	    if(deff<2.8){
	      Float_t tw1tot=tempo[0]- exp_time_pi[0];
	      ht1_texptot->Fill(tw1tot);
	      ht1_texptot->GetXaxis()->SetTitle("t_{1} - t_{exp #pi} (ps)");
	      
	      Float_t tw1corrtw=tw1tot-(offset_tw + p0 *TMath::Min(deff,Float_t(1.75)));//Min(x,1.75) , l'1.75 è stato scelto guardando il fit, come meglio sembrava interploare. Per ora, per quanto ne so, è solo una coincidenza che coincida con pad z
	      htcorrtw->Fill(tw1corrtw);
	      htcorrtw->GetXaxis()->SetTitle("t_corr_tw (ps)");
              
	      h2t1_texp_deff_tw->Fill(deff,tw1corrtw);
	      h2t1_texp_deff_tw->GetXaxis()->SetTitle("deff (cm)");
	      h2t1_texp_deff_tw->GetYaxis()->SetTitle("delay time corr_tw (ps)");
	      
	    }
	  }
	}
      }
    }
  }

    
  TF1 *fstw = new TF1("fstw","gaus",-400.,400.);
  h2t1_texp_deff_tw->FitSlicesY(fstw);
  TH1D *h2t1_texp_deff_tw_1 = (TH1D *) gDirectory->FindObject("h2t1_texp_deff_tw_1");
  TH1D *h2t1_texp_deff_tw_2 = (TH1D *) gDirectory->FindObject("h2t1_texp_deff_tw_2");
  
  h2t1_texp_deff_tw_1->Draw("same");
  h2t1_texp_deff_tw_2->Draw("same");
  
  //    TF1 *f1 = new TF1("f1","pol1",-1.25,0.5);
  //    //hprofx->GetYaxis()->SetRangeUser(-40.,100.);
  //    hprofx->Fit("f1","R");
  //    Double_t offset_p1,x1_p1;
  //    offset_p1= f1->GetParameter(0);
  //    x1_p1= f1->GetParameter(1);
  
  
  //    for(Int_t i=0;i<nentries;i++)
  //    {
  //        T->GetEntry(i);
  //    for(Int_t ip=0;ip<ncluster;ip++)
  //      tempo[ip] -= StartTime;
  //        if(ncluster == 1)
  //        {
  //        if(impulso_trasv>0.8 && impulso_trasv<1.2) // serve per gli exp time
  //        {
  //                //if( TMath::Abs(DeltaZ[0])<1.75)
  //                //{
  //        if(TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigmca che sia un pi*/ )
  //        {
  //
  //            //Float_t posx = (DeltaX[0])* dch;
  //            Float_t tw1=tempo[0]- exp_time_pi[0];
  //            
  //            Float_t tw1corr=tw1-(offset_p1 + x1_p1 *DeltaX[0]);
  //            
  //            hprofxcorr->Fill(DeltaX[0],tw1corr,1);
  //            hprofxcorr->GetXaxis()->SetTitle("Dx1 (cm)");
  //            hprofxcorr->GetYaxis()->SetTitle("t1_corr=t1-t_exp_pi corr(ps)");
  //            
  //            htbest->Fill(tw1corr);
  //            htbest->GetXaxis()->SetTitle("t_best=t1_corr(ps)");
  //        
  //        //}
  //        }
  //        }
  //        }
  //    }
  //    
  //    TF1 *f1c = new TF1("f1c","pol1",-1.25,0.5);
  //    hprofxcorr->Fit("f1c","R");
  
  TFile *fo2 = new TFile("output_ist_tree_1hit.root","RECREATE");
  
  hch->Write();
  ht1->Write();
  h2resxz1->Write();
  hresx1->Write();
  hresz1->Write();
  hresdist1->Write();
  //hdeltax->Write();
  //pxt->Write();
  //pzt->Write();
  //hdeltach->Write();
  ht1_texp->Write();
  ht1_texptot->Write();
  //hprofx->Write();
  //hprofz->Write();
  //hprofxcorr->Write();
  //htbest->Write();
  hprofd->Write();
  h2dxdzt1_texp->Write();
  //h2dxdzt1_texp_dummy->Write();
  h2t1_texp_deff->Write();
  h2t1_texp_deff_1->Write();
  //h2t1_texp_deff_2->Write();
  hprofdeff->Write();
  
  htcorrtw->Write();
  
  h2t1_texp_deff_tw->Write();
  h2t1_texp_deff_tw_1->Write();
  h2t1_texp_deff_tw_2->Write();
  h2t1_texp_TOT->Write();
  hproft1_texp_TOT->Write();
  h2t1_deff_TOT->Write();
  h2t1_TOT_deff->Write();
  fo2->Close();
  
  system("say Ehi you, I have done");
}
