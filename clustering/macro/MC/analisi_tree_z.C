{ //faccio gli istogrammi dal Tree T creato nel file CheckESD.C
    gROOT->Reset();
    gStyle->SetOptStat(0012);
    gStyle->SetOptFit(0111);
    
// definire istogrammi (ricordarsi di fare il write nel file successivamente)
//TH1F *hdeltat = new TH1F("hdeltat","inside the pad (cl_{1}) - cluster along x;t_{1} - t_{2} (ps)",400,-500,500);
TH1F *hdeltaz = new TH1F("hdeltaz","inside the pad (cl_{1}) - cluster along z;#Deltaz_{1} - #Deltaz_{2} (cm)",100,-10,10);
TH1F *hdeltach = new TH1F("hdeltach","inside the pad (cl_{1}) - cluster along x;ch_{1} - ch_{2}",100,-50,50);
TH2F *hdeltaxz = new TH2F("hdeltaxz","Cluster relative distance - cluster along x;#Deltax_{1} - #Deltax_{2} (cm);#Deltaz_{1} - #Deltaz_{2} (cm)",1000,-10,10,1000,-10,10);
TH2F *pxdeltat = new TH2F("pxdeltat","inside the pad (cl_{1}) - cluster along x ;t_{1} - t_{2} (ps);cl_{1} #Deltax (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4); // 24.4 ps quantizzazione TDC
TH2F *pxdeltat2 = new TH2F("pxdeltat2","inside the pad (cl_{1}) - cluster along x;t_{1} - t_{2} (ps);cl_{2} #Deltax (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4);
TH2F *pzdeltat = new TH2F("pzdeltat","inside the pad (cl_{1}) - cluster along z ;t_{1} - t_{2} (ps);cl_{1} #Deltaz (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4); // 24.4 ps quantizzazione TDC
TH2F *pzdeltat2 = new TH2F("pzdeltat2","inside the pad (cl_{1}) - cluster along z;t_{1} - t_{2} (ps);cl_{2} #Deltaz (cm)",41,-20.5*24.4,20.5*24.4,100,-4,4);

//CFC:Istogramma numero cluster(qua coincidono con le hit)
TH1F *hnc = new TH1F("hnc","Number of Hits",12,0.,12.);
//Istogramma differenza tempi TOF
TH1F *hdiff = new TH1F("hdiff","Time difference between hits",41,-510.,510.);

//Istogramma 1D per i residui
TH1F *hresx1 = new TH1F("hresx1","Residui x1",100,-10.,10.);
TH1F *hresx2 = new TH1F("hresx2","Residui x2",100,-10.,10.);
TH1F *hresz1 = new TH1F("hresz1","Residui z1",100,-10.,10.);
TH1F *hresz2 = new TH1F("hresz2","Residui z2",100,-10.,10.);
TH1F *hresdist1 = new TH1F("hresdist1","Residui sqrt(x1^2+z1^2)",100,-10.,10.);
    
//Istogramma 1D per i residui con centroide modificato(a metà tra i due pad)
TH1F *hresz1ricen = new TH1F("hresz1ricen","Residui z1 ricentrato",100,-10.,10.);
TH1F *hresdist1ricen = new TH1F("hresdist1ricen","Residui  sqrt(x1^2+z1ricentrato^2)",100,-10.,10.);
TH1F *hresz2ricen = new TH1F("hresz2ricen","Residui z2 ricentrato",100,-10.,10.);
TH1F *hresdist2ricen = new TH1F("hresdist2ricen","Residui  sqrt(x2^2+z2ricentrato^2)",100,-10.,10.);


//Istogramma 2D per i residui
TH2F *h2resx = new TH2F("h2resx" , "Residui dx1 e dx2", 100, -10. , 10. , 100 , -10. , 10.);
TH2F *h2resz = new TH2F("h2resz" , "Residui dz1 e dz2", 100, -10. , 10. , 100 , -10. , 10.);
TH2F *h2resxz1 = new TH2F("h2resxz1" , "Residui dx1 e dz1", 100, -4. , 4. , 100 , -4. , 4.);
TH2F *h2resxz2 = new TH2F("h2resxz2" , "Residui dx2 e dz2", 100, -4. , 4. , 100 , -4. , 4.);

    
//Istogramma tempi meno tempi attesi
TH1F *ht1_texp = new TH1F("ht1_texp","ht1_texp",100,-2000.,2000.);
TH1F *ht2_texp = new TH1F("ht2_texp","ht2_texp",100,-2000.,2000.);
    
TH2F *h2tw1tw2 = new TH2F("h2tw1tw2" , "h2tw1tw2", 50, -500. , 500. , 50 , -500. , 500.);


//TH1F *hexp_time_pi = new TH1F("hexp_time_pi","hexp_time_pi",1000,0.,30000.);
    
hprof1  = new TProfile("hprof1","Profile t1-t_exp_pi vs dx1",26, -3.,3.);
hprof1corr  = new TProfile("hprof1corr","Profile t1-t_exp_pi corr vs dx1",26, -3.,3.);
hprof2  = new TProfile("hprof2","Profile t2-t_exp_pi vs dx2",26,-3.,3.);
hprof2corr  = new TProfile("hprof2corr","Profile t2-t_exp_pi corr vs dx2",26, -3.,3..);
    
TH1F *htbest = new TH1F("htbest","htbest",100,-2000.,2000.);
    
TH2F *h2bestvst1corr_t2corr = new TH2F("h2bestvst1corr_t2corr" , "Correzione residua", 100, -500. , 500. , 100 , -1000. , 1000.);
hprofbest  = new TProfile("hprofbest","Profile tbest vs t1corr-t2corr",100, -500. , 500.);
hprofbestcorr  = new TProfile("hprofbestcorr","Profile hprofbestcorr vs t1corr-t2corr",100, -500. , 500.);
TH2F *h2bescorrtvst1corr_t2corr = new TH2F("h2bescorrtvst1corr_t2corr" , "Correzione residua", 100, -500. , 500. , 100 , -1000. , 1000.);
    
TH1F *htbestcorr = new TH1F("htbestcorr","htbestcorr",100,-2000.,2000.);

    TFile *f = new TFile("output.root");
    TTree *T = (TTree*)f->Get("T"); //in generale . (e non freccia) se Tfile è un oggetto e NON un puntatore(*)


    //Varibili tree "T"
    //Int_t nevento;
    //Int_t ntracks;
    Int_t ncluster;
    Float_t tempo[100];//con start time sottratto
    Float_t DeltaX[100];
    Float_t DeltaZ[100];
    Int_t ChannelTOF[100];
    Float_t exp_time_el[100];
    Float_t exp_time_mu[100];
    Float_t exp_time_pi[100];
    Float_t exp_time_ka[100];
    Float_t exp_time_pr[100],rTOFused;
    Float_t L[100];
    Float_t impulso_trasv;
    Float_t res[3];
    Int_t charge;
    Float_t phi,eta;
    Float_t secAngle;
    Float_t cval[5],gtime;
    Float_t dxt,dzt,xgl,ygl,zgl;
    Int_t mism;
    Float_t dedx,StartTime,StartTimeRes;
    Float_t interactiontime;
    
    //T->Branch("nevento",&nevento,"nevento/I");
    //T->SetBranchAddress("ntracks",&ntracks);
    T->SetBranchAddress("ncluster",&ncluster);
    T->SetBranchAddress("tempo",tempo);
    T->SetBranchAddress("gtime",&gtime);
    T->SetBranchAddress("DXtrue",&dxt);
    T->SetBranchAddress("DZtrue",&dzt);
    T->SetBranchAddress("DeltaX",DeltaX);
    T->SetBranchAddress("DeltaZ",DeltaZ);
    T->SetBranchAddress("exp_time_el",exp_time_el);
    T->SetBranchAddress("exp_time_mu",exp_time_mu);
    T->SetBranchAddress("exp_time_pi",exp_time_pi);
    T->SetBranchAddress("exp_time_ka",exp_time_ka);
    T->SetBranchAddress("exp_time_pr",exp_time_pr);
    T->SetBranchAddress("L",L);
    T->SetBranchAddress("ChannelTOF",ChannelTOF);
    T->SetBranchAddress("impulso_trasv",&impulso_trasv);
    T->SetBranchAddress("res",res);
    T->SetBranchAddress("charge",&charge);
    T->SetBranchAddress("phi",&phi);
    T->SetBranchAddress("eta",&eta);
    T->SetBranchAddress("secPhi",&secAngle);
    T->SetBranchAddress("cval",cval);
    T->SetBranchAddress("mism",&mism);
    T->SetBranchAddress("dedx",&dedx);
    T->SetBranchAddress("StartTime",&StartTime);
    T->SetBranchAddress("StartTimeRes",&StartTimeRes);
    T->SetBranchAddress("rTOF",&rTOFused);
    T->SetBranchAddress("InteractionTime",&interactiontime);

    
    Int_t nentries = (Int_t)T->GetEntries();

for(Int_t i=0;i<nentries;i++)
{
    T->GetEntry(i);
    

    if(ncluster == 2)
      {
          if(impulso_trasv>0.8 && impulso_trasv<1.2) // serve per gli exp time
          {

          //if(exp_time_pi[0] > 0. && exp_time_pi[1] > 0.)
        
        
          //Int_t dch = TMath::Abs(ChannelTOF[0]-ChannelTOF[1]);
          Int_t dch = ChannelTOF[0]-ChannelTOF[1];
          if((ChannelTOF[0]/96)==(ChannelTOF[1]/96)) // così sono nella stessa strip
          {
          if( TMath::Abs(DeltaX[0])<1.25)
          {
              
          
          if(TMath::Abs(dch) == 48/*&& TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.)// poi dovrei farlo anche per 1
          {
              // hexp_time_pi->Fill(exp_time_pi[0]);
              
              Float_t diff=tempo[0]-tempo[1];
              hdiff->Fill(diff);
              hdiff->GetXaxis()->SetTitle("t1-t2 (ps)");
              
              
              h2resx->Fill(DeltaX[0],DeltaX[1]);
              h2resx->GetXaxis()->SetTitle("Dx1 (cm)");
              h2resx->GetYaxis()->SetTitle("Dx2 (cm)");


              
              h2resxz1->Fill(DeltaX[0],DeltaZ[0]);
              h2resxz1->GetXaxis()->SetTitle("Dx1 (cm)");
              h2resxz1->GetYaxis()->SetTitle("Dz1 (cm)");
              
              h2resxz2->Fill(DeltaX[1],DeltaZ[1]);
              h2resxz2->GetXaxis()->SetTitle("Dx2 (cm)");
              h2resxz2->GetYaxis()->SetTitle("Dz2 (cm)");
              
              hresx1->Fill(DeltaX[0]);
              hresx1->GetXaxis()->SetTitle("Dx1 (cm)");

              hresx2->Fill(DeltaX[1]);
              hresx2->GetXaxis()->SetTitle("Dx2 (cm)");

              hresz1->Fill(DeltaZ[0]);
              hresz1->GetXaxis()->SetTitle("Dz1 (cm)");

              hresz2->Fill(DeltaZ[1]);
              hresz2->GetXaxis()->SetTitle("Dz2 (cm)");

              hresdist1->Fill(sqrt(DeltaX[0]*DeltaX[0]+DeltaZ[0]*DeltaZ[0]));
              hresdist1->GetXaxis()->SetTitle("sqrt(Dx1^2+Dz1^2) (cm)");
              
              h2resz->Fill(DeltaZ[0],DeltaZ[1]);
              h2resz->GetXaxis()->SetTitle("Dz1 (cm)");
              h2resz->GetYaxis()->SetTitle("Dz2 (cm)");
              
              hdeltaxz->Fill(DeltaX[0] - DeltaX[1],DeltaZ[0] - DeltaZ[1]);
              hdeltaxz->GetXaxis()->SetTitle("Dx1-Dx2 (cm)");
              hdeltaxz->GetYaxis()->SetTitle("Dz1-Dz2 (cm)");
              
              Float_t Dz1ricen;
              if((dch)== +48){Dz1ricen=DeltaZ[0]+1.75;}
              if((dch)== -48){Dz1ricen=DeltaZ[0]-1.75;}
              hresz1ricen->Fill(Dz1ricen);
              hresz1ricen->GetXaxis()->SetTitle("Dz1 ricentrato (cm)");
              hresdist1ricen->Fill(sqrt(DeltaX[0]*DeltaX[0]+Dz1ricen[0]*Dz1ricen[0]));
              hresdist1ricen->GetXaxis()->SetTitle("sqrt(Dx1^2+Dz1ricentrato^2) (cm)");
              
              Float_t Dz2ricen;
              if((dch)== +48){Dz2ricen=DeltaZ[1]-1.75;}
              if((dch)== -48){Dz2ricen=DeltaZ[1]+1.75;}
              hresz2ricen->Fill(Dz2ricen);
              hresz2ricen->GetXaxis()->SetTitle("Dz1 ricentrato (cm)");
              hresdist2ricen->Fill(sqrt(DeltaX[1]*DeltaX[1]+Dz1ricen*Dz1ricen));
              hresdist2ricen->GetXaxis()->SetTitle("sqrt(Dx2^2+Dz2ricentrato^2) (cm)");
              
              Float_t res1 = DeltaX[0]*DeltaX[0] + DeltaZ[0]*DeltaZ[0];
              Float_t res2 = DeltaX[1]*DeltaX[1] + DeltaZ[1]*DeltaZ[1];

              Float_t posz = (DeltaZ[0])*(ChannelTOF[0]-ChannelTOF[1])/48;
              Float_t  posz2 = (-3.5*(ChannelTOF[0]-ChannelTOF[1])/48 + (DeltaZ[1]))*(ChannelTOF[0]-ChannelTOF[1])/48;
              
              }
                  
              hdeltaz->Fill(posz);
              hdeltaz->GetXaxis()->SetTitle("posz (cm)");

              pxdeltat->Fill(tempo[0] - tempo[1], DeltaX[0]);
              pxdeltat->GetXaxis()->SetTitle("t1-t2 (cm)");
              pxdeltat->GetYaxis()->SetTitle(" Dx1(cm)");
              
              pxdeltat2->Fill(tempo[0]- tempo[1], DeltaX[1]);
              pxdeltat2->GetXaxis()->SetTitle("t1-t2 (cm)");
              pxdeltat2->GetYaxis()->SetTitle(" Dx2(cm)");

              
              pzdeltat->Fill(tempo[0] - tempo[1], posz);
              pzdeltat->GetXaxis()->SetTitle("t1-t2 (cm)");
              pzdeltat->GetYaxis()->SetTitle("posz1 (cm)");
              
              pzdeltat2->Fill(tempo[0]- tempo[1], posz);
              pzdeltat2->GetXaxis()->SetTitle("t1-t2 (cm)");
              pzdeltat2->GetYaxis()->SetTitle("posz  (cm)");

              

              hdeltach->Fill(ChannelTOF[0]-ChannelTOF[1]);
              hdeltach->GetXaxis()->SetTitle("Dch");

              
             
              //if(TMath::Abs(DeltaZ[0]-DeltaZ[1])<=0.5) // perchè così dovrebbero provenire(credo) da stesso cluster in quanto, se così è , essendo loro adiacenti sulle x, dovrebbero avere circe dtesso residuo in z. NON sono sicura sia necessario fare ciò. // non dovrebbe più servire perchè ho modificato Check ESD spostando i residui dentro l'if
                  
              //h2resz->Fill(DeltaZ[0],DeltaZ[1]);
              Float_t tw1=tempo[0]- gtime[0];
              Float_t tw2=tempo[1]- gtime[1]/*/L[0]*L[1]*/; //non metto exp_time_pi[1] poichè momentaneamente non va
             
           

               /////////////////////PROFILE  1
               
               
               
              
              //riempio istogrammma, poi lo fitto fuori dal loop
               hprof1->Fill(posz,tw1 ,1);
               hprof1->GetXaxis()->SetTitle("posz (cm)");
               hprof1->GetYaxis()->SetTitle("t1-t_gtime (ps)");
              
              
               hprof2->Fill(posz,tw2 ,1);
               hprof2->GetXaxis()->SetTitle("posz (cm)");
               hprof2->GetYaxis()->SetTitle("t2-t_geant (ps)");
              
              

          }
          
          }
          }

    }

}
    TF1 *f1 = new TF1("f1","pol1",-1.75,0.5);
    
    //hprof1->GetYaxis()->SetRangeUser(-40.,100.);
    hprof1->Fit("f1","R");
    
    
    Double_t offset_p1,x1_p1;
    offset_p1= f1->GetParameter(0);
    x1_p1= f1->GetParameter(1);
    
    TF1 *f2 = new TF1("f2","pol1",-1.75,0.5);
    
   // hprof2->GetYaxis()->SetRangeUser(-40.,150.);
    hprof2->Fit("f2","R");
    
    
    Double_t offset_p2,x2_p2;
    offset_p2= f2->GetParameter(0);
    x2_p2= f2->GetParameter(1);
    
    
    for (Int_t nu=0 ; nu<nentries ; nu++)
    {
        T->GetEntry(nu);
        
        if(ncluster == 2)
        {
        if(impulso_trasv>0.8 && impulso_trasv<1.2) // serve per gli exp time
        {
        if((ChannelTOF[0]/96)==(ChannelTOF[1]/96))
        {
        if( TMath::Abs(DeltaX[0])<1.25)
        {
            
        Int_t dch = ChannelTOF[0]-ChannelTOF[1];
            
        if(TMath::Abs(dch) == 48 /*&& TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.)
        {
            
            Float_t posz = (DeltaZ[0])*(ChannelTOF[0]-ChannelTOF[1])/48;
            Float_t  posz2 = (-3.5*(ChannelTOF[0]-ChannelTOF[1])/48 + (DeltaZ[1]))*(ChannelTOF[0]-ChannelTOF[1])/48;
            
            Float_t tw1=tempo[0]- gtime[0];
            Float_t tw2=tempo[1]- gtime[1];
            
            Float_t tw1corr=tw1-(offset_p1 + x1_p1 *posz);
            Float_t tw2corr=tw2-(offset_p2 + x2_p2 *posz);
            
            hprof1corr->Fill(posz,tw1corr,1);
            hprof1corr->GetXaxis()->SetTitle("Dz1 (cm)");
            hprof1corr->GetYaxis()->SetTitle("t1_corr=t1-t_geant corr(ps)");
            
            hprof2corr->Fill(posz2,tw2corr,1);
            hprof2corr->GetXaxis()->SetTitle("Dz2 (cm)");
            hprof2corr->GetYaxis()->SetTitle("t2_corr=t2-t_geant corr(ps)");
            
            Float_t tbest=(tw1corr+tw2corr)/2.;
            htbest->Fill(tbest);
            htbest->GetXaxis()->SetTitle("t_best=(t1_corr +t2_corr)/2 (ps)");
            
            
            h2bestvst1corr_t2corr->Fill((tw1corr-tw2corr), tbest);
            h2bestvst1corr_t2corr->GetXaxis()->SetTitle("t1_corr - t2_corr (ps)");
            h2bestvst1corr_t2corr->GetYaxis()->SetTitle("t_best(ps)");
            

            
            ht1_texp->Fill(tw1);
            ht1_texp->GetXaxis()->SetTitle("t1_corr (ps)");

            ht2_texp->Fill(tw2);
            ht2_texp->GetXaxis()->SetTitle("t2_corr (ps)");

            h2tw1tw2->Fill(tw1, tw2);
            h2tw1tw2->GetXaxis()->SetTitle("t1_corr (ps)");
            h2tw1tw2->GetYaxis()->SetTitle("t2_corr (ps)");

            
            hprofbest->Fill((tw1corr-tw2corr),tbest ,1);
            hprofbest->GetXaxis()->SetTitle("t1corr-t2corr (ps)");
            hprofbest->GetYaxis()->SetTitle("tbest (ps)");

        }
        }
        }
        }
        }
    }
    
    TF1 *f1c = new TF1("f1c","pol1",-1.25,0.5);
    hprof1corr->Fit("f1c","R");
    TF1 *f2c = new TF1("f2c","pol1",-1.25,0.5);
    hprof2corr->Fit("f2c","R");
    
    TF1 *fb = new TF1("fb","pol2",-500. , 500.);
    // hprofbest->GetYaxis()->SetRangeUser(-40.,100.);
    hprofbest->Fit("fb","R");

    Double_t offset_pb,x1_pb,x2_pb;
    offset_pb= fb->GetParameter(0);
    x1_pb= fb->GetParameter(1);
    x2_pb= fb->GetParameter(2);
   
    
    for (Int_t nu=0 ; nu<nentries ; nu++)
    {
        T->GetEntry(nu);
        
        if(ncluster == 2)
        {
        if(impulso_trasv>0.8 && impulso_trasv<1.2) // serve per gli exp time
        {
        if((ChannelTOF[0]/96)==(ChannelTOF[1]/96))
        {
        if( TMath::Abs(DeltaX[0])<1.25)
        {
                    
            Int_t dch = ChannelTOF[0]-ChannelTOF[1];
                    
            if(TMath::Abs(dch) == 48 /*&& TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.)
            {
                        
                        
                        Float_t posz = (DeltaZ[0])*(ChannelTOF[0]-ChannelTOF[1])/48;
                        Float_t  posz2 = (-3.5*(ChannelTOF[0]-ChannelTOF[1])/48 + (DeltaZ[1]))*(ChannelTOF[0]-ChannelTOF[1])/48;
                        
                        Float_t tw1=tempo[0]- gtime[0];
                        Float_t tw2=tempo[1]- gtime[1];
                        
                        Float_t tw1corr=tw1-(offset_p1 + x1_p1 *posz);
                        Float_t tw2corr=tw2-(offset_p2 + x2_p2 *posz);
                        
                        Float_t tbest=(tw1corr+tw2corr)/2.;
                        
                        Float_t tbestcorr=tbest-(offset_pb + x1_pb *(tw1corr-tw2corr) + x2_pb *pow((tw1corr-tw2corr),2));
                        
                        
                        hprofbestcorr->Fill((tw1corr-tw2corr),tbestcorr ,1);
                        hprofbestcorr->GetXaxis()->SetTitle("t1corr-t2corr (ps)");
                        hprofbestcorr->GetYaxis()->SetTitle("tbest_corr (ps)");
                        
                        h2bescorrtvst1corr_t2corr->Fill((tw1corr-tw2corr), tbestcorr);
                        h2bescorrtvst1corr_t2corr->GetXaxis()->SetTitle("t1corr-t2corr (ps)");
                        h2bescorrtvst1corr_t2corr->GetYaxis()->SetTitle("tbest_corr (ps)");
                        
                        
                        htbestcorr->Fill(tbestcorr);
                        htbestcorr->GetXaxis()->SetTitle("t_best_corr");


                    }
                }
                }
            }
        }
    }
    
    TF1 *fbc = new TF1("fbc","pol2",-500. , 500.);
    hprofbestcorr->Fit("fbc","R");

    
  

  
   TFile *fo2 = new TFile("output_z.root","RECREATE"); //DATI VERI
    //TFile *fo2 = new TFile("output_MC_ist_tree_z.root","RECREATE"); //DATI SIMULATI
    
    hdiff->Write();
    h2resx->Write();
    h2resz->Write();
    h2resxz1->Write();
    h2resxz2->Write();
    hresx1->Write();
    hresx2->Write();
    hresz1->Write();
    hresz2->Write();
    hresdist1->Write();
    hresz1ricen->Write();
    hresdist1ricen->Write();
    hresz2ricen->Write();
    hresdist2ricen->Write();
   // hdeltat->Write();
    hdeltaxz->Write();
    hdeltaz->Write();
    pxdeltat->Write();
    pxdeltat2->Write();
    pzdeltat->Write();
    pzdeltat2->Write();
    hdeltach->Write();
    ht1_texp->Write();
    ht2_texp->Write();
    h2tw1tw2->Write();
    hprof1->Write();
    hprof1corr->Write();
    hprof2->Write();
    hprof2corr->Write();
    htbest->Write();
    h2bestvst1corr_t2corr->Write();
    hprofbest->Write();
    hprofbestcorr->Write();
    h2bescorrtvst1corr_t2corr->Write();
    htbestcorr->Write();
    //hexp_time_pi->Write();
      fo2->Close();
}
