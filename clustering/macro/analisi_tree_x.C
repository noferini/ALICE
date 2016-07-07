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

Float_t clusterize(Float_t dx,Float_t dz,Float_t time1, Float_t time2,Float_t tot[2],Int_t chan[2],Float_t &timecomb,Float_t &dxcomb,Float_t &dzcomb); // t1_corr - t2_corr

void newResiduals(Int_t ch[2],Float_t DeltaX[2],Float_t DeltaZ[2],Float_t time[2],Float_t tot[2]);


void analisi_tree_x(){ //faccio gli istogrammi dal Tree T creato nel file CheckESD.C
    gROOT->Reset();
    gStyle->SetOptStat(0012);
    gStyle->SetOptFit(0111);
    
    Bool_t kCal=kFALSE;
    
    TFile *fcal = TFile::Open("calibration.root");
    TProfile *hCalX;
    TProfile *hCalZ;
    if(fcal){
        kCal=kTRUE;
        hCalX = (TProfile *) fcal->Get("hCalX");
        hCalZ = (TProfile *) fcal->Get("hCalZ");
    }
    
    // definire istogrammi (ricordarsi di fare il write nel file successivamente)
    //TH1F *hdeltat = new TH1F("hdeltat","inside the pad (cl_{1}) - cluster along x;t_{1} - t_{2} (ps)",400,-500,500);
    TH1F *hdeltax = new TH1F("hdeltax","inside the pad (cl_{1}) - cluster along x;#Deltax_{1} - #Deltax_{2} (cm)",100,-10,10);
    TH1F *hdeltach = new TH1F("hdeltach","inside the pad (cl_{1}) - cluster along x;ch_{1} - ch_{2}",30,-15,15);
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
    TH1F *hresx1ricen = new TH1F("hresx1ricen","Residui x1 ricentrato",100,-10.,10.);
    TH1F *hresdist1ricen = new TH1F("hresdist1ricen","Residui  sqrt(x1ricentrato^2+z1^2)",100,0.,10.);
    TH1F *hresx2ricen = new TH1F("hresx2ricen","Residui x1 ricentrato",100,-10.,10.);
    TH1F *hresdist2ricen = new TH1F("hresdist2ricen","Residui  sqrt(x2ricentrato^2+z2^2)",100,0.,10.);
    
    //Istogramma 2D per i residui
    TH2F *h2resx = new TH2F("h2resx" , "Residui dx1 e dx2", 100, -10. , 10. , 100 , -10. , 10.);
    TH2F *h2resz = new TH2F("h2resz" , "Residui dz1 e dz2", 100, -10. , 10. , 100 , -10. , 10.);
    TH2F *h2resxz1 = new TH2F("h2resxz1" , "Residui dx1 e dz1", 100, -4. , 4. , 100 , -4. , 4.);
    TH2F *h2resxz2 = new TH2F("h2resxz2" , "Residui dx2 e dz2", 100, -4. , 4. , 100 , -4. , 4.);
    
    
    //Istogramma tempi meno tempi attesi
    TH1F *ht1_texp = new TH1F("ht1_texp","",100,-2000.,2000.);
    TH1F *ht2_texp = new TH1F("ht2_texp","",100,-2000.,2000.);
    
    TH2F *h2tw1tw2 = new TH2F("h2tw1tw2" , "", 50, -500. , 500. , 50 , -500. , 500.);
    
    
    //TH1F *hexp_time_pi = new TH1F("hexp_time_pi","hexp_time_pi",1000,0.,30000.);
    
    TProfile *hprof1  = new TProfile("hprof1","Profile t1-t_exp_pi vs dx1",26, -2.,2.);
    TProfile *hprof1corr  = new TProfile("hprof1corr","Profile t1-t_exp_pi corr vs dx1",26, -2.,2.);
    TProfile *hprof2  = new TProfile("hprof2","Profile t2-t_exp_pi vs dx2",26,-2.,2.);
    TProfile *hprof2corr  = new TProfile("hprof2corr","Profile t2-t_exp_pi corr vs dx2",26, -2.,2.);
    
    TH1F *htbest = new TH1F("htbest","htbest",100,-2000.,2000.);
    
    TH2F *h2bestvst1corr_t2corr = new TH2F("h2bestvst1corr_t2corr" , "Correzione residua", 100, -500. , 500. , 100 , -1000. , 1000.);
    TProfile *hprofbest  = new TProfile("hprofbest","Profile tbest vs t1corr-t2corr",100, -500. , 500.);
    TProfile *hprofbestcorr  = new TProfile("hprofbestcorr","Profile hprofbestcorr vs t1corr-t2corr",100, -500. , 500.);
    TH2F *h2bescorrtvst1corr_t2corr = new TH2F("h2bescorrtvst1corr_t2corr" , "Correzione residua", 100, -500. , 500. , 100 , -1000. , 1000.);
    
    TH1F *htbestcorr = new TH1F("htbestcorr","htbestcorr",100,-2000.,2000.);
    
    //  // Utilizzo TOT.
    //  TH2F *h2t1_texp_TOT= new TH2F("h2t1_texp_TOT" , "Delay time vs TOT", 50, 0. , 100.,100,-400.,400.);
    //  TProfile *hproft1_texp_TOT = new TProfile("hproft1_texp_TOT","Profile t1-t_exp_pi vs TOT",50, 0. , 100.);
    //
    //  TH1F *hdiffTOT = new TH1F("hdiffTOT","TOT difference between hits",25, -40. , 40.);
    //  TH1F *hdifflogTOT = new TH1F("hdifflogTOT","TOT log difference between hits",5, -4. , 4.);
    //
    //  //TH2F *h2distric1_TOT = new TH2F("h2distric1_TOT" , " sqrt(x1ricentrato^2+z1^2) vs TOT",50, 0. , 100. , 100,0.,10.);
    //  //TH2F *h2distric1_diffTOT = new TH2F("h2distric1_diffTOT" , " sqrt(x1ricentrato^2+z1^2) vs difTOT",25, -40. , 40. , 100,0.,10.);
    //  //TH2F *h2distric1_difflogTOT= new TH2F("h2distric1_difflogTOT" , " sqrt(x1ricentrato^2+z1^2) vs difflogTOT",25, -4. , 4. , 100,0.,10.);
    //
    //  //
    //  //TH2F *h2dx1_TOT= new TH2F("h2dx1_TOT" , " dx1 vs TOT1",50, 0. , 100. , 10,-1.25,1.25);
    //  //TH2F *h2dx2_TOT= new TH2F("h2dx2_TOT" , " dx1 vs TOT2",50, 0. , 100. , 10,-1.25,1.25);
    //  //TH2F *h2dx1_difflogTOT= new TH2F("h2dx1_difflogTOT" , " dx1 vs  difflogTOT",25, -4. , 4. , 10,-1.25,1.25);
    //
    //
    //  TH2F *h2DxM_TOTM= new TH2F("h2DxM_TOTM" , " dxM vs TOTMm",25, 0. , 40. , 20,-2.5,2.5);
    //  TH2F *h2DxM_difflogTOTMm= new TH2F("h2DxM_difflogTOTMm" , " dxM vs  difflogTOTMm",25, 0. , 3. , 20,-2.5,2.5);
    //  TH2F *h2DxM_diffnTOTMm= new TH2F("h2DxM_diffnTOTMm" , " dxM vs diffnTOTMm",25, 0. , 1. , 20,-2.5,2.5);
    
    TH2F *h2DxM_TOTM= new TH2F("h2DxM_TOTM" , " dxM vs TOTMm",25, 0. , 40. , 20,-2.5,2.5);
    TH2F *h2DxM_difflogTOTMm= new TH2F("h2DxM_difflogTOTMm" , " dxM vs  difflogTOTMm",100, 0. , 1. , 20,-2.5,2.5);
    TH2F *h2DxM_diffTOTMm= new TH2F("h2DxM_diffTOTMm" , " dxM vs diffTOTMm",25, 0. , 1. , 20,-2.5,2.5);
    
    TH1F *htbestcorrCL = new TH1F("htbestcorrCL","htbestcorrCL",100,-2000.,2000.);
    TH1F *htbestcheck = new TH1F("htbestcheck","htbestcheck",100,-50.,50.);
    
    TH1F *htbestcorrCL_TOTM_pol1_824 = new TH1F("htbestcorrCL_TOTM_pol1_824","htbestcorrCL_TOTM_pol1_824",100,-2000.,2000.);
    TH1F *htbestcorrCL_TOTM_pol3_824 = new TH1F("htbestcorrCL_TOTM_pol3_824","htbestcorrCL_TOTM_pol3_824",100,-2000.,2000.);
    // TH1F *htbestcorrCL_difflog_pol1_048 = new TH1F("htbestcorrCL_difflog_pol1_048","htbestcorrCL_difflog_pol1_048",100,-2000.,2000.);
    TH1F *htbestcorrCL_difflog_pol2 = new TH1F("htbestcorrCL_difflog_pol2","htbestcorrCL_difflog_pol2",100,-2000.,2000.);
    
    TH1F *htbestcorrCL_dTOTsub_pol1 = new TH1F("htbestcorrCL_dTOTsub_pol1","htbestcorrCL_dTOTsub_pol1",100,-2000.,2000.);
    
    
    
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
    //T->SetBranchAddress("thetay",&thetay);
    T->SetBranchAddress("StartTime",&StartTime);
    T->SetBranchAddress("StartTimeRes",&StartTimeRes);
    
    Int_t nentries = (Int_t)T->GetEntries();
    
    Int_t ntotcl=0;
    Int_t n2cl=0;
    
    for(Int_t i=0;i<nentries;i++){
        
        T->GetEntry(i);
        //     if(StartTimeRes > 10) continue;
        
        for(Int_t ip=0;ip<ncluster;ip++){
            tempo[ip] -= StartTime;
            Int_t strip=ChannelTOF[0]/96;
            if(kCal){
                DeltaX[ip] -= hCalX->GetBinContent(strip+1);
                DeltaZ[ip] -= hCalZ->GetBinContent(strip+1);
            }
        }
        
        if(kCal && ncluster==2){ // check again the residual
            if(DeltaX[0]*DeltaX[0] > DeltaX[1]*DeltaX[1]){
                ChannelTOF[99] = ChannelTOF[0];
                tempo[99] = tempo[0];
                DeltaX[99] = DeltaX[0];
                DeltaZ[99] = DeltaZ[0];
                TOT[99] = TOT[0];
                
                ChannelTOF[0] = ChannelTOF[1];
                tempo[0] = tempo[1];
                DeltaX[0] = DeltaX[1];
                DeltaZ[0] = DeltaZ[1];
                TOT[0] = TOT[1];
                
                ChannelTOF[1] = ChannelTOF[9];
                tempo[1] = tempo[9];
                DeltaX[1] = DeltaX[9];
                DeltaZ[1] = DeltaZ[9];
                TOT[1] = TOT[9];
            }
        }
        
        ntotcl++;
        
        if(ncluster == 2){
            n2cl++;
            
            if(impulso_trasv>0.8 && impulso_trasv<1.2){ // serve per gli exp time
                
                //if(exp_time_pi[0] > 0. && exp_time_pi[1] > 0.)
                
                //Int_t dch = TMath::Abs(ChannelTOF[0]-ChannelTOF[1]);
                
                if((ChannelTOF[0]/96)==(ChannelTOF[1]/96) /* così sono nella stessa strip*/ /*&& (ChannelTOF[0]/8)==(ChannelTOF[1]/8) /*così prendo stesso NINO, per vedere cross talk*/ ){


		  newResiduals(ChannelTOF,DeltaX,DeltaZ,tempo,TOT);

		  if( TMath::Abs(DeltaZ[0])<1.75 ){ //prendo che il pad machato sia dentro lungo le z
                        
                        Int_t dch = ChannelTOF[0]-ChannelTOF[1];
                        
                        if(TMath::Abs(dch) == 1 /* seleziono x adiacenti*/ && TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.){// poi dovrei farlo anche per 1
                            // hexp_time_pi->Fill(exp_time_pi[0]);
                            
                            
                            Float_t diff=tempo[0]-tempo[1];
                            hdiff->Fill(diff);
                            hdiff->GetXaxis()->SetTitle("t1-t2 (ps)");
                            
                            //
                            //              Float_t diffTOT=TOT[0]-TOT[1];
                            //              hdiffTOT->Fill(diffTOT);
                            //              hdiffTOT->GetXaxis()->SetTitle("TOT1-TOT2 ");
                            //
                            //              Float_t difflogTOT=log(TOT[0])-log(TOT[1]);
                            //              hdifflogTOT->Fill(difflogTOT);
                            //              hdifflogTOT->GetXaxis()->SetTitle("logTOT1-logTOT2 ");
                            
                            
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
                            
                            
                            Float_t Dx1ricen;
                            if((dch)== +1){Dx1ricen=DeltaX[0]+1.25;}//ricentro il residuo al centro del pad
                            if((dch)== -1){Dx1ricen=DeltaX[0]-1.25;}
                            hresx1ricen->Fill(Dx1ricen);
                            hresx1ricen->GetXaxis()->SetTitle("Dx1 ricentrato (cm)");
                            hresdist1ricen->Fill(sqrt(Dx1ricen*Dx1ricen+DeltaZ[0]*DeltaZ[0]));
                            hresdist1ricen->GetXaxis()->SetTitle("sqrt(Dx1ricentrato^2+Dz1^2) (cm)");
                            
                            //              Float_t Dx2ricen;
                            //              if((dch)== +1){Dx2ricen=DeltaX[1]-1.25;}//ricentro il residuo al centro del pad
                            //              if((dch)== -1){Dx2ricen=DeltaX[1]+1.25;}
                            //              hresx2ricen->Fill(Dx2ricen);
                            //              hresx2ricen->GetXaxis()->SetTitle("Dx2 ricentrato (cm)");
                            //              hresdist2ricen->Fill(sqrt(Dx2ricen*Dx2ricen+DeltaZ[1]*DeltaZ[1]));
                            //              hresdist2ricen->GetXaxis()->SetTitle("sqrt(Dx2ricentrato^2+Dz2^2) (cm)");
                            
                            Float_t res1 = DeltaX[0]*DeltaX[0] + DeltaZ[0]*DeltaZ[0];
                            Float_t res2 = DeltaX[1]*DeltaX[1] + DeltaZ[1]*DeltaZ[1];
                            
                            Float_t posx = (DeltaX[0])*(ChannelTOF[0]-ChannelTOF[1]);
                            Float_t  posx2 = (-2.5*(ChannelTOF[0]-ChannelTOF[1]) + (DeltaX[1]))*(ChannelTOF[0]-ChannelTOF[1]);
                            
                            
                            //cerco correlazione tra TOT e residui vari
                            //              h2dx1_TOT->Fill(TOT[0],DeltaX[0]);
                            //              h2dx1_TOT->GetXaxis()->SetTitle("TOT ");
                            //              h2dx1_TOT->GetYaxis()->SetTitle("Dx1 (cm)");
                            //
                            //              h2dx2_TOT->Fill(TOT[1],DeltaX[0]);
                            //              h2dx2_TOT->GetXaxis()->SetTitle("TOT ");
                            //              h2dx2_TOT->GetYaxis()->SetTitle("Dx2 (cm)");
                            //
                            //              h2dx1_difflogTOT->Fill(difflogTOT,DeltaX[0]);
                            //              h2dx1_difflogTOT->GetXaxis()->SetTitle("logTOT1-logTOT2 ");
                            //              h2dx1_difflogTOT->GetYaxis()->SetTitle("Dx1 (cm)");
                            
                            //cerco correlazione tra TOT e residui vari
                            
                            //              h2distric1_TOT->Fill(TOT[0],sqrt(Dx1ricen*Dx1ricen+DeltaZ[0]*DeltaZ[0]));
                            //              h2distric1_TOT->GetXaxis()->SetTitle("TOT ");
                            //              h2distric1_TOT->GetYaxis()->SetTitle("sqrt(Dx1ricentrato^2+Dz1^2) (cm)");
                            //
                            //              h2distric1_diffTOT->Fill(diffTOT,sqrt(Dx1ricen*Dx1ricen+DeltaZ[0]*DeltaZ[0]));
                            //              h2distric1_diffTOT->GetXaxis()->SetTitle("TOT1-TOT2 ");
                            //              h2distric1_diffTOT->GetYaxis()->SetTitle("sqrt(Dx1ricentrato^2+Dz1^2) (cm)");
                            //
                            //              h2distric1_difflogTOT->Fill(difflogTOT,sqrt(Dx1ricen*Dx1ricen+DeltaZ[0]*DeltaZ[0]));
                            //              h2distric1_difflogTOT->GetXaxis()->SetTitle("logTOT1-logTOT2 ");
                            //              h2distric1_difflogTOT->GetYaxis()->SetTitle("sqrt(Dx1ricentrato^2+Dz1^2) (cm)");
                            
                            Float_t DxM; // quello con TOT magiore e quindi dove plausibilmente è passata la particella(diverso da dx1!!!)
                            Float_t TOT_M;
                            Float_t TOT_m;
                            
                            if(TOT[0]>0 && TOT[1]>0){
                                if( TOT[0]>TOT[1]){
                                    TOT_M=TOT[0];
                                    TOT_m=TOT[1];
                                    
                                    if(ChannelTOF[0]<ChannelTOF[1]){ DxM=DeltaX[0]-1.25;}
                                    if(ChannelTOF[0]>ChannelTOF[1]){ DxM=-DeltaX[0]-1.25;}
                                }
                                
                                if( TOT[0]<TOT[1]){
                                    TOT_M=TOT[1];
                                    TOT_m=TOT[0];
                                    
                                    //if(ChannelTOF[0]<ChannelTOF[1]){ DxM=-DeltaX[1]-1.25;}
                                    //if(ChannelTOF[0]>ChannelTOF[1]){ DxM=DeltaX[1]-1.25;}
                                    if(ChannelTOF[0]<ChannelTOF[1]){ DxM=-DeltaX[0]+1.25;}
                                    if(ChannelTOF[0]>ChannelTOF[1]){ DxM=1.25+DeltaX[0];}
                                }
                                
                                h2DxM_TOTM->Fill(TOT_M,DxM);
                                h2DxM_TOTM->GetXaxis()->SetTitle("TOT_M");
                                h2DxM_TOTM->GetYaxis()->SetTitle("DxM (cm)");
                                
                                h2DxM_difflogTOTMm->Fill((log(TOT_M)-log(TOT_m))/(log(TOT_M)+log(TOT_m)),DxM);
                                h2DxM_difflogTOTMm->GetXaxis()->SetTitle("(logTOTM-logTOTm) /(logTOTM+logTOTm)");
                                h2DxM_difflogTOTMm->GetYaxis()->SetTitle("DxM (cm)");
                                
                                h2DxM_diffTOTMm->Fill((TOT_M-TOT_m)/(TOT_M+TOT_m+0.001),DxM);
                                h2DxM_diffTOTMm->GetXaxis()->SetTitle("(TOT_M-TOT_m)/(TOT_M+TOT_m)");
                                h2DxM_diffTOTMm->GetYaxis()->SetTitle("DxM (cm)");
                                
                                
                            }
                            
                            
                            hdeltax->Fill(posx);
                            hdeltax->GetXaxis()->SetTitle("posx (cm)");
                            
                            pxdeltat->Fill(tempo[0] - tempo[1], posx);
                            pxdeltat->GetXaxis()->SetTitle("t1-t2 (cm)");
                            pxdeltat->GetYaxis()->SetTitle("posx1 (cm)");
                            
                            pxdeltat2->Fill(tempo[0]- tempo[1], posx2);
                            pxdeltat2->GetXaxis()->SetTitle("t1-t2 (cm)");
                            pxdeltat2->GetYaxis()->SetTitle("posx2 (cm)");
                            
                            
                            pzdeltat->Fill(tempo[0] - tempo[1], DeltaZ[0]);
                            pzdeltat->GetXaxis()->SetTitle("t1-t2 (cm)");
                            pzdeltat->GetYaxis()->SetTitle("Dz1 (cm)");
                            
                            pzdeltat2->Fill(tempo[0]- tempo[1], DeltaZ[1]);
                            pzdeltat2->GetXaxis()->SetTitle("t1-t2 (cm)");
                            pzdeltat2->GetYaxis()->SetTitle("Dz2 (cm)");
                            
                            hdeltach->Fill(ChannelTOF[0]-ChannelTOF[1]);
                            hdeltach->GetXaxis()->SetTitle("Dch");
                            
                            //if(TMath::Abs(DeltaZ[0]-DeltaZ[1])<=0.5) // perchè così dovrebbero provenire(credo) da stesso cluster in quanto, se così è , essendo loro adiacenti sulle x, dovrebbero avere circe dtesso residuo in z. NON sono sicura sia necessario fare ciò. // non dovrebbe più servire perchè ho modificato Check ESD spostando i residui dentro l'if
                            
                            //h2resz->Fill(DeltaZ[0],DeltaZ[1]);
                            Float_t tw1=tempo[0]- exp_time_pi[0];
                            Float_t tw2=tempo[1]- exp_time_pi[0]/*/L[0]*L[1]*/; //non metto exp_time_pi[1] poichè momentaneamente non va
                            //              //cerco correlazione tra TOT e delay time
                            //              h2t1_texp_TOT->Fill(TOT[0],tw1);
                            //              h2t1_texp_TOT->GetXaxis()->SetTitle("TOT ");
                            //              h2t1_texp_TOT->GetYaxis()->SetTitle("t_{0}-t_{exp #pi} (ps)");
                            //
                            //              hproft1_texp_TOT->Fill(TOT[0],tw1 ,1);
                            //              hproft1_texp_TOT->GetXaxis()->SetTitle("TOT");
                            //              hproft1_texp_TOT->GetYaxis()->SetTitle("t_{0}-t_{exp #pi} (ps)");
                            
                            ht1_texp->Fill(tw1);
                            ht1_texp->GetXaxis()->SetTitle("t_{1} - t_{exp #pi} (ps)");
                            
                            ht2_texp->Fill(tw2);
                            ht2_texp->GetXaxis()->SetTitle("t_{2} - t_{exp #pi} (ps)");
                            
                            h2tw1tw2->Fill(tw1, tw2);
                            h2tw1tw2->GetXaxis()->SetTitle("t_{1} - t_{exp #pi}  (ps)");
                            h2tw1tw2->GetYaxis()->SetTitle("t_{2} - t_{exp #pi}  (ps)");
                            
                            /////////////////////PROFILE  1
                            
                            //riempio istogrammma, poi lo fitto fuori dal loop
                            hprof1->Fill(posx,tw1 ,1);
                            hprof1->GetXaxis()->SetTitle("posx1 (cm)");
                            hprof1->GetYaxis()->SetTitle("t_{1}-t_{exp #pi} (ps)");
                            
                            
                            hprof2->Fill(posx,tw2 ,1);
                            hprof2->GetXaxis()->SetTitle("posx1 (cm)");
                            hprof2->GetYaxis()->SetTitle("t_{2}-t_t_{exp #pi} (ps)");
                            
                        }
                        
                    }
                }
            }
        }
    }
    TF1 *f1 = new TF1("f1","pol1",-1.25,0.5);
    
    //hprof1->GetYaxis()->SetRangeUser(-40.,100.);
    hprof1->Fit("f1","R");
    
    Double_t offset_p1,x1_p1;
    offset_p1= f1->GetParameter(0);
    x1_p1= f1->GetParameter(1);
    
    TF1 *f2 = new TF1("f2","pol1",-1.25,0.5);
    
    //hprof2->GetYaxis()->SetRangeUser(-40.,150.);
    hprof2->Fit("f2","R");
    
    Double_t offset_p2,x2_p2;
    offset_p2= f2->GetParameter(0);
    x2_p2= f2->GetParameter(1);
    
    
    for (Int_t nu=0 ; nu<nentries ; nu++){
        T->GetEntry(nu);
        //if(StartTimeRes > 10) continue;
        
        for(Int_t ip=0;ip<ncluster;ip++){
            tempo[ip] -= StartTime;
            Int_t strip=ChannelTOF[0]/96;
            if(kCal){
                DeltaX[ip] -= hCalX->GetBinContent(strip+1);
                DeltaZ[ip] -= hCalZ->GetBinContent(strip+1);
            }
        }
        
        if(kCal && ncluster==2){ // check again the residual
            if(DeltaX[0]*DeltaX[0] > DeltaX[1]*DeltaX[1]){
                ChannelTOF[99] = ChannelTOF[0];
                tempo[99] = tempo[0];
                DeltaX[99] = DeltaX[0];
                DeltaZ[99] = DeltaZ[0];
                TOT[99] = TOT[0];
                
                ChannelTOF[0] = ChannelTOF[1];
                tempo[0] = tempo[1];
                DeltaX[0] = DeltaX[1];
                DeltaZ[0] = DeltaZ[1];
                TOT[0] = TOT[1];
                
                ChannelTOF[1] = ChannelTOF[9];
                tempo[1] = tempo[9];
                DeltaX[1] = DeltaX[9];
                DeltaZ[1] = DeltaZ[9];
                TOT[1] = TOT[9];
            }
        }
        
        if(ncluster == 2){
            if(impulso_trasv>0.8 && impulso_trasv<1.2){ // serve per gli exp time
                if((ChannelTOF[0]/96)==(ChannelTOF[1]/96) /*&& (ChannelTOF[0]/8)==(ChannelTOF[1])*/){ // così sono nella stessa strip
                    if( TMath::Abs(DeltaZ[0])<1.75){
                        
		      newResiduals(ChannelTOF,DeltaX,DeltaZ,tempo,TOT);

                        Int_t dch = ChannelTOF[0]-ChannelTOF[1];
                        
                        if(TMath::Abs(dch) == 1 && TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.){
                            Float_t posx = (DeltaX[0])* dch;
                            Float_t  posx2 = (-2.5*(ChannelTOF[0]-ChannelTOF[1]) + (DeltaX[1]))*(ChannelTOF[0]-ChannelTOF[1]);
                            
                            Float_t tw1=tempo[0]- exp_time_pi[0];
                            Float_t tw2=tempo[1]- exp_time_pi[0];
                            
                            Float_t tw1corr=tw1-(offset_p1 + x1_p1 *posx);
                            Float_t tw2corr=tw2-(offset_p2 + x2_p2 *posx);
                            
                            hprof1corr->Fill(posx,tw1corr,1);
                            hprof1corr->GetXaxis()->SetTitle("Dx1 (cm)");
                            hprof1corr->GetYaxis()->SetTitle("t1_corr=t1-t_exp_pi corr(ps)");
                            
                            hprof2corr->Fill(posx2,tw2corr,1);
                            hprof2corr->GetXaxis()->SetTitle("Dx2 (cm)");
                            hprof2corr->GetYaxis()->SetTitle("t2_corr=t2-t_exp_pi corr(ps)");
                            
                            Float_t tbest=(tw1corr+tw2corr)/2.;
                            htbest->Fill(tbest);
                            htbest->GetXaxis()->SetTitle("t_best=(t1_corr +t2_corr)/2 (ps)");
                            
                            h2bestvst1corr_t2corr->Fill((tw1corr-tw2corr), tbest);
                            h2bestvst1corr_t2corr->GetXaxis()->SetTitle("t1_corr - t2_corr (ps)");
                            h2bestvst1corr_t2corr->GetYaxis()->SetTitle("t_best(ps)");
                            
                            
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
    
    
    for (Int_t nu=0 ; nu<nentries ; nu++){
        T->GetEntry(nu);
        //        if(StartTimeRes > 10) continue;
        
        for(Int_t ip=0;ip<ncluster;ip++){
            tempo[ip] -= StartTime;
            Int_t strip=ChannelTOF[0]/96;
            if(kCal){
                DeltaX[ip] -= hCalX->GetBinContent(strip+1);
                DeltaZ[ip] -= hCalZ->GetBinContent(strip+1);
                
            }
        }
        
        if(kCal && ncluster==2){ // check again the residual
            if(DeltaX[0]*DeltaX[0] > DeltaX[1]*DeltaX[1]){
                ChannelTOF[99] = ChannelTOF[0];
                tempo[99] = tempo[0];
                DeltaX[99] = DeltaX[0];
                DeltaZ[99] = DeltaZ[0];
                TOT[99] = TOT[0];
                
                ChannelTOF[0] = ChannelTOF[1];
                tempo[0] = tempo[1];
                DeltaX[0] = DeltaX[1];
                DeltaZ[0] = DeltaZ[1];
                TOT[0] = TOT[1];
                
                ChannelTOF[1] = ChannelTOF[9];
                tempo[1] = tempo[9];
                DeltaX[1] = DeltaX[9];
                DeltaZ[1] = DeltaZ[9];
                TOT[1] = TOT[9];
            }
        }
        
        if(ncluster == 2) {
            if(impulso_trasv>0.8 && impulso_trasv<1.2){ // serve per gli exp time
                if((ChannelTOF[0]/96)==(ChannelTOF[1]/96)/* && (ChannelTOF[0]/8)==(ChannelTOF[1])*/){ // così sono nella stessa strip
                    if( TMath::Abs(DeltaZ[0])<1.75) {
		      newResiduals(ChannelTOF,DeltaX,DeltaZ,tempo,TOT);


		      Int_t dch = ChannelTOF[0]-ChannelTOF[1];
                        
                        if(TMath::Abs(dch) == 1 && TMath::Abs(tempo[0]-exp_time_pi[0])<800./* per avere circa 3 sigma che sia un pi*/ && TMath::Abs(tempo[0] - tempo[1])<470.){
                            Float_t posx = (DeltaX[0])* dch;
                            Float_t  posx2 = (-2.5*(ChannelTOF[0]-ChannelTOF[1]) + (DeltaX[1]))*(ChannelTOF[0]-ChannelTOF[1]);
                            
                            Float_t tw1=tempo[0]- exp_time_pi[0];
                            Float_t tw2=tempo[1]- exp_time_pi[0];
                            
                            Float_t tw1corr=tw1-(offset_p1 + x1_p1 *posx);
                            Float_t tw2corr=tw2-(offset_p2 + x2_p2 *posx);
                            
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
                            
                            Float_t timecomb,dxcomb,dzcomb;
                            clusterize( DeltaX[0],0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            htbestcorrCL->Fill(timecomb);
                            htbestcorrCL->GetXaxis()->SetTitle("timecomb");
                            htbestcheck->Fill(timecomb-tbestcorr);
                            htbestcheck->GetXaxis()->SetTitle("timecomb-tbestcorr");
                            
                            
                            Float_t TOT_M;
                            Float_t TOT_m;
                            if(TOT[0]>0 && TOT[1]>0)
                            {
                                if( TOT[0]>TOT[1])
                                {
                                    TOT_M=TOT[0];
                                    TOT_m=TOT[1];
                                }
                                if( TOT[0]<TOT[1])
                                {
                                    TOT_M=TOT[1];
                                    TOT_m=TOT[0];
                                }
                            }
                            
                            if(TOT_M<8.) TOT_M=8;
                            if(TOT_M>24.) TOT_M=24.;
                            
                            Double_t TOTM_pol1_824=0.309374-TOT_M*0.0299359;
                            clusterize(TOTM_pol1_824,0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            htbestcorrCL_TOTM_pol1_824->Fill(timecomb);
                            htbestcorrCL_TOTM_pol1_824->GetXaxis()->SetTitle("timecomb");
                            
                            Double_t TOTM_pol3_824=-0.811295+0.183438*TOT_M-0.0131374*pow(TOT_M,2)+0.000260201*pow(TOT_M,3);
                            clusterize(TOTM_pol3_824,0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            htbestcorrCL_TOTM_pol3_824->Fill(timecomb);
                            htbestcorrCL_TOTM_pol3_824->GetXaxis()->SetTitle("timecomb");
                            
                            
                            Double_t difflog=(log(TOT_M)-log(TOT_m))/(log(TOT_M)+log(TOT_m));
                            
                            //                            if(difflog>0.48) difflog=0.48;
                            //
                            //                            Double_t difflog_pol1_048=-1.69085-5.19765*difflog;
                            //                            clusterize(difflog_pol1_048,0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            //                            htbestcorrCL_difflog_pol1_048->Fill(timecomb);
                            //                            htbestcorrCL_difflog_pol1_048->GetXaxis()->SetTitle("timecomb");
                            
                            if(difflog>0.07) difflog=0.07;
                            
                            Double_t difflog_pol2=-0.-4.59919*difflog+26.5336*difflog*difflog;
                            clusterize(difflog_pol2,0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            htbestcorrCL_difflog_pol2->Fill(timecomb);
                            htbestcorrCL_difflog_pol2->GetXaxis()->SetTitle("timecomb");
                            
                            
                            Double_t TOTsub=(TOT_M-TOT_m)/(TOT_M+TOT_m+0.001);
                            if(TOTsub>0.24) TOTsub=0.24;
                            Double_t dTOTsub_pol1=0.+-1.22486*TOTsub; //il p0 l ho fissato a 0 per sto fit
                            clusterize(dTOTsub_pol1,0., tw1,tw2, TOT, ChannelTOF,timecomb, dxcomb,dzcomb);
                            htbestcorrCL_dTOTsub_pol1->Fill(timecomb);
                            htbestcorrCL_dTOTsub_pol1->GetXaxis()->SetTitle("timecomb");
                            
                            
                        }
                    }
                }
            }
        }
    }
    Float_t rapp2tot;
    rapp2tot=1.*n2cl/ntotcl;
    cout<<"Il numero di cluster totale è: "<< ntotcl<<endl<< " Il numero di eventi con 2 cluster è: "<<n2cl<<endl<< "Il rapporto è quindi "<< rapp2tot <<endl;

    
    TF1 *fbc = new TF1("fbc","pol2",-500. , 500.);
    hprofbestcorr->Fit("fbc","R");
    
    TFile *fo2 = new TFile("output_ist_tree_x.root","RECREATE");// DATI VERI
    //TFile *fo2 = new TFile("output_MC_ist_tree_x.root","RECREATE");// DATI SIMULATI
    
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
    hresx1ricen->Write();
    hresdist1ricen->Write();
    hresx2ricen->Write();
    hresdist2ricen->Write();
    // hdeltat->Write();
    hdeltaxz->Write();
    hdeltax->Write();
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
    htbestcorrCL->Write();
    htbestcheck->Write();
    htbestcorrCL_TOTM_pol1_824->Write();
    htbestcorrCL_TOTM_pol3_824->Write();
    //htbestcorrCL_difflog_pol1_048->Write();
    htbestcorrCL_difflog_pol2->Write();
    htbestcorrCL_dTOTsub_pol1->Write();
    //hexp_time_pi->Write();
    //h2t1_texp_TOT->Write();
    //hproft1_texp_TOT->Write();
    //hdiffTOT->Write();
    //hdifflogTOT->Write();
    //    h2distric1_TOT->Write();
    //    h2distric1_diffTOT->Write();
    //    h2distric1_difflogTOT->Write();
    //    h2dx1_TOT->Write();
    //    h2dx2_TOT->Write();
    //    h2dx1_difflogTOT->Write();
    h2DxM_TOTM->Write();
    h2DxM_difflogTOTMm->Write();
    h2DxM_diffTOTMm->Write();
    fo2->Close();
    system("say Ehi you, I have done");
}


Float_t clusterize(Float_t dx,Float_t dz,Float_t time1, Float_t time2,Float_t tot[2],Int_t chan[2],Float_t &timecomb,Float_t &dxcomb,Float_t &dzcomb){
    timecomb = time1;
    dxcomb = dx;
    dzcomb = dz;
    Float_t dtime = 9999999;
    
    // check if the are in the same strip.
    if(Int_t(chan[0]/96) != Int_t(chan[1]/96)) return dtime;
    
    // preselection on t1 - t2
    Float_t timewindow = 500;
    if(TMath::Abs(time1 - time2) > timewindow) return dtime;
    
    // deltax limit
    Float_t dxMin = -0.5;
    Float_t dxMax = 1.25;
    // deltaz limit
    Float_t dzMin = -0.5;
    Float_t dzMax = 1.75;
    
    // t1 corr corr = axt1 * dx + bxt1;
    Float_t axt1 = -34.05;
    Float_t bxt1 = 31.67;
    
    // t2 corr corr = axt2 * dx + bxt2;
    Float_t axt2 = -1.02;
    Float_t bxt2 = 79.64;
    
    // t1 corr corr = azt1 * dz + bzt1;
    Float_t azt1 = -13.01;
    Float_t bzt1 = 79.94;
    
    // t2 corr corr = azt2 * dx + bzt2;
    Float_t azt2 = 7.265;
    Float_t bzt2 = 124.5;
    
    
    Int_t mode = 0; // mode 1 = x, 2 = z
    
    Int_t dchan = chan[0] - chan[1];
    if(dchan == 1){ //ATENTA  ho invertito is segni dchan
        mode = 1;
    }
    else if(dchan == -1){
        mode = 1;
        dx *= -1;
    }
    else if(dchan == 48){
        mode = 2;
    }
    else if(dchan == -48){
        mode = 2;
        dz *= -1;
    }
    
    //Momentaneamente sta aprte è commentata perchè voglio vedere se mitorna esattaemten il tempo che mi viene da corresioni
    //if(dx < dxMin) dx = dxMin;
    //else if(dx > dxMax) dx = dxMax;
    //if(dz < dzMin) dz = dzMin;
    //else if(dz > dzMax) dz = dzMax;
    
    // equalization corr along dx= offsetx + par1x * dx + par2x * dx^2
    Float_t offsetx = -24.2;
    Float_t par1x = -0.05069;
    Float_t par2x = 0.001075;
    // equalization corr along dz= offsetz + par1z * dz + par2z * dz^2
    Float_t offsetz = -21.6;
    Float_t par1z = -0.04014;
    Float_t par2z = 0.0009758;
    
    
    switch(mode){
        case 1:
            time1 -= axt1*dx + bxt1;
            time2 -= axt2*dx + bxt2;
            dtime = time1 - time2;
            timecomb = (time1+time2)*0.5 - offsetx - par1x*dtime - par2x*dtime*dtime;
            if(dchan == -1){
                dxcomb -= 1.25;
            }
            else if(dchan == 1){
                dxcomb += 1.25;
            }
            break;
        case 2:
            time1 -= azt1*dz + bzt1;
            time2 -= azt2*dz + bzt2;
            dtime = time1 - time2;
            timecomb = (time1+time2)*0.5 - offsetz - par1z*dtime - par2z*dtime*dtime;
            if(dchan == -48){
                dzcomb -= 1.75;
            }
            else if(dchan == 48){
                dzcomb += 1.75;
            }
            break;
    }
    
    return dtime;
}

void newResiduals(Int_t ch[2],Float_t DeltaX[2],Float_t DeltaZ[2],Float_t time[2],Float_t tot[2]){
  if(tot[1] > tot[0]){
    Float_t tt,tott;
    Int_t cht;
    tt = time[0];
    tott = tot[0];
    cht = ch[0];

    time[0] = time[1];
    tot[0] = tot[1]+0.001;
    ch[0] = ch[1];
    time[1] = tt;
    tot[1] = tott+0.001;
    ch[1] = cht;
  }

  Double_t difflog=(log(tot[0])-log(tot[1]))/(log(tot[0])+log(tot[1]));
                            
  if(difflog>0.07) difflog=0.07;
                            
  Double_t difflog_pol2=-0.-4.59919*difflog+26.5336*difflog*difflog;
  
  if((ch[1]%48)-(ch[0]%48) ==1){
    DeltaX[0]=1.25+difflog_pol2;
    DeltaX[1]=-1.25+difflog_pol2;
  }     
  else if((ch[1]%48)-(ch[0]%48) ==-1){
    DeltaX[0]=-1.25-difflog_pol2;
    DeltaX[1]=1.25-difflog_pol2;
  }
  else{
    DeltaX[0]=0;
    DeltaX[1]=0;
    }
  
  if(((ch[1]/48)%2)-((ch[0]/48)%2) ==1){
    DeltaZ[0]=1.75+difflog_pol2;
    DeltaZ[1]=-1.75+difflog_pol2;
  }     
  else if(((ch[1]/48)%2)-((ch[0]/48)%2) ==-1){
    DeltaZ[0]=-1.75-difflog_pol2;
    DeltaZ[1]=1.75-difflog_pol2;
  }
  else{
    DeltaZ[0]=0;
    DeltaZ[1]=0;
  }
}


