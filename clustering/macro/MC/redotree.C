// definire istogrammi (ricordarsi di fare il write nel file successivamente)
TH1F *hdt = new TH1F("hdtPi","pions 0.9 < p_{T} < 1.5 GeV/c;t - t_{exp}^{#pi}",100,-1000,1000);
TH1F *hdtKa = new TH1F("hdtKa","kaons 0.9 < p_{T} < 1.5 GeV/c;t - t_{exp}^{K}",100,-1000,1000);
TH1F *hdtPr = new TH1F("hdtPr","protons 0.9 < p_{T} < 1.5 GeV/c;t - t_{exp}^{p}",100,-1000,1000);
TH1F *hdeltat = new TH1F("hdeltat","inside the pad (cl_{1}) - cluster along x;t_{1} - t_{2} (ps)",100,-2000,2000);
TH1F *hdeltax = new TH1F("hdeltax","inside the pad (cl_{1}) - cluster along x;#Deltax_{1} - #Deltax_{2} (cm)",100,-10,10);
TH1F *hdeltach = new TH1F("hdeltach","inside the pad (cl_{1}) - cluster along x;ch_{1} - ch_{2}",30,-15,15);
TH2F *hdeltaxz = new TH2F("hdeltaxz","Cluster relative distance - cluster along x;#Deltax_{1} - #Deltax_{2} (cm);#Deltaz_{1} - #Deltaz_{2} (cm)",1000,-10,10,1000,-10,10);
TH2F *pDeltat = new TH2F("pDeltat","inside the pad (cl_{1}) - cluster along x;cl_{1} #Deltax (cm);t_{1} - t_{2} (ps)",100,-10,10,200,-1000,1000);
TH2F *pDeltat2 = new TH2F("pDeltat2","inside the pad (cl_{1}) - cluster along x;cl_{2} #Deltax (cm);t_{1} - t_{2} (ps)",100,-10,10,200,-1000,1000);

//CFC:Istogramma numero cluster(qua coincidono con le hit)
TH1F *hnc = new TH1F("hnc","Number of Hits",12,0.,12.);
//Istogramma differenza tempi TOF
TH1F *hdiff = new TH1F("hdiff","Time difference between hits",40,-1000.,1000.);

//Istogramma 1D per i residui
TH1F *hresx1 = new TH1F("hresx1","Residui x1",10,-10.,10.);
TH1F *hresx2 = new TH1F("hresx2","Residui x2",10,-10.,10.);
TH1F *hresz1 = new TH1F("hresz1","Residui z1",10,-10.,10.);
TH1F *hresz2 = new TH1F("hresz2","Residui z2",10,-10.,10.);


//Istogramma 2D per i residui
TH2F *h2resx = new TH2F("h2resx" , "Residui dx1 e dx2", 40, -10. , 10. , 40 , -10. , 10.);
TH2F *h2resz = new TH2F("h2resz" , "Residui dz1 e dz2", 40, -10. , 10. , 40 , -10. , 10.);
TH2F *h2resxz1 = new TH2F("h2resxz1" , "Residui dx1 e dz1", 40, -10. , 10. , 40 , -10. , 10.);
TH2F *h2resxz2 = new TH2F("h2resxz2" , "Residui dx2 e dz2", 40, -10. , 10. , 40 , -10. , 10.);



//File per il Tree T
TFile *fotree = new TFile("AnalysisResults.root","RECREATE");
//Tree di prova
TTree *T = new TTree("T","test");


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

Float_t smearX = 0.2;
Float_t smearZ = 1.5;

void SimulateTime();
void ReMatch();
void GetTrueCoord();
void CalibrationGeant();
TProfile *hpx,*hpz;

void redotree(){

  TFile *fcal = new TFile("calibration.root");
  hpx = (TProfile *) fcal->Get("hgeantXcal");
  hpz = (TProfile *) fcal->Get("hgeantZcal");

  TChain *t = new TChain("T");
  FILE *fin = fopen("lista","r");
  char name[100];
  Int_t nfile = 0;
  while(fscanf(fin,"%s",name)==1 && nfile < 1000){
    t->AddFile(name);
    nfile++;
  }


  TProfile *hx = new TProfile("hx","x alignement per strip;# strip;#DeltaX (cm)",1700,0,1700);
  TProfile *hz = new TProfile("hz","z alignement per strip;# strip;#DeltaZ (cm)",1700,0,1700);


  t->SetBranchAddress("ncluster",&ncluster);
  t->SetBranchAddress("tempo",tempo);
  t->SetBranchAddress("gtime",&gtime);
  t->SetBranchAddress("DXtrue",&dxt);
  t->SetBranchAddress("DZtrue",&dzt);
  t->SetBranchAddress("DeltaX",DeltaX);
  t->SetBranchAddress("DeltaZ",DeltaZ);
  t->SetBranchAddress("exp_time_el",exp_time_el);
  t->SetBranchAddress("exp_time_mu",exp_time_mu);
  t->SetBranchAddress("exp_time_pi",exp_time_pi);
  t->SetBranchAddress("exp_time_ka",exp_time_ka);
  t->SetBranchAddress("exp_time_pr",exp_time_pr);
  t->SetBranchAddress("L",L);
  t->SetBranchAddress("ChannelTOF",ChannelTOF);
  t->SetBranchAddress("impulso_trasv",&impulso_trasv);
  t->SetBranchAddress("res",res);
  t->SetBranchAddress("charge",&charge);
  t->SetBranchAddress("phi",&phi);
  t->SetBranchAddress("eta",&eta);
  t->SetBranchAddress("secPhi",&secAngle);
  t->SetBranchAddress("cval",cval);
  t->SetBranchAddress("mism",&mism);
  t->SetBranchAddress("dedx",&dedx);
  t->SetBranchAddress("StartTime",&StartTime);
  t->SetBranchAddress("StartTimeRes",&StartTimeRes);
  t->SetBranchAddress("rTOF",&rTOFused);
  t->SetBranchAddress("InteractionTime",&interactiontime);


  TFile *fo = new TFile("output.root","RECREATE");

  TTree *T = new TTree("T","T");
  T->Branch("ncluster",&ncluster,"ncluster/I");
  T->Branch("tempo",tempo,"tempo[ncluster]/F");
  T->Branch("gtime",&gtime,"gtime/F");
  T->Branch("DXtrue",&dxt,"DXtrue/F");
  T->Branch("DZtrue",&dzt,"DZtrue/F");
  T->Branch("DeltaX",DeltaX,"DeltaX[ncluster]/F");
  T->Branch("DeltaZ",DeltaZ,"DeltaZ[ncluster]/F");
  T->Branch("exp_time_el",exp_time_el,"exp_time_el[ncluster]/F");
  T->Branch("exp_time_mu",exp_time_mu,"exp_time_mu[ncluster]/F");
  T->Branch("exp_time_pi",exp_time_pi,"exp_time_pi[ncluster]/F");
  T->Branch("exp_time_ka",exp_time_ka,"exp_time_ka[ncluster]/F");
  T->Branch("exp_time_pr",exp_time_pr,"exp_time_pr[ncluster]/F");
  T->Branch("L",L,"L[ncluster]/F");
  T->Branch("ChannelTOF",ChannelTOF,"ChannelTOF[ncluster]/I");
  T->Branch("impulso_trasv",&impulso_trasv,"impulso_trasv/F");
  T->Branch("res",res,"res[3]/F");
  T->Branch("charge",&charge,"charge/I");
  T->Branch("phi",&phi,"phi/F");
  T->Branch("eta",&eta,"eta/F");
  T->Branch("secPhi",&secAngle,"secPhi/F");
  T->Branch("cval",cval,"cval[5]/F");
  T->Branch("mism",&mism,"mism/I");
  T->Branch("dedx",&dedx,"dedx/F");
  T->Branch("StartTime",&StartTime,"StartTime/F");
  T->Branch("StartTimeRes",&StartTimeRes,"StartTimeRes/F");
  T->Branch("rTOF",&rTOFused,"rTOF/F");
  T->Branch("InteractionTime",&interactiontime,"InteractionTime/F");

  Int_t nev = t->GetEntries();

  for(Int_t i=0;i < nev;i++){
    t->GetEvent(i);

    if((i%1000)==0)printf("%i/%i\n",i,nev);

    GetTrueCoord();

    ReMatch();

    if(dxt > -5 && dzt > -8){
      CalibrationGeant();
      SimulateTime();

      hx->Fill(int(ChannelTOF[0]/96),dxt);
      hz->Fill(int(ChannelTOF[0]/96),dzt);
      T->Fill();

    }
  }

  fo->cd();
  T->Write();
  hx->Write();
  hz->Write();
  fo->Close();
  
}
void CalibrationGeant(){
  Int_t strip=ChannelTOF[0]/96;
  dxt -= hpx->GetBinContent(strip+1);
  dzt -= hpz->GetBinContent(strip+1);
}

void GetTrueCoord(){
  Int_t detId[5];
  AliTOFGeometry::GetVolumeIndices(ChannelTOF[0],detId);
  Float_t pos[3];
  pos[0] = AliTOFGeometry::GetX(detId);
  pos[1] = AliTOFGeometry::GetY(detId);
  pos[2] = AliTOFGeometry::GetZ(detId);

  
  Int_t detId2[5];
  Float_t vx[3];
  if((ChannelTOF[0]%48) < 47){
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[0]+1,detId2);
    vx[0] = (AliTOFGeometry::GetX(detId2)-pos[0])/2.5;
    vx[1] = (AliTOFGeometry::GetY(detId2)-pos[1])/2.5;
    vx[2] = (AliTOFGeometry::GetZ(detId2)-pos[2])/2.5;
  }
  else{
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[0]-1,detId2);
    vx[0] = -(AliTOFGeometry::GetX(detId2)-pos[0])/2.5;
    vx[1] = -(AliTOFGeometry::GetY(detId2)-pos[1])/2.5;
    vx[2] = -(AliTOFGeometry::GetZ(detId2)-pos[2])/2.5;
  }
  Float_t vz[3];
  if(((ChannelTOF[0]/48)%2) < 1){
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[0]+48,detId2);
    vz[0] = -(AliTOFGeometry::GetX(detId2)-pos[0])/3.5;
    vz[1] = -(AliTOFGeometry::GetY(detId2)-pos[1])/3.5;
    vz[2] = -(AliTOFGeometry::GetZ(detId2)-pos[2])/3.5;
  }
  else{
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[0]-48,detId2);
    vz[0] = -(AliTOFGeometry::GetX(detId2)-pos[0])/3.5;
    vz[1] = -(AliTOFGeometry::GetY(detId2)-pos[1])/3.5;
    vz[2] = -(AliTOFGeometry::GetZ(detId2)-pos[2])/3.5;
  }
  
  
  xgl = pos[0] + vx[0]*dxt + vz[0]*dzt;
  ygl = pos[1] + vx[1]*dxt + vz[1]*dzt;
  zgl = pos[2] + vx[2]*dxt + vz[2]*dzt;
}


void ReMatch(){
  Float_t addX = gRandom->Gaus(0,smearX);
  Float_t addZ = gRandom->Gaus(0,smearZ);
  
  Float_t R[100];
  for(Int_t i=0;i < ncluster;i++){
    DeltaX[i] +=addX;
    DeltaZ[i] +=addZ;

    R[i] = DeltaX[i]*DeltaX[i]+DeltaZ[i]*DeltaZ[i];
  }

  Int_t index[100];
  TMath::Sort(ncluster,R,index,kFALSE);

  // adjust true info
  dxt += DeltaX[index[0]] - DeltaX[0];
  dzt += DeltaZ[index[0]] - DeltaZ[0];

  Float_t t[100],dx[100],dz[100];
  Int_t ch[100];
  for(Int_t i=0;i < ncluster;i++){
    t[i] = tempo[index[i]];
    dx[i] = DeltaX[index[i]];
    dz[i] = DeltaZ[index[i]];
    ch[i] = ChannelTOF[index[i]];
  }

  for(Int_t i=0;i < ncluster;i++){
    tempo[i] = t[i];
    DeltaX[i] = dx[i];
    DeltaZ[i] = dz[i];
    ChannelTOF[i] = ch[i];
  }

}

void SimulateTime(){
  Int_t detId[5];
  Int_t detId0[5];
  AliTOFGeometry::GetVolumeIndices(ChannelTOF[0],detId0);
  Float_t dxtt,dztt,tw;
  Float_t dxttFB,dzttFB;
  Float_t delayx,delayz,delay;
  Float_t resX,resZ,res;
  Float_t avalanche = gRandom->Gaus(0,67);

  for(Int_t i=0;i < ncluster;i++){
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[i],detId);

    dxtt = dxt;
    dztt = dzt;

    if(detId0[4] != detId[4]) dxtt += 2.5*(detId[4] - detId0[4]);
    if(detId0[3] != detId[3]) dztt += 3.5*(detId[3] - detId0[3]);

    if(detId[3]==1) dztt*=-1;

    if(dxtt > 1.25 || dxtt < -1.25) dxtt=1.25;
    if(dztt > 1.75) dztt=1.75;
    else if(dztt < -1.75) dztt=-1.75;


    // time walk
    tw = TMath::Sqrt(dxtt*dxtt + (dztt+1.75)*(dztt+1.75))*41 - 71;

    dxttFB = 1.25 - TMath::Abs(dxtt);
    dzttFB = 1.75 - TMath::Abs(dztt);

    if(dxttFB > 1) dxttFB = 1;
    if(dzttFB > 1) dzttFB = 1;

    // delay
    delayx = 90*(1-dxttFB)*(1-dxttFB)*(1-dxttFB);
    delayz = 90*(1-dzttFB)*(1-dzttFB)*(1-dzttFB);
    delay = TMath::Max(delayx,delayz)-30;

    // resolution degradation
    resX = (1-dxttFB)*80;
    resZ = (1-dzttFB)*80;
    if(resX > resZ) res=resX;
    else res=resZ;

    tempo[i] = gtime + tw + delay + gRandom->Gaus(0,res)+ avalanche;

  }
}
