// #include "stdio.h"
// #include "TH1F.h"
// #include "TH2F.h"
// #include "TTree.h"
// #include "TFile.h"
// #include "TLeaf.h"
// #include "AliRunLoader.h"
// #include "AliESDtrack.h"
// #include "AliESDEvent.h"
// #include "AliHeader.h"
// #include "AliESDTOFCluster.h"
// #include "AliTOFGeometry.h"
// #include "AliPIDResponse.h"
// #include "TAlienCollection.h"
// #include "TGrid.h"
// #include "TRandom.h"
// #include "AliRun.h"
// #include "AliStack.h"
// #include "TParticle.h"
// #include "AliGenEventHeader.h"

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
Float_t tot[100];
Float_t impulso_trasv;
Float_t impulso;
Float_t res[3];
Int_t charge;
Float_t phi,eta;
Float_t secAngle;
Float_t cval[5],gtime;
Float_t dxt,dzt,xgl,ygl,zgl;
Int_t mism;
Int_t ntofcl;
Int_t TOFout;

Float_t dedx,StartTime,StartTimeRes;
Float_t interactiontime;
Int_t itrig;
Float_t timetrig;
Float_t smearX = 0.;
Float_t smearZ = 0.;
Int_t pdg;

Bool_t isMC = kTRUE;

// funzione che riempi gli istogrammi
Bool_t CheckSingle(const char* esdFileName = "AliESDs.root",Bool_t kGRID=1); // il vecchio CheckESD

void GetResolutionAtTOF(AliESDtrack *t,Float_t mag,Int_t tofc,Float_t res[3]);

void MakeTrueRes();
void AddDelay();
void ReMatch();

// macro principale che fa il loop sugli eventi e scrive il file
Bool_t CheckESD(const char *lista="wn.xml",Bool_t kGRID=1) //per prendere dalla grid;
//Bool_t CheckESD(const char *lista="lista",Bool_t kGRID=0) // da locale
{
  char name[300];
  Int_t ifile = 0;
  Int_t nmaxfile = 200000; // to limit the number of ESD files
  
  //T->Branch("nevento",&nevento,"nevento/I");
  T->Branch("kTOFout",&TOFout,"kTOFout/I");
  T->Branch("ncluster",&ncluster,"ncluster/I");
  T->Branch("tempo",tempo,"tempo[ncluster]/F");
  if(isMC) T->Branch("pdg",&pdg,"pdg/I");
  if(isMC) T->Branch("gtime",&gtime,"gtime/F");
  if(isMC) T->Branch("DXtrue",&dxt,"DXtrue/F");
  if(isMC) T->Branch("DZtrue",&dzt,"DZtrue/F");
  T->Branch("DeltaX",DeltaX,"DeltaX[ncluster]/F");
  T->Branch("DeltaZ",DeltaZ,"DeltaZ[ncluster]/F");
  T->Branch("exp_time_el",exp_time_el,"exp_time_el[ncluster]/F");
  T->Branch("exp_time_mu",exp_time_mu,"exp_time_mu[ncluster]/F");
  T->Branch("exp_time_pi",exp_time_pi,"exp_time_pi[ncluster]/F");
  T->Branch("exp_time_ka",exp_time_ka,"exp_time_ka[ncluster]/F");
  T->Branch("exp_time_pr",exp_time_pr,"exp_time_pr[ncluster]/F");
  T->Branch("L",L,"L[ncluster]/F");
  T->Branch("TOT",tot,"TOT[ncluster]/F");
  T->Branch("ChannelTOF",ChannelTOF,"ChannelTOF[ncluster]/I");
  T->Branch("impulso",&impulso,"impulso/F");
  T->Branch("impulso_trasv",&impulso_trasv,"impulso_trasv/F");
  T->Branch("res",res,"res[3]/F");
  T->Branch("charge",&charge,"charge/I");
  T->Branch("phi",&phi,"phi/F");
  T->Branch("eta",&eta,"eta/F");
  T->Branch("secPhi",&secAngle,"secPhi/F");
  T->Branch("cval",cval,"cval[5]/F");
  T->Branch("mism",&mism,"mism/I");
  T->Branch("ntofcl",&ntofcl,"ntofcl/I");
  T->Branch("dedx",&dedx,"dedx/F");
  T->Branch("StartTime",&StartTime,"StartTime/F");
  T->Branch("StartTimeRes",&StartTimeRes,"StartTimeRes/F");
  T->Branch("rTOF",&rTOFused,"rTOF/F");
  if(isMC) T->Branch("InteractionTime",&interactiontime,"InteractionTime/F");
  T->Branch("tempoTrig",&timetrig,"tempoTrig/F");
  
  if(! kGRID){
    FILE *fin = fopen (lista,"r");
    
    while(fscanf(fin,"%s",name)==1 && ifile < nmaxfile){
      CheckSingle(name,kGRID),ifile++;
      printf("file %i done\n",ifile);
      system("date");
    }
    
    fclose(fin);
  }
  else{
    TGrid::Connect("alien://");
    
    TAlienCollection *myCollection = (TAlienCollection *) TAlienCollection::Open(lista);
    if (!myCollection)
      {
	Error("CheckESD.C", Form("Cannot create an AliEn collection from %s", lista));
	return 1;
      }
    myCollection->Reset();    
    
    while (myCollection->Next() && ifile < nmaxfile){
      CheckSingle(myCollection->GetTURL("")),ifile++;
      printf("file %i done\n",ifile);
      system("date");
    }
  }

  fotree->cd();
  T->Write(); //scrivo tree
  fotree->Close();
}

// analisi
Bool_t CheckSingle(const char* esdFileName,Bool_t kGRID){
  //inizializzo a zero ncluster (di tree T)
  //for (int ifc=0;ifc<10000;ifc++) ncluster[ifc]=0;
  
  // check the content of the ESD
  
  AliPIDResponse *pidr = new AliPIDResponse();
  
  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()){
    Error("CheckESD", "opening ESD file %s failed", esdFileName);
    return kFALSE;
  }
  
  TString mctrkref(esdFileName);
  mctrkref.ReplaceAll("AliESDs.root","TrackRefs.root");
  TString fgal(esdFileName);
  fgal.ReplaceAll("AliESDs.root","galice.root");
  
  if(kGRID){
    fgal.Insert(0,"alien://");
    mctrkref.Insert(0,"alien://");
  }
  
  TTree *trkref;
  
  printf("ESD = %s\n",esdFileName);
  
  TFile *ftrkref; 
  if(isMC) ftrkref = TFile::Open(mctrkref.Data());
  
  AliHeader *h = new AliHeader();
  
  TFile *fgalice;
  if(isMC) fgalice = TFile::Open(fgal.Data());
  TTree *tgalice;
  if(isMC){
    tgalice = (TTree *) fgalice->Get("TE");
    tgalice->SetBranchAddress("Header",&h);
  }
  
  AliRunLoader* runLoader = NULL;
  
  AliRun *gAlice;
  if(isMC) runLoader = AliRunLoader::Open(fgal.Data());
  if(runLoader){
    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();
    if (!gAlice) {
      Error("CheckESD", "no galice object found");
      return kFALSE;
    }
    runLoader->LoadKinematics();
    runLoader->LoadHeader();
  }
  
  AliESDEvent * esd = new AliESDEvent;
  //  printf("esd object = %x\n",esd);
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree){
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  esd->ReadFromTree(tree); // crea link tra esd e tree
  
  TClonesArray* tofcl;  // array dinamico
  TClonesArray* tofhit;
  TClonesArray* tofmatch;
  
  Int_t nev = tree->GetEntries(); //ogni entries evento
  Float_t mag;
  
  printf("nev = %i\n",nev);
  
  //azzero il contatore delle tracce del TTree T
  //ntracks=0;
  AliStack* stack=NULL;
  
  Int_t trkassociation[1000000];
  
  for(Int_t ie=0;ie < nev;ie++){
    if(runLoader){
      runLoader->GetEvent(ie);
      
      // select simulated primary particles, V0s and cascades
      stack = runLoader->Stack();
    }
    
    if(isMC) trkref = (TTree *) ftrkref->Get(Form("Event%i/TreeTR",ie));
    tree->GetEvent(ie);
    if(isMC) tgalice->GetEvent(ie);
    
    if(isMC) interactiontime = h->GenEventHeader()->InteractionTime()*1E+12;
    
    mag = esd->GetMagneticField();
    
    AliTOFHeader *tofh = esd->GetTOFHeader();
    ntofcl = tofh->GetNumberOfTOFclusters();
    
    esd->ConnectTracks(); // Deve essere sempre chiamato dopo aver letto l'evento (non troverebbe l'ESDevent). Scrivo in tutte le tracce l origine dell evento così poi da arrivare ovunque(tipo al cluster e al tempo quindi).
    
    
    //Riempio variabile del tree "T"
    //nevento=ie;
    
    if(! esd->GetVertex()){
      esd->ResetStdContent();
      continue;// una volta fatto il connect manda un flag ; siccome qua c'era un continue(non si arriva in fondo al ciclo) bisogna resettarlo altrimenti lo trova già attivo.
    }
    
    tofcl = esd->GetESDTOFClusters(); // AliESDTOFCluster *cltof = tofcl->At(i);
    if(tofcl->GetEntries() == 0){
      esd->ResetStdContent();
      continue;
    }
    tofhit = esd->GetESDTOFHits(); // AliESDTOFHit *hittof = tofhit->At(i);
    tofmatch = esd->GetESDTOFMatches(); // AliESDTOFHit *mathctof = tofmatch->At(i);
    
    // loop over tracks
    
    pidr->SetTOFResponse(esd,AliPIDResponse::kTOF_T0); //per recuperare lo start time ("esd", "tipo start time"), tipo cioè o il TOF stesso o il T0 o il best, ovvero la combinazione dei 2
    
    Int_t ntrk = esd->GetNumberOfTracks();
    
    //printf("%i) TPC tracks = %i -- TOF cluster = %i - TOF hit = %i -- matchable info = %i\n",ie,ntrk,tofcl->GetEntries(),tofhit->GetEntries(),tofmatch->GetEntries());
    
    Double_t time[AliPID::kSPECIESC];
    
    
    if(isMC && stack){// create association trackref
      printf("nMC track = %i\n",stack->GetNtrack());
      for(Int_t ist=0;ist < stack->GetNtrack();ist++){
	trkassociation[ist]=-1;
      }
      for(Int_t iref=0;iref < trkref->GetEntries();iref++){
	trkref->GetEvent(iref);
	Int_t trkreference = trkref->GetLeaf("TrackReferences.fTrack")->GetValue();
	if(trkreference > -1 && trkreference < 1000000){
	  trkassociation[trkreference] = iref;
	}
      }
    }
    
    for (Int_t iTrack = 0; iTrack < ntrk; iTrack++){
      AliESDtrack* track = esd->GetTrack(iTrack);
      
      // select tracks of selected particles
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;//almeno un hit nell ITS
      if (track->GetConstrainedChi2() > 1e9) continue; //se brutto X^2
      if ((track->GetStatus() & AliESDtrack::kTOFout) == 0) continue; //se traccia matchata con tof
      
      Float_t p =track->P();
      
      itrig = 0;
      timetrig = 0;
      
      if(p > 0.9 && p < 1.1){
 	track->GetIntegratedTimes(time);
	
	itrig = iTrack;
	timetrig = track->GetTOFsignal() - time[2];
	iTrack = ntrk;
      }
    }
    
    printf("real loop, ntrk = %i\n",ntrk);
  
    for (Int_t iTrack = 0; iTrack < ntrk; iTrack++){
      AliESDtrack* track = esd->GetTrack(iTrack);
      
      // select tracks of selected particles
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;//almeno un hit nell ITS
      if (track->GetConstrainedChi2() > 1e9) continue; //se brutto X^2
      
      //      if ((track->GetStatus() & AliESDtrack::kTOFout) == 0) continue; //se traccia matchata con tof
      
      TOFout = (track->GetStatus() & AliESDtrack::kTOFout) > 0;

      track->GetIntegratedTimes(time);
      
      Float_t dx = track->GetTOFsignalDx(); //leggo i residui tra traccia e canale tof acceso
      Float_t dz = track->GetTOFsignalDz();
      
      mism = 0;
      
      dedx = track->GetTPCsignal();
      
      Int_t label = TMath::Abs(track->GetLabel());
      if(stack){
	TParticle *part=stack->Particle(label);
	pdg = part->GetPdgCode();
      }
      
      Int_t TOFlabel[3];
      track->GetTOFLabel(TOFlabel);
      
      //  printf("%i %i %i %i\n",label,TOFlabel[0],TOFlabel[1],TOFlabel[2]);
      
      ChannelTOF[0] = track->GetTOFCalChannel();
      //      printf("geant time = %f\n",gtime);
      //getchar();
      // if(TMath::Abs(dx) > 1.25 || TMath::Abs(dz) > 1.75) continue; // is inside the pad
      
      //riempio il numro di cludter e impulso trasverso per traccia del TTree T
      ncluster=track->GetNTOFclusters();
      impulso_trasv=track->Pt();
      impulso=track->P();
      
      StartTime = pidr->GetTOFResponse().GetStartTime(track->P());
      StartTimeRes = pidr->GetTOFResponse().GetStartTimeRes(track->P());
      
      if(track->Pt() > 0.9 && track->Pt() < 1.5){  //impulso non troppo alto per separazione tra particelle
	Float_t dt = track->GetTOFsignal() - time[2] - pidr->GetTOFResponse().GetStartTime(track->P());//tempo TOF(è lo stesso di Gettime, solo che lo prendo dale tracce)(già calibrato) -ip del PI (posizione 0 e, posizione 1 mu, pos 2 PI, pos 3 K,pos 4 p) -start time
	Float_t dtKa = track->GetTOFsignal() - time[3] - pidr->GetTOFResponse().GetStartTime(track->P());
	Float_t dtPr = track->GetTOFsignal() - time[4] - pidr->GetTOFResponse().GetStartTime(track->P());
	hdt->Fill(dt);
	hdtKa->Fill(dtKa);
	hdtPr->Fill(dtPr);
      }
      
      for (int i=0;i<(track->GetNTOFclusters());i++){
	int idummy=track->GetTOFclusterArray()[i];
        
	AliESDTOFCluster *cl = (AliESDTOFCluster *) tofcl->At(idummy);
        
	tempo[i]=cl->GetTime();
	tot[i]=cl->GetTOT();
        
	ChannelTOF[i]=cl->GetTOFchannel();
	
	charge = track->Charge();
	phi = track->Phi();
	eta = track->Eta();
	
	if(i==0){
	  GetResolutionAtTOF(track,mag,ChannelTOF[i],res);
	}
	
	for(int im=cl->GetNMatchableTracks();im--;){ //o così o da n-1 a 0 //for(int im=cl->GetNMatchableTracks();im>0;im--) non andava bene perchè non prendeva mai lo 0  
	  
	  //	    if(track->GetNTOFclusters()==2) printf("-- %i) %f %f\n",im,cl->GetLength(im),cl->GetIntegratedTime(2,im));
	  
	  if(cl->GetTrackIndex(im) == track->GetID()){
	    exp_time_el[i] = cl->GetIntegratedTime(0,im); // pi = 2
	    exp_time_mu[i] = cl->GetIntegratedTime(1,im); // pi = 2
	    exp_time_pi[i] = cl->GetIntegratedTime(2,im); // pi = 2
	    exp_time_ka[i] = cl->GetIntegratedTime(3,im); // pi = 2
	    exp_time_pr[i] = cl->GetIntegratedTime(4,im); // pi = 2
	    L[i] = cl->GetLength(im);
	    //		  if(track->GetNTOFclusters()==2)printf("%i) %f %f\n",i,L[i],exp_time_pi[i]);
	    DeltaX[i]=cl->GetDx(im); // mettendolo dentro questo if dovrei prendere i residui di una stessa traccia
	    DeltaZ[i]=cl->GetDz(im);
	  }
	}
      }
         
      //ReMatch();
      
      Int_t jref=0;
      if(isMC){
	if(TOFlabel[0] > -1 && TOFlabel[0] < 1000000){
	  trkref->GetEvent(trkassociation[TOFlabel[0]]);
	  if(TOFlabel[0] == trkref->GetLeaf("TrackReferences.fTrack")->GetValue()){
	    //  printf("trk -> %i (%i)\n",trkref->GetLeaf("TrackReferences.fTrack")->GetValue(),trkref->GetLeaf("TrackReferences.fTrack")->GetValue(jref));	
	    while(jref > -1 && trkref->GetLeaf("TrackReferences.fTrack")->GetValue(jref) != 0){
	      //printf("det = %i\n",trkref->GetLeaf("TrackReferences.fDetectorId")->GetValue(jref));
	      if(trkref->GetLeaf("TrackReferences.fDetectorId")->GetValue(jref) == 4){
		gtime=trkref->GetLeaf("TrackReferences.fTime")->GetValue(jref)*1E+12;
		xgl = trkref->GetLeaf("TrackReferences.fX")->GetValue(jref);
		ygl = trkref->GetLeaf("TrackReferences.fY")->GetValue(jref);
		zgl = trkref->GetLeaf("TrackReferences.fZ")->GetValue(jref);
		MakeTrueRes();
		jref =  100;
	      }
	      jref++;
	    }
	  }
	}
      }
      
      
      if(TMath::Abs(label) != TOFlabel[0] && stack){
	mism=2;
	
	while(TOFlabel[0] != -1 && TOFlabel[0] != label){
	  TOFlabel[0] = stack->Particle(TOFlabel[0])->GetMother(0);
	}
	
	if(label == TOFlabel[0])
	  mism=1;	
	
      }
      
      //AddDelay();
      T->Fill(); //cout<<"riempio il tree  "<<endl; //Riempio tree "T"
    
      
      
      //incremento il contatore delle tracce del TTree T matchate e che superano i tagli
      //ntracks++;
      
    }//end of for(tracks)
      
      
    
    esd->ResetStdContent();
    
    
    
  } //end of for(events)

  if(runLoader){
    runLoader->UnloadHeader();
    runLoader->UnloadKinematics();
    delete runLoader;
  }
  
  esdFile->Close();
  if(isMC) ftrkref->Close();
  if(isMC) fgalice->Close();
}

void GetResolutionAtTOF(AliESDtrack *t,Float_t mag,Int_t tofc,Float_t res[3]){
  Int_t detId[5];
  AliTOFGeometry::GetVolumeIndices(tofc,detId);
  
  Float_t xTOFch[3];
  xTOFch[0] = AliTOFGeometry::GetX(detId);
  xTOFch[1] = AliTOFGeometry::GetY(detId);
  xTOFch[2] = AliTOFGeometry::GetZ(detId);
  
  Float_t rTOF = TMath::Sqrt(xTOFch[0]*xTOFch[0] + xTOFch[1]*xTOFch[1]);
  
  rTOFused = rTOF;
  
  secAngle = (Int_t (TMath::ATan2(xTOFch[1],xTOFch[0])/TMath::Pi()*9)+ 0.5) * TMath::Pi() / 9;
  Float_t thetaAngle = AliTOFGeometry::GetAngles(detId[1]/*iplate*/, detId[2]/*istrip*/)/180*TMath::Pi();

  Float_t versX[3] = {TMath::Sin(secAngle),-TMath::Cos(secAngle),0};
  Float_t versY[3] = {TMath::Cos(thetaAngle)*TMath::Cos(secAngle),TMath::Cos(thetaAngle)*TMath::Sin(secAngle),-TMath::Sin(thetaAngle)};
  Float_t versZ[3] = {TMath::Sin(thetaAngle)*TMath::Cos(secAngle),TMath::Sin(thetaAngle)*TMath::Sin(secAngle),TMath::Cos(thetaAngle)};

  Double_t cv[21];
  const AliExternalTrackParam* trkExt = t->GetOuterParam();
  
  AliExternalTrackParam trk(*trkExt);
  Int_t value = trk.PropagateTo(rTOF,mag);
  //   printf("propto = %i\n",value);
  while(!value && rTOFused > 0){
    rTOFused-=1;
    value = trk.PropagateTo(rTOFused,mag);
    //     printf("re-propto = %i (%f)\n",value,rTOF);
  }
  
  trk.GetCovarianceXYZPxPyPz(cv);
  
  Double_t pos[3];
  trk.GetXYZ(pos);
  Float_t pos2[] = {pos[0],pos[1],pos[2]};
  //  printf("residuals recomputed at %f: %f %f (%f %f)\n",sqrt(pos[0]*pos[0]+pos[1]*pos[1]),AliTOFGeometry::GetPadDx(pos2),AliTOFGeometry::GetPadDz(pos2),t->GetTOFsignalDx(),t->GetTOFsignalDz());
  
  for(Int_t i=0;i<5;i++) cval[i] = cv[i];
  
  Float_t mat[3][3];
  mat[0][0] = cv[0];
  mat[0][1] = cv[1];
  mat[0][2] = cv[3];
  mat[1][0] = cv[1];
  mat[1][1] = cv[2];
  mat[1][2] = cv[4];
  mat[2][0] = cv[3];
  mat[2][1] = cv[4];
  mat[2][2] = cv[5];
  
  res[0] = 0;
  res[1] = 0;
  res[2] = 0;
  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<3;j++){
      res[0] += mat[i][j]*versX[i]*versX[j];
      res[1] += mat[i][j]*versY[i]*versY[j];
      res[2] += mat[i][j]*versZ[i]*versZ[j];
    }
  }
  
  res[0] = TMath::Sqrt(TMath::Abs(res[0]));//cv[0]);
  res[1] = TMath::Sqrt(TMath::Abs(res[1]));//cv[5]);
  res[2] = TMath::Sqrt(TMath::Abs(res[2]));//cv[5]);
  
  res[0] *= rTOF/rTOFused;
  res[1] *= rTOF/rTOFused;
  res[2] *= rTOF/rTOFused;
}

void AddDelay(){
  // xgl,ygl,zgl
  // dxt,dzt
  // ChannelTOF[0]
  Int_t detId[5];
  Int_t detId2[5];
  
  for(Int_t i=0;i < ncluster;i++){
    AliTOFGeometry::GetVolumeIndices(ChannelTOF[i],detId);
    
    //  printf("A) det: %i %i %i %i %i\n",detId[0],detId[1],detId[2],detId[3],detId[4]);
    Float_t pos[3] = {xgl,ygl,zgl};
    
    detId2[0] = AliTOFGeometry::GetSector(pos);
    detId2[1] = AliTOFGeometry::GetPlate(pos);
    detId2[2] = AliTOFGeometry::GetStrip(pos);
    detId2[3] = AliTOFGeometry::GetPadZ(pos);
    detId2[4] = AliTOFGeometry::GetPadX(pos);
    
    //  printf("B) det: %i %i %i %i %i\n",detId[0],detId[1],detId[2],detId[3],detId[4]);
    
    //  printf("truth %f) %f %f %f\n",sqrt(xgl*xgl+ygl*ygl),dxt,dyt,dzt);
    Float_t dxtt = AliTOFGeometry::GetPadDx(pos);
    Float_t dytt = AliTOFGeometry::GetPadDy(pos);
    Float_t dztt = AliTOFGeometry::GetPadDz(pos);
    
    //  printf("truth %f) %f %f %f\n",sqrt(xgl*xgl+ygl*ygl),dxt,dyt,dzt);
    
    if(detId2[4] != detId[4]) dxtt += 2.5*(detId2[4] - detId[4]);
    if(detId2[3] != detId[3]) dztt += 3.5*(detId2[3] - detId[3]);

    if(TMath::Abs(dxtt) > 1.25 || TMath::Abs(dztt) > 1.75) tempo[i] += 60;
    
    if(dxtt < -1.25) dxtt = -1.25;
    else if (dxtt > 1.25) dxtt = 1.25;
    if(dztt < -1.75) dztt = -1.75;
    else if (dztt > 1.75) dztt = 1.75;
    
    dxtt = 1.25 - TMath::Abs(dxtt);
    dztt = 1.75 - TMath::Abs(dztt);
    
    if(dxtt > 1) dxtt = 1;
    if(dztt > 1) dztt = 1;
    
    Float_t delayx = 100*(1-dxtt)*(1-dxtt)*(1-dxtt);
    Float_t delayz = 100*(1-dztt)*(1-dztt)*(1-dztt);
    Float_t delay = TMath::Max(delayx,delayz);
    
    // tempo[i] += delay;
    
  }
}

void MakeTrueRes(){
  // xgl,ygl,zgl
  // dxt,dzt
  // ChannelTOF[0]
  Int_t detId[5];
  Int_t detId2[5];
  AliTOFGeometry::GetVolumeIndices(ChannelTOF[0],detId);
  
  //  printf("A) det: %i %i %i %i %i\n",detId[0],detId[1],detId[2],detId[3],detId[4]);
  Float_t pos[3] = {xgl,ygl,zgl};
  
  detId2[0] = AliTOFGeometry::GetSector(pos);
  detId2[1] = AliTOFGeometry::GetPlate(pos);
  detId2[2] = AliTOFGeometry::GetStrip(pos);
  detId2[3] = AliTOFGeometry::GetPadZ(pos);
  detId2[4] = AliTOFGeometry::GetPadX(pos);
  
  //  printf("B) det: %i %i %i %i %i\n",detId[0],detId[1],detId[2],detId[3],detId[4]);

  dxt = AliTOFGeometry::GetPadDx(pos);
  Float_t dyt = AliTOFGeometry::GetPadDy(pos);
  dzt = AliTOFGeometry::GetPadDz(pos);
  
  //  printf("truth %f) %f %f %f\n",sqrt(xgl*xgl+ygl*ygl),dxt,dyt,dzt);
  
  if(detId2[4] != detId[4]) dxt += 2.5*(detId2[4] - detId[4]);
  if(detId2[3] != detId[3]) dzt += 3.5*(detId2[3] - detId[3]);
  
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
