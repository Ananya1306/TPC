#if 0
#include "/osx_sphenix/coresoftware/simulation/g4simulation/g4main/PHG4HitContainer.h"
#include "/osx_sphenix/coresoftware/simulation/g4simulation/g4main/PHG4Hit.h"
#else
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#endif
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libphg4hit.so)

bool IsOverFrame(double r, double phi){
  //these parameters are taken from Feb 12 drawings of frames.
  double tpc_frame_side_gap=0.8;//mm //space between radial line and start of frame
  double tpc_frame_side_width=2.6;//mm //thickness of frame
  double tpc_margin=0.0;//mm // extra gap between edge of frame and start of GEM holes
  
  double tpc_frame_r3_outer=758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner=583.5;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r2_outer=574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner=411.4;//mm outer edge of smaller-r frame of r3
 
  double tpc_frame_r1_outer=402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner=221.0;//mm outer edge of smaller-r frame of r3
 
  double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  //if the coordinate is in the radial spaces of the frames, return true:
  if (r<tpc_frame_r1_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r1_outer-tpc_margin  && r<tpc_frame_r2_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r2_outer-tpc_margin  && r<tpc_frame_r3_inner+tpc_margin)
    return true;
  if (r>tpc_frame_r3_outer-tpc_margin)
    return true;

  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:
  float sectorangle=(TMath::Pi()/6);
  float nsectors=phi/sectorangle;
  int nsec=floor(nsectors);
  float reduced_phi=phi-nsec*sectorangle; //between zero and sixty degrees.
  float dist_to_previous=r*sin(reduced_phi);
  float dist_to_next=r*sin(sectorangle-reduced_phi);
  if (dist_to_previous<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  if (dist_to_next<tpc_frame_side_gap+tpc_frame_side_width+tpc_margin)
    return true;
  
  return false;
}

void CreateSpacechargeHist(const char *dirname, const char *filename, int istart=0, int maxend=0, int freqKhz=22, int BeamStarting=0, int beamSubset = 40000, bool saveTree=false){
  printf("are you running with the sphenix env?  This probably doesn't work without that!\n");

  gSystem->Load("libg4testbench.so");
  gSystem->Load("libphg4hit.so");

  
  // TFile *file=TFile::Open("/osx_sphenix/G4Hits_sHijing_9-11fm_09900_10000.root","R");
  TFile *file=TFile::Open(Form("%s%s",dirname,filename),"R");
  assert(file.IsOpen());
  TTree *T=(TTree*)file->Get("T");
  int neve=T->GetEntries();

  
  PHG4HitContainer *eleHits=new PHG4HitContainer();
  printf("eleHits pointer is %p\n",(void *)eleHits);
  //older data:  T->SetBranchAddress("DST.TPC.G4HIT_TPC",&eleHits);
  //newer data:
  T->SetBranchAddress("DST#TPC#G4HIT_TPC",&eleHits);
  printf("after set-address, eleHits pointer is %p\n",(void *)eleHits);




  float us=1.0,ms=us*1e3,s=ms*1e3;
  float um=1e-4, mm=um*1e3, cm=mm*10,m=mm*1e3; //changed to make 'cm' 1.0, for convenience.
  float Hz=1/s, kHz=1/ms, MHz=1/us;
  float V=1;
  //used two ways:  1) to apply units to variables when defined
  //                2) to divide by certain units so that those variables are expressed in those units.

  float ionMobility=3.37*cm*cm/V/s;
  float vIon=ionMobility*400*V/cm;
  //float vIon=16.0*um/us;
  float ampGain=2e3;
  float ampIBFfrac=0.02;
  float ionsPerEle=ampGain*ampIBFfrac;
  float mbRate=freqKhz*kHz;
  float xingRate = 9.383*MHz;
  float mean = mbRate/xingRate;
  float z_rdo=105.5*cm;
  float rmin=20*cm;
  float rmax=78*cm;

  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm
  double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6; //electrons per gev.
 
   
  TFile *outfile = TFile::Open(Form("%s_%dkHz_%d.rcc_sc.hist.root",filename,freqKhz,BeamStarting),"RECREATE");
  int nr=159;
  int nphi=360;
  int nz=62*2;
  TH3D *hCharge=new TH3D("sphenix_minbias_charge","SC (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,nz,0,z_rdo/cm);
  TH3D *hPrimary=new TH3D("sphenix_minbias_primary","Primary (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,nz,0,z_rdo/cm);
  TH3D *hIBF=new TH3D("sphenix_minbias_IBF","IBF (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,nz,0,z_rdo/cm);
  TH3D *hPrimaryNoDrift=new TH3D("sphenix_minbias_raw","Undrifted Primary (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,nz,0,z_rdo/cm);

  
  //TH3D *hChargeBack=new TH3D("sphenix_minbias_charge_backward","SC (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,110,0,z_rdo/cm);
  double hrstep=(rmax-rmin)/cm/nr;
  double hphistep=6.28319/nphi;
  double hzstep=z_rdo/cm/nz;
  TTree *rawHits;


  //to fully populate the detector half, we need events covering t0=0 to t0=z_rdo/vIon.  These occur at mbRate, so:
  printf("vIon=%f cm/s\tz=%f cm\t rate=%f kHz\n ==> need z/v*r=%f events to cover the detector. T has %d entries.\n",vIon/(cm/s),z_rdo/(cm),mbRate/(kHz),(z_rdo/vIon*mbRate), neve);
  printf("IBF per measured e: %f , e per GeV deposited: %f\n", ionsPerEle,Tpc_ElectronsPerGeV);
  
  float x,y,z,zibf,zprim;
  float r,phi;
  float ne;

  float driftedZ;//distance these particular particles have drifted
  //float driftedZtile=tilesize*vIon/mbRate;  //drift distance between tile repetitions.

  int testi=5;
  int i;
  if (saveTree){
    rawHits=new TTree("hTree","tpc hit tree for ionization");
    rawHits->Branch("x",&x);
    rawHits->Branch("y",&y);
    rawHits->Branch("zorig",&z);
    rawHits->Branch("zibf",&zibf);
    rawHits->Branch("zprim",&zprim);
    rawHits->Branch("r",&r);
    rawHits->Branch("phi",&phi);
    rawHits->Branch("ev",&i);
    rawHits->Branch("ne",&ne);
  }

  
  T->GetEntry(testi);
  printf("after GetEntry(%d), eleHits pointer is %p\n",testi,(void *)eleHits);

  //return;
  
  /* for (i=1;i<neve;i++){
    T->GetEntry(i);
    //eventually I ought to roll a random number for each bunch crossing and use that to detemrine the number of events, but for now I just take one per minbiasRate.
    //must be sure to convert the incoming positions into local units.
    //load prim
    printf("loading eve=%d\n",i);
    float t0=(i+istart)/mbRate;
  */
  //float poisson = TRandom3::Poisson(mean);

  cout<<"mean = mbRate/xingRate = "<<mean<<endl;

  TRandom3 *rand = new TRandom3();
  //rand->Poisson(mean);
  cout<<"Get Seed = "<<rand->GetSeed()<<endl;

  int eventnum = 0;

  int startingEve = 0;
  int stoppingEve = 0;

  int ourStart = istart;
  int ourEnd = ourStart + neve;

  int beamEnd = BeamStarting + beamSubset;

  int nBeams = z_rdo/(vIon/xingRate); 

  PHG4HitContainer::ConstRange range;

   float f=0.5;//for now, just pick the middle of the hit.  Do better later.
   int bin;   double hr, vol;
   bool overFrame;

   //testing something
   unsigned long now, start;
   int nticks = 0;
   int tock = 100;
   start = gSystem->Now();

   //for(int i=0; i<subsets; i++){

   //rand->SetSeed(beamXing[i]);
   rand->SetSeed(BeamStarting);
   cout<<"Seed = "<<rand->GetSeed()<<endl;

   //TFile *outfile = TFile::Open(Form("%s_%dkHz.rcc_sc_beamxingtest%s.hist.root",filename,freqKhz,freqName[i].c_str()),"RECREATE");

      for(int beamxing = BeamStarting; beamxing<beamEnd&&beamxing<nBeams; beamxing++){
   // for(int beamxing = beamXing[i]; beamxing<beamXing[i+1]; beamxing++){
   //   for(int beamxing = 730000; beamxing<nBeams; beamxing++){
    //    printf("loading beamxing, eventnum = %d , %d\n ", beamxing, eventnum);
   

    float t0 = (beamxing/xingRate); //units in microseconds
    driftedZ=t0*vIon;//drift position in local units

    nticks++;
    if(nticks>tock){
      nticks = 0;
      now = gSystem->Now();
      printf("time since last tock is %lu\n",(unsigned long)(now-start));
      start = now;
    }


    int nEventsHere = rand->Poisson(mean); //This generates a double, not an int
  
    //    cout<<"nEventsHere = "<<nEventsHere<<endl;
    //cout<<"maxend = "<<maxend<<endl;
    startingEve = stoppingEve%maxend;
    stoppingEve = startingEve + nEventsHere;
    // eventnum++;
    printf("loading beamxing %d of %d xings,  eventnum = %d, startingEve = %d, stoppingEve = %d, nEventsHere = %d\n",beamxing, nBeams,eventnum,startingEve,stoppingEve,nEventsHere);
    //printf("loading beamxing %d of %d xings\n", beamxing, nBeams);

    if(startingEve>=ourStart && startingEve<ourEnd){    
  
  //place the events starting at the first unused, and continuing until we have placed nEventsHere here.
      //   cout<<"if statement satisfied "<<" ourStart = "<<ourStart<<" ourEnd = "<<ourEnd<<endl;
    for(eventnum=startingEve; eventnum<stoppingEve; eventnum++){
      // printf("Reading eve num from %d from file name %s\n",eventnum%neve, filename);
      //printf("Reading evt num %d",eventnum);

      //  T->GetEntry(23);

     
       T->GetEntry(eventnum%neve); //the '%' guarantees we wrap around instead of asking for an event the file doesnt have.

      range=eleHits->getHits();
      //cout<<"range"<<endl;
      assert(range);
      //cout<<"assert range works!"<<endl;


      //continue;
      
      /*
      range=NULL;
      range=eleHits->getHits();
      if(range == NULL){
	cout<<"No range object for this event!!!"<<endl;
	assert(false);
      }


    assert(range);
   
      */

    //continue;
      /*
      cout<<"Checking hiter for range.first"<<endl;
      cout<<"ne = "<<range.first->second->get_eion()<<endl;
      cout<<"x = "<<range.first->second->get_x(0)<<endl;
      cout<<"y = "<<range.first->second->get_y(0)<<endl;
      cout<<"z = "<<range.first->second->get_z(0)<<endl;
      */
    for (PHG4HitContainer::ConstIterator hiter=range.first;hiter!=range.second;hiter++) {
      ne=hiter->second->get_eion()*Tpc_ElectronsPerGeV;
      //load the three coordinates in with units of cm.
      x = (hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0)))*(cm);
      y = (hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0)))*(cm);
      z = (hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0)))*(cm);

      if (z<0) continue;
      r=sqrt(x*x+y*y);
      phi=atan2(x,y);
      zprim=z-driftedZ;
      zibf=z_rdo-driftedZ;
      
      if (phi<0) phi+=6.28319;
      //compute the bin volume:
      bin=hCharge->GetYaxis()->FindBin(r/(cm));
      hr=hCharge->GetYaxis()->GetBinLowEdge(bin);
      vol=(hzstep*hphistep*(hr+hrstep*0.5)*hrstep)/cm/cm/cm;
      overFrame=IsOverFrame(r/(mm),phi);
      //cout<<"values in iterator loop filled"<<endl;
      hCharge->Fill(phi,r/(cm),zprim/(cm),ne/vol); //primary ion, drifted by t0, in cm
      // printf("filling in the histograms %d %d");
      //cout<<"Filling in hCharge histogram = "<<phi<<"\t"<<r/(cm)<<"\t"<<zprim/(cm)<<"\t"<<ne/vol<<endl;
      if (!overFrame) {
	hCharge->Fill(phi,r/(cm),zibf/(cm),ne*ionsPerEle/vol); //amp ion, drifted by t0, in cm
	hIBF->Fill(phi,r/(cm),zibf/(cm),ne*ionsPerEle/vol);
      }
      hPrimary->Fill(phi,r/(cm),zprim/(cm),ne/vol);
      hPrimaryNoDrift->Fill(phi,r/(cm),z/(cm),ne/vol);

      
      if (saveTree){
	rawHits->Fill();
	
      }
      // cout<<"histograms filles"<<endl;
      //     cout<<"for loop for iterator works 2"<<endl;
    }
    //cout<<"for eventnum loop works"<<endl;
    }
    //cout<<"if startingEve loop works"<<endl;
    }
    // cout<<"beamxing completed"<<endl;
      }
    
      cout<<"Writing out the files"<<endl;
      

outfile->cd();
  
//     outfile1->cd();
  hCharge->Write();
  hPrimary->Write();
  hIBF->Write();
  hPrimaryNoDrift->Write();
  if (saveTree){
    rawHits->Write();
  }
  outfile->Close();
   
   //outfile->Close();
  // outfile1->Close();
  return;
}
