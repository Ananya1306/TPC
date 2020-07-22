#include <iostream>
#include <fstream>
#include <string>
#include <math.h>


using namespace std;

void MinBiasCharge() 

{

  const string freq[10] = {"00000_00050","00050_00100","00100_00150","00150_00200","00200_00250","00250_00300","00300_00350","00350_00400","00400_00450","00450_00500"};

  TH3D *hCharge[10];
  TFile *infile[10];
  char fileName[100];
  char histName[20];
  char histTitle[20];

  for(int i=0; i<10; i++){
    //sprintf(histName, "hCharge_%d",i);
    //sprintf(histTitle, "hChargeHist_%d",i);
    sprintf(fileName,"G4Hits_sHijing_0-12fm_%s.root_200kHz.rcc_sc.hist.root",freq[i].c_str());
    //hCharge[i] = new TH3D(histName, histTitle, 360,0,6.28319,159,20,78,124,0,105.5);
    infile[i] = new TFile(fileName);
    hCharge[i] = (TH3D*)infile[i]->Get("sphenix_minbias_charge");
    hCharge[i]->Project3D("xy")->Draw("colz");
  }
  /*
  TCanvas *c1 = new TCanvas("c1","ProjectXY");
  c1->Divide(2,2);

  TH3D *hm1 = new TH3D("hm1","hm1",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile1 = new TFile("G4Hits_sHijing_0-12fm_00000_00050.root_200kHz.rcc_sc.hist.root"); 
  hm1 = (TH3D*)hfile1->Get("sphenix_minbias_charge");
  c1->cd(1);
  hm1->Project3D("xy")->Draw("colz");
  
TH3D *hm2 = new TH3D("hm2","hm2",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile2 = new TFile("G4Hits_sHijing_0-12fm_00050_00100.root_200kHz.rcc_sc.hist.root"); 
  hm2 = (TH3D*)hfile2->Get("sphenix_minbias_charge");
  c1->cd(2);
  hm2->Project3D("xy")->Draw("colz");
  
TH3D *hm3 = new TH3D("hm3","hm3",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile3 = new TFile("G4Hits_sHijing_0-12fm_00100_00150.root_200kHz.rcc_sc.hist.root"); 
  hm3 = (TH3D*)hfile3->Get("sphenix_minbias_charge");
  c1->cd(3);
  hm3->Project3D("xy")->Draw("colz");

TH3D *hm4 = new TH3D("hm4","hm4",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile4 = new TFile("G4Hits_sHijing_0-12fm_00150_00200.root_200kHz.rcc_sc.hist.root"); 
  hm4 = (TH3D*)hfile4->Get("sphenix_minbias_charge");
  c1->cd(4);
  hm4->Project3D("xy")->Draw("colz");
  
TH3D *hm5 = new TH3D("hm5","hm5",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile5 = new TFile("G4Hits_sHijing_0-12fm_00200_00250.root_200kHz.rcc_sc.hist.root"); 
  hm5 = (TH3D*)hfile5->Get("sphenix_minbias_charge");
  c1->cd(5);
  hm5->Project3D("xy")->Draw("colz");

TH3D *hm6 = new TH3D("hm6","hm6",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile6 = new TFile("G4Hits_sHijing_0-12fm_00250_00300.root_200kHz.rcc_sc.hist.root"); 
  hm6 = (TH3D*)hfile6->Get("sphenix_minbias_charge");
  c1->cd(6);
  hm6->Project3D("xy")->Draw("colz");

TH3D *hm7 = new TH3D("hm7","hm7",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile7 = new TFile("G4Hits_sHijing_0-12fm_00300_00350.root_200kHz.rcc_sc.hist.root"); 
  hm7 = (TH3D*)hfile7->Get("sphenix_minbias_charge");
  c1->cd(7);
  hm7->Project3D("xy")->Draw("colz");

TH3D *hm8 = new TH3D("hm8","hm8",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile8 = new TFile("G4Hits_sHijing_0-12fm_00350_00400.root_200kHz.rcc_sc.hist.root"); 
  hm8 = (TH3D*)hfile8->Get("sphenix_minbias_charge");
  c1->cd(8);
  hm8->Project3D("xy")->Draw("colz");

TH3D *hm9 = new TH3D("hm9","hm9",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile9 = new TFile("G4Hits_sHijing_0-12fm_00400_00450.root_200kHz.rcc_sc.hist.root"); 
  hm9 = (TH3D*)hfile9->Get("sphenix_minbias_charge");
  c1->cd(9);
  hm9->Project3D("xy")->Draw("colz");

TH3D *hm10 = new TH3D("hm10","hm10",360,0,6.28319,159,20,78,124,0,105.5);
  TFile *hfile10 = new TFile("G4Hits_sHijing_0-12fm_00450_00500.root_200kHz.rcc_sc.hist.root"); 
  hm10 = (TH3D*)hfile10->Get("sphenix_minbias_charge");
  c1->cd(10);
  hm10->Project3D("xy")->Draw("colz");

  */

}
