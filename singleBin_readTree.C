// Read Pythia produced c->e and b->e events from MinBias Pythia p+p 200 GeV events
// Updated Oct. 7, 2015 - Z. Miller

#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TVector2.h"

const float PI = TMath::Pi();
const float PP0[6]={0.849029,0.854208,0.861378,0.860819,0.865031,0.868127};   //PP0,PP1,PP2  3 parameters for efficiency calculation, 6 array is for Cu centrality
const float PP1[6]={0.079329,0.0689266,0.0588677,0.0675269,0.0543948,0.0501144};
const float PP2[6]={0.850396,0.8306,0.798866,0.884421,0.809765,0.784555};
const float pt_trig_up = 3.5;
const float pt_trig_lo = 2.5;             
const float pt_asso_up = 20.5;
const float pt_asso_lo = 0.3;//0.15;
const float deletacut = 0.5;
const float EtaCut = 0.7;
const float hEtaCut = 1.05;
const Float_t Vz_offset = .0;
const Float_t Vz_cut = 30;
const Float_t refmult_binWidth = 25;
const Int_t refmult_binNumber = 5;   //Vz bin 20
const int Phibin = 200;

// Prepare output file and Hists
TFile fout("readTreeOut.root","RECREATE");

//defining histograms
TH1D* hEventTallyce = new TH1D("ceEventTally","ceEvent Tally",10,0,1);
TH1D* hEventTallybe = new TH1D("beEventTally","beEvent Tally",10,0,1);
TH1D* hEventTallybce = new TH1D("bceEventTally","bceEvent Tally",10,0,1);
TH1F *hRefMult = new TH1F("hRefMult","hRefMult",800,0,800);
TH1D *hdPhiRawce = new TH1D("hdPhiRawce","hdPhiRawce",Phibin, -10,10);
TH1D *hdPhiRawceN;
TH1D *hdPhiRawbe = new TH1D("hdPhiRawbe","hdPhiRawbe",Phibin, -10,10);
TH1D *hdPhiRawbeN;
TH1D *hdPhiRawbce = new TH1D("hdPhiRawbce","hdPhiRawbce",Phibin, -10,10);
//        TH1D *hdPhiBg = new TH1D("hdPhiBg","hdPhiBg", Phibin, -PI,PI);
//        TH1D *hdPhiMix = new TH1D("hdPhiMix","hdPhiMix",Phibin, -PI,PI);
TH1D *hdEtaRawce = new TH1D("hdEtaRawce","hdEtaRawce",100, -5,5);
TH1D *hdEtaRawbe = new TH1D("hdEtaRawbe","hdEtaRawbe",100, -5,5);
TH1D *hrefmult = new TH1D("hrefmult","hrefmult",1000,0,500);
TH1D *hept = new TH1D("hept","hept",100,0.,20.);

void Loop()
{

  TH1F::SetDefaultSumw2(); 
  TH1D::SetDefaultSumw2(); 

  // Make Chain
  TChain* chain = new TChain("tree");
  int nfile = 0;
  nfile += chain->Add("pythia_tree_Oct6_1.root");
  //nfile += chain->Add("liweiTemplate_part2.root");
  cout <<"Added "<<nfile<<" files"<<endl;
  cout<<"# entries in chain: "<<chain->GetEntries()<<endl;
  if (chain == 0) return;

  int ceNtrigger=0;
  int beNtrigger=0;
  int bceNtrigger=0;

  //define variables
  Int_t   Event, numberofcElectrons, numberofbElectrons,numberofbcElectrons, numberofHadrons, noTracks;   //
  Float_t celectron_id,celectron_status,celectron_pt,celectron_pz,celectron_phi,celectron_eta,celectron_y;                //track info
  Float_t belectron_id,belectron_status,belectron_pt,belectron_pz,belectron_phi,belectron_eta,belectron_y;
  Float_t bcelectron_id,bcelectron_status,bcelectron_pt,bcelectron_pz,bcelectron_phi,bcelectron_eta,bcelectron_y;
  Float_t assoh_id,assoh_status,assoh_pt,assoh_pz,assoh_phi,assoh_eta,assoh_y;

  int goodEvent = 0;


  Long64_t nentries = chain->GetEntriesFast();
  int x = 0; int n = nentries; int w = 25;

  for(int i = 0; i < nentries; i++){
    Long64_t ientry = chain->LoadTree(i);
    if (ientry < 0) break;
    TBranch *b_destep = chain->GetBranch("hf2eDecay");

    TLeaf* leaf_Event_id            = b_destep->GetLeaf("Event_id");
    TLeaf* leaf_refMult             = b_destep->GetLeaf("refMult");
    TLeaf* leaf_numberofcElectrons  = b_destep->GetLeaf("numberofcElectrons");
    TLeaf* leaf_numberofbElectrons  = b_destep->GetLeaf("numberofbElectrons");
    TLeaf* leaf_numberofbcElectrons = b_destep->GetLeaf("numberofbcElectrons");
    TLeaf* leaf_numberofHadrons     = b_destep->GetLeaf("numberofHadrons");
    TLeaf* leaf_noTracks            = b_destep->GetLeaf("noTracks");

    TLeaf* leaf_ce_id               = b_destep->GetLeaf("ce_id");
    TLeaf* leaf_ce_status           = b_destep->GetLeaf("ce_status");
    TLeaf* leaf_ce_pt               = b_destep->GetLeaf("ce_pt");
    TLeaf* leaf_ce_pz               = b_destep->GetLeaf("ce_pz");
    TLeaf* leaf_ce_phi              = b_destep->GetLeaf("ce_phi");
    TLeaf* leaf_ce_eta              = b_destep->GetLeaf("ce_eta");
    TLeaf* leaf_ce_y                = b_destep->GetLeaf("ce_y");

    TLeaf* leaf_be_id               = b_destep->GetLeaf("be_id");
    TLeaf* leaf_be_status           = b_destep->GetLeaf("be_status");
    TLeaf* leaf_be_pt               = b_destep->GetLeaf("be_pt");
    TLeaf* leaf_be_pz               = b_destep->GetLeaf("be_pz");
    TLeaf* leaf_be_phi              = b_destep->GetLeaf("be_phi");
    TLeaf* leaf_be_eta              = b_destep->GetLeaf("be_eta");
    TLeaf* leaf_be_y                = b_destep->GetLeaf("be_y");

    TLeaf* leaf_bce_id              = b_destep->GetLeaf("bce_id");
    TLeaf* leaf_bce_status          = b_destep->GetLeaf("bce_status");
    TLeaf* leaf_bce_pt              = b_destep->GetLeaf("bce_pt");
    TLeaf* leaf_bce_pz              = b_destep->GetLeaf("bce_pz");
    TLeaf* leaf_bce_phi             = b_destep->GetLeaf("bce_phi");
    TLeaf* leaf_bce_eta             = b_destep->GetLeaf("bce_eta");
    TLeaf* leaf_bce_y               = b_destep->GetLeaf("bce_y");

    TLeaf* leaf_hadron_id           = b_destep->GetLeaf("hadron_id");
    TLeaf* leaf_hadron_status       = b_destep->GetLeaf("hadron_status");
    TLeaf* leaf_hadron_pt           = b_destep->GetLeaf("hadron_pt");
    TLeaf* leaf_hadron_pz           = b_destep->GetLeaf("hadron_pz");
    TLeaf* leaf_hadron_phi          = b_destep->GetLeaf("hadron_phi");
    TLeaf* leaf_hadron_eta          = b_destep->GetLeaf("hadron_eta");
    TLeaf* leaf_hadron_y            = b_destep->GetLeaf("hadron_y");

    x = i+1;
    // Process Completion bar
    if ( (x != n) && (x % (n/100) == 0) )
    {
      float ratio  =  (float)x/(float)n;
      int   c      =  ratio * w;

      if(ratio < 1){
        cout << setw(3) << (int)(ratio*100) << "% [";
        for (int x=0; x<c; x++) cout << "=";
        for (int x=c; x<w; x++) cout << " ";
        cout << "]\r" << flush ;
      }
    }

    chain->GetEntry(i);

    Event   = (int)leaf_Event_id->GetValue(0);
    numberofcElectrons = (int)leaf_numberofcElectrons->GetValue(0);
    //   cout << numberofcElectrons << " ";
    numberofbElectrons = (int)leaf_numberofbElectrons->GetValue(0);
    numberofbcElectrons = (int)leaf_numberofbcElectrons->GetValue(0);
    numberofHadrons = (int)leaf_numberofHadrons->GetValue(0);
    // cout << "numHad: " << numberofHadrons << endl;
    noTracks        = leaf_noTracks->GetValue(0);
    hrefmult -> Fill(noTracks);

    //loop through matched primary tracks electron find c decayed electron
    int ceNtrigcount=0;
    for(int trki = 0; trki < numberofcElectrons; trki++){
      celectron_id      = (int)leaf_ce_id->GetValue(trki);
      celectron_status  = (int)leaf_ce_status->GetValue(trki);
      celectron_pt      = leaf_ce_pt->GetValue(trki);
      //  cout << "pt: " << celectron_pt << endl;
      celectron_pz      = leaf_ce_pz->GetValue(trki);
      celectron_phi     = leaf_ce_phi->GetValue(trki);
      celectron_eta     = leaf_ce_eta->GetValue(trki);
      celectron_y       = leaf_ce_y->GetValue(trki);
      if(celectron_pt < pt_trig_lo) continue;              
      if(celectron_eta > EtaCut || celectron_eta < -EtaCut) continue;
      hept->Fill(celectron_pt);

      for(int trkj = 0; trkj < numberofHadrons; trkj++){

        assoh_id = (int)leaf_hadron_id->GetValue(trkj);
        assoh_status  = (int)leaf_hadron_status->GetValue(trkj);
        assoh_pt      = leaf_hadron_pt->GetValue(trkj);
        assoh_pz      = leaf_hadron_pz->GetValue(trkj);
        assoh_phi     = leaf_hadron_phi->GetValue(trkj);
        assoh_eta     = leaf_hadron_eta->GetValue(trkj);
        assoh_y       = leaf_hadron_y->GetValue(trkj);
        if(assoh_eta > hEtaCut || assoh_eta < -hEtaCut) continue;

        float deltPhi = assoh_phi - celectron_phi;
        float deltEta = assoh_eta - celectron_eta;
        if(deltPhi < -PI)  deltPhi += 2*PI;
        if(deltPhi >  PI) deltPhi -= 2*PI;
        if(abs(deltEta) > deletacut) continue;
        if(celectron_pt>pt_trig_lo && celectron_pt<pt_trig_up && assoh_pt>pt_asso_lo)
        {
          hdPhiRawce->Fill(deltPhi); 
          hdEtaRawce->Fill(deltEta);
          ceNtrigcount++;
        }
      }
    }

    if(ceNtrigcount>0)hEventTallyce->Fill("ce non photonic electron",1);
    if(ceNtrigcount>0)ceNtrigger++;

    //b decayed electron-------------------------------------------------- 
    int beNtrigcount=0;
    for(int trki = 0; trki < numberofbElectrons; trki++){
      belectron_id = (int)leaf_be_id->GetValue(trki);
      belectron_status  = (int)leaf_be_status->GetValue(trki);
      belectron_pt      = leaf_be_pt->GetValue(trki);
      belectron_pz      = leaf_be_pz->GetValue(trki);
      belectron_phi     = leaf_be_phi->GetValue(trki);
      belectron_eta     = leaf_be_eta->GetValue(trki);
      belectron_y        = leaf_be_y->GetValue(trki);
      if(belectron_pt < pt_trig_lo) continue;                   //here right to set pt cut? YES,it is to electron cut------------------------------------       
      if(belectron_eta > EtaCut || belectron_eta < -EtaCut) continue;
      hept->Fill(belectron_pt);
      for(int trkj = 0; trkj < numberofHadrons; trkj++){

        assoh_id = (int)leaf_hadron_id->GetValue(trkj);
        assoh_status  = (int)leaf_hadron_status->GetValue(trkj);
        assoh_pt      = leaf_hadron_pt->GetValue(trkj);
        assoh_pz      = leaf_hadron_pz->GetValue(trkj);
        assoh_phi     = leaf_hadron_phi->GetValue(trkj);
        assoh_eta     = leaf_hadron_eta->GetValue(trkj);
        assoh_y       = leaf_hadron_y->GetValue(trkj);
        if(assoh_eta > hEtaCut || assoh_eta < -hEtaCut) continue;

        float deltPhi = assoh_phi - belectron_phi;
        float deltEta = assoh_eta - belectron_eta;
        if(deltPhi < -PI)  deltPhi += 2*PI;
        if(deltPhi >  PI) deltPhi -= 2*PI;
        if(abs(deltEta) > deletacut) continue;
        if(belectron_pt>pt_trig_lo && belectron_pt<pt_trig_up && assoh_pt>pt_asso_lo)
        {
          hdPhiRawbe->Fill(deltPhi);
          hdEtaRawbe->Fill(deltEta);  
          beNtrigcount++;
        }
      }
    }
    if(beNtrigcount>0)hEventTallybe->Fill("be non photonic electron",1);
    if(beNtrigcount>0)beNtrigger++;

    //bce decayed electron-----------------------------------------------------------------                                                                                           
    int bceNtrigcount=0;
    for(int trki = 0; trki < numberofbcElectrons; trki++){
      bcelectron_id = (int)leaf_bce_id->GetValue(trki);
      bcelectron_status  = (int)leaf_bce_status->GetValue(trki);
      bcelectron_pt      = leaf_bce_pt->GetValue(trki);
      bcelectron_pz      = leaf_bce_pz->GetValue(trki);
      bcelectron_phi     = leaf_bce_phi->GetValue(trki);
      bcelectron_eta     = leaf_bce_eta->GetValue(trki);
      bcelectron_y        = leaf_bce_y->GetValue(trki);
      if(bcelectron_pt < pt_trig_lo) continue;                   //here right to set pt cut? YES,it is to electron cut------------------------------------      
      //  if(bcelectron_eta > EtaCut || bcelectron_eta < -EtaCut) continue;
      hept->Fill(bcelectron_pt);
      for(int trkj = 0; trkj < numberofHadrons; trkj++){

        assoh_id = (int)leaf_hadron_id->GetValue(trkj);
        assoh_status  = (int)leaf_hadron_status->GetValue(trkj);
        assoh_pt      = leaf_hadron_pt->GetValue(trkj);
        assoh_pz      = leaf_hadron_pz->GetValue(trkj);
        assoh_phi     = leaf_hadron_phi->GetValue(trkj);
        assoh_eta     = leaf_hadron_eta->GetValue(trkj);
        assoh_y       = leaf_hadron_y->GetValue(trkj);
        if(assoh_eta > hEtaCut || assoh_eta < -hEtaCut) continue;

        float deltPhi = assoh_phi - bcelectron_phi;
        float deltEta = assoh_eta - celectron_eta;
        if(deltPhi < -PI)  deltPhi += 2*PI;
        if(deltPhi >  PI) deltPhi -= 2*PI;
        if(abs(deltEta) > deletacut) continue;
        if(bcelectron_pt>pt_trig_lo && bcelectron_pt<pt_trig_up && assoh_pt>pt_asso_lo)
        {
          hdPhiRawbe->Fill(deltPhi);
          hdEtaRawbe->Fill(deltEta);  
          bceNtrigcount++;
        }
      }
    }
    if(bceNtrigcount>0)hEventTallybe->Fill("be non photonic electron",1);
    if(bceNtrigcount>0)beNtrigger++;


    /* if(ceNtrigger+bceNtrigger+beNtrigger > 0){
       cout<<"ce Trigger electron number = "<<ceNtrigger<<endl;
       cout<<"be Trigger electron number = "<<beNtrigger<<endl;
       cout<<"bce Trigger electron number = "<<bceNtrigger<<endl;
       }*/
  }

  // After Fill Manipulations

  hdPhiRawceN = (TH1D*)hdPhiRawce->Clone();
  hdPhiRawceN -> SetName("hdPhiRawceN");
  hdPhiRawceN -> Sumw2();
  hdPhiRawceN -> Scale(1./(Double_t)hEventTallyce->GetBinContent(1));
  hdPhiRawbeN = (TH1D*)hdPhiRawbe->Clone();
  hdPhiRawbeN -> SetName("hdPhiRawbeN");
  hdPhiRawbeN -> Sumw2();
  hdPhiRawbeN -> Scale(1./(Double_t)hEventTallybe->GetBinContent(1));


  cout << "100% [";
  for (int x=0; x<w; x++) cout << "=";
  cout << "]" << endl;
  fout.Write();
  fout.Close();
  delete chain;
}
