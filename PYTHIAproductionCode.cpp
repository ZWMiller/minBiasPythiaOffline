//==============================================================================
//  pmainHF2e.cpp
//    
//  This is an example program to study c -> e and b -> e decays
//  in 200 GeV pp collisions with Pythia8.
//
//  The decays are stored in a ROOT tree and written to file.
//
//  Once written most things can be controlled through the runcard,
//  so there's no need to recompile.
//
//  Usage: pmainHF2e  runcard  rootfile
// 
//  Author: Thomas Ullrich  
//  Last update: September 9, 2008
//==============================================================================
#include <cmath>
#include "Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <fstream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#define PR(x) std::cout << #x << " = " << (x) << std::endl;
using namespace Pythia8; 

//
//  Forward declarations
//
//bool isInAcceptance(int, const Event&);  // acceptance filter
//int myEvent(Pythia&, double,int);            // event handler (analyze event and

// stores data in tuple)
//
// This structure contains all the info we 
// collect. This info is later stored in a tree.
// This is our own business and has nothing to do
// with Pythia directly.
//
struct hf2eDecay_t {
  int Event_id;
  int refMult;
  int numberofcElectrons;
  int numberofbElectrons;
  int numberofbcElectrons;
  int numberofHadrons;   
  int noTracks; 
  /*    int hf_id;         // mother (c/b hadron)
	int hf_status;
	float hf_pt;
	float hf_pz;
	float hf_phi;
	float hf_eta;
	float hf_y;
  */    
  int ce_id[3];            // c to electron
  int ce_status[3];
  float ce_pt[3];
  float ce_pz[3];
  float ce_phi[3];
  float ce_eta[3];
  float ce_y[3];
   

  int be_id[3];            // b  to electron
  int be_status[3];
  float be_pt[3];
  float be_pz[3];
  float be_phi[3];
  float be_eta[3];
  float be_y[3];


  int bce_id[3];            // b to c to electron
  int bce_status[3];
  float bce_pt[3];
  float bce_pz[3];
  float bce_phi[3];
  float bce_eta[3];
  float bce_y[3];

	
  int hadron_id[800];            // hadron
  int hadron_status[800];
  float hadron_pt[800];
  float hadron_pz[800];
  float hadron_phi[800];
  float hadron_eta[800];
  float hadron_y[800];

  /*
    int   q1_id;          // q1 and q2 stand for what
    float q1_x;
    int   q2_id;
    float q2_x;
    float Q2fac;
    float alphas;
    float ptHat;
    int   nFinal;
    float pdf1;
    float pdf2;
    int   code;
    float sigmaGen;
    float weight;   // useful for normalization/x-section
  */
};

hf2eDecay_t hf2eDecay;

int main(int argc, char* argv[]) {
    
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " runcard  rootfile" << endl;
    return 2;
  }
  char* runcard  = argv[1];
  char* rootfile = argv[2];
  const char* xmlDB    = "/global/u1/l/liwei/minepythia/PYTHIA/pythia8135/xmldoc";
    
  //--------------------------------------------------------------
  //  Initialization
  //--------------------------------------------------------------
    
  //
  //  ROOT
  //
  TFile *hfile  = new TFile(rootfile,"RECREATE");
  TTree tree("tree","b->e c->e and b->c->e decays pp at 500 GeV");
  tree.Branch("hf2eDecay",&hf2eDecay.Event_id,
	      "Event_id/I:refMult/I:numberofcElectrons/I:numberofbElectrons/I:numberofbcElectrons/I:numberofHadrons/I:noTracks/I:"
	      "ce_id[3]/I:ce_status[3]/I:ce_pt[3]/F:ce_pz[3]/F:ce_phi[3]/F:ce_eta[3]/F:ce_y[3]/F:"
	      "be_id[3]/I:be_status[3]/I:be_pt[3]/F:be_pz[3]/F:be_phi[3]/F:be_eta[3]/F:be_y[3]/F:"
	      "bce_id[3]/I:bce_status[3]/I:bce_pt[3]/F:bce_pz[3]/F:bce_phi[3]/F:bce_eta[3]/F:bce_y[3]/F:"
	      "hadron_id[800]/I:hadron_status[800]/I:hadron_pt[800]/F:hadron_pz[800]/F:hadron_phi[800]/F:hadron_eta[800]/F:hadron_y[800]/F");
  //  Create instance of Pythia 
  //
  Pythia pythia(xmlDB); // the init parameters are read from xml files
  // stored in the xmldoc directory. This includes
  // particle data and decay definitions.
    
  //
  // Shorthand for (static) settings
  //
  Settings& settings = pythia.settings;
    
  //
  //  Read in runcard
  //
  pythia.readFile(runcard);  
  cout << "Runcard '" << runcard << "' loaded." << endl;
    
  //
  //  Retrieve number of events and other parameters from the runcard.
  //  We need to deal with those settings ourself. Getting
  //  them through the runcard just avoids recompiling.
  //
  int  maxNumberOfEvents = settings.mode("Main:numberOfEvents");
  int  nList     = settings.mode("Main:numberToList");
  int  nShow     = settings.mode("Main:timesToShow");
  int  maxErrors = settings.mode("Main:timesAllowErrors");
  bool showCS    = settings.flag("Main:showChangedSettings");
  bool showAS    = settings.flag("Main:showAllSettings");
  int  pace = maxNumberOfEvents/nShow;
	
  //
  //  Remark: in this example we do NOT alter the
  //  BRs since they are different for the various charm
  //  hadrons making the normalization at the end rather
  //  cumbersome. In a production version this is what
  //  one probably would implement to save processing time.
  //
    
  //
  //  Initialize Pythia, ready to go
  //
  pythia.init();
    
  //
  // List changed or all data
  //
  if (showCS) settings.listChanged();
  if (showAS) settings.listAll();
    
  //--------------------------------------------------------------
  //  Event loop
  //--------------------------------------------------------------
  int ievent = 0;
  int numberOfElectrons = 0;
  int iErrors = 0;
  while (ievent < maxNumberOfEvents) {
        
    if (!pythia.next()) {
      if (++iErrors < maxErrors) continue;
      cout << "Error: too many errors in event generation - check your settings & code" << endl;
      break;
    }

    Event &event = pythia.event;
    int eventId=ievent+maxNumberOfEvents;
    hf2eDecay.Event_id   = eventId;
    //    if(!event){ cout << " Error : No Event" << endl;return 0; }           //how to ignore wrong event------------------
    //    int nelectrons = 0;
    int total_tracks = 0;
    int ic =0;
    int e_number=0;
    int nce=0,nbe=0,nbce=0;
    int hadron_number=0;
    int nh=0;
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------
    for (int i = 0; i < event.size(); i++) {                                      //event.size means total track number
      //       if(abs(event[i].id())>100)
      if(event[i].isHadron()&&event[i].isFinal())                     //select hadron
	{
	  hf2eDecay.hadron_id[nh] = event[i].id();
	  hf2eDecay.hadron_status[nh] = event[i].status();
	  hf2eDecay.hadron_pt[nh] = event[i].pT();
	  hf2eDecay.hadron_pz[nh] = event[i].pz();
	  hf2eDecay.hadron_phi[nh] = event[i].phi();
	  hf2eDecay.hadron_eta[nh] = event[i].eta();
	  hf2eDecay.hadron_y[nh] = event[i].y();

	  nh++;
	}
      if (abs(event[i].id()) == 11) {

	//
	//  Check if mother is a c/b hadron
	//
	vector<int> mothers = event.motherList(i);
	if (mothers.size() != 1) {
	  cout << "Error: electron has more than one mother. Stop." << endl;
	  //abort();
	  //  return 0;
	  continue;
	}
	ic = mothers[0];
	int ic_id = abs(event[ic].id());
	int flavor = static_cast<int>(ic_id/pow(10.,static_cast<int>(log10(ic_id)))); //this will miss excited states like D0*
	//          if (flavor != 4 && flavor != 5) continue; // c (b) hadrons start with 4(5)  
	if (flavor != 4 && flavor != 5) continue;
	//  Acceptance filter
	//            if (!(isInAcceptance(i, event))) continue;                    //--------------------------------------------------------------
	if(flavor == 4)
	  {
            vector<int> grandmothers = event.motherList(ic);
	    /*            if(grandmothers.size()==0)
			  {   hf2eDecay.ce_id[nce] = event[i].id();
			  hf2eDecay.ce_status[nce] = event[i].status();
			  hf2eDecay.ce_pt[nce] = event[i].pT();
			  hf2eDecay.ce_pz[nce] = event[i].pz();
			  hf2eDecay.ce_phi[nce] = event[i].phi();
			  hf2eDecay.ce_eta[nce] = event[i].eta();
			  hf2eDecay.ce_y[nce] = event[i].y();
			  nce++;
			  }
	    */
	    if(grandmothers.size()==1)
	      {
                int orig = grandmothers[0];
                int orig_id = abs(event[orig].id());     //need absolute value here? -----------------------------------------
                int flavor_grandmon = static_cast<int>(orig_id/pow(10.,static_cast<int>(log10(orig_id)))); //this will miss excited states like D0*
                if(flavor_grandmon == 5)
		  {
		    hf2eDecay.bce_id[nbce] = event[i].id();
		    hf2eDecay.bce_status[nbce] = event[i].status();
		    hf2eDecay.bce_pt[nbce] = event[i].pT();
		    hf2eDecay.bce_pz[nbce] = event[i].pz();
		    hf2eDecay.bce_phi[nbce] = event[i].phi();
		    hf2eDecay.bce_eta[nbce] = event[i].eta();
		    hf2eDecay.bce_y[nbce] = event[i].y();
		    nbce++;
		  }
		else
		  {
		    hf2eDecay.ce_id[nce] = event[i].id();
		    hf2eDecay.ce_status[nce] = event[i].status();
		    hf2eDecay.ce_pt[nce] = event[i].pT();
		    hf2eDecay.ce_pz[nce] = event[i].pz();
		    hf2eDecay.ce_phi[nce] = event[i].phi();
		    hf2eDecay.ce_eta[nce] = event[i].eta();
		    hf2eDecay.ce_y[nce] = event[i].y();
		    nce++;
		  }


	      }

	  }
	if(flavor == 5)
	  {
	    hf2eDecay.be_id[nbe] = event[i].id();
	    hf2eDecay.be_status[nbe] = event[i].status();
	    hf2eDecay.be_pt[nbe] = event[i].pT();
	    hf2eDecay.be_pz[nbe] = event[i].pz();
	    hf2eDecay.be_phi[nbe] = event[i].phi();
	    hf2eDecay.be_eta[nbe] = event[i].eta();
	    hf2eDecay.be_y[nbe] = event[i].y();
	    nbe++;
	  }

	//            ne++;                     //counts have some problems
	//            nelectrons++;

	//
	// Get grandmother (origin of c/b hadron)
	//
	//  Store in tuple
	//

	// If no origin or more than 1 than id == 0
	// and the status identifies what happens: -1 no mother, 
	// -2 more than 1. This should not happen at all, but ...
      }

    }              //end of this track loop                            
   
    e_number=nbe+nce+nbce;
    hadron_number=nh;
    total_tracks = event.size();
    hf2eDecay.numberofcElectrons =nce;
    hf2eDecay.numberofbElectrons =nbe;
    hf2eDecay.numberofbcElectrons =nbce;
    hf2eDecay.numberofHadrons = hadron_number;
    hf2eDecay.noTracks = e_number+hadron_number;    //event.size();
    if(e_number > 0) tree.Fill();
    for(int i=0;i<nh;i++)
      {
	hf2eDecay.hadron_id[i] = 0;
	hf2eDecay.hadron_status[i] = 0;
	hf2eDecay.hadron_pt[i] = 0.;
	hf2eDecay.hadron_pz[i] = 0.;
	hf2eDecay.hadron_phi[i] = 0.;
	hf2eDecay.hadron_eta[i] = 0.;
	hf2eDecay.hadron_y[i] = 0.;
      }
    for(int i=0;i<nce;i++)
      {
	hf2eDecay.ce_id[i] = 0;
	hf2eDecay.ce_status[i] = 0;
	hf2eDecay.ce_pt[i] = 0.;
	hf2eDecay.ce_pz[i] = 0.;
	hf2eDecay.ce_phi[i] = 0.;
	hf2eDecay.ce_eta[i] = 0.;
	hf2eDecay.ce_y[i] = 0.;

      }
    for(int i=0;i<nbe;i++)
      {
	hf2eDecay.be_id[i] = 0;
	hf2eDecay.be_status[i] = 0;
	hf2eDecay.be_pt[i] = 0.;
	hf2eDecay.be_pz[i] = 0.;
	hf2eDecay.be_phi[i] = 0.;
	hf2eDecay.be_eta[i] = 0.;
	hf2eDecay.be_y[i] = 0.;

      }
    for(int i=0;i<nbce;i++)
      {
	hf2eDecay.bce_id[i] = 0;
	hf2eDecay.bce_status[i] = 0;
	hf2eDecay.bce_pt[i] = 0.;
	hf2eDecay.bce_pz[i] = 0.;
	hf2eDecay.bce_phi[i] = 0.;
	hf2eDecay.bce_eta[i] = 0.;
	hf2eDecay.bce_y[i] = 0.;

      }

    //        n= myEvent(pythia, maxNumberOfEvents,ievent);  // in myEvent we deal with the whole event and return
    //        numberOfElectrons += n; 
    //        if (n) tree.Fill();                 //here only fill once no matter how many electrons have been found?     
    ievent++;
    numberOfElectrons = numberOfElectrons+e_number;
    if (ievent%1000 == 0) {
      cout << "# of events generated = " << ievent 
	   << ", # of electrons from c or b hadron decays generated so far = " <<  numberOfElectrons << endl;
    }
        
    // List first few events.
    if (ievent < nList) {
      pythia.info.list(); 
      pythia.process.list();
      pythia.event.list();
    }
  }
    
  //--------------------------------------------------------------
  //  Finish up
  //--------------------------------------------------------------
  pythia.statistics();
  cout << "Writing File" << endl;
  hfile->Write();
  return 0;
}                                 //end of main function-------------------------------------------------------------------

//
//  Event analysis
//
/*int myEvent(Pythia& pythia, double nMaxEvt, int eventnumber)                        // for current event
  {
  Event &event = pythia.event;
  int eventId=eventnumber; 
  hf2eDecay.Event_id   = eventId;
  //    if(!event){ cout << " Error : No Event" << endl;return 0; }           //how to ignore wrong event------------------
  int nelectrons = 0;
  int ic = 0;
  int e_number=0,ne=0;
  int hadron_number=0,nh=0;
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------
  for (int i = 0; i < event.size(); i++) {                                      //event.size means total track number
  if(abs(event[i].id())>100)
  {           
  hf2eDecay.hadron_id[nh] = event[i].id();     
  hf2eDecay.hadron_status[nh] = event[i].status();
  hf2eDecay.hadron_pt[nh] = event[i].pT();
  hf2eDecay.hadron_pz[nh] = event[i].pz();
  hf2eDecay.hadron_phi[nh] = event[i].phi();
  hf2eDecay.hadron_eta[nh] = event[i].eta();
  hf2eDecay.hadron_y[nh] = event[i].y();

  nh++; 
  }
  if (abs(event[i].id()) == 11) {
			
  //
  //  Check if mother is a c/b hadron
  //
  vector<int> mothers = event.motherList(i);
  if (mothers.size() != 1) {
  cout << "Error: electron has more than one mother. Stop." << endl;
  //abort();
  return 0;
  }
  ic = mothers[0];
  int ic_id = abs(event[ic].id());
  int flavor = static_cast<int>(ic_id/pow(10.,static_cast<int>(log10(ic_id)))); //this will miss excited states like D0*
  if (flavor != 4 && flavor != 5) continue; // c (b) hadrons start with 4(5)  
  // if (flavor != 5) continue;
  //
  //  Acceptance filter
  //    
  if (!(isInAcceptance(i, event))) continue;
  hf2eDecay.e_id[ne] = event[i].id();
  hf2eDecay.e_status[ne] = event[i].status();
  hf2eDecay.e_pt[ne] = event[i].pT();
  hf2eDecay.e_pz[ne] = event[i].pz();
  hf2eDecay.e_phi[ne] = event[i].phi();
  hf2eDecay.e_eta[ne] = event[i].eta();
  hf2eDecay.e_y[ne] = event[i].y();
             
  ne++;
  nelectrons++;
            
  //
  // Get grandmother (origin of c/b hadron)
  //
  //  Store in tuple
  //
            
  // If no origin or more than 1 than id == 0
  // and the status identifies what happens: -1 no mother, 
  // -2 more than 1. This should not happen at all, but ...
  }                          
       
  }              //end of this track loop                             
  e_number=ne;
  hadron_number=nh;  
  hf2eDecay.numberofElectrons = e_number;
  hf2eDecay.numberofHadrons = hadron_number;
  if(e_number > 0) tree.Fill();

  return nelectrons;
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------
  }
*/
//
//  Acceptance filter
//
/*
  bool isInAcceptance(int i, const Event& event)
  {
  // accept all (useful for many studies)
  return true;
    
  // limit to STAR TPC/BEMC/ToF acceptance
  //double eta = event[i].eta();
  //if (fabs(eta) < 1)
  //    return true;
  //else
  //    return false;
  }
*/
