//copied from test33.cc    This code is using the gcc compilor
// This code is using for the coin time purpose
//JULY 01, 2019
// To run this code do  make and then ./bin/coin_time
// This code  is using the gcc compilor
//Should also have a make file to run the gcc compilor
// In order to run this code you must have a separate makefile
// the code for the  makefile is saved  in the bottom of this code
// do emacs makefile and paste the text from the bottom of this code.

#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
using namespace std;

#include "TApplication.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TRandom.h"

int main(int argc, char** argv){
  
  TApplication *theApp = new TApplication("App", &argc, argv);
  TChain *T = new TChain("T");
  Double_t tdcTime = 56.23e-12;    // in ns:   F1 TDCs, 56.23 ps / Ch
  //Double_t c_light = 0.3;// ns
  // Double_t tdcTime1 = 58.0e-12; //after 111368 runs for left arm
  gStyle->SetOptStat(11111);
 
  //TH1F *h1 = new TH1F("h1","a1<40 and a2>1900,cut<7200.0", 500,-40,40); //-0.1e-6,0.8e-6);
  TH1F *h1 = new TH1F("h1","histogram1", 500, -40.0, 40.0);
  
 
  TString filename;
  
  // for(int irun = 111171;irun<111172;irun++){ //100677;irun<100679
  // if(irun ==111193) continue;
  //    T->Add(Form("./Rootfiles/tritium_%d*.root",irun));
    //  }
   T->Add(Form("./Rootfiles/tritium_111171*.root"));
  //   T->Add(Form("./Rootfiles/out.root"));	 
  cout<<"irun done"<<endl;

  /////////////////////////////////////////////////////////
  ///Defining variables for Left HRS 
  /////////////////////////////////////////////////////////
  Double_t coin_trig;
  /// variables to be used for Left HRAS
  Double_t RF_s2L_mean;

  // TString nLHRS_trig = "DR.evtypebits";
  TString nCOINCOIN = "DR.T5";
  TString nLs2_pad = "L.s2.t_pads";
  TString nLs2_nthit = "L.s2.nthit";
  TString nLs2_tdchit = "LTDC.F1FirstHit";
  TString nL_trx = "L.tr.x";
  TString nL_trth = "L.tr.th";
  TString nLs2_ladc = "L.s2.la_c";
  TString nL_trpath = "L.tr.pathl";
  TString nLs2_trpath = "L.s2.trpath";
  // Define TTree leaf variables to hold the values
  // Double_t LHRS_trig;
  Double_t COINCOIN;
  Double_t  Ls2_pad[100];
  Double_t Ls2_nthit;
  Double_t Ls2_tdchit[100];
  Double_t L_trx[100];
  Double_t  L_trth[100];
  Double_t Ls2_ladc[100];
  Double_t L_trpath;
  Double_t Ls2_trpath;
  

 
  
  // Now set the branch address LHRS
  T->SetBranchStatus("*",0);
  // T->SetBranchStatus(nLHRS_trig, 1); T->SetBranchAddress(nLHRS_trig, &LHRS_trig);
  // T->SetBranchAddress("DR.T5", &COINCOIN);
  T->SetBranchStatus(nCOINCOIN,1);   T->SetBranchAddress(nCOINCOIN, &COINCOIN); 
  T->SetBranchStatus(nLs2_pad,1);    T->SetBranchAddress(nLs2_pad, &Ls2_pad);
  T->SetBranchStatus(nLs2_nthit,1);  T->SetBranchAddress(nLs2_nthit, &Ls2_nthit);
  T->SetBranchStatus(nLs2_tdchit,1); T->SetBranchAddress(nLs2_tdchit, &Ls2_tdchit);
  T->SetBranchStatus(nL_trx, 1);     T->SetBranchAddress(nL_trx, &L_trx);
  T->SetBranchStatus(nL_trth, 1);    T->SetBranchAddress(nL_trth, &L_trth);
  T->SetBranchStatus(nLs2_ladc,1);   T->SetBranchAddress(nLs2_ladc, &Ls2_ladc);
  T->SetBranchStatus(nL_trpath,1);   T->SetBranchAddress(nL_trpath, &L_trpath);
  T->SetBranchStatus(nLs2_trpath,1); T->SetBranchAddress(nLs2_trpath, &Ls2_trpath);
  
  // Defining array LHRS
  Double_t corr_L_x[14] = {9.5982e-09, 2.39686e-09, 5.50452e-09, 8.67284e-09, 7.88134e-09, 
			   9.39930e-09,   9.09441e-09, 8.13305e-09, 8.36477e-09, 8.74297e-09, 
			   7.745e-09,  5.94972e-09, 6.22836e-09, 5.52765e-09};
  
  Double_t cL =-9.15734e-10;   // previous one  4.87486E-11;                -3.9094e-09;
  for(int l=0;l<14;l++) {
    corr_L_x[l] = corr_L_x[l] + cL;
  }
  
  Double_t corr_L_th[14] = {-5.3783e-08,  - 3.32898e-08, -4.14532e-08, -4.08767e-08, 
			    -4.07972e-08, -3.63437e-08,  -3.67840e-08, -3.54952e-08, 
			    -3.63706e-08,-3.39145e-08, -3.43925e-08,  -3.05521e-08,
			    -3.07010e-08, -3.79624e-08};
  Double_t cL1 = 1.75759e-9;   // previous one  4.87486E-11;                -3.9094e-09;
  for(int m=0;m<14;m++) {
    corr_L_th[m] = corr_L_th[m] + cL1;
  }
  
  Double_t corr_L_adc[14] = {- 1.592e-12, -1.24122e-12, -1.18518e-12, -1.16133e-12, 
			     -1.24632e-12, -1.22617e-12, -1.02470e-12, -6.57058e-13, 
			     -1.14584e-12, -1.3259e-12, -1.816135e-12, -1.15547e-12,  
			     -1.23475e-12, -1.50406e-12};
  Double_t alignment_L[14] = {1.0319760e-9, -1.0e-9, -0.35e-9, 9.985e-10, 9.835e-10,  
			      4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 
			      9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10 };
  
  //Double_t alignment_L[14] = {7.319760e-9, -1.0e-9, -0.35e-9, 9.985e-10, 9.835e-10,  4.748e-10, 1.257e-10, 0, -1.785e-10, -7.9345e-10, 9.985e-10, 9.975e-10,  1.485e-10,7.9375e-10 };
  
  
  
  /////////////////////////////////////////////
  //now defining same for the RHRS system
  ////////////////////////////////////////
  
  
  // Defining variables for RHRS 
  
  // Double_t S0_mean_time;
  Double_t RF_s2R_mean;
  // Define Tree Leaf name for the RHRS that will be used.
  
  TString nRs2_pads = "R.s2.t_pads";
  TString nRs2_nthit = "R.s2.nthit";
  TString nRs2_tdchit = "RTDC.F1FirstHit";
  TString nR_trx = "R.tr.x";
  TString nR_trth = "R.tr.th";
  TString nRs2_lac = "R.s2.la_c";
  TString nR_trbeta="R.tr.beta";
  TString nR_a1 ="R.a1.asum_c";
  TString nR_a2 ="R.a2.asum_c";
  TString nR_trpath = "R.tr.pathl";
  TString nRs2_trpath = "R.s2.trpath";
  TString nR_sh="R.sh.asum_c";
  TString nR_ps="R.ps.asum_c";
 
  //  TString nRs2_trpad = "R.s2.trpad";
  // NOw define the tree leaf variable
  // Double_t HRS_trig;
  Double_t Rs2_pads[100];
  Double_t Rs2_nthit;
  Double_t Rs2_tdchit[100];
  Double_t R_trx[100];
  Double_t R_trth[100];
  Double_t Rs2_lac[100];
  Double_t R_trbeta;
  Double_t R_a1;
  Double_t R_a2;
  Double_t R_trpath;
  Double_t Rs2_trpath;
  Double_t R_sh;
  Double_t R_ps;

  Double_t sum_cut;
  // Double_t Rs2_trpad;
  // Set the branch address
  T->SetBranchAddress(nRs2_pads, &Rs2_pads);
  T->SetBranchAddress(nRs2_nthit, &Rs2_nthit);
  T->SetBranchAddress(nRs2_tdchit, &Rs2_tdchit);
  T->SetBranchAddress(nR_trx, &R_trx);
  T->SetBranchAddress(nR_trth, &R_trth);
  T->SetBranchAddress(nRs2_lac, &Rs2_lac);
  T->SetBranchAddress(nR_trbeta, &R_trbeta);
  T->SetBranchAddress(nR_a1, &R_a1);
  T->SetBranchAddress(nR_a2, &R_a2);
  T->SetBranchAddress(nR_trpath, &R_trpath);
  T->SetBranchAddress(nRs2_trpath, &Rs2_trpath);
  T->SetBranchAddress(nR_sh, &R_sh);
  T->SetBranchAddress(nR_ps, &R_ps);

  // T->SetBranchAddress(nRs2_trpad, &Rs2_trpad);
  // Defining some array variables
  
  Double_t corr_R_x[14] = { 4.91939e-09, 5.41034e-09, 6.87688e-09, 9.22121e-09, 0.795016e-8, 0.802636e-08, 0.787479e-08, 7.54862e-09, 7.5127e-9, 8.48865e-9, 5.42156e-9, 6.27864e-9, 5.78027e-9, 8.605971e-10};
  
  
  Double_t cx =4.87486E-11;   // previous one  4.87486E-11;                -3.9094e-09;
  for(int l=0;l<14;l++)  {
    corr_R_x[l] = corr_R_x[l] + cx;
  }
  
  
  
  Double_t corr_R_th[14] = {-4.46911e-08, -4.79507e-08, -3.59540e-08, -3.04303e-08, -3.25733e-08, 
			    -3.10881e-08, -3.18107e-08, -3.5057e-08,  -3.59703e-08, -3.76206e-08, 
			    -3.68166e-08,  -3.51979e-08, -3.54868e-08, -1.56998e-08};
  //-3.06204e-09; 
  
  Double_t c_th =-3.06204E-09 ;   // previous one-3.06204E-09;     -3.9094e-09;
  for(int m=0;m<14;m++) {
    corr_R_th[m] = corr_R_th[m] + c_th;
  }
  
  
  Double_t corr_R_adc[14] = {-8.66369e-13, -3.84824e-13, -1.45016e-12, -12.08217e-13, -6.95305e-13, 
			     -5.37148e-13, -5.95287e-13, -1.01789e-12, -1.02612e-12, -1.87664e-12, 
			     -3.19282e-12, -3.37812e-12, -7.80362e-13,-7.65658e-13};
  
  
   
  Double_t alignment_R[14] = {-1.915e-9, -1.917e-9, 0.85e-9, 1.90e-9,2.0e-10, 
			      6.200e-10, 0, 9.50e-10, 1.0e-10, 2.2e-10, 
			      2.20e-9, 2.0e-09, 1.6e-09, 1.615e-9};
 
  

 
  // Loop over entries
  Long64_t nentries = T->GetEntries();

  for(Long64_t k=0;k<nentries;k++)  {
  
    T->GetEntry(k);    
    if(k%10000 ==0){cout<<k<<" "<<nentries<<endl;}
    ////////////////////////////////////////
    ///// loop for LHRS
    /////////////////////////////////////////
  
    for (int i=1;i<15;i++){
     
      // up to here is ok
      // cout<<"my name is bishnu"<<endl;
      // cout<< "Rs2_nthit = " <<COINCOIN <<endl;
      if(COINCOIN>0 && Ls2_nthit==1 &&  Ls2_pad[0]==i){
	     
	RF_s2L_mean = (((Ls2_tdchit[30] - Ls2_tdchit[i])*tdcTime 
			+ (Ls2_tdchit[37] - Ls2_tdchit[i+48])*tdcTime )/2.0 
		       + corr_L_x[i-1]*L_trx[0] + corr_L_th[i-1]*L_trth[0] 
		       + corr_L_adc[i-1]*Ls2_ladc[i] + alignment_L[i-1]);	
	

	//	cout<<"my name is bishnu"<<endl;
	///////////////////////////////////////////////////////////////
	//Loop for RHRS paddles
	//////////////////////////////////////////////////////////////
	
	for( int j=2;j<16;j++){
	
	  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  sum_cut = (R_sh+R_ps+R_a2);
	  // cout<< "sum_cut = " <<   R_a2   <<endl;
	 
	    if( Rs2_nthit==1 && Rs2_pads[0]==j && 
		R_trbeta>0.76 && R_trbeta< 2 && R_a1<100.0  && R_a2>1900.0  && sum_cut<7400.0 ){ 
	    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	   
	       
	    
	    RF_s2R_mean =  (((Rs2_tdchit[9]-Rs2_tdchit[j+16])*tdcTime 
			     + (Rs2_tdchit[46]-Rs2_tdchit[j+48])*tdcTime)/2.0 
			    + corr_R_x[j-2]*R_trx[0] + corr_R_th[j-2]*R_trth[0] 
			    + corr_R_adc[j-2]*Rs2_lac[j] + alignment_R[j-2]);
	    
	   
	   
	    coin_trig = (RF_s2L_mean -RF_s2R_mean);
	    
	      
	    coin_trig = coin_trig*1e+9  + 1.4375*R_trx[0] -249.05;      /// changing in tonano second scale
	    //   cout<<"the value of the coin trig = "<<coin_trig<<endl;
	    h1->Fill(coin_trig);
	  
	  }
	 
	} // if ends
	
	
	 
	
      } // end of j<16 loop
      
      
      
    }/// end of i<15 loop
    
    
    
  }// end of entry loop
  
  
  
  
  
 
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
 
  
    h1->GetXaxis()->SetTitle("Coin_Time(ns)");
    h1->GetXaxis()->CenterTitle();
    h1->GetXaxis()->SetTitleSize(0.035);
    h1->GetXaxis()->SetTitleOffset(1.02);
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetYaxis()->CenterTitle();
    h1->GetYaxis()->SetTitleSize(0.035);
    h1->GetYaxis()->SetTitleOffset(1.02);
  
  
   gPad->SetLogz(1);
 
 
  h1->Draw();
  //TFile* fout = new TFile("hist/output.root", "RECREATE");// uncoment when save a root file
  // fout->cd(); // uncoment when save a root file
  //  h1->Write(); // uncomen when save a root file 
  // h2->Write();
  // && h2->Write();
  // &&h3->Write();
  // && h4->Write();
  //c1->Write();// uncoment when save a root file
  /////////TObjArray h(0);
  /////////h.Add(c1);
  /////////h.Add(h2);
  /////////h.Add(h3);
  /////////h.Write();
  //fout->Close();// uncoment when save a root file
  
  theApp->Run();
  return 0;
}

// makefile   just make file no.cc or .txt just makefile 
// CC = gcc 
// CXX = g++ 
// CFLAGS  = -O2

// BINDIR = ./bin
// LIBDIR = ./lib

// ROOTFLAGS = $(shell root-config --cflags)
// ROOTLIBS = $(shell root-config --libs)
// ROOTGLIBS = $(shell root-config --glibs)
// CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
// CXXLIBS = $(ROOTLIBS)



// TARGET1= coin_time
// OBJS1=   coin_time.o



// #aall: $(TARGET1)  $(TARGET2)  $(TARGET3)\
// aall: $(TARGET1)

// $(LIBDIR)/%.o : %.cc
// 	$(CXX) $(CFLAGS) -c -o $@ $< $(CXXFLAGS)

// $(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
// 	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)


// #$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
// #	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)


// .PHONY: clean
// clean:
// 	rm -f $(LIBDIR)/*.o core $(BINDIR)/*


