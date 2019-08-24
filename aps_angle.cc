// March 12, 2019.. copied from angle_recon.cc to plot the angle theta vs phi for aps 2019 presentation
// On March 13, code seems working for theta and phi with 4th order matrix and Zt built in it. total 126 parameter
// Tmarker is included by March 14 2019 and position of the marker is ok now. The grid or tmarker are landed on the correct location
// March 19, 2019 I am going to plot SS_X vs SS_Y by using some additional information. The holes looks ok in X_ss vs Y_ss histogram
// This code can be used to optimize the HRS angle that is theta and phi angle by using the seive slit data.
// can also be used to optimize the Z vertex and so on
extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_ph(double* P, 
			 double xf, double xpf,
			 double yf, double ypf, double); // this part is ok

extern double calcf2t_zt(double* P, 
			 double xf, double xpf,
			 double yf, double ypf);
const double  XFPm=-0.7,  XpFPm=-0.15; 
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74;
const double  XFPr=1.3,   XpFPr=0.27;
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; 
const double  PLm = 25.4, PLr=0.7; 
const double  Ztm = -0.15,Ztr=0.35;
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);

const int nfoil = 10;

double fcent[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
		       0.00, 0.025, 0.05, 0.10, 0.125}; 

double fcent_real[nfoil] = {-0.125, -0.100, -0.075, -0.050, -0.025,
			    0.000, 0.025, 0.050, 0.100, 0.125}; 

double selection_width = 0.0078; 
double row_width = 0.0032;
const int nParamT = 35;
const int nmax = 10000;
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;
double Pzt_opt[nmax]; 


// ========================================
const int nParamT2 = 4; 
double parRaster[nParamT2]; 

extern double calcRasterCor(double a, double b, double c){
  return a*b + c; 
}
double Opt_Par[nParamT2];
double RasterCor; // for x position calculation
double Ras_curx[nmax];
//========================================

const int nrow = 11; //was 9 beefore 
const int ncol = 13; // was 10 before
const int nsshole = nrow*ncol;
double refx[nsshole];// nominal values for each holes across the column
double refy[nsshole];

const int Angle_Par =126;
double theta_opt[nmax];
double phi_opt[nmax];
double theta_recon[nmax];
double phi_recon[nmax];
double theta_real[nrow][nfoil];
double theta[nrow][nfoil];
double phi_real[ncol][nfoil];
double phi[ncol][nfoil];
double  theta_flag_;
double  phi_flag_;
int theta_flag[nmax];


const double hrs_ang = 13.2 * 3.14159/180.; 
double ztR_wRC_[nmax];
//=========== additional variable to plot the SS_X vs SS_Y=======================
const double R = 1.0034; // in meter
double DR_z; // delta_R_z = -z
double DX_sy; // Delta_X_sy = -Beam_Y
double X_ss; // seive slit x 
// for  ss Y pos===================
double DY_sx; 
double DY_sz;
double DR_x;
double DR_z1; // radius change due to different z
double Y_ss;
// ==================================


void aps_angle(){ // Probably the nmain function
  // ======================================== //
  // ======= Opening a ROOT file ============ //
  // ======================================== //
  TChain * t1 = new TChain("T");
  t1->Add("../Rootfiles/tritium_111721_5.root");
  // t1->Add("../Rootfiles/tritium_111718*.root");
  
  Double_t trig1;
  double ent = t1->GetEntries();
  cout<<"entry in the t1=="<<ent<<endl;
  
  int evshift = 30000;
  
  const int max = 100;
  double ltime_s0[max]; 
  double ltime_s2[max];
  
  
  
  double mom2[max]; 
  const int f1n = 64;
  
  double lvz[max]; 
  double th2[max], ph2[max];
  double th3[max], ph3[max]; 
  Int_t runnum; 
  double hallap;
  
  char tempc1[500]; 
  int temp1;
  
  ifstream *ifs = new ifstream("./dat_file/theta_real.dat");
  for (int foil = 0; foil<nfoil;foil++){
    
    for(int row =0;row<nrow;row++){
      *ifs >> theta_real[row][foil];
    }
  }
  /*
  for(int foil = 0;foil<nfoil;foil++){
    cout<< "theta_real "<< foil << ": ";
    for(int row = 0;row<nrow;row++){
      cout<<theta_real[row][foil] << " ";
    }
    cout<< endl;
  }
  */
  char tempc2[500]; 
  int temp2;
  
  ifstream *ifs1 = new ifstream("./dat_file/theta.dat");
  for (int foil = 0; foil<nfoil;foil++){
    for(int row =0;row<nrow;row++){
      *ifs1 >> theta[row][foil];
    }
  }
  /*
  for(int foil = 0;foil<nfoil;foil++){
    cout<< "theta "<< foil << ": ";
    for(int row = 0;row<nrow;row++){
      cout<<theta[row][foil] << " ";
    }
    cout<< endl;
  }
  */
 char tempc3[500]; 
 int temp3;
 ifstream *ifs2 = new ifstream("./dat_file/phi_real.dat");
    for (int foil = 0; foil<nfoil;foil++){
     
      for(int col =0;col<ncol;col++){
	*ifs2 >> phi_real[col][foil]; 
      }
 }

    //================================

    char tempc4[500]; 
    int temp4;
  
    ifstream *ifs3 = new ifstream("./dat_file/phi.dat");
    for (int foil = 0; foil<nfoil;foil++){
      
      for(int col =0;col<ncol;col++){
	*ifs3 >> phi[col][foil];
      }
    }
    //==================================
    TMarker* mark[nfoil][nsshole];// red mark on the histogram
  
    for(int foil = 0;foil<nfoil;foil++){
      int nhole = 0;
      for(int i = 0;i<nrow;i++){
	for(int j = 0; j<ncol;j++){
	  if(foil ==5){ // plotting all z together for X_SS vs Y_ss
	  refx[nhole] = theta_real[i][foil]; // y before
	  refy[nhole] = phi_real[j][foil];
	  }
	  else{
	    refx[nhole] = 0; // =+++++++
	    refy[nhole] = 0; //++++++++++++
	  }
	  //  cout<<refx[nhole]<<endl;
	  mark[foil][nhole] = new TMarker(refy[nhole],refx[nhole],28);
	  nhole++;
	}
      }
    }

 
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];
  const int n = 16;
  
  double lbeta[max];
  double ctime[max]; 
  
  double Lrb; 
  double Lrb_y;
  double TGT_x;
  double rpr;
  
  
  double Ras_x; 
  double Ras_y; 
  double cer_asum;
  
  t1->SetBranchAddress("fEvtHdr.fRun", &runnum   );
  t1->SetBranchAddress("HALLA_p", &hallap   );
  t1->SetBranchAddress("DR.T1", &trig1  );
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("Lrb.x", &Lrb);
  t1->SetBranchAddress("Lrb.y", &Lrb_y); // +++++++++++

  t1->SetBranchAddress("rpr.z", &rpr); 
  
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &Ras_x); 
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &Ras_y);
  
  t1->SetBranchAddress("L.tr.vz", &lvz);
  t1->SetBranchAddress("L.tr.tg_th", &th3);
  t1->SetBranchAddress("L.tr.tg_ph", &ph3);
  t1->SetBranchAddress("L.s0.time", &ltime_s0);
  t1->SetBranchAddress("L.s2.time", &ltime_s2);
  
  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  t1->SetBranchAddress("L.cer.asum_c",  &cer_asum); 
  
  TFile* fnew = new TFile("./output_root/angle_lhrs.root","recreate"); 
  TTree* tnew = new TTree("tree","For z calibration (LHRS)");
  double ztR[max]; 
  double ztR_wRC[max]; 
  
  tnew->Branch("fEvtHdr.fRun", &runnum,"fEvtHdr.fRun/D");
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D" );
  tnew->Branch("DR.T1", &trig1, "DR.T1/D"    );
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D"  );
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D"  );
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D" );
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D" );
  tnew->Branch("L.tr.vz_TG",  &ztR,   "L.tr.vz_TG[100]/D" );
  tnew->Branch("L.tr.vz_TG2",  &ztR_wRC, "L.tr.vz_TG2[100]/D" ); 
  tnew->Branch("Lrb.x", &Lrb, "Lrb.x/D"    ); 
  tnew->Branch("Lrb.y", &Lrb_y, "Lrb.y/D"); // +++++
  tnew->Branch("rpr.z", &rpr, "rpr.z/D"    ); 
  tnew->Branch("L.cer.asum_c", &cer_asum, "L.cer.asum_c/D"    ); 
  tnew->Branch("L.tr.tg_th_TH2", &th2, "L.tr.tg_th_TH2[100]/D");
  tnew->Branch("L.tr.tg_ph_PH2", &ph2, "L.tr.tg_ph_PH2[100]/D");
  
  tnew->Branch("Lrb.Raster2.rawcur.x", &Ras_x, "Lrb.Raster2.rawcur.x/D"); 
  tnew->Branch("Lrb.Raster2.rawcur.y", &Ras_y, "Lrb.Raster2.rawcur.y/D"); 
  tnew->Branch("RasterCor", &RasterCor, "RasterCor/D"); 
  
  double XFP, XpFP; 
  double YFP, YpFP;
  
  ntune_event = 0; 
  for(int i=0 ; i<nParamT ; i++){
    Pzt_opt[i] = -2222.0;
  }
  
  char name_Mzt_L[500];
  // this is the third order mareix
  sprintf(name_Mzt_L,"testmat/opt_par.dat"); // copied the opt_par.dat from NonOpt_49(optimized for non rastered beam) Jan26,2019
  // sprintf(name_Mzt_L,"matrices/zt_LHRS_2.dat"); //  original matrix here 
  ifstream Mzt_L(name_Mzt_L);
  double Pzt_L[nParamT];
  
  for (int i=0;i<nParamT;i++){
    double par=0.;
    int p=0;
    Mzt_L >> par >> p >> p >> p >> p;
    
    Pzt_L[i]=par;
    Pzt_opt[i] =  Pzt_L[i]; 
    
  }
  Mzt_L.close();
  
  
  
  //================================================================================
  ntune_event = 0;
  for(int j=0 ; j<nParamT2 ; j++){ // 
    Opt_Par[j] = -2222.0;
  }
  
  char name_Raster_L[500];
  sprintf(name_Raster_L,"./rastered_mat/raster_opt.dat"); //  copied from newpar_19.dat optimized by toshi code
  //  sprintf(name_Raster_L,"./rastered_mat/raster_halla_opt.dat"); // using original hall A optimized p0 and p1.
  ifstream Raster_L(name_Raster_L);
  for (int i = 0;i<nParamT2; i++){
    Raster_L >> Opt_Par[i];
    
  }
  Raster_L.close();
  // =============================================================
  ntune_event = 0; 
  for(int i=0 ; i<Angle_Par; i++){ 
    theta_opt[i] = -2222.0;
  }
  char name_Angle_L[500];
  sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");
  ifstream Angle_L(name_Angle_L);
  double Theta_L[Angle_Par];
  
  for(int i =0; i<Angle_Par;i++){
    double par1 =0.0;
    int p1 =0;
    Angle_L>>par1>>p1>>p1>>p1>>p1>>p1;
    Theta_L[i]=par1;
    theta_opt[i] = Theta_L[i];
  }
  Angle_L.close();
  //===================================
  // for phi information input
  ntune_event = 0;
  for(int i =0;i<Angle_Par;i++){
    phi_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  sprintf(name_angle_phi,"./phi_matrices/phi_LHRS_3rd_Opt_9.dat");  /// +++++++++++++ this part is ok
  ifstream angle_phi(name_angle_phi);
  double PHI_L[Angle_Par];
  for(int i =0; i<Angle_Par;i++){
    double par2 =0.0;
    int p2 =0;
    angle_phi>>par2>>p2>>p2>>p2>>p2>>p2;
    PHI_L[i] = par2;
    phi_opt[i]= PHI_L[i];
  }
  angle_phi.close();
  
  
  /* // for one dimensional histogram for example, counts vs theta or counts vs  phi
  // TH2F *h = new TH2F("h",";L.tr.tg_ph;L.tr.tg_th",200,-0.27,-0.19,200,-0.1,0.09); // for theta phi calculation
  //TH2F *h = new TH2F("h",";L.tr.tg_ph;L.tr.tg_th",200,-0.035,0.04,200,-0.08,0.08); // 2 D histogram
  TH1F *h = new TH1F("h",";L.tr.tg_th;Counts",500,-0.1,0.1); // 200 March 3
  // TH1F* h = new TH1F("h","Run # 111721;L.tr.tg_ph;counts",100,-0.04,0.04);// 1 D histogram
  TH1F* h1 = new TH1F("h1","Multi-C-foil target",750,-0.5,0.5); // -0.5, 0.5
  gStyle->SetOptStat(111111);
  h1->GetXaxis()->SetTitle("z_LHRS(m)");
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetRangeUser(-0.2,0.2); // range is  -0.2,0.2
  h1->GetYaxis()->SetTitle("number of counts/2.5mm");
  h1->GetYaxis()->CenterTitle();
  TH1F* h2[nfoil];
  TH1F *h3[nfoil];
  TH1F *h3_new[nrow][nfoil]; 
  char tempc[500];
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH1F(tempc,tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax());
    h2[i]->GetXaxis()->SetTitle("L.tr.tg_th(rad)");
    h2[i]->GetXaxis()->SetRangeUser(-0.2,0.2);
    
    sprintf(tempc,"h3_new_%d",i);
    h3[i] = (TH1F*)h2[i]->Clone(tempc);
    for(int row =0;row<nrow;row++){
      sprintf(tempc,"h3_new_%d",i);
      h3_new[row][i] = (TH1F*)h2[i]->Clone(tempc);
    }
    
  }
 */ 
  // TH2F* h = new TH2F("h","Multi-C-foil target",2000,-0.003,0,500,-0.1, 0.1); // new range
  TH2F *h = new TH2F("h",";Y_ss(m);X_ss(m)",1000,-0.06,0.06,400,-0.1,0.1); // for X_s vs Y_ss
  //  TH2F *h = new TH2F("h",";L.tr.tg_ph;L.tr.tg_th",1000,-0.06,0.06,500,-0.2,0.2); // 2 D histogram
  //  TH2F *h = new TH2F("h",";L.tr.tg_th;Counts",500,-0.003,0,500,-0.5,0.5); // TGT_X vs ph2[0]
 //  TH1F* h = new TH1F("h","Run # 111721;L.tr.tg_ph;counts",100,-0.1,0.1);// 1 D histogram
  TH2F* h1 = new TH2F("h1","Multi-C-foil target",2000,-0.5,0.5,500,-0.1, 0.1); // -0.5, 0.5 original
  gStyle->SetOptStat(11111);
  h1->GetXaxis()->SetTitle("phi(rad)");
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetRangeUser(-0.05,0.05); // range is  -0.2,0.2
  h1->GetYaxis()->SetTitle("theta(rad)");
  h1->GetYaxis()->CenterTitle();
  h1->GetYaxis()->SetRangeUser(-0.1,0.1);
  TH2F* h2[nfoil];
  TH2F *h3[nfoil];
  TH2F *h3_new[nrow][nfoil]; 
  char tempc[500];
  for(int i=0 ; i<nfoil ; i++){
    sprintf(tempc,"h2_%d",i);
    h2[i] = new TH2F(tempc,tempc,
		     h1->GetXaxis()->GetNbins(),
		     h1->GetXaxis()->GetXmin(),
		     h1->GetXaxis()->GetXmax(),
		     h1->GetYaxis()->GetNbins(),
		     h1->GetYaxis()->GetXmin(),
		     h1->GetYaxis()->GetXmax());
  h2[i]->GetXaxis()->SetTitle("L.tr.tg_ph(rad)");
  h2[i]->GetXaxis()->SetRangeUser(-0.05,0.05);
  h2[i]->GetYaxis()->SetTitle("L.tr.tg_th(rad)");
  h2[i]->GetYaxis()->SetRangeUser(-0.1,0.1);
  
  sprintf(tempc,"h3_new_%d",i);
  h3[i] = (TH2F*)h2[i]->Clone(tempc);
  for(int row =0;row<nrow;row++){
    sprintf(tempc,"h3_new_%d",i);
    h3_new[row][i] = (TH2F*)h2[i]->Clone(tempc);
  }
  
}

// up to here is the Histograms




  bool ltrig = false;
  
  
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    foil_flag[i] = -1; 
  }
  
  for (int i=0 ; i< ent ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp[j]  = -2222.0;
      l_th_fp[j] = -2222.0; 
      l_y_fp[j]  = -2222.0;
      l_ph_fp[j] = -2222.0;
    }
    trig1 = 0.0;
    ltrig = false;
    
    if(i+evshift<ent) t1->GetEntry(i+evshift);
    else t1->GetEntry(i-evshift);
    if(trig1>1.0) ltrig = true;
    else ltrig = false;
    
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    if(ltrig==true
       && fabs(XFP)  <2.0 
       && fabs(XpFP) <0.1
       && fabs(YFP)  <0.5
       && fabs(YpFP) <0.1
       && cer_asum>500){
      
      XFP  = (XFP-XFPm)/XFPr; 
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      ztR[0] =calcf2t_zt(Pzt_L, XFP, XpFP, YFP, YpFP); 
      
      ztR[0] = ztR[0] * Ztr + Ztm; 
      
      
      
      RasterCor  = calcRasterCor(Ras_x, Opt_Par[2], Opt_Par[0]); 
      
      TGT_x = -Lrb_y;
      
     
      ztR_wRC[0] = ztR[0] + RasterCor/tan(hrs_ang); // using raster parameter p0 and p1
     
      //ztR_wRC[0]  =  ztR[0] + Lrb/tan(hrs_ang); /// directly using the Lrb.x
      
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  ztR_wRC[0]);
      th2[0] = th2[0]*Xptr + Xptm;
      
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP,ztR_wRC[0]);
      ph2[0] = ph2[0]*Yptr + Yptm;
      
      
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      // =======additional variables for SS_X vs SS_Y
      DR_z = - ztR_wRC[0];
      DX_sy = -Lrb_y;
      X_ss = (R + DR_z)*th3[0] + DX_sy; // th2 our calculated th3 hall A calculated
      // fo ss y  position
      DY_sx = Lrb*cos(hrs_ang);
      DY_sz = - ztR_wRC[0]*sin(hrs_ang);
      DR_x = Lrb*sin(hrs_ang);
      DR_z1 = - ztR_wRC[0]*cos(hrs_ang);
      Y_ss = (R+ DR_z1 -  DR_x)*ph3[0] +  DY_sx +  DY_sz; // ph2 for our matrix and ph3 for hall A matrix


      //=====================================
      
      
      // h->Fill(Y_ss);
      h->Fill(Y_ss,X_ss); // for two  dimensionlal
      //  h->Fill(ph2[0],th2[0]); // for one dimensional
      
      tnew->Fill();
      
      
      
      bool foilfoilflag=false;
      bool row_flag = false;
      int foil_with_hit= -1; 
      for(int j=0 ; j<nfoil ; j++){
	if(fcent[j]-selection_width<ztR_wRC[0]
	   && ztR_wRC[0]< fcent[j]+selection_width){
	  h3[j]->Fill(ph2[0],th2[0]);// this part can be changed as th2 and ph2 
	  foilfoilflag=true; 
	  foil_with_hit=j;
	}
	else foilfoilflag=false;
	for(int row =0;row<nrow;row++){
	  if( foilfoilflag==true && theta[row][j]-row_width<th2[0] 
	      && th2[0]< theta[row][j]+row_width){
	    row_flag=true;
	    theta_flag_ = row;	 
	    
	    // h3_new[row][j]->Fill(th2[0]);
	    //  //  if(j+2!=10) h3_new[row][j]->SetLineColor(row+2);
	    //  else h3_new[row][j]->SetLineColor(2);
	    // h3_new[row][j]->SetLineStyle(9);
	  }
	}
	
	if(ntune_event<nmax && row_flag==true){
	  foil_flag[ntune_event] = foil_with_hit;// TG Feb 27
	  theta_flag[ntune_event] = theta_flag_;
	  
	  
	  x[ntune_event]  = XFP;
	  y[ntune_event]  = YFP;
	  xp[ntune_event] = XpFP;
	  yp[ntune_event] = YpFP;
	  theta_recon[ntune_event] = th2[0];
	  z_recon[ntune_event] = ztR[0];
	  ztR_wRC_[ntune_event] = ztR[0] + RasterCor; // z position with raster correction
	  Ras_curx[ntune_event] = Ras_x;
	  ntune_event++;
	  
	}
	
      }//int j
      
      
    }
    
  }
  
  
  
  tnew->Write(); //
  // fnew->Close(); 
  
  // =================================== 
  // ======== Draw histograms ========== 
  // =================================== 
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  h->Draw("colz");

  for(int i =0; i<nfoil;i++){
    for(int j = 0;j<nsshole;j++){
   
      mark[i][j]->SetMarkerColor(2);
      mark[i][j]->Draw("same");
    }
  }


  /*
   // to see the event selection plot
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h1->Draw();
  
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  c3->Divide(2,5);
  for(int i=0 ; i<nfoil ; i++){
    
    c3->cd(i+1);
    h3[i]->Draw();
    for(int row =0;row<nrow;row++){
      h3_new[row][i]->Draw("same");
    }
  }
  
  */
  
  // one can vas one histogram
  TCanvas *c4[nfoil];
  
  for(int i=0 ; i<nfoil ; i++){
    c4[i]= new TCanvas(Form("c4_%d",i),Form("c4_%d",i));
    c4[i]->cd(i+1);
    h3[i]->Draw("colz");
    for(int j =0;j<nsshole;j++){
      mark[i][j]->SetMarkerColor(2);
      mark[i][j]->Draw("same");
      // cout<< "The number of ss hole is" <<nsshole<<endl;
      //  cout<<" value of mark [i] is" <<mark[i]<<endl;
     
    }
  }
 
  //============================================
  // tuning starts here
  //============================================   
  
  const int nite =0;
  double temp[nite];
  double x[nite];
  
  if (nite>0) cout << " Tuning started: " << endl;
  
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    temp[i] = tune(theta_opt,i);
    // cout<<"This is the first location test. Up to here looks  fine"<<endl;
    sprintf(tempc, "./matrices/theta_3rd_LHRS_Opt_%d.dat",i); 
    ofstream * ofs = new ofstream(tempc); 
    int nppp = 0;
    cout<< "This is the second location for test"<<endl; 
    const int nn = 4; 
    for(int i=0; i<nn+1; i++){
      for(int e=0; e<nn+1; e++){
	for(int d=0; d<nn+1; d++){
	  for(int c=0; c<nn+1; c++){
	    for(int b=0; b<nn+1; b++){
	      for(int a=0; a<nn+1; a++){  
		if(a+b+c+d+e==i){
		  *ofs << theta_opt[nppp]
		       << " " << a 
		       << " " << b
		       << " " << c
		       << " " << d
		       << " " << e << endl;
		  nppp++; 
		  
		}
	      }
	    }
	  }
	}
      } // int e = 0
    }
    ofs->close();
    ofs->clear();
    
    
    
    
    
    cout << temp[i]<<endl; 
  }
  
  
  if(nite>0){
    TGraph * gr = new TGraph(nite,x,temp);  //+++++++++= closed on 21 Decemebr 2018
    TCanvas * c4 = new TCanvas("c4","",600,600); // ++++++++
    gr->Draw("*la"); //+++
  }
} //end of  main function

//////////////////////////////////////////////////
double calcf2t_zt(double* P, double xf, double xpf, 
		  double yf, double ypf)
//////////////////////////////////////////////////
{
  // -----3rd order ----- 
  // This is the third order claculation byb  using 35 parameter
  const int nMatT=3;  
  const int nXf=3;
  const int nXpf=3;
  const int nYf=3;
  const int nYpf=3;
  
  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0;
  
  for (int n=0;n<nMatT+1;n++){
    for (d=0;d<n+1;d++){
      for (c=0;c<n+1;c++){
	for (b=0;b<n+1;b++){
	  for (a=0;a<n+1;a++){ 
	    
	    if (a+b+c+d==n){
	      if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf){
		x = pow(xf,double(a))*pow(xpf,double(b))*
		  pow(yf,double(c))*pow(ypf,double(d));
	      }
	      else{
		x = 0.;
	      }
	      Y += x*P[npar]; 
	      npar++;
	     
	    }
	    
	  }
	}
      }    
    }
  }
   
  return Y; 
}
//_______________________________________________________________________________________________________________________________________________

//////////////////////////////////////////////////
double calcf2t_th(double* P, double xf, double xpf, 
		  double yf, double ypf,double zt)
//////////////////////////////////////////////////
{
  // -----4th order ----- 
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){ 
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
 
  return Y; 
}
  // =========================================================


//////////////////////////////////////////////////
  double calcf2t_ph(double* P, double xf, double xpf, 
		    double yf, double ypf,double zt) // This par is ok
//////////////////////////////////////////////////
{
  // -----4th order ----- 
  const int nMatT=4;  
  const int nXf=4;
  const int nXpf=4;
  const int nYf=4;
  const int nYpf=4;
  const int nZt=4;

  double Y=0.;
  double x=1.; 
  int npar=0;
  int a=0,b=0,c=0,d=0,e=0;
  
  for (int n=0;n<nMatT+1;n++){ 
    for(e=0;e<n+1;e++){
      for (d=0;d<n+1;d++){
	for (c=0;c<n+1;c++){ 
	  for (b=0;b<n+1;b++){
	    for (a=0;a<n+1;a++){ 
	      
	      if (a+b+c+d+e==n){
		if (a<=nXf && b<=nXpf && c<=nYf && d<=nYpf && e<=nZt){
		  x = pow(xf,double(a))*pow(xpf,double(b))*
		    pow(yf,double(c))*pow(ypf,double(d))*pow(zt,double(e));
		}
		else{
		  x = 0.;
		}
		Y += x*P[npar]; 
		npar++;
	      }
	      
	    }
	  }
	}    
      }
    }
  }// n = 
  
  return Y; 
}
// #############################################################
double tune(double* pa, int j) // tune fun defn
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10]; 
  int ierflg = 0;
  int allparam = Angle_Par; 
  
  TMinuit* minuit = new TMinuit(allparam);
  minuit->SetFCN(fcn);
  
  // ~~~ Chi-square ~~~~
  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);
  
  minuit -> SetPrintLevel(-1);
  double start[allparam];
  double step[allparam];
  double LLim[allparam];
  double ULim[allparam];
  char pname[500];
 
  for(int i=0 ; i<allparam ; i++){
    sprintf(pname,"param_%d",i+1);
    start[i] = pa[i];
    step[i] = 1.0e-3; 
    
    LLim[i] = pa[i] - pa[i]*0.8;
    ULim[i] = pa[i] + pa[i]*0.8; 
    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  arglist[0] = 2.0;
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ one of the way to get optimized parameter
  arglist[0] = 20000;
  arglist[1] = 0.01;
  minuit -> mnexcm("MINImize",arglist,2,ierflg); 
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double e;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(0,amin); 
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){
   
 
    minuit -> GetParameter(i,theta_opt[i],e);
  }
  
  return chi2; 
}

// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
  
{
  double chi2 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  const double sigma = 0.0045;
  
  double ztR;
  double th2;
  double refz = 0.0; 
  double residual = 0.0;
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;
    ztR = 0.0;
    th2 = 0.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    refz  = theta_real[theta_flag[i]][foil_flag[i]];// peak position of each foil +++++++++++++++++++++++++++++++++++
    
  
    XFP   =(XFP -XFPm)/XFPr; 
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    
    double ztR_wRC; 
    ztR_wRC = ztR_wRC_[i];
    
    th2 = calcf2t_th(param, XFP, XpFP, YFP, YpFP, ztR_wRC); // need the exact value
    th2 = th2*Xptr + Xptm;
    residual = th2-refz;
    chi2 = chi2 + pow(residual,2.0);
   
  }
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}

