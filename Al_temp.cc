// 11/20/2019
// copied from ht_optics.cc
// purpose is to do the Al analysis suggested by Dr. Tang on November 20, 2019
extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_ph(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_mom(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);



const double  XFPm=-0.7,  XpFPm=-0.15; // m is the mean from the old definition
const double  YFPm=-0.05, YpFPm=-0.18;
const double  Xptm=-0.07, Yptm=-0.2, Momm=1.74; // tm = target offset.. MOmm is the momentum offset
const double  XFPr=1.3,   XpFPr=0.27; // r is the scaling factor or range
const double  YFPr=0.1,   YpFPr=0.10; 
const double  Xptr=0.15,  Yptr=0.08, Momr=0.18; // tr is the target range
const double  PLm = 25.4, PLr=0.7; // m is the offset and PLr is the path laegth range
const double  Ztm = -0.15,Ztr=0.35; //Ztm  z position at target  point offset
extern void fcn(int &nPar, double* /*grad*/, 
		double &fval, double* param, int /*iflag*/);
extern double tune(double* pa, int j);

const int nfoil = 10;

const int nmax = 3000; 
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;
double Pzt_opt[nmax];

const int npeak = 2;
double Lambda_width[npeak] = {16.0,16.0}; //6.5,8.8}
double Lambda_cent[npeak] ={1117.48,1188.65};
double Lambda_real[npeak] ={1115.683,1192.642}; // Me
//========================================

//========================================


const int Angle_Par =126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];

const int Mom_Par = 252;
double momL_opt[nmax];
double momR_opt[nmax];
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
double ztR_wRC_[nmax];

const double me = 0.000511;
const double mk = 0.493677;
//////const double mp = 2.808944; // for tritium by DR. Tang
//////const double mp =25.1267; // Al target mass 
const double mp = 0.938272; // H target mass
//const double mp = 2.808921; // for tritium by Gogami Tritium target mass
const double mL = 1.115683;
extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);

void Al_temp(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
  
  TChain * t1 = new TChain("T"); 
  
  // t1->Add("./Rootfiles/NOV11_Rootfiles/T_221_830.root"); // T kinematics updated on November 13,2019. replayed by  by Ole 
  t1->Add("./Rootfiles/NOV11_Rootfiles/HT_552_716.root");// H data with tritium kinematics updated on November 13,2019. replayed by Ole
  //  t1->Add("./Rootfiles/NOV11_Rootfiles/H_149_542.root"); // H/H kinematics updated on November 13,2019. replayed by  by Ole

 
  double ent = t1->GetEntries();
  cout<<"entry in the t1=="<<ent<<endl;
  
  int evshift = 10; 
  
  const int max = 100;
  Double_t trig5[max]; 
  double ltime_s0[max];
  double ltime_s2[max];
 
  double mom2[max];
  double mom1[max];
  double momL[max];
  double momR[max];
  double DL[max];
  double DR[max];
  double pep_L[max];
  double pk_R[max];
  double theta_epk[max];
  const int f1n = 64;
  
  double lvz[max],rvz[max];// raster corrected
  double LVZ[max], RVZ[max]; // with out raster
  double th1[max], ph1[max];// RHRS angle 
  double th2[max], ph2[max]; 
  double th3[max], ph3[max];
  double th4[max], ph4[max];
  double Th_lhrs[max]; // to store the value of angle July 16
  double Ph_lhrs[max]; 
  double TH_L[max];// The theta angle in the spectrometer geometry
  double PH_L[max];
 
  double Th_rhrs[max]; // to store the value of angle July 16
  double Ph_rhrs[max]; 
  double TH_R[max];// The theta angle in the spectrometer geometry
  double PH_R[max];
  double delta_pep[max];     // target straggling
  double pep_real[max]; 
  double delta_pk[max];
  double pk_real[max];
  double par_ep[3];
  double par_k[3];
  double mm,mm_1st_cor;
  double Dpe_off =-0.0489648;
  double Dpk_off =-0.0836116;   
  

  double delta[max];
  double a1, a2;
  Int_t runnum; 
  double hallap;
 
  char tempc1[500]; 
  int temp1;
 
  // ++++++++++++++++++++++++++++++++++++++++++++++
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];

  double r_th_fp[max];
  double r_ph_fp[max];
  double r_x_fp[max];
  double r_y_fp[max];
  const int n = 16;
  
  double lbeta[max];
  double ctime; 
  double Lrb;
  double rpr;     
  double Ras_x; 
  double Ras_y;
  double cer_asum;
  double z_av[max];
  double z_av_1[max];
 
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("R.tr.p", &mom1);
  t1->SetBranchAddress("L.tr.tg_dp", &DL);
  t1->SetBranchAddress("R.tr.tg_dp", &DR);
  t1->SetBranchAddress("Lrb.x", &Lrb);

  
  t1->SetBranchAddress("Lrb.Raster2.rawcur.x", &Ras_x); 
  t1->SetBranchAddress("Lrb.Raster2.rawcur.y", &Ras_y);
  t1->SetBranchAddress("R.a1.asum_c", &a1);
  t1->SetBranchAddress("R.a2.asum_c", &a2);
  
  t1->SetBranchAddress("L.tr.vz", &LVZ); // no raster
  t1->SetBranchAddress("R.tr.vz", &RVZ);
  t1->SetBranchAddress("L.tr.tg_th", &th3);
  t1->SetBranchAddress("L.tr.tg_ph", &ph3);
  t1->SetBranchAddress("R.tr.tg_th", &th4);
  t1->SetBranchAddress("R.tr.tg_ph", &ph4);
  t1->SetBranchAddress("L.s0.time", &ltime_s0);
  t1->SetBranchAddress("L.s2.time", &ltime_s2);
  
  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);
   
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp);

  t1->SetBranchAddress("L.tr.beta",  &lbeta);
  t1->SetBranchAddress("L.cer.asum_c",  &cer_asum); 
  t1->SetBranchAddress("coin_time",  &ctime); 
  t1->SetBranchAddress("ztR_wRC",  &rvz);
  t1->SetBranchAddress("ztL_wRC",  &lvz);

  double ztL[max];
  double zL_RC[max];  
  double ztR[max]; // July 01, 2019
  double ZR_RC[max];   
  
  TFile* fnew = new TFile("./output_root/angle_lhrs.root","recreate"); 
  TTree* tnew = new TTree("tree","For z calibration (LHRS)");
 

  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D");
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D");
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D");
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D");
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D");
  tnew->Branch("L.tr.vz_TG",  &ztL,   "L.tr.vz_TG[100]/D");
  tnew->Branch("L.tr.vz_TG2",  &zL_RC, "L.tr.vz_TG2[100]/D"); 
  tnew->Branch("Lrb.x", &Lrb, "Lrb.x/D"); 
  tnew->Branch("L.cer.asum_c", &cer_asum, "L.cer.asum_c/D"); 
  tnew->Branch("L.tr.tg_th_TH2", &th2, "L.tr.tg_th_TH2[100]/D");
  tnew->Branch("L.tr.tg_ph_PH2", &ph2, "L.tr.tg_ph_PH2[100]/D");
  
  tnew->Branch("Lrb.Raster2.rawcur.x", &Ras_x, "Lrb.Raster2.rawcur.x/D"); 
  tnew->Branch("Lrb.Raster2.rawcur.y", &Ras_y, "Lrb.Raster2.rawcur.y/D"); 
  // tnew->Branch("RasterCor", &RasterCor, "RasterCor/D"); 
  
  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;

 
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Angle_Par; i++){
    thetaL_opt[i] = -2222.0;
  }
  
  char name_Angle_L[500]; 
  //  sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");
  sprintf(name_Angle_L,"./matrices/theta_L4th_4th_6.dat");// optimized on OCT 23, 2019

  ifstream Angle_L(name_Angle_L);
  double Theta_L[Angle_Par];    
  for(int i =0; i<Angle_Par;i++){
    double par1 =0.0;
    int p1 =0;
    Angle_L>>par1>>p1>>p1>>p1>>p1>>p1;
    Theta_L[i]=par1;
    thetaL_opt[i] = Theta_L[i];
    
  }
  Angle_L.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // for LHRS  phi information input
  //--------------------------------------------------------
  
  ntune_event = 0;
  for(int i =0;i<Angle_Par;i++){
    phiL_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  //  sprintf(name_angle_phi,"./matrices/phi_LHRS_3rd_Opt_9.dat");
  sprintf(name_angle_phi,"./matrices/phi_L4th_5th_5.dat");// optimized on OCT 23, 2019
  ifstream angle_phi(name_angle_phi);
  double PHI_L[Angle_Par];
  for(int i =0; i<Angle_Par;i++){
    double par2 =0.0;
    int p2 =0;
    angle_phi>>par2>>p2>>p2>>p2>>p2>>p2;
    PHI_L[i] = par2;
    phiL_opt[i]= PHI_L[i];
  }
  angle_phi.close();
  // LHRS momentum information========================July 20, 2019
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momL_opt[i] = -2222.0;
  }
  char name_Mom_lhrs[500];
   // sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_LHRS_5_5th_2.dat");
  sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_L5_4th_2.dat"); // MATRIX PRODUCED ON nOV, 19, 2019
  ifstream Mom_lhrs(name_Mom_lhrs);
  double mom_L[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par5 = 0.0;
    int p5 =0;
    Mom_lhrs>>par5>>p5>>p5>>p5>>p5>>p5;
    mom_L[i]= par5;
    momL_opt[i] = mom_L[i];
  }
  Mom_lhrs.close();
  // up to here thelhrs momentum matrix======================
  
  
  // =======RHRS theta input information
  ntune_event =0;
  for(int i =0;i<Angle_Par;i++){
    thetaR_opt[i] = -2222.0;
  }
  char name_Angle_R[500]; 
  sprintf(name_Angle_R,"./All_Matrices/xpt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream Angle_R(name_Angle_R);
  double Theta_R[Angle_Par];
  for(int i =0; i<Angle_Par;i++){
    double par3 =0.0;
    int p3 = 0;
    Angle_R>>par3>>p3>>p3>>p3>>p3>>p3;
    Theta_R[i]=par3;
    thetaR_opt[i] = Theta_R[i];
  }
  Angle_R.close();
  //====================================================
  //=======RHRS phi input information===============
  ntune_event = 0;
  for(int i =0;i<Angle_Par;i++){
    phiR_opt[i] = -2222.0;
  }
  char name_phi_Rhrs[500];
  sprintf(name_phi_Rhrs,"./All_Matrices/ypt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
   ifstream phi_Rhrs(name_phi_Rhrs);
  double PHI_R[Angle_Par];
  for(int i =0; i<Angle_Par;i++){
    double par4 =0.0;
    int p4 =0;
    phi_Rhrs>>par4>>p4>>p4>>p4>>p4>>p4;
    PHI_R[i] = par4;
    phiR_opt[i]= PHI_R[i];
  }
  phi_Rhrs.close();
  //==================================================
  // =====RHRS momentum recon==========================
 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];
  sprintf(name_Mom_rhrs,"./MOM_MATRICES/mom_R5_8th_2.dat"); // matrix prodeced on the Nov 19, 2019
  ifstream Mom_rhrs(name_Mom_rhrs);
  double mom_R[Mom_Par];
  for(int i = 0; i<Mom_Par;i++){
    double par6 = 0.0;
    int p6 =0;
    Mom_rhrs>>par6>>p6>>p6>>p6>>p6>>p6;
    mom_R[i]= par6;
    momR_opt[i] = mom_R[i];
  }
  Mom_rhrs.close();
  // =====RHRS momentum recon up to here=============== 

  // new replay histograms
  TH1F *h20 = new TH1F("h20","Background from accidental ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200 ); //H/H 250 and H/T 150  and 130 for T data Aug 19, 2019
  TH1F *h21 = new TH1F("h21","nnL spectrum ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",210,-150,200); // 232 
  TH1F *h24 = new TH1F("h24"," ; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);
  TH1F *h26 = new TH1F("h26","Background from Al real;-B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);
  TH1F *h28 = new TH1F("h28","; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);
  TH1F *h29 = new TH1F("h29","; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);
  TH1F *hal = new TH1F("hal","; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);
  TH1F *h50 = new TH1F("h50","; -B_{#Lambda}(MeV);Counts/2.0 MeV ",175,-150,200);

  TH1F *h36 = new TH1F("h36"," H in T kinematics ; missing mass(MeV/c^{2});Counts/2.5 MeV ",160,900,1300); // to see H in T kinematics
  TH1F *h30 = new TH1F("h30"," ;",100,-0.3,0.3); 

  /// old replay
  TH1F *hz = new TH1F("hz"," ;Z(LHRS) ",150, -0.3,0.3); // HT data
  // TH1F *h7 = new TH1F("h7"," ;;Counts ",200,-0.25,0.25); 
  TH1F *H1 = new TH1F("H1"," ;Missing Mass(MeV/c^{2});Counts/2MeV ",125, 1025,1275); // HT data
  TH1F *H = new TH1F("H"," ;Missing Mass(MeV/c^{2});Counts/MeV ",250, 1025,1275); //H/H 250 and H/T 150  and 130 for T data Aug 19, 2019
  TH1F *ht = new TH1F("ht","; -B_{#Lambda}(MeV);Counts/1.5 MeV ",200,-100,200); // for tritium
  TH1F *hb = new TH1F("hb","; -B_{#Lambda}(MeV);Counts/1.5 MeV ",200,-100,200); // for tritium
  TH1F *h = new TH1F("h","Z average < -0.115(m) OR Z average > 0.115(m); -B_{#Lambda}(MeV);Counts/1.5 MeV ",200,-100,200); // for tritium
  TH1F *h2 = new TH1F("h2"," ;",200,-100,200); 
  TH1F *h4 = new TH1F("h4"," ;",200,-100,200);
  TH1F *hf = new TH1F("hf","; -B_{#Lambda}(MeV);Counts/1.5 MeV ",200,-100,200); // for tritium
 
  H->GetXaxis()->CenterTitle(); 
  H->GetYaxis()->CenterTitle(); 
  gStyle->SetOptStat(11111);
  TH1F * h3 = new TH1F("h3", "Z average < -0.115(m) OR Z average > 0.115(m); -B_{#Lambda}(MeV);Counts/1.5 MeV ", 200,-100,200);
  h3->GetXaxis()->CenterTitle(); 
  h3->GetYaxis()->CenterTitle(); 
  gStyle->SetOptStat(111111);
 
  
  char tempc[500];
 
  //===================================================================  
 
  bool ltrig = false;
  bool rtrig = false; 
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
      th1[j] = -2222.0;
      th2[j] = -2222.0;
      ph1[j] =-2222.0;
      ph2[j] =-2222.0;
      Th_lhrs[j] = -2222.0; // July 16, 2019
      TH_L[j] = -2222.0;
      Th_rhrs[j] = -2222.0;
      TH_R[j] = -2222.0;
      delta[j] = -2222.0;
      theta_epk[j]=-2222.0;
      delta_pep[j]= -2222.0;
      pep_real[j] =-2222.0;
      delta_pk[j]= -2222.0;
      pk_real[j] = -2222.0;
     
      r_x_fp[j]  = -2222.0;
      r_th_fp[j] = -2222.0;
      r_y_fp[j]  = -2222.0;
      r_ph_fp[j] = -2222.0;
    
      
      
      trig5[j] = 0.0;
      rtrig = false;
    }
   
    trig5[0] = 0.0;
    rtrig = false;
   
    t1->GetEntry(i); 
   
   
    if(trig5[0]>1.0) rtrig = true; //JUly 01, 2019
    else rtrig = false;

    z_av[0] = (lvz[0] + rvz[0])/2.0;
    z_av_1[0] =  z_av[0];
   

    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];

 
    //  z_av[0] =(z_av[0]- Ztm)/Ztr;	
    // use z_av[0] in stead of individual z
    //  if(rtrig==true &&  fabs(ctime)<1.0 && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.10){
   
  
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;

      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;
   
      z_av[0] =(z_av[0]- Ztm)/Ztr;	
      //  cout<<"the value of a_av is " << z_av[0] << "and that of z-az_a[0] is  "<< z_av_1[0] <<endl;

      // LHRS angle and momentum calculation
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; // reconstructed theta LHRS	
    
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;
 
     
     
      momL[0] =  calcf2t_mom(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]);
      momL[0] =momL[0]*Momr + Momm; 
    
      //   momL[0] = momL[0]*2.218/2.10; 
     

      // Target struggling LHRS step #7 
      if( z_av_1[0]>0.08)
	{delta_pep[0] = 6.23409e-3*ph2[0] + 4.03363e-1;} 
      else	
	{delta_pep[0] = -1.35758*sin(-4.59571* ph2[0]) + 2.09093;} 
      
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; 
      
     
        // RHRS angle and momentum calculation
      th1[0] = calcf2t_th(Theta_R, R_XFP, R_XpFP, R_YFP, R_YpFP,   z_av[0]);
      th1[0] = th1[0]*Xptr + Xptm;
   
      ph1[0] = calcf2t_ph(PHI_R, R_XFP, R_XpFP, R_YFP, R_YpFP, z_av[0]);
      ph1[0] = ph1[0]*Yptr + Yptm;
        
    
      momR[0] =  calcf2t_mom(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;
        
      // target struggling step #11	
      if(z_av_1[0]>0.08) 
	{delta_pk[0] = 3.158e-2*ph1[0] + 4.05819e-1;} 
      else 
	{delta_pk[0] =-1.31749*sin(-4.61513* ph1[0]) + 2.03687;}
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; 
     
      // ==========================================================

      // missing mass calculation==============================
      par_ep[0] = pep_real[0]; 
      par_ep[1] = th2[0]; 
      par_ep[2] = ph2[0];

      par_k[0] = pk_real[0];  
      par_k[1] = th1[0];; 
      par_k[2] = ph1[0];
     
      hallap = hallap-0.1843;
      hallap = hallap/1000.0; 
    
 
      mm = CalcMM(hallap, par_ep, par_k, mp); 
      
      mm = (mm)*1000.;
      if(rtrig==true &&  fabs(ctime)<1.0 && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.10){
    
      H->Fill(mm);
      H1->Fill(mm);
      mm = mm -2994.814;
      h21->Fill(mm); 
     //  if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0]) < 0.10 ){
    //   h36->Fill(mm); //H in T kinematics
     
    // }

   
     
      /////// To see the background undr the real Lambda be spectrum on Tritium data
      /////   mm = mm -25312.4; // for Al study when I consider real Al mass
       // mm = mm -2994.814; // for tritium data only By TOSHI when consider the tritium mass
      // mm = mm -1115.683; // for H/H and tritium kinematics
      ///// new bg analysis
     //  if(rtrig==true && ((ctime>-49.445 && ctime< -7.063) ||(ctime>15.135 && ctime< 48.70)) /*&& fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.10 */){ 
     //   h20->Fill(mm); // for back ground analysis
     // }

      //if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0]) < 0.10 ){
       // h21->Fill(mm); 
     
      //  }

     //  if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.045 && (z_av[0]<-0.10 || z_av[0]>0.10) ){
     // 	h26->Fill(mm);  // for the Al Analysis Nov, 21, 2019
     
     //  }

     //  if(rtrig==true && ((ctime>-49.445 && ctime< -7.063) ||(ctime>15.135 && ctime< 48.70)) && fabs(lvz[0]-rvz[0])<0.045  &&(z_av[0]<-0.10 || z_av[0]>0.10)){ 
     //   h50->Fill(mm); // for back ground analysis
     // }



  

      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      
      R_XFP  = R_XFP*XFPr +XFPm ; 
      R_XpFP = R_XpFP*XpFPr+XpFPm;
      R_YFP  = R_YFP*YFPr+ YFPm;
      R_YpFP = R_YpFP*YpFPr +YpFPm;
     
      z_av[0] =z_av[0]*Ztr + Ztm;
   
    
      tnew->Fill();
      
       bool lambdaflag=false;  
      int peak_with_hit= -1; 
      for(int j=0 ; j<npeak ; j++){
	if(Lambda_cent[j]-Lambda_width[j]<mm
	   &&mm < Lambda_cent[j]+Lambda_width[j]){	 
	  
	  lambdaflag=true;
	  peak_with_hit=j; 
	 

	}
	else lambdaflag=false;  
	
	if(ntune_event<nmax && lambdaflag==true){
	  foil_flag[ntune_event] = peak_with_hit;
	  
	  
	  
	  x[ntune_event]  = XFP;
	  y[ntune_event]  = YFP;
	  xp[ntune_event] = XpFP;
	  yp[ntune_event] = YpFP;
	  
	  ntune_event++;
	  
	}
	
      }//int j	
  }
    
}
  
  tnew->Write(); 
  // fnew->Close(); 
  
  // =================================== 
  // ======== Draw histograms ========== 
  //  ===================================
  
  // TCanvas* cz = new TCanvas("cz","cz",600,600);
  // cz->cd();
  // hz->Draw();

  //// TCanvas* c7 = new TCanvas("c7","c7",600,600);
  // h7->Draw();

  //H/H data
  // TF1 *f1 = new TF1("f1","gaus",1111.4,1119.97);
  // TF1 *f2 = new TF1("f2","gaus",1188.38,1196.39);
  
  // TCanvas* c2 = new TCanvas("c2","c2",600,600);
  // c2->cd();
  // H->Draw();
  // f1->SetLineWidth(1);
  // f2->SetLineWidth(1);
  // H->Fit("f1","","",1111.4,1119.97);
  // H->Fit("f2","","",1188.38,1196.39);
  // f1->Draw("same"); 
    

  // TLatex l;
  // l.SetTextSize(0.025);
  // l.DrawLatex(1130,60,Form("#Lambda"));
 
  // l.DrawLatex(1130,65,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  // l.DrawLatex(1130,70,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  // l.DrawLatex(1200,20,Form("#Sigma^{0}"));
  // l.DrawLatex(1200,25,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  // l.DrawLatex(1200,30,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));


  //  // For H data with T kinematics
  TF1 *f1_2 = new TF1("f1_2","gaus",1110.99,1120.62);
  TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  c2_2->cd();
  H1->Draw();
  
  f1_2->SetLineWidth(1);
  H1->Fit("f1_2","","",1110.99,1120.62);
 
  TLatex l2;
  l2.SetTextSize(0.025);
  l2.DrawLatex(1130,32,Form("#Lambda"));
  
  l2.DrawLatex(1130,36,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  l2.DrawLatex(1130,40,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));
  

 
  
  // h->Fit("f1","","",-0.1362,-0.11509);
  // h->Fit("f2","","",0.1120,0.1334);
  /// f1->Draw("same");
  // TCanvas* c3 = new TCanvas("c3","c3",600,600);
  //  c3->cd();
  // h5->Draw();
  // TH1F *h6 = (TH1F*) h->Clone();
  //  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  // c4->cd();
  //  h6->Add(h5);
  // h6->Scale(0.5);
  //  h6->Draw();
  //  h6->Fit("f1","","",-0.1362,-0.11509);
  // h6->Fit("f2","","",0.11386,0.13541);
  // f1->Draw("same");


 
  //////++++++++++++++++++for new  background analysis on November 14, 2019
  // TCanvas* c20 =  new TCanvas("c20","c20", 600,600); //background from accidental
  // c20->cd();
  // h20->Draw();
  // h20->Scale(1.0/38.0);
  // TH1F *h22 = (TH1F*)h20->Clone();
  // TCanvas* c21 =  new TCanvas("c21","c21", 600,600); // real tritium spectrum
  // c21->cd();
  // h21->Draw();
  // TH1F* h23 = (TH1F*)h21->Clone();
  // h20->Draw("E2 same");
  // h20->SetFillStyle(3002);
  // h20->SetMarkerStyle(28);
  // h20->SetMarkerColor(kGreen);

  // TCanvas* c26 =  new TCanvas("c26","c26", 600,600); // for Al analysis before scaling
  // c26->cd();
  // h26->Draw();
  // //// TH1F *hal = (TH1F*)h26->Clone();
  // h50->Draw("E2 same");
  // h50->Scale(1.0/38.0);
  // h50->SetFillStyle(3002);
  // h50->SetMarkerStyle(28);
  // h50->SetMarkerColor(kGreen);

  // TCanvas* c16 =  new TCanvas("c16","c16", 600,600); // for Al analysis before scaling
  // c16->cd(); 
  // hal->Add(h26,h50,1.0,-1.0);
  // TH1F *h16 =(TH1F*)hal->Clone();
  // hal->Draw();




  // TCanvas* cal =  new TCanvas("cal","cal", 600,600); // for Al analysis scaled
  // cal->cd();  
  // h16->Draw();
  // h16->Scale(72.0/797.0);
  // // hal->Scale(12.0/585.0);
  // TH1F *h27 = (TH1F*)h16->Clone();


  // TCanvas *c27 = new TCanvas("c27", "c27", 600,600);
  // c27->cd();
  // h21->Draw();
  // h22->Draw("E2 same");
  // h22->SetFillStyle(3002);
  // h22->SetMarkerStyle(28);
  // h22->SetMarkerColor(kGreen);

  // h27->Draw("E2 same");
  // h27->SetFillStyle(3002);
  // h27->SetMarkerStyle(28);
  // h27->SetMarkerColor(kRed);

  // TCanvas *c28 = new TCanvas("c28", "c28", 600,600); // total bg grom Al + accidental
  // c28->cd();
  // h21->Draw();
  // h28->Add(h22,h27,1.0,1.0);
  // h28->Draw("E2 same");
  // h28->SetFillStyle(3002);
  // h28->SetMarkerStyle(28);
  // h28->SetMarkerColor(kGreen);
 
  
  // TCanvas *c29 = new TCanvas("c29", "c29", 600,600); // total bg grom Al + accidental
  // c29->cd();
  // h29->Add(h21,h28,1.0,-1.0);
  // h29->Draw();
  // TH1F *hff = (TH1F*)h29->Clone();

 
  // TCanvas *cff = new TCanvas("cff", "cff", 600,600); //to print the error bars
  // cff->cd();
  // hff->Draw("e");

 


 // TCanvas* c22 =  new TCanvas("c22","c22", 600,600);
  // c22->cd();
  // h24->Add(h23,h22,1.0,-1.0);
  // h24->Draw();
  // TH1F *h25 = (TH1F*)h24->Clone();
  // TCanvas* c23 =  new TCanvas("c23","c22", 600,600);
  // c23->cd();
  // h25->Draw("e");


  // //// H in T kinematics
  // TCanvas* c36 =  new TCanvas("c36","c36", 600,600);
  // c36->cd();
  // h36->Draw();

  // TCanvas* c30 =  new TCanvas("c30","c30", 600,600);
  // c30->cd();
  // h30->Draw();





  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const int nite =0;
  double temp[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    temp[i] = tune(thetaL_opt,i);

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
		  *ofs << thetaL_opt[nppp] 
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
    
   
    cout << temp[i]<<endl; // can see the chi2 values after the tunung by this way +++++++++++++++ 12/12/18  
  }    
  
  if(nite>0){
    TGraph * gr = new TGraph(nite,x,temp);  
    TCanvas * c4 = new TCanvas("c4","",600,600); 
    gr->Draw("*la"); 
  }
} //end of  main function

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
// ////////////////////////////////////////////////
//////////////////////////////////////////////////
double calcf2t_ph(double* P, double xf, double xpf, 
		  double yf, double ypf, double zt)
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

//////////////////////////////////////////////////
double calcf2t_mom(double* P, double xf, double xpf, 
		  double yf, double ypf, double zt)
//////////////////////////////////////////////////
{
  // -----4th order -----   
  const int nMatT=5;  
  const int nXf=5;
  const int nXpf=5;
  const int nYf=5;
  const int nYpf=5;
  const int nZt=5;

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


// missing mass function definition====================
double CalcMM(double ee, double* par_ep, double* par_k, double mt){
  double pe = ee; 
  double Ee = sqrt(me*me + pe*pe); 
  //  Ee = Ee +0.01188;// GeV
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = par_ep[0]; 
  double xpep = -par_ep[1];
  double ypep = -par_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  vec_ep.RotateX(hrs_ang);
  double Eep = sqrt(pep*pep + me*me);  

  double pk  = par_k[0];
  double xpk = -par_k[1];
  double ypk = -par_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
  vec_k.RotateX(-hrs_ang);
  double Ek = sqrt(pk*pk + mk*mk);
 
  double missingE2 = 0.0, missingP2 = 0.0, missingM2 = 0.0;
  missingE2 = pow(Ee + mt - Ek - Eep, 2.0); 
  missingP2 = (vec_e - vec_ep - vec_k) * (vec_e - vec_ep - vec_k);
  missingM2 = missingE2 - missingP2;
  
  double MissingMass = 0.0;
  MissingMass = sqrt(missingM2);
  return MissingMass;  
}


//############### up to hear missing mass #####################

// #############################################################
double tune(double* pa, int j) // tune fun defn
// #############################################################
{
  double chi2 = 0.0;
  double arglist[10];
  int ierflg = 0;
  int allparam = Mom_Par; 
  //cout << allparam << endl;
  TMinuit* minuit = new TMinuit(allparam); 
  minuit->SetFCN(fcn); // very imp function setying for chi square
  
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
    minuit -> GetParameter(i,thetaL_opt[i],e);
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
  
  double ztL;
  double th2;
  double refz = 0.0; 
  double residual = 0.0;
  
  for(int i=0 ; i<ntune_event ; i++){
    residual = 0.0;
    refz = 0.0;
    ztL = 0.0;
    th2 = 0.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    //  refz  = theta_real[theta_flag[i]][foil_flag[i]];    
   
    XFP   =(XFP -XFPm)/XFPr;  
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
   
    // double zL_RC; 
    // zL_RC = ztR_wRC_[i];
    // th2 = calcf2t_th(param, XFP, XpFP, YFP, YpFP, zL_RC); 
    // th2 = th2*Xptr + Xptm;
    // residual = th2-refz;
    chi2 = chi2 + pow(residual,2.0);
   
  }
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}
