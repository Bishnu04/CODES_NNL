// COpied from COMBINED.cc sept 30, 2019
// The purpose is to use the fifth ordr momentum matrix for tune
// This code is for the optimization of the matrix parameters
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


const int nmax = 3000; // can go up to 3000
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;

const int npeak = 2;
double Lambda_width[npeak] = {2.93, 2.74}; //6.5,8.8}
double Lambda_cent[npeak] ={1115.71,1192.37};
double Lambda_real[npeak] ={1115.683,1192.642}; // Mev



double p10[nmax],p11[nmax],p12[nmax];
double p13[nmax],p14[nmax],p15[nmax];
double p16[nmax],p17[nmax],p18[nmax],p19[nmax];
double phir[nmax];
double phil[nmax];
// ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
const int nmax_2 = 2400;
double x_2[nmax_2], y_2[nmax_2]; 
double xp_2[nmax_2], yp_2[nmax_2];
double z_recon_2[nmax_2];
int foil_flag_2[nmax_2];

const int npeak_2 = 1;
double Lambda_width_2[npeak_2] = {2.86}; //6.5,8.8}
double Lambda_cent_2[npeak_2] ={1115.79};
double Lambda_real_2[npeak_2] ={1115.683}; // Mev



double p10_2[nmax_2],p11_2[nmax_2],p12_2[nmax_2];
double p13_2[nmax_2],p14_2[nmax_2],p15_2[nmax_2];
double p16_2[nmax_2];
double phir_2[nmax_2];
double phil_2[nmax_2];
int ntune_event_2 = 0;

//))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))



//========================================

const int Total_Par = 126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];
double momL_opt[nmax];
double momR_opt[nmax];
const int Mom_Par = 252;
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;
extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);

void HT_optics(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
 
  TChain * t1 = new TChain("T");  
  TChain * t2 = new TChain("T"); 

  // t1->Add("./Rootfiles/July22_Rootfiles/H2_HAug06.root"); // for tuning optics use this data updated on Aug 06, 2019 old replY
  // t2->Add("./Rootfiles/July22_Rootfiles/H22_T_Aug05.root");// H data with tritium kinematics OLD REPLAY
  
  t1->Add("./Rootfiles/NOV11_Rootfiles/H_149_542.root"); // H/H kinematics updated on November 13,2019. replay by Ole
  t2->Add("./Rootfiles/NOV11_Rootfiles/HT_552_716.root");// H data with tritium kinematics updated on November 13,2019. replay by Ole
  
  double ent = t1->GetEntries();
  double ent_2 = t2->GetEntries();
  cout<<"entry in the t1=="<<ent<<endl;
  cout<<"entry in the t2=="<<ent_2<<endl;
  
  
  const int max = 100;
  Double_t trig5[max]; // JUly 01, 2019 
  double momL[max];
  double momR[max]; 
  
  double lvz[max],rvz[max];// raster corrected 
  double th1[max], ph1[max];// RHRS angle 
  double th2[max], ph2[max];  
  double delta_pep[max];     // target straggling
  double pep_real[max]; 
  double delta_pk[max];
  double pk_real[max];
  double par_ep[3];
  double par_k[3];
  double mm;  
  double hallap;
  
  double l_th_fp[max];
  double l_ph_fp[max];
  double l_x_fp[max];
  double l_y_fp[max];

  double r_th_fp[max];
  double r_ph_fp[max];
  double r_x_fp[max];
  double r_y_fp[max];
  const int n = 16; 
  double ctime; 
 
 
  double z_av[nmax];
  double z_av_1[nmax];
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  Double_t trig5_2[max]; // JUly 01, 2019 
  double momL_2[max];
  double momR_2[max]; 
  
  double lvz_2[max],rvz_2[max];// raster corrected 
  double th1_2[max], ph1_2[max];// RHRS angle 
  double th2_2[max], ph2_2[max];  
  double delta_pep_2[max];     // target straggling
  double pep_real_2[max]; 
  double delta_pk_2[max];
  double pk_real_2[max];
  double par_ep_2[3];
  double par_k_2[3];
  double mm_2;  
  double hallap_2;
  
  double l_th_fp_2[max];
  double l_ph_fp_2[max];
  double l_x_fp_2[max];
  double l_y_fp_2[max];

  double r_th_fp_2[max];
  double r_ph_fp_2[max];
  double r_x_fp_2[max];
  double r_y_fp_2[max];
  double ctime_2; 
 
  double z_av_2[nmax];
  double z_av_1_2[nmax];

  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
 
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);  
  t1->SetBranchAddress("L.tr.x",   &l_x_fp);
  t1->SetBranchAddress("L.tr.y",   &l_y_fp);
  t1->SetBranchAddress("L.tr.th",  &l_th_fp);
  t1->SetBranchAddress("L.tr.ph",  &l_ph_fp);   
  t1->SetBranchAddress("R.tr.x",   &r_x_fp);
  t1->SetBranchAddress("R.tr.y",   &r_y_fp);
  t1->SetBranchAddress("R.tr.th",  &r_th_fp);
  t1->SetBranchAddress("R.tr.ph",  &r_ph_fp); 
  t1->SetBranchAddress("coin_time",  &ctime); 
  t1->SetBranchAddress("ztR_wRC",  &rvz);
  t1->SetBranchAddress("ztL_wRC",  &lvz);
  // ((((((((((((((((((((((((((((((((   t2   (((((((((((((((((((((((((((((((((((((((((((((((((((
  t2->SetBranchAddress("HALLA_p", &hallap_2);  
  t2->SetBranchAddress("DR.T5", &trig5_2);  
  t2->SetBranchAddress("L.tr.x",   &l_x_fp_2);
  t2->SetBranchAddress("L.tr.y",   &l_y_fp_2);
  t2->SetBranchAddress("L.tr.th",  &l_th_fp_2);
  t2->SetBranchAddress("L.tr.ph",  &l_ph_fp_2);   
  t2->SetBranchAddress("R.tr.x",   &r_x_fp_2);
  t2->SetBranchAddress("R.tr.y",   &r_y_fp_2);
  t2->SetBranchAddress("R.tr.th",  &r_th_fp_2);
  t2->SetBranchAddress("R.tr.ph",  &r_ph_fp_2); 
  t2->SetBranchAddress("coin_time",  &ctime_2); 
  t2->SetBranchAddress("ztR_wRC",  &rvz_2);
  t2->SetBranchAddress("ztL_wRC",  &lvz_2);



  //))))))))))))))))))))))))))))))))   t2   ))))))))))))))))))))))))))))))))))))))))))))))))))))))
   
  TFile* fnew = new TFile("./output_root/angle_lhrs.root","recreate"); 
  TTree* tnew = new TTree("tree","For z calibration (LHRS)");
  tnew->Branch("HALLA_p", &hallap,"HALLA_p/D");
  tnew->Branch("L.tr.vz", &lvz, "L.tr.vz[100]/D");
  tnew->Branch("L.tr.x",   &l_x_fp, "L.tr.x[100]/D");
  tnew->Branch("L.tr.y",   &l_y_fp, "L.tr.y[100]/D");
  tnew->Branch("L.tr.th",  &l_th_fp,"L.tr.th[100]/D");
  tnew->Branch("L.tr.ph",  &l_ph_fp,"L.tr.ph[100]/D");  
  tnew->Branch("L.tr.tg_th_TH2", &th2, "L.tr.tg_th_TH2[100]/D");
  tnew->Branch("L.tr.tg_ph_PH2", &ph2, "L.tr.tg_ph_PH2[100]/D");
  
 
  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;

  // (((((((((((((((((((((((((((((((((((((((((((
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double R_XFP_2, R_XpFP_2; 
  double R_YFP_2, R_YpFP_2;

  //)))))))))))))))))))))))))))))))))))))))))))))

 
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Total_Par; i++){
    thetaL_opt[i] = -2222.0;
  }
  
  char name_Angle_L[500]; 
  // sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");//theta_3rd_LHRS_Opt_7.dat
  sprintf(name_Angle_L,"./matrices/theta_L4th_4th_6.dat");// optimized on OCT 23, 2019
  ifstream Angle_L(name_Angle_L);
  double Theta_L[Total_Par];    
  for(int i =0; i<Total_Par;i++){
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
  for(int i =0;i<Total_Par;i++){
    phiL_opt[i] = -2222.0;
  }
  char name_angle_phi[500];
  //  sprintf(name_angle_phi,"./matrices/phi_LHRS_3rd_Opt_9.dat");
  sprintf(name_angle_phi,"./matrices/phi_L4th_5th_5.dat");// optimized on OCT 23, 2019  
   ifstream angle_phi(name_angle_phi);
  double PHI_L[Total_Par];
  for(int i =0; i<Total_Par;i++){
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
  
 
  // sprintf(name_Mom_lhrs,"./All_Matrices/mom_LHRS_5_upto3.dat"); 
  // sprintf(name_Mom_lhrs,"./All_Matrices/mom30_L_6th_1.dat");
  // sprintf(name_Mom_lhrs,"./All_Matrices/momLH_3rd_1.dat");
  // sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_LHRS_5_0.dat");// orignal first 5th order matrix Nov 15, 2019
  sprintf(name_Mom_lhrs,"./MOM_MATRICES/mom_L5_4th_2.dat");
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
  for(int i =0;i<Total_Par;i++){
    thetaR_opt[i] = -2222.0;
  }
  char name_Angle_R[500]; 
  sprintf(name_Angle_R,"./All_Matrices/xpt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
 
  ifstream Angle_R(name_Angle_R);
  double Theta_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
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
  for(int i =0;i<Total_Par;i++){
    phiR_opt[i] = -2222.0;
  }
  char name_phi_Rhrs[500];
 
  sprintf(name_phi_Rhrs,"./All_Matrices/ypt_RHRS_4_upto2.dat"); //This is the RHRS Phi matrix optimized by Gogami with ss data
  ifstream phi_Rhrs(name_phi_Rhrs);
  double PHI_R[Total_Par];
  for(int i =0; i<Total_Par;i++){
    double par4 =0.0;
    int p4 =0;
    phi_Rhrs>>par4>>p4>>p4>>p4>>p4>>p4;
    PHI_R[i] = par4;
    phiR_opt[i]= PHI_R[i];
  }
  phi_Rhrs.close();
  //==================================================
  // =====RHRS momentum recon==========================6
 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];
 
  //sprintf(name_Mom_rhrs,"./All_Matrices/MOMR_30th_5.dat");
  //  sprintf(name_Mom_rhrs,"./All_Matrices/mom_RHRS_5_upto3.dat");
  //  sprintf(name_Mom_rhrs,"./All_Matrices/mom30_R_2nd_0.dat");
  // sprintf(name_Mom_rhrs,"./All_Matrices/momRH_3rd_4.dat");
  //  sprintf(name_Mom_rhrs,"./MOM_MATRICES/mom_RHRS_5_0.dat"); // matrix prodeced on the Nov 15, 2019
  sprintf(name_Mom_rhrs,"./MOM_MATRICES/mom_R5_8th_2.dat"); // matrix prodeced on the Nov 15, 2019
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
 
 
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts/MeV ",250,1025,1275); //H/H 245 and H/T 225
  TH1F *h_2 = new TH1F("h_2"," ;Missing Mass(MeV/c^{2});Counts/2MeV ",125,1025,1275); //H/H 245 and H/T 225
  TH1F *h_1 = new TH1F("h_1"," ",120,1000, 1300);
  TH2F * h10 = new TH2F("h10","H data with H kinematics;Missing Mass(MeV/c^{2});L.tr.th",200,1000,1300, 200,1.2, 2.5);
  h->GetXaxis()->CenterTitle(); 
  h->GetYaxis()->CenterTitle(); 
  h_2->GetXaxis()->CenterTitle();
  h_2->GetYaxis()->CenterTitle();
  TH2F *h5 = new TH2F("h5",";Z(LHRS);Z(RHRS)",500,-0.8,0.8,500,-20,20); 
  TH1F *h6 = new TH1F("h6",";Coin_Time(ns);Counts",150,-0.2,0.2); 
  char tempc[500];
  // ======================================================
  
  //===================================================================  
 
  
  bool rtrig = false; 
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    z_av[i] = -2222.0;
    z_av_1[i] = -2222.0;
    phir[i] = -2222.0;
    phil[i] = -2222.0;

    foil_flag[i] = -1;
  }

  // ((((((((((((((((((((((((((((((((((((((((((((

 bool rtrig_2 = false; 
  for(int i=0 ; i<nmax_2 ; i++){
    x_2[i]    = -2222.0; 
    y_2[i]    = -2222.0; 
    xp_2[i]   = -2222.0;
    yp_2[i]   = -2222.0;
    z_av_2[i] = -2222.0;
    z_av_1_2[i] = -2222.0;
    phir_2[i] = -2222.0;
    phil_2[i] = -2222.0;
    foil_flag_2[i] = -1;
  }

  //)))))))))))))))))))))))))))))))))))))))))))))

  // +++++++++++++++++++++++++ for t1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    
    if(rtrig==true &&  fabs(ctime)<1.0  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.10){     
      
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;
      
      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;

      z_av[0] =(z_av[0]- Ztm)/Ztr;
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; 
     
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;
          
     
      momL[0] =  calcf2t_mom(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]); 
      momL[0] =momL[0]*Momr + Momm;
    
      
      // Target struggling LHRS step #7
      if( z_av_1[0]>0.08)
	{delta_pep[0] = 6.23409e-3*ph2[0] + 4.03363e-1;} 
      else	
	{delta_pep[0] = -1.35758*sin(-4.59571* ph2[0]) + 2.09093;}       
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV
      
      
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
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; // kaon momentum at the reaction point
      
      // missing mass calculation==============================
      par_ep[0] = pep_real[0] ; 
      par_ep[1] = th2[0];
      par_ep[2] = ph2[0];
      
      par_k[0] = pk_real[0];     
      par_k[1] = th1[0]; 
      par_k[2] = ph1[0];
      
      hallap = hallap-0.1843 ;// must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
     
      mm = CalcMM(hallap, par_ep, par_k, mp);     
      mm = (mm)*1000.; // MeV--->GeV
      
      h->Fill(mm);
     
      
      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;
      
      
      R_XFP  = R_XFP*XFPr +XFPm ; 
      R_XpFP = R_XpFP*XpFPr+XpFPm;
      R_YFP  = R_YFP*YFPr+ YFPm;
      R_YpFP = R_YpFP*YpFPr +YpFPm; 
      
      z_av[0] =z_av[0]*Ztr + Ztm;
      cout<<"the value of z-av is "<<z_av[0] << " and that of z_av_1 is "<<z_av_1[0]<<endl; 
      tnew->Fill();
      
      bool lambdaflag=false;  
      int peak_with_hit= -1; 
      for(int j=0 ; j<npeak ; j++){
	if(Lambda_cent[j]-Lambda_width[j]<mm
	   &&mm < Lambda_cent[j]+Lambda_width[j]){	 
	  
	  lambdaflag=true;
	  peak_with_hit=j; 
	  h_1 ->Fill(mm);
	  h_1 ->SetLineColor(j+2);
	  
	}
	else lambdaflag=false;  

	if(ntune_event<nmax && lambdaflag==true){
	  foil_flag[ntune_event] = peak_with_hit;
	  
	  p10[ntune_event]  = par_ep[0];
	  p11[ntune_event]  = par_ep[1];
	  p12[ntune_event]  = par_ep[2];
	  p13[ntune_event]  = par_k[0];
	  p14[ntune_event]  = par_k[1];
	  p15[ntune_event]  = par_k[2];
	  p16[ntune_event]  = hallap;
	    
	  x[ntune_event]  = R_XFP; //RHRS
	  y[ntune_event]  = R_YFP;
	  xp[ntune_event] = R_XpFP;
	  yp[ntune_event] = R_YpFP;

	  // x[ntune_event]  = XFP; //// LHRS
	  // y[ntune_event]  = YFP;
	  // xp[ntune_event] = XpFP;
	  // yp[ntune_event] = YpFP;

	  z_recon[ntune_event] = z_av_1[0];
	  phir[ntune_event] =ph1[0];
	  phil[ntune_event] =ph2[0];
	 
	  ntune_event++;
	    
	}
	  
      }//int j	
    }
    
  }
  
  tnew->Write();
  // ((((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
  for (int i=0 ; i< ent_2 ; i++){
    for(int j=0 ; j<max ; j++){
      l_x_fp_2[j]  = -2222.0;
      l_th_fp_2[j] = -2222.0; 
      l_y_fp_2[j]  = -2222.0;
      l_ph_fp_2[j] = -2222.0;
      th1_2[j] = -2222.0;
      th2_2[j] = -2222.0;
      ph1_2[j] =-2222.0;
      ph2_2[j] =-2222.0;    
     
      delta_pep_2[j]= -2222.0;
      pep_real_2[j] =-2222.0;
      delta_pk_2[j]= -2222.0;
      pk_real_2[j] = -2222.0;
     
      r_x_fp_2[j]  = -2222.0;
      r_th_fp_2[j] = -2222.0;
      r_y_fp_2[j]  = -2222.0;
      r_ph_fp_2[j] = -2222.0;
    
      
      
      trig5_2[j] = 0.0;
      rtrig_2 = false;
    }
   
    trig5_2[0] = 0.0;
    rtrig_2 = false;
   
    t2->GetEntry(i);
   
   
    if(trig5_2[0]>1.0) rtrig_2 = true; //JUly 01, 2019
    else rtrig_2 = false;

    z_av_2[0] = (lvz_2[0] + rvz_2[0])/2.0;
    z_av_1_2[0] =  z_av_2[0];
   
    XFP_2   = l_x_fp_2[0];
    XpFP_2  = l_th_fp_2[0];
    YFP_2   = l_y_fp_2[0];
    YpFP_2  = l_ph_fp_2[0];
    
    R_XFP_2   = r_x_fp_2[0];
    R_XpFP_2  = r_th_fp_2[0];
    R_YFP_2   = r_y_fp_2[0];
    R_YpFP_2  = r_ph_fp_2[0];
    
    if(rtrig_2==true &&  fabs(ctime_2)<1.0  && fabs(lvz_2[0]-rvz_2[0])<0.045  && fabs(z_av_2[0])<0.10 ){     
      
      XFP_2  = (XFP_2-XFPm)/XFPr;
      XpFP_2 = (XpFP_2-XpFPm)/XpFPr;
      YFP_2  = (YFP_2-YFPm)/YFPr;
      YpFP_2 = (YpFP_2-YpFPm)/YpFPr;
      
      R_XFP_2  = (R_XFP_2-XFPm)/XFPr; 
      R_XpFP_2 = (R_XpFP_2-XpFPm)/XpFPr;
      R_YFP_2  = (R_YFP_2-YFPm)/YFPr;
      R_YpFP_2 = (R_YpFP_2-YpFPm)/YpFPr;

      z_av_2[0] =(z_av_2[0]- Ztm)/Ztr;
      
      th2_2[0] = calcf2t_th(Theta_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_1_2[0]);
      th2_2[0] = th2_2[0]*Xptr + Xptm; 
    
      ph2_2[0] = calcf2t_ph(PHI_L, XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_1_2[0] );
      ph2_2[0] = ph2_2[0]*Yptr + Yptm;   

     
     
      momL_2[0] =  calcf2t_mom(mom_L, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2[0]); 
      momL_2[0] =momL_2[0]*Momr + Momm; 
    
      momL_2[0] = momL_2[0]*2.218/2.10; //for H/T and tritium only 
     
      // Target struggling LHRS step #7
      if( z_av_1_2[0]>0.08)
	{delta_pep_2[0] = 6.23409e-3*ph2_2[0] + 4.03363e-1;} 
      else	
	{delta_pep_2[0] = -1.35758*sin(-4.59571* ph2_2[0]) + 2.09093;} 
      
      pep_real_2[0] = momL_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV
      
      
      // RHRS angle and momentum calculation      
      th1_2[0] = calcf2t_th(Theta_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,   z_av_2[0]);
      th1_2[0] = th1_2[0]*Xptr + Xptm;
    
      ph1_2[0] = calcf2t_ph(PHI_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2, z_av_2[0]);
      ph1_2[0] = ph1_2[0]*Yptr + Yptm;
       
     
      
      
      momR_2[0] =  calcf2t_mom(mom_R, R_XFP_2, R_XpFP_2, R_YFP_2, R_YpFP_2,  z_av_2[0]);
      momR_2[0] = momR_2[0]*Momr+Momm;
      
      // target struggling step #11	
      if(z_av_1_2[0]>0.08) 
	{delta_pk_2[0] = 3.158e-2*ph1_2[0] + 4.05819e-1;} 
      else 
	{delta_pk_2[0] =-1.31749*sin(-4.61513* ph1_2[0]) + 2.03687;}
      pk_real_2[0] = momR_2[0] + delta_pk_2[0]/1000.0; // kaon momentum at the reaction point
      
      // missing mass calculation==============================
      par_ep_2[0] = pep_real_2[0] ; 
      par_ep_2[1] = th2_2[0];
      par_ep_2[2] = ph2_2[0];
      
      par_k_2[0] = pk_real_2[0];      
      par_k_2[1] = th1_2[0]; 
      par_k_2[2] = ph1_2[0];
      
      hallap_2 = hallap_2-0.1843 ;// must be -ve
      hallap_2 = hallap_2/1000.0; // MeV-->GeV
     
      mm_2 = CalcMM(hallap_2, par_ep_2, par_k_2, mp);     
      mm_2 = (mm_2)*1000.; // MeV--->GeV
      
      h_2->Fill(mm_2);
     
      
      XFP_2 = XFP_2 * XFPr + XFPm;
      XpFP_2 = XpFP_2 * XpFPr + XpFPm;
      YFP_2 = YFP_2 * YFPr + YFPm;
      YpFP_2 = YpFP_2 * YpFPr + YpFPm;
      
      
      R_XFP_2  = R_XFP_2*XFPr +XFPm ; 
      R_XpFP_2 = R_XpFP_2*XpFPr+XpFPm;
      R_YFP_2  = R_YFP_2*YFPr+ YFPm;
      R_YpFP_2 = R_YpFP_2*YpFPr +YpFPm; 
      
      z_av_2[0] =z_av_2[0]*Ztr + Ztm;
      
      tnew->Fill();
      
      bool lambdaflag_2=false;  
      int peak_with_hit_2= -1; 
      for(int j=0 ; j<npeak_2 ; j++){
	if(Lambda_cent_2[j]-Lambda_width_2[j]<mm_2
	   &&mm_2 < Lambda_cent_2[j]+Lambda_width_2[j]){	 
	  
	  lambdaflag_2=true;
	  peak_with_hit_2=j; 
	 
	  
	}
	else lambdaflag_2=false;  

	if(ntune_event_2<nmax_2 && lambdaflag_2==true){
	  foil_flag_2[ntune_event_2] = peak_with_hit_2;
	  
	  p10_2[ntune_event_2]  = par_ep_2[0];
	  p11_2[ntune_event_2]  = par_ep_2[1];
	  p12_2[ntune_event_2]  = par_ep_2[2];
	  p13_2[ntune_event_2]  = par_k_2[0];
	  p14_2[ntune_event_2]  = par_k_2[1];
	  p15_2[ntune_event_2]  = par_k_2[2];
	  p16_2[ntune_event_2]  = hallap_2;
	    
	  x_2[ntune_event_2]  = R_XFP_2;// RHRS
	  y_2[ntune_event_2]  = R_YFP_2;
	  xp_2[ntune_event_2] = R_XpFP_2;
	  yp_2[ntune_event_2] = R_YpFP_2;

	  // x_2[ntune_event_2]  = XFP_2; //LHRS
	  // y_2[ntune_event_2]  = YFP_2;
	  // xp_2[ntune_event_2] = XpFP_2;
	  // yp_2[ntune_event_2] = YpFP_2;

	  z_recon_2[ntune_event_2] = z_av_1_2[0];
	  phir_2[ntune_event_2] =ph1_2[0];
	  phil_2[ntune_event_2] =ph2_2[0];
	 
	  ntune_event_2++;
	    
	}
	  
      }//int j	
    }
    
  }
  
  tnew->Write();

 // ))))))))))))))))))))))))))))))))))))))))) t2 ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
  
  // fnew->Close(); 
  // gStyle->SetOptFit(111);
  // =================================== 
  // ======== Draw histograms ========== 
  // ===================================
  
  TF1 *f1 = new TF1("f1","gaus",1111.4,1119.97);
  TF1 *f2 = new TF1("f2","gaus",1188.38,1196.39);
  
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  h->Fit("f1","","",1111.4,1119.97);
  h->Fit("f2","","",1188.38,1196.39);
  f1->Draw("same"); 
    
  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(1130,60,Form("#Lambda"));
 
  l.DrawLatex(1130,65,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  l.DrawLatex(1130,70,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  l.DrawLatex(1200,20,Form("#Sigma^{0}"));
  l.DrawLatex(1200,25,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  l.DrawLatex(1200,30,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));


 
  
  // For H data with T kinematics
 
  TF1 *f1_2 = new TF1("f1_2","gaus",1110.99,1120.62);
  TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  c2_2->cd();
  h_2->Draw();
  
  f1_2->SetLineWidth(1);
  h_2->Fit("f1_2","","",1110.99,1120.62);
 
  TLatex l2;
  l2.SetTextSize(0.025);
  l2.DrawLatex(1130,32,Form("#Lambda"));
  
  l2.DrawLatex(1130,36,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  l2.DrawLatex(1130,40,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));
  


  // TLatex l;
  // l.SetTextSize(0.025);
  // l.DrawLatex(1125,80,Form("#Lambda"));
 
  // l.DrawLatex(1125,85,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  // l.DrawLatex(1125,90,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  // l.DrawLatex(1200,32,Form("#Sigma^{0}"));
  // l.DrawLatex(1200,37,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  // l.DrawLatex(1200,42,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));


 

  
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const int nite =0;
  double temp[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    temp[i] = tune(momR_opt,i);
    sprintf(tempc, "./MOM_MATRICES/mom_R5_8th_%d.dat",i);
  
    ofstream * ofs = new ofstream(tempc); 
    int nppp = 0;
    const int nn = 5;  // 5 for 5th order 4 for 4th order
    for(int i=0; i<nn+1; i++){
      for(int e=0; e<nn+1; e++){
	for(int d=0; d<nn+1; d++){ 
	  for(int c=0; c<nn+1; c++){
	    for(int b=0; b<nn+1; b++){
	      for(int a=0; a<nn+1; a++){  
		if(a+b+c+d+e==i){
		  *ofs <<momR_opt[nppp] 
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
      }
    }
    ofs->close();
    ofs->clear();
    
    cout << temp[i]<<endl; 
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
  // Ee = Ee + 0.0116;// GeV
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
  int allparam =Mom_Par; 

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
    step[i] = 1.0e-3; /// Original 
    // step[i] = 2.0*0.5;  // for rough matrix, when matrix is far from reality
    
    // LLim[i] = pa[i] -10; // pa[i]*0.8; // KI
    // ULim[i] = pa[i] + 10; //pa[i]*0.8; // KI
    LLim[i] = pa[i] - pa[i]*0.8; // KI
    ULim[i] = pa[i] + pa[i]*0.8; // KI


    minuit -> mnparm(i,pname,start[i],step[i],LLim[i],ULim[i],ierflg);
  }
  // ~~~~ Strategy ~~~~
  //  arglist[0] = 2.0; // was active before
  arglist[0] = 1.0;  // KI
  minuit->mnexcm("SET STR",arglist,1,ierflg);
  
  
  // ~~~~ Migrad + Simplex  ~~~~ one of the way to get optimized parameter
  arglist[0] = 20000;
  arglist[1] = 0.01; // To make more presise
  minuit -> mnexcm("MINImize",arglist,2,ierflg); 
  
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  double er;
  
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat); 
  minuit -> mnprin(0,amin);
  if(amin>0) chi2=amin;
  
  for(int i=0 ; i<allparam ; i++){  
     
    minuit -> GetParameter(i,momR_opt[i],er); 
  }
  
  return chi2; 
}


// #############################################################
void fcn(int &nPar, double* /*grad*/, double &fval, double* param, int /*iflag*/)
// #############################################################
  
{
  double chi2 = 0.0;
  double chi_12 = 0.0;
  double XFP, XpFP;
  double YFP, YpFP;
  const double sigma = 0.0045; 
  double ref_mm = 0.0; 
  double residual = 0.0;

  double par_ep[3];
  double par_k[3];
  double halla_p;
  double momr[100];
  double moml[100];
  double z_av;
  double z_av_sc;
  double rvz;
  double ph1;
  double th1;
  double ph2;
  double th2;
  double  delta_pk[100];
  double delta_pep[100];
  double pk_real[100];
  double pep_real[100];
  double MM;
  double THL;
  double PHL;
  double THR;
  double PHR;
  // (((((((((((((((((((((((((((((((((((((((( t2 (((((((((((((((((((((((((((((((((
 
  double chi_22 = 0.0;
  double XFP_2, XpFP_2;
  double YFP_2, YpFP_2;
  double ref_mm_2 = 0.0; 
  double residual_2 = 0.0;

  double par_ep_2[3];
  double par_k_2[3];
  double halla_p_2;
  double momr_2[100];
  double moml_2[100];
  double z_av_2;
  double z_av_sc_2;
  double ph1_2;
  double th1_2;
  double ph2_2;
  double th2_2;
  double delta_pk_2[100];
  double delta_pep_2[100];
  double pk_real_2[100];
  double pep_real_2[100];
  double MM_2;
  double THL_2;
  double PHL_2;
  double THR_2;
  double PHR_2;


  //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))
  
  for(int i=0 ; i<ntune_event ; i++){ 
    residual = 0.0;
    ref_mm = 0.0; 
    ref_mm  = Lambda_real[foil_flag[i]];    
    ref_mm = ref_mm/1000.0;
    
    XFP   = x[i];
    XpFP  = xp[i];
    YFP   = y[i];
    YpFP  = yp[i];
    z_av = z_recon[i]; 
    ph1 = phir[i];  // open when calibrate the Momentum 
    ph2 = phil[i];    
   
    XFP   =(XFP -XFPm)/XFPr;  
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    z_av_sc = (z_av - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    moml[0] = moml[0]*Momr+Momm; 
    if( z_av>0.08)
      {delta_pep[0] = 6.23409e-3*ph2 + 4.03363e-1;} 
    else	
      {delta_pep[0] = -1.35758*sin(-4.59571* ph2) + 2.09093;} 
    
    pep_real[0] = moml[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point in GeV    
   
    // THL =  calcf2t_th(param, XFP, XpFP, YFP, YpFP,  z_av);
    // THL = THL*Xptr + Xptm;


    // PHL =  calcf2t_ph(param, XFP, XpFP, YFP, YpFP,  z_av);
    // PHL = PHL*Yptr + Yptm;
   

 
    
    //  par_ep[0] = pep_real[0];// Wwhen LHRS momentum  tuned
    par_ep[0] = p10[i];
    par_ep[1] = p11[i];
    par_ep[2] = p12[i];  
    
   
    //  for the RHRS momentum tunning
    momr[0] =  calcf2t_mom(param, XFP, XpFP, YFP, YpFP,  z_av_sc);
    momr[0] = momr[0]*Momr+Momm;    
    if(z_av>0.08) 
      {delta_pk[0] = 3.158e-2*ph1 + 4.05819e-1;} 
    else 
      {delta_pk[0] =-1.31749*sin(-4.61513* ph1) + 2.03687;}
    pk_real[0] = momr[0] + delta_pk[0]/1000.0; 


    // THR =  calcf2t_th(param, XFP, XpFP, YFP, YpFP, z_av_sc);
    // THR = THR*Xptr + Xptm;
    
    // PHR =  calcf2t_th(param, XFP, XpFP, YFP, YpFP, z_av_sc);
    // PHR = PHR*Yptr + Yptm;
   
    
    par_k[0] = pk_real[0];// when RHRS matrix tuned
    //  par_k[0] = p13[i];    
    par_k[1] = p14[i]; 
    par_k[2] = p15[i];     
    
    halla_p = p16[i];
    MM = CalcMM(halla_p, par_ep, par_k, mp);    
    residual = MM-ref_mm;


    //   chi_12 = chi_12 + pow(residual,2.0);
    ///////  if need to use the sigma statistical weigh

    if(foil_flag[i] ==0)
      {chi_12 = chi_12 + pow(residual,2.0);}
    else
      {chi_12 = chi_12 +3*pow(residual,2.0);}
    
  }
    // (((((((((((((((((((((((((((((((((((((((( t2 ((((((((((((((((((((((((((((((((( 
  for(int i=0 ; i<ntune_event_2; i++){ 
    residual_2 = 0.0;
    ref_mm_2 = 0.0; 
    ref_mm_2  = Lambda_real_2[foil_flag_2[i]];    
    ref_mm_2 = ref_mm_2/1000.0;
    
    XFP_2   = x_2[i];
    XpFP_2  = xp_2[i];
    YFP_2   = y_2[i];
    YpFP_2  = yp_2[i];
    z_av_2 = z_recon_2[i]; 
    ph1_2 = phir_2[i];  // open when calibrate the Momentum 
    ph2_2 = phil_2[i];    
   
    XFP_2   =(XFP_2 -XFPm)/XFPr;  
    XpFP_2  =(XpFP_2-XpFPm)/XpFPr;
    YFP_2   =(YFP_2 -YFPm)/YFPr;
    YpFP_2  =(YpFP_2-YpFPm)/YpFPr;
    z_av_sc_2 = (z_av_2 - Ztm)/Ztr;       
  
    // For the LHRS momentum tunning
    moml_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    moml_2[0] = moml_2[0]*Momr+Momm;

    moml_2[0] = moml_2[0]*2.218/2.1; //for H/T and tritium only
   
    if( z_av_2>0.08)
      {delta_pep_2[0] = 6.23409e-3*ph2_2 + 4.03363e-1;} 
    else	
      {delta_pep_2[0] = -1.35758*sin(-4.59571* ph2_2) + 2.09093;}    
    pep_real_2[0] = moml_2[0] + delta_pep_2[0]/1000.0; //LHRS  momentum at the reaction point in GeV 
    
    // THL_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2);
    // THL_2 = THL_2*Xptr + Xptm;


    // PHL_2 =  calcf2t_ph(param,XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_2);
    // PHL_2 = PHL_2*Yptr + Yptm;
   




    
    //  par_ep_2[0] = pep_real_2[0];// Wwhen LHRS momentum  tuned
    par_ep_2[0] = p10_2[i];
    par_ep_2[1] = p11_2[i];
    par_ep_2[2] = p12_2[i];  
    
   
    //  for the RHRS momentum tunning
    momr_2[0] =  calcf2t_mom(param, XFP_2, XpFP_2, YFP_2, YpFP_2,  z_av_sc_2);
    momr_2[0] = momr_2[0]*Momr+Momm;    
    if(z_av_2>0.08) 
      {delta_pk_2[0] = 3.158e-2*ph1_2 + 4.05819e-1;} 
    else 
      {delta_pk_2[0] =-1.31749*sin(-4.61513* ph1_2) + 2.03687;}
    pk_real_2[0] = momr_2[0] + delta_pk_2[0]/1000.0;

    // THR_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_sc_2);
    // THR_2 = THR_2*Xptr + Xptm;
   
    // PHR_2 =  calcf2t_th(param,XFP_2, XpFP_2, YFP_2, YpFP_2, z_av_sc_2);
    // PHR_2 = PHR_2*Yptr + Yptm;
    
    
    par_k_2[0] = pk_real_2[0];// when RHRS matrix tuned
    //par_k_2[0] = p13_2[i];    
    par_k_2[1] = p14_2[i]; 
    par_k_2[2] = p15_2[i];     
    
    halla_p_2 = p16_2[i];
    MM_2 = CalcMM(halla_p_2, par_ep_2, par_k_2, mp);    
    residual_2 = MM_2-ref_mm_2;
  
    chi_22 = chi_22 +11*pow(residual_2,2.0);
    // chi_22 = chi_22 +pow(residual_2,2.0);
    //)))))))))))))))))))))))))))))))))))))))))))))   t2  ))))))))))))))))))))))))))))
    
  }
  chi2 = chi_12 +chi_22;
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}
