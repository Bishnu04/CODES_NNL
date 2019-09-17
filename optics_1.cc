// This code is just checking the histogram after the COMBINED.cc optimize the matrix. Thsi code can also be used for the matrix tune( but one kinematics one time)
// SEPT 17, 2019.. code is working

extern double calcf2t_th(double* P, 
			 double xf, double xpf,
			 double yf, double ypf,double);
extern double calcf2t_ph(double* P, 
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

const int Mom_Par = 126;
double momL_opt[nmax];
double momR_opt[nmax];
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
//const double hrs_ang1 = 13.1 * 3.14159/180.;
double ztR_wRC_[nmax];

const double me = 0.000511;
const double mk = 0.493677;
//////const double mp = 2.808944; // for tritium by DR. Tang
const double mp = 0.938272;
//const double mp = 2.808921; // for tritium by Gogami
const double mL = 1.115683;
extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);

void optics(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
  
  TChain * t1 = new TChain("T"); 
  // t1->Add("./Rootfiles/July22_Rootfiles/Aug19_all_T.root"); // tritium data 
  // //// t1->Add("./Rootfiles/July22_Rootfiles/final_H2_July22.root"); // a1<70 && a2>1800 && a2<600 
  // //// t1->Add("./Rootfiles/July22_Rootfiles/H22_T_Aug11.root"); // Aug 11,H data with T kinematics
  // t1->Add("./Rootfiles/July22_Rootfiles/H22_T_Aug05.root");  // a1<70 && a2>1800 && a2<6000 H/T kinematics
  t1->Add("./Rootfiles/July22_Rootfiles/H2_HAug06.root");  // a1<120 && a2>1600 && a2<7000 H/H kinematics
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
  sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");//theta_3rd_LHRS_Opt_7.dat
  //  sprintf(name_Angle_L,"./Angle_matrices/THL_3rd_4.dat");
  // sprintf(name_Angle_L,"./All_Matrices/LTH_3rd_2.dat");
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
  sprintf(name_angle_phi,"./matrices/phi_LHRS_3rd_Opt_9.dat");  
  //sprintf(name_angle_phi,"./Angle_matrices/PHL_4th_4.dat");
  //  sprintf(name_angle_phi,"./All_Matrices/PHL_5th_0.dat"); 
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
  //  sprintf(name_Mom_lhrs,"./matrices/mom_LHRS_4.dat");
  //sprintf(name_Mom_lhrs,"./OPT_matrix/momL_19th_3.dat");
  sprintf(name_Mom_lhrs,"./All_Matrices/MOML_28th_3.dat");//MOML_2nd_2.dat
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
  //  sprintf(name_Angle_R,"./matrices/xpt_RHRS_4.dat");//theta_rhrs_2nd_11.dat
  // sprintf(name_Angle_R,"./Angle_matrices/thetaR_3rd_3.dat");
  // sprintf(name_Angle_R,"./Angle_matrices/THR_6th_4.dat");
  //  sprintf(name_Angle_R,"./All_Matrices/RTH_8th_2.dat");
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
  //  sprintf(name_phi_Rhrs,"./matrices/ypt_RHRS_4.dat");//phi_LHRS_3rd_Opt_9.dat  initial matrix
  //  sprintf(name_phi_Rhrs,"./Angle_matrices/phiR_2nd_2.dat");
  //  // sprintf(name_phi_Rhrs,"./Angle_matrices/PHR_4th_2.dat");
  //  sprintf(name_phi_Rhrs,"./All_Matrices/RPH_12th_0.dat");
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
  // sprintf(name_Mom_rhrs,"./matrices/mom_RHRS_4.dat");
  //  sprintf(name_Mom_rhrs,"./OPT_matrix/momR_20th_1.dat");
  //sprintf(name_Mom_rhrs,"./All_Matrices/MOMR_17th_1.dat");
  sprintf(name_Mom_rhrs,"./All_Matrices/MOMR_27th_4.dat");
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
 
  // TH1F *h7 = new TH1F("h7"," ;;Counts ",200,-0.25,0.25); 
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts ",250, 1025,1275); //H/H 250 and H/T 150  and 130 fot T data Aug 19, 2019
  // TH1F *h = new TH1F("h"," ; -B_{#Lambda}(MeV);Counts ",130,-200,200); // for tritium
  TH1F *h2 = new TH1F("h2"," ;",130,-200,200); 
  // TH2F *h2 = new TH2F("h2"," ; ",150,-1.5,1.5,150,1.2,2.4); // 130 bin -40 to 40 for back ground analysis
  h->GetXaxis()->CenterTitle(); 
  h->GetYaxis()->CenterTitle(); 
  gStyle->SetOptStat(0);
  TH2F * h10 = new TH2F("h10","",100,1025,1275, 200,-0.5,1);
 
  TH2F *h5 = new TH2F("h5",";Z(LHRS);Z(RHRS)",500,-0.8,0.8,500,-20,20); 
  TH1F *h6 = new TH1F("h6",";Coin_Time(ns);Counts",150,-0.2,0.2);
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
      // Dpe_off[j] =0.0;
      // Dpk_off[j] = 0.0;
      
      
      trig5[j] = 0.0;
      rtrig = false;
    }
   
    trig5[0] = 0.0;
    rtrig = false;
   
    t1->GetEntry(i); //  may the following 2 lines are open this line will be closed
    //   if(i+evshift<ent) t1->GetEntry(i+evshift); 
    //  else t1->GetEntry(i-evshift); 
   
    if(trig5[0]>1.0) rtrig = true; //JUly 01, 2019
    else rtrig = false;

    z_av[0] = (lvz[0] + rvz[0])/2.0;
    z_av_1[0] =  z_av[0];
    // z_av[0] =(z_av[0]- Ztm)/Ztr;

    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];


      
    // use z_av[0] in stead of individual z
    if(rtrig==true &&  fabs(ctime)<1  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.115 ){
    //  if(rtrig==true && ((ctime>-39.3765 && ctime< -7.28) ||(ctime>15.0499 && ctime< 39.0487))   && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.115 ){
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;

      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;
    
      z_av[0] =(z_av[0]- Ztm)/Ztr;	
      // z_av[0] =z_av[0] *Ztr+ Ztm;

      // LHRS angle and momentum calculation
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av_1[0]);
      th2[0] = th2[0]*Xptr + Xptm; // reconstructed theta LHRS	
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av_1[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;// reconstructed phi LHRS
      ph2[0] = ph2[0]-hrs_ang;
     
     
      momL[0] =  calcf2t_th(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]);
      momL[0] =momL[0]*Momr + Momm; 
    
      //  momL[0] = momL[0]*2.21748/2.10; //for H/T and tritium only 
    

      // Target struggling LHRS step #7 ( momentum loss  when the particle hits the wallod the gas cell)
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
      ph1[0] = ph1[0] + hrs_ang;  
     

   
    
      momR[0] =  calcf2t_th(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;
        
      // target struggling step #11	
      if(z_av_1[0]>0.08) 
	{delta_pk[0] = 3.158e-2*ph1[0] + 4.05819e-1;} 
      else 
	{delta_pk[0] =-1.31749*sin(-4.61513* ph1[0]) + 2.03687;}
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; // kaon momentum at the reaction point
     
      // ==========================================================

      // missing mass calculation==============================
      par_ep[0] = pep_real[0]; // offset included
      par_ep[1] = th2[0]; 
      par_ep[2] = ph2[0];

      par_k[0] = pk_real[0];  // offset included    
      par_k[1] = th1[0];; 
      par_k[2] = ph1[0];
     
      hallap = hallap-0.148;//loss when elctorn entres the target cell must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
      // hallap = hallap + 0.000439299 + 0.00169121*(z_av[0]+0.115);  // offset included
 
      mm = CalcMM(hallap, par_ep, par_k, mp); 
      
      mm = (mm)*1000.;
    
      /////////// mm = mm -2994.761; // for tritium data only By Dr. Tang
      //   mm = mm -2994.814; // for tritium data only By TOSHI

      h->Fill(mm);
   
      //   h2->Fill(ctime,R_XFP);
      // To see the background undr the real Lambda be spectrum on Tritium data
      // if(rtrig==true &&  fabs(ctime)<1  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.115){ //===============mm calculation up to here=============
      // 	h->Fill(mm);
      // }
      
      // if(rtrig==true && ((ctime>-39.3765 && ctime< -7.28) ||(ctime>15.0499 && ctime< 39.0487))   && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.115){    
      // 	h2->Fill(mm); // for back ground analysis
      // }
     


      //  h2->Fill(ctime,R_YpFP);
    


  
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
	  //  h_1 ->Fill(mm); // see  how much is the svent selection from the peak
	  // h_1 ->SetLineColor(j+2);

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
  
  //// TCanvas* c7 = new TCanvas("c7","c7",600,600);
  // h7->Draw();
  
  TF1 *f1 = new TF1("f1","gaus",1109.15,1122.13);
  TF1 *f2 = new TF1("f2","gaus",1186.09,1197.4);
  
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  h->Fit("f1","","",1109.15,1122.13);
  h->Fit("f2","","",1186.09,1197.4);
  f1->Draw("same"); 
    
  

  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(1130,60,Form("#Lambda"));
 
  l.DrawLatex(1130,65,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  l.DrawLatex(1130,70,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  l.DrawLatex(1200,20,Form("#Sigma^{0}"));
  l.DrawLatex(1200,25,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  l.DrawLatex(1200,30,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));


  //  // For H data with T kinematics
  // TF1 *f1_2 = new TF1("f1_2","gaus",1108.28,1123.16);
  // TCanvas* c2_2 = new TCanvas("c2_2","c2_2",600,600);
  // c2_2->cd();
  // h->Draw();
  
  // f1_2->SetLineWidth(1);
  // h->Fit("f1_2","","",1108.28,1123.16);
 
  // TLatex l2;
  // l2.SetTextSize(0.025);
  // l2.DrawLatex(1130,32,Form("#Lambda"));
  
  // l2.DrawLatex(1130,36,Form("#color[2]{#sigma = %.6g}",f1_2->GetParameter(2)));
  // l2.DrawLatex(1130,40,Form("#color[2]{mean = %.6g}",f1_2->GetParameter(1)));
  

 
  
  // h->Fit("f1","","",-0.1362,-0.11509);
  // h->Fit("f2","","",0.1120,0.1334);
  //  f1->Draw("same");
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

  //for background analysis 
  // TCanvas* c2 = new TCanvas("c2","c2",600,600);
  // c2->cd();
  // h->Draw("e");
  // //  h2->Draw("same");
  // h2->SetFillStyle(3002);
  // // h2->SetFillColor(51);
  // // h2->DrawCopy("hist same");
  // h2->Draw("E2 same");
  // h2->Scale(1.0/28.0);
  // h2->SetMarkerStyle(28);
  // h2->SetMarkerColor(2);
  // //  h2->SetFillColor(2);
  





  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const int nite =0;
  double temp[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;
  for(int i=0 ; i<nite ; i++){
    x[i] = i+1;
    temp[i] = tune(thetaL_opt,i);
    cout<<"This is the first location test. Up to here looks  fine"<<endl;
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
// missing mass function definition====================
double CalcMM(double ee, double* par_ep, double* par_k, double mt){
  double pe = ee; 
  double Ee = sqrt(me*me + pe*pe); 
  Ee = Ee -0.0003;//+0.01107254 original value in GEV
  TVector3 vec_e (0.0, 0.0, pe);
  
  double pep  = par_ep[0]; 
  double xpep = par_ep[1];
  double ypep = par_ep[2];
  double px_ep, py_ep, pz_ep;
  pz_ep = pep / sqrt(1.0 + xpep*xpep + ypep*ypep);
  px_ep = xpep * pz_ep;
  py_ep = ypep * pz_ep;
  TVector3 vec_ep (px_ep, py_ep, pz_ep);
  double Eep = sqrt(pep*pep + me*me);  

  double pk  = par_k[0];
  double xpk = par_k[1];
  double ypk = par_k[2];
  double px_k, py_k, pz_k;
  pz_k = pk / sqrt(1.0 + xpk*xpk + ypk*ypk);
  px_k = xpk * pz_k;
  py_k = ypk * pz_k;
  TVector3 vec_k (px_k, py_k, pz_k);
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
  int allparam = Angle_Par; 
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
