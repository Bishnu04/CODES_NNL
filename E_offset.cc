//AUG 26 2019
// copied from par_opt.cc
//purpose to optimize the Beam  energy and the momentum offsets just like the raster correction, 12 parameters are optimized by using this code which are mostly optimized by Kazuki by using geant4. KAzuki optimized energy and momentum loss parameters are optimized
// Aug27, 2019, code looks working
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


const int nmax = 3000; // can go up to 3000
double x[nmax], y[nmax]; 
double xp[nmax], yp[nmax];
double z_recon[nmax];
int foil_flag[nmax];
int ntune_event = 0;

const int npeak = 2;
double Lambda_width[npeak] = {2.5, 2.45}; //6.5,8.8}
double Lambda_cent[npeak] ={1115.68,1192.51};
double Lambda_real[npeak] ={1115.683,1192.642}; // Mev 
// double Lambda_width[npeak] = {2.5}; // to select only sigma peak
// double Lambda_cent[npeak] ={1192.73};
// double Lambda_real[npeak] ={1192.642}; // Mev 


double p10[nmax],p11[nmax],p12[nmax];
double p13[nmax],p14[nmax],p15[nmax];
double p16[nmax],p17[nmax],p18[nmax],p19[nmax];
double ZR[nmax];
double phir[nmax];
double phil[nmax];
double l_mom[nmax];
double r_mom[nmax];
//========================================
double del_1 = 2.55664e-2;
double del_2 = -9.95197e-1;
// double del_3 = 2.55664e-2;
// double del_4 = -9.95197e-1;

const int Angle_Par =126;
double thetaL_opt[nmax];
double phiL_opt[nmax];
double thetaR_opt[nmax];
double phiR_opt[nmax];

const int Mom_Par = 126;
double momL_opt[nmax];
double momR_opt[nmax];
const int nParamL = 12;
double parLambda[nParamL];
//++++++++++++++++++++++++++++++++++++++++++
const double hrs_ang = 13.2 * 3.14159/180.; 
const double me = 0.000511;
const double mk = 0.493677;
const double mp = 0.938272;
const double mL = 1.115683;
extern double CalcMM(double ee, double* par_ep, double* par_k, double mt);

void E_offset(){
  // ========================================
  // ======= Opening a ROOT file ============ 
  // ======================================== 
 
  TChain * t1 = new TChain("T"); 
  // t1->Add("./Rootfiles/July22_Rootfiles/final_H2_July22.root"); // root file updated on July 22, 2019(/final_H2_July2.root)
  t1->Add("./Rootfiles/July22_Rootfiles/H2_HAug06.root"); // for tuning optics use this data updated on Aug 06, 2019
  //  t1->Add("./Rootfiles/July22_Rootfiles/H22_T_Aug05.root");// H data with tritium kinematics
  double ent = t1->GetEntries();
  cout<<"entry in the t1=="<<ent<<endl;
  
  int evshift = 0; 
  
  const int max = 100;
  Double_t trig5[max]; // JUly 01, 2019
  
  double mom2[max];
  double mom1[max];
  double momL[max];
  double momR[max];
  // double DL[max];
  // double DR[max];
  // double pep_L[max];
  // double pk_R[max];
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
  // double Dpe_off =-0.0489648; //-0.00683418;
  // double Dpk_off =-0.0836116;    //-0.0515421;
  

  double delta[max];
  double a1, a2;
  Int_t runnum; 
  double hallap;
 
  char tempc1[500]; 
  int temp1;
  
  // ifstream *ifs = new ifstream("./dat_file/theta_real.dat");
  // for (int foil = 0; foil<nfoil;foil++){
  //   for(int row =0;row<nrow;row++){
  //     *ifs >> theta_real[row][foil];
  //   }
  // }
  
  // //==============================================================
  // char tempc2[500]; 
  // int temp2;
  
  // ifstream *ifs1 = new ifstream("./dat_file/theta.dat");
  // for (int foil = 0; foil<nfoil;foil++){
  //   for(int row =0;row<nrow;row++){
  //     *ifs1 >> theta[row][foil]; 
  //   }
  // }
 
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
  double Ras_x; 
  double Ras_y;
  double cer_asum;
  double z_av[nmax];
 
  t1->SetBranchAddress("HALLA_p", &hallap);  
  t1->SetBranchAddress("DR.T5", &trig5);
  t1->SetBranchAddress("L.tr.p", &mom2);
  t1->SetBranchAddress("R.tr.p", &mom1);
  // t1->SetBranchAddress("L.tr.tg_dp", &DL);
  // t1->SetBranchAddress("R.tr.tg_dp", &DR);
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
  //tnew->Branch("RasterCor", &RasterCor, "RasterCor/D"); 
  
  double XFP, XpFP;
  double YFP, YpFP;
  double R_XFP, R_XpFP; 
  double R_YFP, R_YpFP;

  /// offset parameters  


  for(int j=0 ; j<nParamL ; j++){ //
    parLambda[j] = -2222.0; 
    //  Opt_ParL[j]   = -2222.0;   // Parameters after tuning
  }
  
  char name_Energy_L[500];
  sprintf(name_Energy_L,"./z_matrix/offset_2nd.dat"); //copied form newpar_19.dat after second tune

  ifstream Energy_L(name_Energy_L);
  for(int j=0;j<nParamL;j++){
    Energy_L >> parLambda[j];
    //Opt_ParL[j] =  parLambda[j];
    
  }
  Energy_L.close();
  
  // --- Check input parameters by cout --- //
  // for(int i = 0; i<nParamL; i++){
  // cout << " Input offset parameters: "
  //      << " " << parLambda[i] << endl; 
  // }
  
   
  // ===============or LHRS  theta information input==========   
  ntune_event = 0;
  for(int i=0 ; i<Angle_Par; i++){
    thetaL_opt[i] = -2222.0;
  }
  
  char name_Angle_L[500]; 
  //  sprintf(name_Angle_L,"./matrices/theta_3rd_LHRS_Opt_7.dat");//theta_3rd_LHRS_Opt_7.dat
  //sprintf(name_Angle_L,"./Angle_matrices/THL_3rd_4.dat");
  sprintf(name_Angle_L,"./All_Matrices/LTH_3rd_2.dat");
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
  // sprintf(name_angle_phi,"./Angle_matrices/PHL_4th_4.dat"); 
  sprintf(name_angle_phi,"./All_Matrices/PHL_5th_0.dat");
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
  //  sprintf(name_Mom_lhrs,"./OPT_matrix/momL_9th_4.dat");
  //  sprintf(name_Mom_lhrs,"./OPT_matrix/momL_19th_3.dat");
 sprintf(name_Mom_lhrs,"./All_Matrices/LMOM_10th_4.dat");
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
  //sprintf(name_Angle_R,"./matrices/xpt_RHRS_4.dat");//theta_rhrs_2nd_11.dat
  // sprintf(name_Angle_R,"./Angle_matrices/THR_6th_4.dat");
  sprintf(name_Angle_R,"./All_Matrices/RTH_8th_2.dat");
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
  // sprintf(name_phi_Rhrs,"./matrices/ypt_RHRS_4.dat");//phi_LHRS_3rd_Opt_9.dat  
  // sprintf(name_phi_Rhrs,"./Angle_matrices/PHR_4th_2.dat"); 
   sprintf(name_phi_Rhrs,"./All_Matrices/RPH_12th_0.dat");
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
  // =====RHRS momentum recon==========================6
 
  ntune_event = 0;
  for(int i =0;i<Mom_Par;i++){
    momR_opt[i] = -2222.0;
  }
  char name_Mom_rhrs[500];
  //sprintf(name_Mom_rhrs,"./matrices/mom_RHRS_4.dat"); 
  // sprintf(name_Mom_rhrs,"./OPT_matrix/momR_11th_3.dat");
  //  sprintf(name_Mom_rhrs,"./OPT_matrix/momR_20th_1.dat");
  sprintf(name_Mom_rhrs,"./All_Matrices/MOMR_10th_1.dat"); 

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
 
 
  TH1F *h = new TH1F("h"," ;Missing Mass(MeV/c^{2});Counts ",250,1025,1275); //H/H 245 and H/T 225
  TH1F *h_1 = new TH1F("h_1"," ",120,1000, 1300);
  TH2F * h10 = new TH2F("h10","H data with H kinematics;Missing Mass(MeV/c^{2});L.tr.th",200,1000,1300, 200,1.2, 2.5);
  h->GetXaxis()->CenterTitle(); 
  h->GetYaxis()->CenterTitle(); 
  TH2F *h5 = new TH2F("h5",";Z(LHRS);Z(RHRS)",500,-0.8,0.8,500,-20,20); 
  TH1F *h6 = new TH1F("h6",";Coin_Time(ns);Counts",150,-0.2,0.2); 
  char tempc[500];
  // ======================================================
  
  //===================================================================  
 
  bool ltrig = false;
  bool rtrig = false; 
  for(int i=0 ; i<nmax ; i++){
    x[i]    = -2222.0; 
    y[i]    = -2222.0; 
    xp[i]   = -2222.0;
    yp[i]   = -2222.0;
    z_av[i] = -2222.0;
    ZR[i] = -2222.0;
    phir[i] = -2222.0;
    phil[i] = -2222.0;
    l_mom[i] = -2222.0;
    r_mom[i] = -2222.0;

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
   
    XFP   = l_x_fp[0];
    XpFP  = l_th_fp[0];
    YFP   = l_y_fp[0];
    YpFP  = l_ph_fp[0];
    
    R_XFP   = r_x_fp[0];
    R_XpFP  = r_th_fp[0];
    R_YFP   = r_y_fp[0];
    R_YpFP  = r_ph_fp[0];
      
    if(rtrig==true &&  fabs(ctime)<1  && fabs(lvz[0]-rvz[0])<0.045  && fabs(z_av[0])<0.115 ){ //&& fabs(lvz[0])<0.112){
      // Note when analyze Z_average, the cut must be fabs(zL-ZR)<0.045 and  z vertex cut will be +/- 0.12 decided on July 17, 2019
	
      XFP  = (XFP-XFPm)/XFPr;
      XpFP = (XpFP-XpFPm)/XpFPr;
      YFP  = (YFP-YFPm)/YFPr;
      YpFP = (YpFP-YpFPm)/YpFPr;

      R_XFP  = (R_XFP-XFPm)/XFPr; 
      R_XpFP = (R_XpFP-XpFPm)/XpFPr;
      R_YFP  = (R_YFP-YFPm)/YFPr;
      R_YpFP = (R_YpFP-YpFPm)/YpFPr;
    
      
      th2[0] = calcf2t_th(Theta_L, XFP, XpFP, YFP, YpFP,  z_av[0]);
      th2[0] = th2[0]*Xptr + Xptm; // reconstructed theta LHRS	
      ph2[0] = calcf2t_ph(PHI_L, XFP, XpFP, YFP, YpFP, z_av[0] );
      ph2[0] = ph2[0]*Yptr + Yptm;// reconstructed phi LHRS
      ph2[0] = ph2[0]-hrs_ang;
      ph3[0] = ph3[0] -hrs_ang;

     
     
      momL[0] =  calcf2t_th(mom_L, XFP, XpFP, YFP, YpFP,  z_av[0]); 
      momL[0] =momL[0]*Momr + Momm; // make sure with TG
    
      //   momL[0] = momL[0]*2.21682/2.10;//2.2160/2.10
     
      // Target struggling LHRS step #7
      if( z_av[0]>0.08) 
	{delta_pep[0] =del_1*ph2[0] + 3.62363;} 
      else	
	{delta_pep[0] =del_2*sin(-4.66838* ph2[0]) + 1.59176;}       
      pep_real[0] = momL[0] + delta_pep[0]/1000.0; //LHRS  momentum at the reaction point     
     
     
   
      th1[0] = calcf2t_th(Theta_R, R_XFP, R_XpFP, R_YFP, R_YpFP,   z_av[0]);
      th1[0] = th1[0]*Xptr + Xptm;
      ph1[0] = calcf2t_ph(PHI_R, R_XFP, R_XpFP, R_YFP, R_YpFP, z_av[0]);
      ph1[0] = ph1[0]*Yptr + Yptm;
      ph1[0] = ph1[0] + hrs_ang;  
      ph4[0] = ph4[0] + hrs_ang; 

       
      momR[0] =  calcf2t_th(mom_R, R_XFP, R_XpFP, R_YFP, R_YpFP,  z_av[0]);
      momR[0] = momR[0]*Momr+Momm;
          
      // target struggling step #11	
      if(z_av[0]>0.08)  
	{delta_pk[0] = del_1*ph1[0] + 3.62363;} 
      else 
	{delta_pk[0] =del_2*sin(-4.66838* ph1[0]) + 1.59176;}
      pk_real[0] = momR[0] + delta_pk[0]/1000.0; // kaon momentum at the reaction point     

      // missing mass calculation==============================
      par_ep[0] = pep_real[0] + 0.00299898; // offset included
      par_ep[1] = th2[0]; 
      par_ep[2] = ph2[0];

      par_k[0] = pk_real[0] - 0.00277373;      
      par_k[1] = th1[0]; 
      par_k[2] = ph1[0];
     
      hallap = hallap - 0.148 ;// must be -ve
      hallap = hallap/1000.0; // MeV-->GeV
      hallap = hallap + 0.000316903; // offset included
      mm = CalcMM(hallap, par_ep, par_k, mp); 
      
      mm = (mm)*1000.; // GeV--->MeV
        
      h->Fill(mm);
      h10->Fill(mm,YpFP);

      XFP = XFP * XFPr + XFPm;
      XpFP = XpFP * XpFPr + XpFPm;
      YFP = YFP * YFPr + YFPm;
      YpFP = YpFP * YpFPr + YpFPm;


      R_XFP  = R_XFP*XFPr +XFPm ; 
      R_XpFP = R_XpFP*XpFPr+XpFPm;
      R_YFP  = R_YFP*YFPr+ YFPm;
      R_YpFP = R_YpFP*YpFPr +YpFPm;     
   
      
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
	    
	  // x[ntune_event]  = R_XFP; // RHRS
	  // y[ntune_event]  = R_YFP;
	  // xp[ntune_event] = R_XpFP;
	  // yp[ntune_event] = R_YpFP;

	  x[ntune_event]  = XFP; //LHRS
	  y[ntune_event]  = YFP;
	  xp[ntune_event] = XpFP;
	  yp[ntune_event] = YpFP;

	  z_recon[ntune_event] = z_av[0];
	  ZR[ntune_event] = rvz[0];
	  phir[ntune_event] =ph1[0];
	  phil[ntune_event] =ph2[0];
	  l_mom[ntune_event] =  momL[0];
	  r_mom[ntune_event] =  momR[0];
	  ntune_event++;
	    
	}
	  
      }//int j	
    }
    
  
  }  
  tnew->Write();
   
  // fnew->Close(); 
  // gStyle->SetOptFit(111);
  // =================================== 
  // ======== Draw histograms ========== 
  // ===================================
  
  TF1 *f1 = new TF1("f1","gaus",1111.82,1119.61);
  TF1 *f2 = new TF1("f2","gaus",1189.35,1195.58);
  
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->cd();
  h->Draw();
  // h_1->Draw("same");
  f1->SetLineWidth(1);
  f2->SetLineWidth(1);
  h->Fit("f1","","",1111.82,1119.61);
  h->Fit("f2","","",1189.35,1195.58);
  f1->Draw("same"); 
    
   
   
  
  TLatex l;
  l.SetTextSize(0.025);
  l.DrawLatex(1130,80,Form("#Lambda"));
 
  l.DrawLatex(1130,85,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  l.DrawLatex(1130,90,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  l.DrawLatex(1205,40,Form("#Sigma^{0}"));
  l.DrawLatex(1205,47,Form("#color[2]{#sigma = %.6g}",f2->GetParameter(2)));
  l.DrawLatex(1205,55,Form("#color[2]{mean = %.6g}",f2->GetParameter(1)));


   
  
  
  // For H data with T kinematics
  // TF1 *f1 = new TF1("f1","gaus",1111.20,1117.83);
  // TCanvas* c2 = new TCanvas("c2","c2",600,600);
  // c2->cd();
  // h->Draw();
  
  // f1->SetLineWidth(1);
  // h->Fit("f1","","",1111.20,1117.83);
 
  // TLatex l;
  // l.SetTextSize(0.025);
  // l.DrawLatex(1130,30,Form("#Lambda"));
  
  // l.DrawLatex(1130,33,Form("#color[2]{#sigma = %.6g}",f1->GetParameter(2)));
  // l.DrawLatex(1130,35,Form("#color[2]{mean = %.6g}",f1->GetParameter(1)));
  


   
  
 // //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
//   const int nite =0;
//   double temp[nite]; 
//   double x[nite];
//   if (nite>0) cout << " Tuning started: " << endl;
//   for(int i=0 ; i<nite ; i++){
//     x[i] = i+1;
//     //  temp[i] = tune(phiL_opt,i); 
//     //  temp[i] = tune(phiL_opt,i); 
//     temp[i] = tune(phiL_opt,i); // for RHRS tune 
//     sprintf(tempc, "./All_Matrices/PHL_6th_%d.dat",i);
//     //  sprintf(tempc, "./Angle_matrices/PHL_5th_%d.dat",i);
//     ofstream * ofs = new ofstream(tempc); 
//     int nppp = 0;
//     const int nn = 4; 
//     for(int i=0; i<nn+1; i++){
//       for(int e=0; e<nn+1; e++){
// 	for(int d=0; d<nn+1; d++){ 
// 	  for(int c=0; c<nn+1; c++){
// 	    for(int b=0; b<nn+1; b++){
// 	      for(int a=0; a<nn+1; a++){  
// 		if(a+b+c+d+e==i){
// 		  *ofs <<phiL_opt[nppp] 
// 		       << " " << a 
// 		       << " " << b
// 		       << " " << c
// 		       << " " << d
// 		       << " " << e << endl;
// 		  nppp++; 
		  
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     ofs->close();
//     ofs->clear();
    
//     cout << temp[i]<<endl; 
//   }    
  
//   if(nite>0){
//     TGraph * gr = new TGraph(nite,x,temp);  
//     TCanvas * c4 = new TCanvas("c4","",600,600); 
//     gr->Draw("*la"); 
//   }


//=================================================================================================
  const int nite = 0; 
  double chi_sq[nite]; 
  double x[nite];
  if (nite>0) cout << " Tuning started: " << endl;  
  for(int i=0 ; i<nite ; i++){    
    x[i] = i+1;
    
    chi_sq[i] = tune(parLambda,i);
    cout << " Tuning# = " << i << ": chisq = " << chi_sq[i] << endl; 
    cout << endl;    
    
    sprintf(tempc, "./E_off_matrix/Eoff_6th_%d.dat",i); 
    ofstream * ofs = new ofstream(tempc);
    *ofs << parLambda[0] << " " << parLambda[1] << " " << parLambda[2]
	 << parLambda[3] << " " << parLambda[4] << " " << parLambda[5]
	 << parLambda[6] << " " << parLambda[7] << " " << parLambda[8] 
	 << parLambda[9] << " " << parLambda[10] << " " << parLambda[11] << endl;


    //  cout << parLambda[0] << " " << parLambda[1] << " " << parLambda[2] << endl;
    // *ofs << Opt_ParL[0] << " " << Opt_ParL[1] << " " << Opt_ParL[2] << endl;
    // cout << Opt_Par[0] << " " << Opt_Par[1] << " " << Opt_Par[2] << endl;
    
    ofs->close();
    ofs->clear();
  }
  
  TGraph * gr = new TGraph(nite,x,chi_sq);
  gr->SetName("gr");
  if(nite>0){
    TCanvas * c2 = new TCanvas("c2","c2",600,600);
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
  // double DE_off = -0.0128793;// -0.0554373;//-0.0543673
  double pe = ee; 
  double Ee = sqrt(me*me + pe*pe); 
  Ee = Ee +0.002835;//+0.00736 original value in GEV
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
  //double Ek = sqrt(vec_k * vec_k);
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
  int allparam = nParamL; 
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
    //  minuit -> GetParameter(i,momL_opt[i],e); // For LHRS tune

    minuit -> GetParameter(i,parLambda[i],e); // For RHRS tune
    // minuit -> GetParameter(i,thetaR_opt[i],e); // For RHRS tune
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
 
  double ref_mm = 0.0; 
  double residual = 0.0;

  double par_ep[3];
  double par_k[3];
  double halla_p;
  double momr[100];
  double moml[100];
  double z_av;
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
    rvz = ZR[i];
    moml[0] = l_mom[i];
    momr[0] = r_mom[i];
    ph1 = phir[i];  // open when calibrate the Momentum 
    ph2 = phil[i];
    
    
   
    XFP   =(XFP -XFPm)/XFPr;  
    XpFP  =(XpFP-XpFPm)/XpFPr;
    YFP   =(YFP -YFPm)/YFPr;
    YpFP  =(YpFP-YpFPm)/YpFPr;
    // for angle optimization
    // th2 = calcf2t_th(param, XFP, XpFP, YFP,YpFP, z_av);
    // th2 = th2*Xptr + Xptm;
    // ph2 = calcf2t_ph(param, XFP, XpFP, YFP, YpFP, z_av);
    // ph2 = ph2*Yptr + Yptm;
    // ph2 = ph2 - hrs_ang; 
    
    // //
    // For the LHRS momentum tunning
    // moml[0] =  calcf2t_th(param, XFP, XpFP, YFP, YpFP,  z_av);
    // moml[0] = moml[0]*Momr+Momm; 
    
    if(z_av>0.08)
      {delta_pep[0] = param[0]*ph2 + param[1];} 
    else 
      {delta_pep[0] = param[2]*sin(-4.66838* ph2) + param[3];}
    pep_real[0] = moml[0] + delta_pep[0]/1000.0;
    
    
    par_ep[0] = pep_real[0] + param[4];// Wwhen LHRS momentum  tuned
    //  par_ep[0] = p10[i] ;
    par_ep[1] = p11[i];
    par_ep[2] = p12[i];
    
    // Note any parameter need to be optimized should use the param[0] form just like I did in this portion of the code
    
    // th1 = calcf2t_th(param, XFP, XpFP, YFP,YpFP, z_av);
    // th1 = th1*Xptr + Xptm;
    // ph1 = calcf2t_ph(param, XFP, XpFP, YFP, YpFP, z_av);
    // ph1 = ph1*Yptr + Yptm;
    // ph1 = ph1 + hrs_ang;  

    //
      //  for the RHRS momentum tunning
    // momr[0] =  calcf2t_th(param, XFP, XpFP, YFP, YpFP,  z_av);
    // momr[0] = momr[0]*Momr+Momm; 
    
    if(z_av>0.08)
      {delta_pk[0] = param[5]*ph1 + param[6];} 
    else 
      {delta_pk[0] = param[7]*sin(-4.66838* ph1) + param[8];}
    pk_real[0] = momr[0] + delta_pk[0]/1000.0;     
   
    
    par_k[0] = pk_real[0] + param[9]; // when RHRS matrix tuned
    //  par_k[0] = p13[i] ;    
    par_k[1] = p14[i]; // theta tune keep close
    par_k[2] = p15[i];
   
    // NOte: All of the momentums are in GeV, so that the calculated offsets must be in GeV
    
    
    halla_p = p16[i] +param[10]+  param[11]*(z_av + 0.115);
    //  cout<< "The value are z_av is  = :"<< z_av<<endl;
    MM = CalcMM(halla_p, par_ep, par_k, mp);    
    residual = MM-ref_mm;
    chi2 = chi2 + pow(residual,2.0);
    // if(foil_flag[i] ==0)
    //   {chi2 = chi2 + pow(residual,2.0);}
    // else
    //   {chi2 = chi2 +4*pow(residual,2.0);}
    
    
    
  }
  chi2 = sqrt(chi2)/(double)ntune_event/sigma;
  fval = chi2;
}





