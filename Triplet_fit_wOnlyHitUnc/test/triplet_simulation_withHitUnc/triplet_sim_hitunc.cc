#include "Libraries.h"
#include <bits/stdc++.h>
#define PI 3.1415926
#define B_avg 2.0   // in Tesla
   
using namespace std;
 
double calculate_x_eff(vector<double> &p_dir, vector<double> &det_dir)                   // function to calculate the effective path length of the particle inside the material
{ 
   double x_eff;
   x_eff = 0.01/ (p_dir[0]*det_dir[0] + p_dir[1]*det_dir[1] + p_dir[2]*det_dir[2]);
    
   //cout<<x_eff<<"\t"<<p_dir[0]<<"\t"<<p_dir[1]<<"\t"<<p_dir[2]<<endl;
     
   return x_eff;    
} 
 
double calculate_hit_uncertainty(double x_coord)
{
  double dx;
  double pixel_width = 0.1/sqrt(12);    // 100 µm = 0.1 mm
  
  TRandom3  *rand = new TRandom3(0);     
  
      dx = rand->Gaus(0, pixel_width);
      x_coord += dx;
       
  return x_coord;
}

vector<double> calculate_det_dir(double x, double y, double z, int det_type)
{
    vector<double> detector_orientation;
    double x_comp, y_comp;
    detector_orientation.clear();
           
    if(det_type == 0)
    {  
       x_comp = x/sqrt(x*x + y*y);
       y_comp = y/sqrt(x*x + y*y); 
       detector_orientation.push_back(x_comp);
       detector_orientation.push_back(y_comp);
       detector_orientation.push_back(0.0);
       //cout<<"The detector orientation"<<endl;
       //cout<<detector_orientation[0]<<"\t"<<detector_orientation[1]<<"\t"<<detector_orientation[2]<<endl;    
    }
    else
    {
       detector_orientation.push_back(0.0);
       detector_orientation.push_back(0.0);
       detector_orientation.push_back(1.0);  
    }
    
    return detector_orientation;
}

double calculate_MS_uncertainty( double p_3d, double theta, double phi, double x, double y, double z, int det_type)
{     
   double theta_rms;
   double X_eff;
   //double p_3d;
   double x_comp, y_comp, z_comp;
   
   vector<double> p_dir, det_dir;
   
   p_dir.clear(); 
   det_dir.clear();
   
   det_dir = calculate_det_dir(x,y,z,det_type);
   
   //if(det_type == 0)
   //{cout<<det_dir[0]<<"\t"<<det_dir[1]<<"\t"<<det_dir[2]<<endl;}
   
   x_comp = sin(theta)*cos(phi);
   y_comp = sin(theta)*sin(phi);
   z_comp = cos(theta);
   
   //if(det_type == 0)
   //cout<<x_comp<<"\t"<<y_comp<<"\t"<<z_comp<<"\t"<<theta<<"\t"<<phi<<endl;
   
   p_dir.push_back(x_comp);
   p_dir.push_back(y_comp);
   p_dir.push_back(z_comp);
   
   //p_3d = pt/sin(theta);
   
   //if(det_type == 0)
   //cout<<p_3d<<endl;
   
   X_eff = calculate_x_eff(p_dir, det_dir);
   
   //if(det_type == 0)
   //cout<<X_eff<<endl;
   
   theta_rms = (13.6/p_3d) * sqrt(X_eff)*(1.0 + 0.038*log(X_eff));      // Highland Formula
   
   //if(det_type == 0)
   //cout<<theta_rms<<endl;
   
   return theta_rms;
} 

vector<double> cross_layer(double x, double y, double z, double r, double theta, double phi, double R, int distantpoint)  // including the hit uncertainties
{
  /* x,y,z and phi are updated */
  double r3d;
  double xc,yc,rc,phic;
  double omega;
  //double alpha;
  double xnew,ynew,phinew, znew;
  double xnew1,ynew1,phinew1;
  double xnew2,ynew2,phinew2;
  //double dist;
  double PHI,PHI1,PHI2;
  double dz;
  double arg;
  //double eps=1e-7;
  double alpha1,alpha2;
  //double beta1,beta2;
  double phi0;
  double phi_wc;  
  double du = 0.0;
  double dx, dy;
  
  vector<double> new_coord;
  new_coord.clear();
  
  // Multiple Scattering
  //  printf("\n*x=%f,*y=%f,*phi=%f\n",x,*y,*phi);

  // azimuthal inclination angle at start point

  r3d = r/sin(theta);

  /* center of circle */            
  xc=x-sin(phi)*r;                   // with respect to the original coordinates                                   
  yc=y+cos(phi)*r;
  rc=sqrt(xc*xc+yc*yc);

  phi0=atan2(y-yc, x-xc);            //  azimuthal angle w.r.t. the shifted origin to the center of the particle's trajectory                  
  
  new_coord.push_back(phi0);

  // azimuthal angle to center
  phic=atan2(yc,xc);

  //  printf("xc=%f  yc=%f    (phic=%f)\n",xc,yc,phic);


  // angle between circle center and next hit point 
  arg=0.5*(R*R+rc*rc-r*r)/R/rc;
  
  if (fabs(arg)>1.0) {                                   
    cout<<"Error detected!!!"<<endl;  
  } else {
    omega=acos(arg);
  }

  // two solutions: phic+omega and phic-omega - find out right one
  alpha1=phic+omega;
  alpha2=phic-omega;

  //  printf("xc=%f  yc=%f   phic=%f omega=%f alpha=%f  phi=%f   theta=%f  y=%f  r=%f\n",xc,yc,phic,omega,alpha,*phi,theta,*y,r);

  // azimuthal angle to next hit point 

  // new coordinates
  xnew1=R*cos(alpha1);
  ynew1=R*sin(alpha1);

  xnew2=R*cos(alpha2);
  ynew2=R*sin(alpha2);

  phinew1=atan2(ynew1-yc,xnew1-xc);           //  the azimuthal angles for the two cross points w.r.t. the center of the particle's trajectory
  phinew2=atan2(ynew2-yc,xnew2-xc);

  PHI1=phinew1-phi0;                          //  the azimuthal distances of the particle from the two crossings 
  PHI2=phinew2-phi0;

  if (r3d<0) {                                
    PHI1=-PHI1;
    PHI2=-PHI2;
  }
  
  while (PHI1<0) PHI1+=2*PI;                           
  while (PHI2<0) PHI2+=2*PI;                      
  while (PHI1>=2*PI) PHI1-=2*PI;                   
  while (PHI2>=2*PI) PHI2-=2*PI;                

  if (PHI1<PHI2) {
    if (distantpoint==0) {
      PHI=PHI1;
      xnew=xnew1;
      ynew=ynew1;
      phi_wc = alpha1; // i have defined for my personal use
    }else{
      PHI=PHI2;
      xnew=xnew2;
      ynew=ynew2;
    }
  }else{
    if (distantpoint==0) {
      PHI=PHI2;
      xnew=xnew2;
      ynew=ynew2;
      phi_wc = alpha2;
    }else{
      PHI=PHI1;
      xnew=xnew1;
      ynew=ynew1;
    }
  }
  
  if (r<0) PHI=-PHI;

  
  // printf("cross_layer: PHI=%f  (PHI1=%f  PHI2=%f)\n",PHI,PHI1,PHI2);

  // new azimuthal angle
  phinew = phi+PHI;        
  
  new_coord.push_back(phi_wc);
  new_coord.push_back(phinew);
  
  dz = r * PHI/tan(theta);                  // it is a relation between the arc length and the longitudinal distance z
  
  //update
  //cout<<"without hit uncertainty : "<<xnew<<"\t"<<ynew<<endl;
  
  //cout<<"with hit uncertainty : "<<xnew<<"\t"<<ynew<<endl;
  //cout<<endl;
  //cout<<endl;
  
  new_coord.push_back(xnew);
  new_coord.push_back(ynew);
  znew = z + dz;
   
  new_coord.push_back(znew);
  new_coord.push_back(PHI);
  new_coord.push_back(xc);
  new_coord.push_back(yc);
  
  //*phi=phinew; 

  //  printf("xnew=%f   ynew=%f   znew=%f   (phinew=%f)\n",xnew,ynew,*z,phinew);

  return new_coord;   // It contains phi0, phi_wc, phinew, xnew, ynew, znew respectively
}

/*--------------------------------------------------Main Function---------------------------------------------------------*/

 
int main()
{   
   gROOT->SetStyle("Plain");
                                          
   TLatex Tl;                         
   Tl.SetTextSize(0.03);              
    
   TMultiGraph  *mg1  = new TMultiGraph();
   TMultiGraph  *mg2  = new TMultiGraph();
   TMultiGraph  *mg3  = new TMultiGraph();
   TMultiGraph  *mg4  = new TMultiGraph();
   TH1D* h1 = new TH1D("h1", "Triplet Fit; R - R_{cal}; Entries;", 100, -0.001, 0.001);
   
   int Ev_Id;
   std::map<uint, ROOT::Math::XYZPoint> coorMap;
   double PHI_initial;
   int particle_id = 0;
   
   TFile *fout_data = new TFile("triplet_sim_output.root","recreate");                   // tree to store the x,y,z positions of the hits 
   
   TTree *treeOut = new TTree("T","Triplet_Info");
   
   treeOut -> Branch("Evt_ID",&Ev_Id);
   treeOut -> Branch("Coord0", &coorMap[0]);
   treeOut -> Branch("Coord1", &coorMap[1]);
   treeOut -> Branch("Coord2", &coorMap[2]);
   treeOut -> Branch("In_Phi", &PHI_initial);
   
      
   double q = 1.0;                                     
   double phi, x,y,z, In_phi;                      
   double p3d = 2000.0;     // in MeV   
   double R3d = p3d/(0.3*q*B_avg);   // in mm  
   //cout<<R3D<<endl;
   double R;   // radius of the barrel layer
   double theta;
   double phi0, phi_1_MS, theta1_MS, phi_wc, PHI;
   int distantpoint = 0;
   double theta_RMS, xc, yc;
   vector<double> new_coordinates;                 // number of particles generated 
   vector<double> cx, cy;
   int n = 10;         // this is here to define the number of steps taken to traverse the bending angle
   double arc_s,s;
   vector<double> arc_length, z_pos;
   int ni;
   double z_endcap_ini = 300.0;   // the position of the first endcap layer in mm
   double y_ini, z_ini, phi_ini;
   int det_type1 = 0;
   int det_type2 = 1;
   double r;
   //vector<double> p3D_vec;
   vector<double> phiMS_vec, s_vec, phi_vec;
   vector<double> thetaMS_vec, Cc_vec, RR1_vec, RR2_vec, RR3_vec;
   vector<double> vec_Xc;
   vector<double> vec_Yc;
   vector<double> vec_Zc;
   
   double Cc, d01, d02, d12;
   double x_hitUnc, y_hitUnc, z_hitUnc;
   double dx, dy, du, dv;
   
   int count_dev = 0;
   
       phiMS_vec.clear(); 
       thetaMS_vec.clear();
     
   TRandom3  *rand1 = new TRandom3(0);          // 0 is there for random seeds 
   TRandom3  *rand2 = new TRandom3(0);   
     
   fstream myFile;
   myFile.open("triplet_sim.txt", ios::in);
   
   ofstream fout;
   fout.open ("triplet_sim_out.txt");
    fout<<"x"<<"\t"<<"\t"<<"y"<<"\t"<<"\t"<<"z"<<"\t"<<"\t"<<"phi"<<"\t"<<"\t"<<"theta"<<"\t"<<"\t"<<"r"<<"\t"<<"R"<<"\t"<<"pt"<<endl;       
   
   ofstream fout2;
   fout2.open ("triplet_sim_out2.txt");
   fout2<<"x"<<"\t"<<"\t"<<"y"<<"\t"<<"\t"<<"z"<<"\t"<<"\t"<<"phi"<<"\t"<<"\t"<<"theta"<<"\t"<<"\t"<<"r"<<"\t"<<"\t"<<"pt"<<endl;       
   
   const char *path="/home/mtalha/Ph.D._work/Triplet_fit_sim/Triplet_fit_wOnlyHitUnc/Method_Trees/Fit_hitUncOnly";    //Triplet_MS_Fit";
   ofstream fout_1;
   fout_1.open ("/home/mtalha/Ph.D._work/Triplet_fit_sim/Triplet_fit_wOnlyHitUnc/Method_Trees/Fit_hitUncOnly/triplet_sim_out_true12.txt");
   //fout_1<<"x"<<"\t"<<"\t"<<"y"<<"\t"<<"\t"<<"z"<<endl;       
   
   ofstream fout_2;
   fout_2.open ("triplet_sim_out2_true.txt");
   fout_2<<"x"<<"\t"<<"\t"<<"y"<<"\t"<<"\t"<<"z"<<endl;       
   
   fout_data -> cd();
   
   
   if(myFile.is_open())
   {  int counter = 0;
      int i = 0;
    while( myFile >> x >> y >> z >> phi && counter < 10000)    // 'counter' is for the number of particle simulations one wants to see.
    {  
       // p3d = 100.0 + i * 0.20;     // in MeV                   // varying the initial momentum
       // R3d = p3d/(0.3*B_avg);   // in mm  
       // i++; 
       y_ini = y;                // this is done for endcap
       z_ini = z;
       In_phi = phi; 
       phi_ini = phi;
       
       counter++;
       arc_s = 0.0;
       theta = M_PI/3.0;
        
       r = R3d * sin(theta);    // radius for the transverse plane
       
       //pt = p3d * sin(theta);
    //cout<<"--------------------------------------------------------------------"<<endl;                   
    for(int j = 1; j <= 3; j++)
    {   
     R = j*50.0;            // the radius of the barrel layer
     
     arc_length.push_back(arc_s);   // this information is needed for s vs z plot
     z_pos.push_back(z);
     
     double zin = z;
              
     new_coordinates = cross_layer( x, y, z, r, theta, In_phi, R, distantpoint);          // r and distantpoint remain unchanged
     phi0 = new_coordinates[0];               
     phi_wc = new_coordinates[1];            
          
     In_phi = new_coordinates[2];    
     x = new_coordinates[3];
     y = new_coordinates[4];
     z = new_coordinates[5];
     PHI = new_coordinates[6];
     xc  = new_coordinates[7];
     yc  = new_coordinates[8];
     
     vec_Xc.push_back(x);
     vec_Yc.push_back(y);
     vec_Zc.push_back(z);
     
    double R_cal = sqrt(x*x + y*y);
    if(R == 50.0)
    {RR1_vec.push_back(R-R_cal); }
    if(R == 100.0)
    {RR2_vec.push_back(R-R_cal); }
    if(R == 150.0)
    {RR3_vec.push_back(R-R_cal); }
    
    h1 -> Fill(R-R_cal);
    if(abs(R-R_cal) > 0.0)
     {
         count_dev++;
         //cout<<abs(R-R_cal)<<endl;
                  
     }
    //cout<<R<<""<<R_cal<<endl;
     
     s = r*PHI;                // PHI is the bending angle  
     if(j > 1) 
     { s_vec.push_back(s);
       //phi_vec.push_back(phi_ini);
     }
     
     arc_s+=s;                 
     arc_length.push_back(arc_s);      // again this is here for s vs z plot 
     z_pos.push_back(z);
     ni = z_pos.size();
     
     new_coordinates.clear();
     
     du = 0.0;                                                                                         // hit uncertainty implementation-------------------------------
     dv = 0.0;                                     
     du = calculate_hit_uncertainty(du);
     dv = calculate_hit_uncertainty(dv);      
     dx = du * sin(phi_wc);   // need to check this 
     dy = du * cos(phi_wc);
     x_hitUnc = x + dx;
     y_hitUnc = y + dy;
     z_hitUnc = z + dv;
     
     Ev_Id = particle_id;
     coorMap[j-1] = ROOT::Math::XYZPoint(x_hitUnc,y_hitUnc,z_hitUnc);
     PHI_initial = (phi_ini * 180.0 / PI);
     
     //treeOut -> Fill();
     
     
     fout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<phi<<"\t"<<theta<<"\t"<<r<<"\t"<<R<<"\t"<<p3d<<endl; 
     fout_1<<x_hitUnc<<"\t"<<y_hitUnc<<"\t"<<z_hitUnc<<"\t"<<phi_ini<<endl; 
   
   // visualizing the hits on the barrel layers 
    
     double h = PHI/n;                  // stepsize to traverse through the bending angle
     for(int k = 0; k <= n; k++)                    
     {
         cx.push_back(r*cos(phi0 + k*h) + xc);  
         cy.push_back(r*sin(phi0 + k*h) + yc); 
         //cout<<cx[k]<<"\t"<<cy[k]<<endl; 
     } 
      
      TGraph *gr1 = new TGraph(n+1, &cx[0], &cy[0]);
      gr1->SetLineColor(kBlue);
      mg1->Add(gr1);
      
     /* TGraph *gr10 = new TGraph(1, &cxc[0], &cyc[0]); 
     // gr10->SetLineColor(kBlack);
      gr10->SetMarkerStyle(2);
      gr10->SetMarkerSize(5);
      gr10->SetLineColor(i+2);
      gr10->SetLineWidth(3);
      gr10->SetMarkerColor(4);
      mg1->Add(gr10);
      */
      TGraph *gr2 = new TGraph(ni, &z_pos[0], &arc_length[0]);
      gr2->SetLineColor(kBlue);
      mg2->Add(gr2);
      
      cx.clear();
      cy.clear();
      arc_length.clear();
      z_pos.clear();              
   }
      
      particle_id++;
   
      d01 = sqrt(pow((vec_Xc[1] -vec_Xc[0]),2) + pow((vec_Yc[1] -vec_Yc[0]),2));
      d12 = sqrt(pow((vec_Xc[2] -vec_Xc[1]),2) + pow((vec_Yc[2] -vec_Yc[1]),2));
      d02 = sqrt(pow((vec_Xc[2] -vec_Xc[0]),2) + pow((vec_Yc[2] -vec_Yc[0]),2));     
 
      Cc =  (2.0 * ((vec_Xc[1] - vec_Xc[0])*(vec_Yc[2]- vec_Yc[1]) - (vec_Yc[1]- vec_Yc[0])*(vec_Xc[2]-vec_Xc[1]))) / (d01 * d12 * d02);         
      Cc_vec.push_back(Cc);
      phi_vec.push_back(phi_ini * 180.0 / PI);
            
      vec_Xc.clear();
      vec_Yc.clear();
      vec_Zc.clear();
   
   //fout<<endl;
   //fout<<endl; 
   
  // cout<<endl;
  // cout<<endl;
   //------------------------------------------------- for endcap define the z positions of the first three layers----------------------------------------------------// 
   
   double theta_endcap = M_PI / 9.0;
   
   double s_endcap, s_ec_vis = 0.0; 
   double alpha_endcap;
   double phi0_endcap;
   double phi_wc_endcap = M_PI / 2.0;
   double theta1_MS_endcap, phi_1_MS_endcap, theta_RMS_endcap;
   
   double phi_endcap = phi_ini;
   
   double x_endcap, y_endcap, z_endcap, z_endcap_co;
   vector<double> y_ec, z_ec, s_ec;
   y_ec.push_back(y_ini);
   z_ec.push_back(z_ini);
   s_ec.push_back(0.0);
    
    //cout<<theta_endcap<<endl;
    //cout<<"init theta"<<endl;
       
   for(int l = 0; l < 3; l++)
   {   
       
       phi0_endcap   = -1.0*theta_endcap;
       z_endcap_co = z_endcap_ini + l* 50.0;
       
       if(l==0)
       {z_endcap = (z_endcap_ini - z_ini) + l* 50.0;}
       else
       { z_endcap = 50.0;}
       
       s_endcap = z_endcap * tan(theta_endcap);
       
       s_ec_vis += s_endcap;
       //cout<<s_ec_vis<<endl;
       
       alpha_endcap = s_endcap / r;  // bending angle
       
       phi_endcap += alpha_endcap;
       
       theta_RMS_endcap = calculate_MS_uncertainty( p3d, theta_endcap, phi_endcap, 0, 0, 1, det_type2);
            
       phi_1_MS_endcap = rand1->Gaus(0, theta_RMS_endcap);
       theta1_MS_endcap = rand2->Gaus(0, theta_RMS_endcap);   
       phi_endcap += phi_1_MS_endcap;
       theta_endcap += theta1_MS_endcap; 
       //cout<<theta_endcap<<endl;
    
       x_endcap = r* cos(phi_endcap);
       y_endcap = r* sin(phi_endcap);
       
       y_ec.push_back(y_endcap);
       z_ec.push_back(z_endcap_co);
       s_ec.push_back(s_ec_vis);
       
       fout2<<x_endcap<<"\t"<<y_endcap<<"\t"<<z_endcap<<"\t"<<"\t"<<phi_endcap<<"\t"<<theta_endcap<<"\t"<<r<<"\t"<<p3d<<endl;
       fout_2<<x<<"\t"<<y<<"\t"<<z<<endl; 
    
    }   
    //cout<<endl;
    //cout<<endl;
   
      //cout<<"end of 1 event"<<endl;
      TGraph *gr4 = new TGraph(4, &z_ec[0], &y_ec[0]);
      gr4->SetLineColor(kBlue);
      mg3->Add(gr4);
      
      TGraph *gr6 = new TGraph(4, &z_ec[0], &s_ec[0]);
      gr6->SetLineColor(kBlue);
      mg4->Add(gr6);
 
    fout2<<endl;
    fout2<<endl; 
    y_ec.clear();
    z_ec.clear();
    s_ec.clear();      
   
     treeOut -> Fill();    
    
   } 
     
     myFile.close();
   } 
 // --------------------------------------------------------------------------------------------Event Loop closed ----------------------------------------------------------------------  
   //treeOut -> Fill();
   fout_data -> cd();
   treeOut -> Write();
   //fout_data -> Write();
   fout_data -> Close();
   
   
   
//------------------------------------ visualization of the simulation of the barrel layers ----------------------------------------------------------------------------------------------                
   int nint = 20;
   double hphi = 2.0*M_PI / nint;
   double w = 0;
   vector<double> Rx, Ry;
   for(int jr = 1; jr <= 3; jr++)
   {     R = 50.0 * jr;      
   for(int p = 0; p <= nint; p++)
   {
      Rx.push_back(R*cos(w + p*hphi)); 
      Ry.push_back(R*sin(w + p*hphi)); 
   }
     TGraph *gr3 = new TGraph(nint+1, &Rx[0], &Ry[0]);
     gr3->SetLineColor(kRed);
     mg1->Add(gr3);
     Rx.clear();
     Ry.clear();
   }
   
/*--------- -----------------------------------------------------------------------------------ENDCAP------------------------------------------------------------------------------------*/   
   
   vector<double> y_ec_vis, z_ec_vis;
   double y_vis, z_vis;
   
   for(int jk = 0; jk < 3; jk++)
   {
      z_vis = z_endcap_ini + jk*50.0;
      
      for(int j = 0; j < 2; j++) 
       {if(j==0){
       y_vis = 4000.0;
       y_ec_vis.push_back(y_vis); z_ec_vis.push_back(z_vis);}
       else
       {y_vis = -4000.0;
       y_ec_vis.push_back(y_vis); z_ec_vis.push_back(z_vis);}
       }
     TGraph *gr5 = new TGraph( 2, &z_ec_vis[0], &y_ec_vis[0]);
     gr5->SetLineColor(kRed);
     mg3->Add(gr5);
     
     y_ec_vis.clear();
     z_ec_vis.clear();
       
   }
   
   //cout<<count_dev<<endl;
   
   TFile fileout("triplet_sim_plots.root","recreate");
   fileout.cd();
   
   TCanvas *c = new TCanvas();
   c->SetCanvasSize(1200, 1200);
   c->SetWindowSize(1000, 1000);
   //c->Divide();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   
   TCanvas *c2 = new TCanvas();
   c2->SetCanvasSize(1800, 1200);
   c2->SetWindowSize(1800, 1000);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(1800, 1200);
   c4->SetWindowSize(1800, 1000);
   c4->Divide(2,2);
   
   TCanvas *c5 = new TCanvas();
   c5->SetCanvasSize(1800, 1200);
   c5->SetWindowSize(1800, 1000);
         
   c->cd();
   //mg1->SetTitle("");
   mg1->GetXaxis()->SetTitle("X [mm]");
   mg1->GetYaxis()->SetTitle("Y [mm]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   mg1->Draw("AC");
   
   //TLegend *leg=new TLegend(0.7,0.7,0.9,0.9);
   //leg->AddEntry(mg1, "measured points", "l");                 //Maybe you can try replacing mg with gr
   //leg->Draw();
    
   c1->cd();
   mg2->GetXaxis()->SetTitle("Z [mm]");
   mg2->GetYaxis()->SetTitle("S [mm]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   mg2->Draw("AL*");
   
   c2->cd();
   mg3->GetXaxis()->SetTitle("Z [mm]");
   mg3->GetYaxis()->SetTitle("Y [mm]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   mg3->Draw("AL");
   
   c3->cd();
   mg4->GetXaxis()->SetTitle("Z [mm]");
   mg4->GetYaxis()->SetTitle("S [mm]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   mg4->Draw("AL*");
   
   /*c4->cd();   
   TGraph *gr8 = new TGraph(10000, &thetaMS_vec[0], &phiMS_vec[0]);
   gr8->SetTitle("Particles simulated with polar angle #frac{#pi}{4} and P_{3D} = 100 MeV");
   gr8->GetYaxis()->SetTitle("#Phi_{MS}");
   gr8->GetXaxis()->SetTitle("#Theta_{MS}");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr8->SetLineColor(kBlue);
   gr8->Draw("AP*");
   */
   
   int sz = phi_vec.size();
   
   c4->cd(1);   
   TGraph *g1 = new TGraph(sz, &phi_vec[0], &RR1_vec[0]);
   g1->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 20GeV");
   g1->GetYaxis()->SetTitle("R1-R1_{cal} [mm]");
   g1->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   g1->SetLineColor(kRed);
   g1->Draw("AP*");
   
   c4->cd(2);   
   TGraph *g2 = new TGraph(sz, &phi_vec[0], &RR2_vec[0]);
   g2->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 20GeV");
   g2->GetYaxis()->SetTitle("R2-R2_{cal} [mm]");
   g2->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   g2->SetLineColor(kRed);
   g2->Draw("AP*");
   
   c4->cd(3);   
   TGraph *g3 = new TGraph(sz, &phi_vec[0], &RR3_vec[0]);
   g3->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 20GeV");
   g3->GetYaxis()->SetTitle("R3-R3_{cal} [mm]");
   g3->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   g3->SetLineColor(kRed);
   g3->Draw("AP*");
    
   c4->cd(4);   
   TGraph *g = new TGraph(sz, &phi_vec[0], &Cc_vec[0]);
   g->SetTitle("Particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 20GeV");
   g->GetYaxis()->SetTitle("C2D");
   g->GetXaxis()->SetTitle("#phi (degrees)");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   g->SetLineColor(kBlue);
   g->Draw("AP*");
   
   
   Cc_vec.clear();
   
   c5->cd();
   h1->Draw("hist PLC");
   c5->Modified(); c5->Update();
   
   c->Write();
   c->Clear();
   
   c1->Write();
   c1->Clear();
   
   c2->Write();
   c2->Clear();
   
   c3->Write();
   c3->Clear();
   
   c4->Write();
   c4->Clear();
   
   c5->Write();
   c5->Clear();
   
   fileout.Close();
       
   return 0;
   
  }










