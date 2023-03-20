#include "Libraries.h"
#include <bits/stdc++.h>
#define PI 3.1415926
#define B_avg 2.0   // in Tesla

using namespace std;

int main()
{  
   gROOT->SetStyle("Plain");
                                         
   TLatex Tl;                         
   Tl.SetTextSize(0.03);              
   
   TH1D* h1 = new TH1D("h1", "Triplet Fit; p(fit) / p(gen); Entries;", 50, 0.0, 2.0);
   TH1D* h2 = new TH1D("h2", "Triplet Fit; sigma_{p}(fit); Entries;", 50, 0.0, 20.0);
   TH1D* h3 = new TH1D("h3", "Triplet Fit; chi2; Entries;", 50, 0.0, 10.0);
   TH1D* h4 = new TH1D("h4", "Triplet Fit; (p(fit) - p(gen)) / sigma_{p}; Entries;", 20, -5.0, 5.0);
   
   double P_gen = 2000.0;  // in MeV
   double Rc;
   double d01;
   double d12;
   double d02;
   double PHI1C;  
   double PHI2C;
   double Z01;
   double Z12;
   double R3D1C;
   double R3D2C;
   double Theta1C;
   double Theta2C;
   double x0,x1,x2;
   double y0,y1,y2;
   double z0,z1,z2;
   double alpha1, alpha2;   // index parameters
   double R3D_min;
   double r, rt;
   double eta, beta;
   double chi2;
   double P3D;
   double sigma_R3D;
   double sigma_P3D;
   double Theta_MS;
   double PHI_MS;
   
   //vector<double> P_t;
   
   int num_events = 10000;
   
   fstream myFile;
   myFile.open("triplet_sim_out_true.txt", ios::in);
   
   if(myFile.is_open())
   {  
      for (int x = 0; x < num_events; x++)
       {
	    myFile >> x0 >> y0 >> z0;
            myFile >> x1 >> y1 >> z1;
            myFile >> x2 >> y2 >> z2;
            
            Z01 = z1 - z0;
            Z12 = z2 - z1;
               
      	    Rc = (d01 * d12 * d02) / (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1)));	
            
            PHI1C = 2.0 * asin(d01 / (2.0 * Rc));
            PHI2C = 2.0 * asin(d12 / (2.0 * Rc));
           
            R3D1C = sqrt(Rc*Rc + (Z01*Z01)/(PHI1C*PHI1C));
            R3D2C = sqrt(Rc*Rc + (Z12*Z12)/(PHI2C*PHI2C));
            
            Theta1C = acos(Z01/(PHI1C * R3D1C));
            Theta2C = acos(Z12/(PHI2C * R3D2C)); 
            
            alpha1 = (Rc*Rc*PHI1C*PHI1C  + Z01*Z01) / (0.5*Rc*Rc*pow(PHI1C,3)*cos(PHI1C/2.0)/sin(PHI1C/2.0) + Z01*Z01); 
              	
            alpha2 = (Rc*Rc*PHI2C*PHI2C  + Z12*Z12) / (0.5*Rc*Rc*pow(PHI2C,3)*cos(PHI2C/2.0)/sin(PHI2C/2.0) + Z12*Z12);
            
            PHI_tilde = -0.5* (PHI1C * alpha1 + PHI2C * alpha2);
            Theta_tilde = Theta2C - Theta1C - ((1.0 - alpha2) * cos(Theta2C)/sin(Theta2C) - (1.0 - alpha1) * cos(Theta1C)/sin(Theta1C));
            
            eta = (PHI1C * alpha1/(2.0*R3D1C)) +  (PHI2C * alpha2/(2.0*R3D2C));
            
            beta = ((1.0 - alpha2) * cos(Theta2C)/sin(Theta2C)/R3D2C) - ((1.0 - alpha1) * cos(Theta1C)/sin(Theta1C)/R3D1C);
            
            R3D_min = -1.0 * (eta * PHI_tilde * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);
            
            chi2 = (1.0 / (sigma_MS * sigma_MS)) * pow(( beta * PHI_tilde - eta * Theta_tilde), 2) / (eta * eta + beta * beta / pow(sin((Theta1C + Theta2C) / 2.0), 2));
            
            P3D = 0.3 * B_avg * R3D_min;
            
            sigma_R3D = sigma_MS * sqrt(1.0 / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta));
              
            sigma_P3D = 0.3 * B_avg * sigma_R3D;
            
            PHI_MS = (beta * beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);
            
            Theta_MS = (-1.0 * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) * beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);
            
            
            r = P3D / P_gen;
            
            h1 -> Fill(r);
              		  
        }                     
        
                  
     myFile.close();
   } 
   
   
/*--------------------------------------------------------------------------------------------ENDCAP------------------------------------------------------------------------------------*/   
   
   
   TFile fileout("triplet_fit_output.root","recreate");
   fileout.cd();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   
   c1->cd();
   h1->Draw("hist");
   c1->Modified(); c1->Update();
  
   c1->Write();
   c1->Clear();
   
   fileout.Close();
       
   return 0;
   
  }










