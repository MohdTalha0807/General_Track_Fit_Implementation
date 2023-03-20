#include "Libraries.h"
#include <bits/stdc++.h>
#include "other_lib.h"

#define PI 3.1415926
#define B_avg 2.0   // in Tesla

using namespace std;

int main()
{  
   gROOT->SetStyle("Plain");
                                           
   TLatex Tl;                          
   Tl.SetTextSize(0.03);              
   
   TH1D* h1 = new TH1D("h1", "Triplet Fit P3D = 100 MeV; p(fit) / p(gen); Entries;", 100, 0.8, 1.2);
   TH1D* h2 = new TH1D("h2", "Triplet Fit; sigma_{p}(fit); Entries;", 100, 2.0, 6.0);
   TH1D* h3 = new TH1D("h3", "Triplet Fit; chi2; Entries;", 100, 0.0, 10.0);
   TH1D* h4 = new TH1D("h4", "Triplet Fit; (p(fit) - p(gen)) / sigma_{p}; Entries;", 80, -5.0, 5.0); 
   TH1D* h5 = new TH1D("h5", "Triplet Fit; #Phi_{MS}; Entries;", 100, -0.001, 0.001);
   TH1D* h6 = new TH1D("h6", "Triplet Fit; #Theta_{MS}; Entries;", 100, -0.1, 0.1); 
   
     
   double P_gen = 200.0;  // in MeV    
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
   double Theta_tilde, PHI_tilde;
   double R3D_min;
   double r, rt;
   double eta, beta;
   double chi2;
   double P3D;
   double sigma_MS;
   double sigma_R3D;
   double sigma_P3D;
   double Theta_MS;
   double PHI_MS;
   double pull_parameter;
   double phi1, phi2, phi01;
   double sigma_MS1, sigma_MS2;
   double sum_chi2 = 0.0;
   double avg_chi2;
   
   int det_type = 0;    // 0 for barrel layers, 1 for the endcaps
   
   vector<double> PhiMS, ThetaMS;
   
   int num_events = 10000; 
   
   fstream myFile;
   myFile.open("triplet_sim_out_true12.txt", ios::in);
   
   if(myFile.is_open())
   {  
      for (int x = 0; x < num_events; x++)
       {
	    myFile >> x0 >> y0 >> z0;
            myFile >> x1 >> y1 >> z1;
            myFile >> x2 >> y2 >> z2;
            
            d01 = sqrt(pow((x1 -x0),2) + pow((y1 -y0),2));
            d12 = sqrt(pow((x2 -x1),2) + pow((y2 -y1),2));
            d02 = sqrt(pow((x2 -x0),2) + pow((y2 -y0),2));
            
            Z01 = z1 - z0;
            Z12 = z2 - z1;
               
      	    Rc = (d01 * d12 * d02) / (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1)));	
            
            //cout<<Rc<<endl; 
             
            PHI1C = 2.0 * asin(d01 / (2.0 * Rc));
            PHI2C = 2.0 * asin(d12 / (2.0 * Rc));
           
            R3D1C = sqrt(Rc*Rc + (Z01*Z01)/(PHI1C*PHI1C));
            R3D2C = sqrt(Rc*Rc + (Z12*Z12)/(PHI2C*PHI2C));
            
            Theta1C = acos(Z01/(PHI1C * R3D1C));
            Theta2C = acos(Z12/(PHI2C * R3D2C)); 
            
            alpha1 = (Rc*Rc*PHI1C*PHI1C  + Z01*Z01) / (0.5*Rc*Rc*pow(PHI1C,3)*(cos(PHI1C/2.0)/sin(PHI1C/2.0)) + Z01*Z01); 
              	
            alpha2 = (Rc*Rc*PHI2C*PHI2C  + Z12*Z12) / (0.5*Rc*Rc*pow(PHI2C,3)*(cos(PHI2C/2.0)/sin(PHI2C/2.0)) + Z12*Z12);
            
            PHI_tilde = -0.5* (PHI1C * alpha1 + PHI2C * alpha2);
            
            Theta_tilde = Theta2C - Theta1C - ((1.0 - alpha2) * (cos(Theta2C)/sin(Theta2C)) - (1.0 - alpha1) * (cos(Theta1C)/sin(Theta1C)));
            
            eta = (PHI1C * alpha1/(2.0*R3D1C)) +  (PHI2C * alpha2/(2.0*R3D2C));
            
            beta = ((1.0 - alpha2) * (cos(Theta2C)/sin(Theta2C))/R3D2C) - ((1.0 - alpha1) * (cos(Theta1C)/sin(Theta1C))/R3D1C);
            
            R3D_min = -1.0 * (eta * PHI_tilde * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);
            
            P3D = 0.3 * B_avg * R3D_min;
            //cout<<P3D<<endl;
            
            phi01 = atan2((y1-y0),(x1-x0));    // verify this
            phi1 = 0.5 * PHI1C + phi01;
            phi2 = PHI1C + phi1;
           
            sigma_MS1 = calculate_MS_uncertainty( P3D, Theta1C, phi1, x1, y1, z1, det_type);
            
            sigma_MS2 = calculate_MS_uncertainty( P3D, Theta2C, phi2, x2, y2, z2, det_type);
            
            sigma_MS = (sigma_MS1 + sigma_MS2) / 2.0;
            
            chi2 = (1.0 / (sigma_MS * sigma_MS)) * pow(( beta * PHI_tilde - eta * Theta_tilde), 2) / (eta * eta + beta * beta / pow(sin((Theta1C + Theta2C) / 2.0), 2));
            
            sum_chi2 += chi2;
            
            sigma_R3D = sigma_MS * sqrt(1.0 / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta));
              
            sigma_P3D = 0.3 * B_avg * sigma_R3D;
            
            PHI_MS = (beta) * (beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta); // there was a mistake here
            
            Theta_MS = (-1.0 * eta * pow(sin((Theta1C + Theta2C) / 2.0),2)) * (beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);  // there was a mistake here as well 
             
            r = P3D / P_gen;
            
            //cout<<r<<endl;
            //cout<<r<<"\t"<<P3D<<"\t"<<sigma_P3D<<endl;
            
            PhiMS.push_back(PHI_MS);
            ThetaMS.push_back(Theta_MS);
            
            h5->Fill(PHI_MS);
            h6->Fill(Theta_MS);
            
            //cout<<Theta_MS<<"\t"<<PHI_MS<<endl;
            
            h1 -> Fill(r);
            h2 -> Fill(sigma_P3D);
            h3 -> Fill(chi2);
            pull_parameter = (P3D - P_gen) / sigma_P3D;
            h4 -> Fill(pull_parameter);              		  
        
        }                     
        
            //PhiMS.clear();
            //ThetaMS.clear();
                  
     myFile.close();
   } 
   
     avg_chi2 = sum_chi2 / num_events;
     
     cout<<avg_chi2<<endl;
     //cout<<sigma_P3D<<endl;
/*--------------------------------------------------------------------------------------------ENDCAP------------------------------------------------------------------------------------*/   
   
   
   TFile fileout("triplet_fit_output.root","recreate");
   fileout.cd();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   c1->Divide(2,2);
   
   TCanvas *c2 = new TCanvas();
   c2->SetCanvasSize(1800, 1200);
   c2->SetWindowSize(1800, 1000);
   //c1->Divide(2,2);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   c3->Divide(1,2);
   
   c1->cd(1);
   h1->Draw("hist");
   c1->Modified(); c1->Update();
   
   c1->cd(2);
   h2->Draw("hist");
   c1->Modified(); c1->Update();
   
   c1->cd(3);
   h3->Draw("hist");
   c1->Modified(); c1->Update();
   
   c1->cd(4);
   h4->Draw("hist");
   c1->Modified(); c1->Update();
    
   c2->cd();   
   TGraph *gr1 = new TGraph(num_events, &ThetaMS[0], &PhiMS[0]);
   gr1->SetTitle("Fit result for particles simulated with polar angle #frac{#pi}{4} and P_{3D} = 100 MeV");
   gr1->GetYaxis()->SetTitle("#Phi_{MS}");
   gr1->GetXaxis()->SetTitle("#Theta_{MS}");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr1->SetLineColor(kRed);
   gr1->Draw("AP*");
  
   c3->cd(1);
   h5->Draw("hist");
   c3->Modified(); c3->Update();
   
   c3->cd(2);
   h6->Draw("hist");
   c3->Modified(); c3->Update();
   
   c1->Write();
   c1->Clear();
   
   c2->Write();
   c2->Clear();
   
   c3->Write();
   c3->Clear();
  
   fileout.Close();
       
   return 0;
   
  }










