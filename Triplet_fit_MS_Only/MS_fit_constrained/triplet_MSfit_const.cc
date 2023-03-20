
// Triplet MS fit for the Constrained case 
 
#include "Libraries.h"
#include <bits/stdc++.h>
#include "other_lib.h"
     
#define PI 3.1415926 
#define B_avg 2.0 // in Tesla 
  
using namespace std; 
   
int main()
{  
   gROOT->SetStyle("Plain");
                                           
   TLatex Tl;                          
   Tl.SetTextSize(0.03);              
   
   TH1D* h1 = new TH1D("h1", "Triplet MS Fit P3D = 20 GeV (constrained); p(fit) / p(gen); Entries;", 100, 0.8, 1.2);
   TH1D* h2 = new TH1D("h2", "Triplet MS Fit; sigma_{p}(fit) [MeV]; Entries;", 150, 600.0, 1200.0);
   TH1D* h3 = new TH1D("h3", "Triplet MS Fit (constrained); chi2_min_con; Entries;", 100, 0.0, 12.0);
   TH1D* h4 = new TH1D("h4", "Triplet MS Fit; (p(fit) - p(gen)) / sigma_{p}; Entries;", 80, -5.0, 5.0); 
   //TH1D* h5 = new TH1D("h5", "Triplet Fit; #Phi_{MS}; Entries;", 100, -0.001, 0.001);
   //TH1D* h6 = new TH1D("h6", "Triplet Fit; #Theta_{MS}; Entries;", 100, -0.1, 0.1); 
    
      
   double P_gen = 20000.0;  // in MeV    
   double Cc;   // transverse curvature for the circle solution 
   double d01;
   double d12;
   double d02;
   double PHI1C;  
   double PHI2C;
   double Z01;
   double Z12;
   double C3D1C;
   double C3D2C;
   double Theta1C;
   double Theta2C;
   double x0,x1,x2;
   double y0,y1,y2;
   double z0,z1,z2;
   double n1C, n2C;   // index parameters
   double Theta0_tilde, PHI0_tilde;
   double C3D_min_con;
   double r, rt;
   double RhoPHI, RhoTheta;
   double chi2_min_con;
   double P3D;
   double sigma_MS;
   double sigma_C3D_con;
   double sigma_P3D;
   //double Theta_MS;
   //double PHI_MS;
   double pull_parameter;
   double phi1, phi2, phi01;
   double sigma_MS1, sigma_MS2;
   double sum_chi2 = 0.0;
   double avg_chi2;
   double sum_sigma, sum_C3D; 
      
   int det_type = 0;    // 0 for barrel layers, 1 for the endcaps
   
  // vector<double> PhiMS, ThetaMS;
               
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
               
      	    Cc =  (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1))) / (d01 * d12 * d02);	
            
            double xx = (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1)));
            //cout<<Cc<<endl; 
             
            PHI1C = 2.0 * asin((d01 * Cc) / 2.0);
            PHI2C = 2.0 * asin((d12 * Cc) / 2.0);
            
            C3D1C = PHI1C / sqrt(Z01 * Z01 + d01 * d01 * (PHI1C * PHI1C)/(4.0 * pow(sin(PHI1C/2.0),2))); 
            C3D2C = PHI2C / sqrt(Z12 * Z12 + d12 * d12 * (PHI2C * PHI2C)/(4.0 * pow(sin(PHI2C/2.0),2)));
            
            //cout<<C3D1C<<endl;
            
            Theta1C = acos((Z01 * C3D1C)/PHI1C);
            Theta2C = acos((Z12 * C3D2C)/PHI2C); 
            
            n1C = 1.0/((PHI1C/2.0) * (cos(PHI1C/2.0)/sin(PHI1C/2.0)) * pow(sin(Theta1C),2) + pow(cos(Theta1C),2));
            n2C = 1.0/((PHI2C/2.0) * (cos(PHI2C/2.0)/sin(PHI2C/2.0)) * pow(sin(Theta2C),2) + pow(cos(Theta2C),2));      
            
              
            PHI0_tilde = 0.5* (PHI1C * n1C + PHI2C * n2C);
            
            Theta0_tilde = Theta2C - Theta1C + ((1.0 - n2C) * (cos(Theta2C)/sin(Theta2C)) - (1.0 - n1C) * (cos(Theta1C)/sin(Theta1C)));
            
            RhoPHI = -0.5 * ((PHI1C * n1C / C3D1C) + (PHI2C * n2C / C3D2C));
            
            RhoTheta = (1.0 - n1C) * (cos(Theta1C)/sin(Theta1C)) / C3D1C  -  (1.0 - n2C) * (cos(Theta2C)/sin(Theta2C)) / C3D2C; 
            
            C3D_min_con =  -(RhoPHI * PHI0_tilde * pow(sin((Theta1C + Theta2C) / 2.0),2) + RhoTheta * Theta0_tilde) / (RhoPHI * RhoPHI * pow(sin((Theta1C + Theta2C) / 2.0),2) + RhoTheta * RhoTheta);
            
            //cout<<C3D_min_con<<endl;
             
            P3D = 0.3 * B_avg / C3D_min_con;
            //cout<<P3D<<endl;
            
            phi01 = atan2((y1-y0),(x1-x0));    // verify this
            phi1 = 0.5 * PHI1C + phi01;
            phi2 = PHI1C + phi1;
           
            sigma_MS1 = calculate_MS_uncertainty( P3D, Theta1C, phi1, x1, y1, z1, det_type);
            
            sigma_MS2 = calculate_MS_uncertainty( P3D, Theta2C, phi2, x2, y2, z2, det_type);
            
            sigma_MS = (sigma_MS1 + sigma_MS2) / 2.0;
            
            chi2_min_con = (1.0 / (sigma_MS * sigma_MS)) * pow(( RhoTheta * PHI0_tilde - RhoPHI * Theta0_tilde), 2) / ( RhoPHI * RhoPHI + RhoTheta * RhoTheta / pow(sin((Theta1C + Theta2C) / 2.0), 2));
            
            sum_C3D += C3D_min_con/ num_events;
            
           // if(isnan(C3D_min_con))
           // { cout<<Cc<<"      "<<xx<<"    "<<x0<<"    "<<x1<<"   "<<x2<<endl;
           //   cout<<y0<<"    "<<y1<<"   "<<y2<<endl; 
           //  }
             
            sum_chi2 += chi2_min_con / num_events;
            
            sigma_C3D_con = sigma_MS * sqrt(1.0 / (RhoPHI * RhoPHI * pow(sin((Theta1C + Theta2C) / 2.0),2) + RhoTheta * RhoTheta));
              
            sigma_P3D = P3D * sigma_C3D_con / C3D_min_con;
            
            sum_sigma += sigma_P3D/ num_events;
            
           // PHI_MS = (beta) * (beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta); // there was a mistake here
            
           // Theta_MS = (-1.0 * eta * pow(sin((Theta1C + Theta2C) / 2.0),2)) * (beta * PHI_tilde - eta * Theta_tilde) / (eta * eta * pow(sin((Theta1C + Theta2C) / 2.0),2) + beta * beta);  // there was a mistake here as well 
             
            r = P3D / P_gen;
            
            //cout<<r<<endl;
            //cout<<r<<"\t"<<P3D<<"\t"<<sigma_P3D<<endl;
            
           // PhiMS.push_back(PHI_MS);
           // ThetaMS.push_back(Theta_MS);
            
           // h5->Fill(PHI_MS);
           // h6->Fill(Theta_MS);
            
            //cout<<Theta_MS<<"\t"<<PHI_MS<<endl;
            
            h1 -> Fill(r);
            h2 -> Fill(sigma_P3D);
            h3 -> Fill(chi2_min_con);
            pull_parameter = (P3D - P_gen) / sigma_P3D;
            h4 -> Fill(pull_parameter);              		  
        
        }                     
                   
            //PhiMS.clear();
            //ThetaMS.clear();
                  
     myFile.close();
   } 
    
     avg_chi2 = sum_chi2; 
      
     cout<<sum_sigma<<endl;
           
     cout<<avg_chi2<<endl;
     //cout<<sigma_P3D<<endl;
/*--------------------------------------------------------------------------------------------ENDCAP------------------------------------------------------------------------------------*/   
   
   
   TFile fileout("triplet_fit_output.root","recreate");
   fileout.cd();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   //c1->Divide(2,2);
   
   TCanvas *c2 = new TCanvas();
   c2->SetCanvasSize(1800, 1200);
   c2->SetWindowSize(1800, 1000);
   //c1->Divide(2,2);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   //c3->Divide(1,2);
   
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(1800, 1200);
   c4->SetWindowSize(1800, 1000);
   
   c1->cd();
   h1->Draw("hist PLC");
   c1->Modified(); c1->Update();
   
   c2->cd();
   h2->Draw("hist PLC");
   c2->Modified(); c2->Update();
   
   c3->cd();
   h3->Draw("hist PLC");
   c3->Modified(); c3->Update();
   
   c4->cd();
   h4->Draw("hist PLC");
   c4->Modified(); c4->Update();
    
  /* c2->cd();   
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
   */
   c1->Write();
   c1->Clear();
   
   c2->Write();
   c2->Clear();
   
   c3->Write();
   c3->Clear();
   
   c4->Write();
   c4->Clear();
  
   fileout.Close();
       
   return 0;
   
  }










