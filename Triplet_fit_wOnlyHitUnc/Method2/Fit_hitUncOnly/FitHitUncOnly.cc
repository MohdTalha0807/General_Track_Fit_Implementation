
// Triplet fit for the hit uncertainty only case 
 
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
    
   TH1D* h1 = new TH1D("h1", "Triplet Fit P3D = 20 GeV (with only hit uncertainties); p(fit) / p(gen); Entries;", 60, 0.96, 1.04);
   TH1D* h2 = new TH1D("h2", "Triplet Fit; sigma_{p}(fit) [MeV]; Entries;", 50, 19.4, 20.6);
   TH1D* h3 = new TH1D("h3", "Triplet Fit (with only hit uncertainties); chi2_min_hit; Entries;", 100, 0.0, 15.0);
   TH1D* h4 = new TH1D("h4", "Triplet Fit; (p(fit) - p(gen)) / sigma_{p}; Entries;", 100, -15.0, 15.0); 
   TH1D* h5 = new TH1D("h5", "Triplet Fit; #Theta_min_hit (degrees); Entries;", 100, 44.8, 45.2);
   TH1D* h6 = new TH1D("h6", "Triplet Fit; #sigma_{#theta, hit} (degrees); Entries;", 100, 0.0403, 0.0407); 
    
      
   double P_gen = 20000.0;  // in MeV    
   double Cc;   // transverse curvature for the circle solution 
   double d01;
   double d12;
   double d02;
   double PHI1C;  
   double PHI2C;
   double Z01;
   double Z12;
   double x0,x1,x2;
   double y0,y1,y2;
   double z0,z1,z2;
   double r, rt;
   double chi2_min_hit;
   double P3D;
   double sigma_MS;
   double sigma_C3D_hit;
   double sigma_P3D;
   double pull_parameter;
   double sum_chi2 = 0.0;
   double avg_chi2;
   double sum_sigma, sum_C3D; 
   double S01, S12;
   
   double sigma_z = 0.1 / sqrt(12);   // 100 µm = 0.1 mm
    
   double Wk = 12.0 / 0.01;   // introduced weights
   
   double zeta0, zeta1, zeta2;
   double sigma_theta_hit;
   double sum_WkSkZk, sum_Wk, sum_WkZk, sum_WkSk, sum_WkSk2;
   double numerator, denominator;
   double Theta1, Theta2, Theta_min_hit;
   double C3D_hit; 
   double th;     
   //int det_type = 0;    // 0 for barrel layers, 1 for the endcaps
   double S0, S1, S2;
   double phi; 
   vector<double> phi_vec; 
   vector<double> Cc_vec, chi2_vec, sigmaTH_vec;
                  
   int num_events = 10000; 
      
   fstream myFile;
   myFile.open("triplet_sim_out_true12.txt", ios::in);
    
      //double xx = (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1)));
            //cout<<Cc<<endl; 
   if(myFile.is_open())
   {  
      for (int x = 0; x < num_events; x++)
       {
	    myFile >> x0 >> y0 >> z0 >> phi;
            myFile >> x1 >> y1 >> z1 >> phi;
            myFile >> x2 >> y2 >> z2 >> phi;
            
            phi_vec.push_back(phi * 180.0/PI);
            
            d01 = sqrt(pow((x1 -x0),2) + pow((y1 -y0),2));
            d12 = sqrt(pow((x2 -x1),2) + pow((y2 -y1),2));
            d02 = sqrt(pow((x2 -x0),2) + pow((y2 -y0),2));
            
            Z01 = z1 - z0;
            Z12 = z2 - z1;
               
      	    Cc =  (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1))) / (d01 * d12 * d02);   // or there is a problem here	
            
            Cc_vec.push_back(Cc);
             
            PHI1C = 2.0 * asin((d01 * Cc) / 2.0);
            PHI2C = 2.0 * asin((d12 * Cc) / 2.0);
           
            S01 = PHI1C / Cc;
            S12 = PHI2C / Cc;
             
            S0 = -S01;
            S1 = 0.0;
            S2 = S12;
          /*  
            sum_WkSkZk =  Wk * S0 * z0  +  Wk * S1 * z1 + Wk * S2 * z2; 
            sum_Wk     =  3.0 * Wk;
            sum_WkZk   =  Wk * (z0 + z1 + z2);
            sum_WkSk   =  Wk * (S0 + S1 + S2);
            sum_WkSk2  =  Wk * (S0 * S0 + S1 * S1 + S2 * S2);            
            numerator   = sum_WkSkZk * sum_Wk  - sum_WkZk * sum_WkSk;
            denominator = sum_WkSk2 * sum_Wk  - sum_WkSk * sum_WkSk;
            */
            
            numerator   =  S12*S12/Wk + S01*S01/Wk + S02*S02/Wk; 
            denominator =  S12*Z12/Wk + S01*Z01/Wk + S02*Z02/Wk;
            
            Theta_min_hit = atan2(numerator, denominator);
            
            //th = Theta_min_hit *180.0 /PI;
            //cout<<Theta_min_hit<<endl;
            
            h5 -> Fill(Theta_min_hit * 180.0/PI);
            
            C3D_hit = Cc * sin(Theta_min_hit);    // either this relationship is incorrect
            
            P3D = 0.3 * B_avg / C3D_hit;
            //cout<<P3D<<endl;  
            
            Theta1 = atan2(S01, Z01); 
            
            Theta2 = atan2(S12, Z12);
            
            zeta0 = - (Cc * PHI1C)/(PHI1C * PHI1C + Z01*Z01 * Cc * Cc); 
            zeta1 =  (Cc * PHI1C)/(PHI1C * PHI1C + Z01*Z01 * Cc * Cc)  + (Cc * PHI2C)/(PHI2C * PHI2C + Z12*Z12 * Cc * Cc);
            zeta2 =  - (Cc * PHI2C)/(PHI2C * PHI2C + Z12*Z12 * Cc * Cc); 
            
            sigma_theta_hit = sigma_z * sqrt(pow(zeta0,2) + pow(zeta1,2) + pow(zeta2,2)); 
            
            sigmaTH_vec.push_back(sigma_theta_hit * 180.0/PI);
            
            h6 -> Fill(sigma_theta_hit * 180.0/PI);
            
            //cout<<sigma_theta_hit<<endl;
            
            chi2_min_hit = pow((Theta2 - Theta1),2) / pow(sigma_theta_hit,2); 
            
            chi2_vec.push_back(chi2_min_hit);
            
            //sum_C3D += C3D_hit/ num_events;
               
            sum_chi2 += chi2_min_hit / num_events;
            
            sigma_C3D_hit = Cc * sigma_theta_hit;
              
            sigma_P3D = P3D * sigma_C3D_hit / C3D_hit;
            
            //cout<<sigma_P3D<<endl;
             
            r = P3D / P_gen;
               
            h1 -> Fill(r);
            h2 -> Fill(sigma_P3D);
            h3 -> Fill(chi2_min_hit);
            pull_parameter = (P3D - P_gen) / sigma_P3D;
            //cout<<pull_parameter<<endl;
            h4 -> Fill(pull_parameter);              		  
        
        }                     
                    
     myFile.close();
   } 
    
     
     avg_chi2 = sum_chi2; 
      
     //cout<<sum_sigma<<"\t"<<"\t"<<sum_C3D<<endl;
           
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
   c2->Divide(1,2);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   //c3->Divide(1,2);
 /*  
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(1800, 1200);
   c4->SetWindowSize(1800, 1000);
   */
   c1->cd(1);
   h1->Draw("hist PLC");
   c1->Modified(); c1->Update();
   
   c1->cd(2);
   h2->Draw("hist PLC");
   c1->Modified(); c1->Update();
   
   c1->cd(3);
   h3->Draw("hist PLC");
   c1->Modified(); c1->Update();
   
   c1->cd(4);
   h4->Draw("hist PLC");
   c1->Modified(); c1->Update();
    
   c2->cd(1);
   h5->Draw("hist");
   c2->Modified(); c2->Update();
   
   c2->cd(2);
   h6->Draw("hist");
   c2->Modified(); c2->Update();
   
   c3->cd();   
   TGraph *gr1 = new TGraph(num_events, &phi_vec[0], &sigmaTH_vec[0]);
   gr1->SetTitle("#sigma_{#theta, hit} values for different #phi for particles simulated with polar angle #frac{#pi}{4} and P_{3D} = 20GeV");
   gr1->GetYaxis()->SetTitle("#sigma_{#theta, hit} (degrees)");
   gr1->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr1->SetLineColor(kRed);
   gr1->Draw("AP*");
  
   Cc_vec.clear();
   phi_vec.clear();
   chi2_vec.clear();
   
 /* 
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
 /*  
   c4->Write();
   c4->Clear();
  */
   fileout.Close();
       
   return 0;
   
  }










