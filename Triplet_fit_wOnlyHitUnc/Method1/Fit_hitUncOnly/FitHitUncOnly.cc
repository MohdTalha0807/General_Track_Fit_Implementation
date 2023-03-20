
// Triplet fit for the hit uncertainty only case 
 
#include "Libraries.h"
#include <bits/stdc++.h>
#include "other_lib.h"
 
#define PI 3.1415926
#define B_avg 2.0 // in Tesla
  
using namespace std;
 
double calculate_w(double Theta_min_hit, double phi1, double PHI_C, double x, double y, double z) 
{  
   double phi0;
   double R;
   double w0;    // angle between the detector direction and the particle's direction
   double w;         
   
   phi0 = phi1 + PHI_C/2.0;
   //cout<<phi0<<endl;
// phi1 = phi01 - PHI1C/2.0;
// phi2 = phi12 - PHI2C/2.0;
    
   R = sqrt(x*x  + y*y);
   
   w0 = acos(sin(Theta_min_hit)*(cos(phi0) * x  +  sin(phi0) * y)/R);
   
   w = PI/2.0  - w0;
    
   return w; 
} 
      
double calculate_sigma_Cc(double PHI, double w0, double w1, double w2, double w3, double d01, double d12, double d02)
{
  double sigma_C;
  double sigma_phi;
  double sigma_t = 0.1/sqrt(12.0);
  
  sigma_phi = sqrt(pow(sin(w0)*sigma_t / d01, 2) + pow(sin(w3)*sigma_t / d12, 2) + pow(((sin(w1)*sigma_t / d01) + (sin(w2)*sigma_t / d12)), 2));
   
   sigma_C = cos(PHI/2.0) * sigma_phi / d02; 

   return sigma_C;
}
      
int main()
{   
   gROOT->SetStyle("Plain");
                                           
   TLatex Tl;                          
   Tl.SetTextSize(0.03);              
    
   TH1D* h1 = new TH1D("h1", "Triplet Fit P3D = 1 GeV (with only hit uncertainties); p(fit) / p(gen); Entries;", 60, 0.75, 1.4);
   TH1D* h2 = new TH1D("h2", "Triplet Fit; sigma_{p}(fit) [MeV]; Entries;", 100, 50.0, 90.0);
   TH1D* h3 = new TH1D("h3", "Triplet Fit (with only hit uncertainties); chi2_min_hit; Entries;", 100, 0.0, 15.0);
   TH1D* h4 = new TH1D("h4", "Triplet Fit; (p(fit) - p(gen)) / #sigma_{p}; Entries;", 120, -16.0, 16.0); 
   TH1D* h5 = new TH1D("h5", "Triplet Fit; #Theta_min_hit (degrees); Entries;", 100, 59.8, 60.2);
   //TH1D* h6 = new TH1D("h6", "Triplet Fit; #sigma_{#theta, hit} (degrees); Entries;", 100, 0.0403, 0.0407); 
   TH1D* h7 = new TH1D("h7", "Triplet Fit; R-R_{cal}; Entries;", 100, -0.001, 0.001); 
   TH1D* h8 = new TH1D("h8", "Triplet Fit; ((cot#theta)_{fit} - (cot#theta)_{gen}) / #sigma(cot#theta); Entries;", 60, -6.0, 6.0); 
    
      
   double P_gen = 1000.0;  // in MeV    
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
   double pull_parameter2;
   double sum_chi2 = 0.0;
   double avg_chi2;
   double sum_sigma, sum_C3D; 
   double S01, S12;
    
   double sigma_z = 0.1 / sqrt(12);   // 100 µm = 0.1 mm
    
   double Wk = 1.0 / (sigma_z * sigma_z);   // introduced weights
 
   double numerator, denominator, numerator2;
   double Theta1, Theta2, Theta_min_hit;
   double C3D_hit; 
   double th;     
   //int det_type = 0;    // 0 for barrel layers, 1 for the endcaps
   double S02, Z02;
   double phi; 
   double w0,w1,w2,w3;
   double PHI;
   double sigma_cotTheta;
   double sigma_Cc; 
   double phi01, phi12;
   double cotTheta_gen = cos(PI/3.0)/sin(PI/3.0);
   double cotTheta_fit;
   
   vector<double> phi_vec, Cc_vec, RR1_vec, RR3_vec, RR2_vec; 
   vector<double> sigmaTH_vec, ThetaMin_vec, S_vec; 
                  
   int num_events = 10000; 
      
   fstream myFile;
   myFile.open("triplet_sim_out_true12.txt", ios::in);
    
      //double xx = (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1)));
                  //cout<<Cc<<endl;
            
 double     R1 = 50.0;
 double     R2 = 100.0;
 double     R3 = 150.0;
      int count_dev = 0.0;  
              
   if(myFile.is_open())  
   {  
      for (int x = 0; x < num_events; x++)
       {
	    myFile >> x0 >> y0 >> z0 >> phi;
            myFile >> x1 >> y1 >> z1 >> phi;
            myFile >> x2 >> y2 >> z2 >> phi;
            
            double R1_cal = sqrt(x0*x0 + y0*y0);
            double R2_cal = sqrt(x1*x1 + y1*y1);
            double R3_cal = sqrt(x2*x2 + y2*y2);
          
            RR1_vec.push_back(R1-R1_cal);
            RR2_vec.push_back(R2-R2_cal);
            RR3_vec.push_back(R3-R3_cal);
          
         //cout<<R1_cal<<"    "<<R1<<endl;
         //cout<<R2_cal<<"    "<<R2<<endl;
         //cout<<R3_cal<<"    "<<R3<<endl;
         
            h7 -> Fill(R1-R1_cal);
            h7 -> Fill(R2-R2_cal);
            h7 -> Fill(R3-R3_cal);
            
            if(abs(R1-R1_cal) > 0.0 or abs(R2-R2_cal) > 0.0 or abs(R3-R3_cal) > 0.0)
           {
              count_dev++;
              //cout<<abs(R1-R1_cal)<<endl;       
           }
     
            phi_vec.push_back(phi * 180.0/PI);
            
            d01 = sqrt(pow((x1 -x0),2) + pow((y1 -y0),2));
            d12 = sqrt(pow((x2 -x1),2) + pow((y2 -y1),2));
            d02 = sqrt(pow((x2 -x0),2) + pow((y2 -y0),2));
            
            Z01 = z1 - z0;
            Z12 = z2 - z1; 
            Z02 = z2 - z0;
               
      	    Cc =  (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1))) / (d01 * d12 * d02);         
            Cc_vec.push_back(Cc);
             
            PHI1C = 2.0 * asin((d01 * Cc) / 2.0);
            PHI2C = 2.0 * asin((d12 * Cc) / 2.0);
            PHI = 2.0 * asin((Cc * d02) / 2.0);
            
            S01 = PHI1C / Cc;
            S12 = PHI2C / Cc;  
            S02 = PHI / Cc;
            
            S_vec.push_back(S02);
             
            Theta1 = atan2(S01, Z01); 
            
            Theta2 = atan2(S12, Z12); 
            
            numerator   =  S12*S12/Wk + S01*S01/Wk + S02*S02/Wk; 
            denominator =  S12*Z12/Wk + S01*Z01/Wk + S02*Z02/Wk;
            
            Theta_min_hit = atan2(numerator, denominator);
            
            cotTheta_fit = cos(Theta_min_hit)/sin(Theta_min_hit);
               
            ThetaMin_vec.push_back(Theta_min_hit * 180.0/PI);
                  
            h5 -> Fill(Theta_min_hit * 180.0/PI);
            
            C3D_hit = Cc * sin(Theta_min_hit); 
            
            P3D = 0.3 * B_avg / C3D_hit;
            //cout<<P3D<<endl;  
               
            //sigmaTH_vec.push_back(sigma_theta_hit * 180.0/PI);
            //h6 -> Fill(sigma_theta_hit * 180.0/PI);      
            //cout<<sigma_theta_hit<<endl;
            
            numerator2 = S01*S01*S12*S12* pow((cos(Theta2)/sin(Theta2) - cos(Theta1)/sin(Theta1)),2); 
            
            chi2_min_hit =  numerator2/ numerator;
            
            //chi2_vec.push_back(chi2_min_hit);
               
            sum_chi2 += chi2_min_hit / num_events; 
            
            sigma_cotTheta = sqrt( (3.0*Wk)  / (Wk*Wk*Wk*numerator));
            
            
            //sigmaTH_vec.push_back(sigma_cotTheta * 180.0/PI);
            //h6 -> Fill(sigma_cotTheta * 180.0/PI);      
            //cout<<sigma_cotTheta<<endl;
            
            phi01 = atan2(y1-y0, x1-x0);
            phi12 = atan2(y2-y1, x2-x1);
            //cout<<phi01<<"  "<<phi12<<endl;
          //double phi0 = phi01 + PHI1C/2.0;
          //double phi2 = phi12 - PHI2C/2.0;
          
          //cout<<phi0<<"    "<<phi2<<endl; 
             
            w0 = calculate_w(Theta_min_hit, phi01, -PHI1C, x0, y0, z0);
            w1 = calculate_w(Theta_min_hit, phi01, PHI1C, x1, y1, z1);
            w2 = w1;
            w3 = calculate_w(Theta_min_hit, phi12, PHI2C, x2, y2, z2);
            
            //cout<<w0<<"\t"<<w1<<"\t"<<w2<<"\t"<<w3<<endl;
            
            sigma_Cc  =  calculate_sigma_Cc(PHI, w0, w1, w2, w3, d01, d12, d02);   // figure out this 
            
            //cout<<sigma_Cc<<endl;
            
             
            //sigma_C3D_hit  = C3D_hit * sqrt(pow(sin(2.0*Theta_min_hit),2) * pow(sin(Theta_min_hit),2) * pow(sigma_cotTheta,2));
              
            sigma_C3D_hit  = C3D_hit * sqrt(pow((sigma_Cc/Cc), 2) + pow(sin(2.0*Theta_min_hit),2)/4.0 * pow(sigma_cotTheta,2));
             
            sigma_P3D = P3D * sigma_C3D_hit / C3D_hit;
            
            //cout<<sigma_C3D_hit/C3D_hit<<endl;
            cout<<sigma_P3D<<endl;
             
            r = P3D / P_gen;
               
            h1 -> Fill(r);
            h2 -> Fill(sigma_P3D);
            h3 -> Fill(chi2_min_hit);
            pull_parameter = (P3D - P_gen) / sigma_P3D;
            //cout<<pull_parameter<<endl;
            h4 -> Fill(pull_parameter); 
            
            pull_parameter2 = (cotTheta_fit - cotTheta_gen) / sigma_cotTheta;
            h8 -> Fill(pull_parameter2);
                         		  
              cout<<"----------------------------------------------------------------------------------------------------"<<endl;    
        }                     
                    
     myFile.close();
   } 
    
      
     avg_chi2 = sum_chi2; 
      
     //cout<<sum_sigma<<"\t"<<"\t"<<sum_C3D<<endl;
     cout<<count_dev<<endl;      
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
   //c2->Divide(1,2);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   c3->Divide(2,2);
   
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(1800, 1200);
   c4->SetWindowSize(1800, 1000);
   
   TCanvas *c5 = new TCanvas();
   c5->SetCanvasSize(1800, 1200);
   c5->SetWindowSize(1800, 1000);
   
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
    
   c2->cd();
   h5->Draw("hist");
   c2->Modified(); c2->Update();
   
  // c2->cd(2);
  // h6->Draw("hist");
  // c2->Modified(); c2->Update();
   
   c3->cd(1);   
   TGraph *gr1 = new TGraph(num_events, &phi_vec[0], &RR1_vec[0]);
   gr1->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 1GeV");
   gr1->GetYaxis()->SetTitle("R1-R1_{cal} [mm]");
   gr1->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr1->SetLineColor(kRed);
   gr1->Draw("AP*");
   
   c3->cd(2);   
   TGraph *gr2 = new TGraph(num_events, &phi_vec[0], &RR2_vec[0]);
   gr2->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 1GeV");
   gr2->GetYaxis()->SetTitle("R2-R2_{cal} [mm]");
   gr2->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr2->SetLineColor(kRed);
   gr2->Draw("AP*");
   
   c3->cd(3);   
   TGraph *gr3 = new TGraph(num_events, &phi_vec[0], &RR3_vec[0]);
   gr3->SetTitle("R-Rcal values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 1GeV");
   gr3->GetYaxis()->SetTitle("R3-R3_{cal} [mm]");
   gr3->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr3->SetLineColor(kRed);
   gr3->Draw("AP*");
   
   c3->cd(4);   
   TGraph *gr4 = new TGraph(num_events, &phi_vec[0], &Cc_vec[0]);
   gr4->SetTitle("C2D values for different #phi for particles simulated with polar angle #frac{#pi}{3} and P_{3D} = 1GeV");
   gr4->GetYaxis()->SetTitle("Cc");
   gr4->GetXaxis()->SetTitle("#phi [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr4->SetLineColor(kRed);
   gr4->Draw("AP*");
   
   ThetaMin_vec.clear();
   Cc_vec.clear();
   //chi2_vec.clear();
   RR1_vec.clear();
   RR2_vec.clear();
   RR3_vec.clear();
   
   c4->cd();
   h7->Draw("hist");
   c4->Modified(); c4->Update();
   
   c5->cd();
   h8->Draw("hist PLC");
   c5->Modified(); c5->Update();
   
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
   
   c4->Write();
   c4->Clear();
   
   c5->Write();
   c5->Clear();
  
   fileout.Close();
       
   return 0;
   
  }










