
// General fit: dominating hit uncertainty case 
 
#include "Libraries.h"
#include <bits/stdc++.h>
#include "other_lib.h"
     
#define PI 3.1415926 
#define B_avg 2.0 // in Tesla 
  
using namespace std; 
  
vector<double> calculate_h_transPlane_t( double d01, double d12, double d02)
{
   double h0;
   double h1;
   double h2;
   vector<double> h_t;
   
   h0 = 1.0 / d01;
   
   h1 =  -1.0/d12  - 1.0/d01;
   
   h2 = 1.0 /d12;
  
   h_t.push_back(h0);
   h_t.push_back(h1);
   h_t.push_back(h2);
   
   return h_t;
}  

vector<double> calculate_h_longPlane_z( double d01, double d12, double d02, double Z01, double Z12)
{
   double h0;
   double h1;
   double h2;
   vector<double> h_z;
   
   h0 = - d01 / (Z01 * Z01 + d01 * d01);
   
   h1 = d01 / (Z01 * Z01 + d01 * d01)  + d12 / (Z12 * Z12 + d12 * d12); 
   
   h2 = - d12 / (Z12 * Z12 + d12 * d12);
  
   h_z.push_back(h0);
   h_z.push_back(h1);
   h_z.push_back(h2);
   
   return h_z;
}  

   
int main()
{  
   gROOT->SetStyle("Plain");
                                           
   TLatex Tl;                          
   Tl.SetTextSize(0.03);              
   
   TH1D* h1 = new TH1D("h1", "Triplet Fit P3D = 1 GeV (0.1mmXYZ uncertainty); p(fit) / p(gen); Entries;", 100, 0.9, 1.1); 
   TH1D* h2 = new TH1D("h2", "Triplet Fit; sigma_{p}(fit) [MeV]; Entries;", 150, 4.0, 14.0);
   TH1D* h3 = new TH1D("h3", "Triplet Fit ; chi2_min_con; Entries;", 100, 0.0, 12.0);
   TH1D* h4 = new TH1D("h4", "Triplet Fit; (p(fit) - p(gen)) / sigma_{p}; Entries;", 80, -5.0, 5.0);  
    
   TFile *fileIn = new TFile("triplet_sim_output.root","READ");
   TTree *treeIn = (TTree*)fileIn->Get("T"); 
   
   vector<int> *Event_ID = 0;
   vector<double> *Cx = 0;
   vector<double> *Cy = 0;
   vector<double> *Cz = 0;
   vector<double> *phi_gen = 0; 
   vector<double> *theta_gen = 0;
   
   treeIn->SetBranchAddress("Evt_ID", &Event_ID); 
   treeIn->SetBranchAddress("CoordX",&Cx);  
   treeIn->SetBranchAddress("CoordY",&Cy);    
   treeIn->SetBranchAddress("CoordZ",&Cz);
   treeIn->SetBranchAddress("In_Phi",&phi_gen);
   treeIn->SetBranchAddress("In_Theta",&theta_gen);
      
         
   double P_gen = 1000.0;  // in MeV    
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
   double C3D_hit_min;
   double r, rt;
   double RhoPHI, RhoTheta;
   double chi2_3D_min_hit;
   double P3D;
   double sigma_C3D_hit;
   double sigma_P3D;
   double pull_parameter;
   double phi1, phi2, phi01;
   double sum_chi2 = 0.0;
   double avg_chi2;
   double sum_sigma, sum_C3D; 
   double phi; 
   double theta;
   vector<double> h_transPlane_t;
   vector<double> h_longPlane_z;
   vector<double> theta_vec, sigmaC3D_vec;
   vector<double> rhoPhi_vec, rhoTheta_vec; 
   vector<double> gamma_long_vec, gamma_trans_vec; 
   
   double gamma_long2, gamma_trans2;
   
   double sigma_z = 0.1 / sqrt(12.0);
   double sigma_t = 0.1 / sqrt(12.0);
      
   int det_type = 0;    // 0 for barrel layers, 1 for the endcaps
               
   int num_events = 10000; 
   
   double     R1 = 100.0;    // radii of the barrel layers all in mm
   double     R2 = 200.0; 
   double     R3 = 300.0;  
   
   int nentries = int(treeIn->GetEntries());
   cout<<"number of entries in TTree = "<<nentries<<endl;
        UInt_t j1 = 0, j2 = 1, j3 = 2;        
      
      for (int l = 0; l < nentries; l++)
     {
	    treeIn -> GetEntry(l);
             
            x0 = Cx->at(j1); y0 = Cy->at(j1); z0 = Cz->at(j1);
            x1 = Cx->at(j2); y1 = Cy->at(j2); z1 = Cz->at(j2);
            x2 = Cx->at(j3); y2 = Cy->at(j3); z2 = Cz->at(j3);
            phi = phi_gen->at(j1); 
            theta = theta_gen->at(j1);
            
            theta_vec.push_back(theta);
            
            d01 = sqrt(pow((x1 -x0),2) + pow((y1 -y0),2));
            d12 = sqrt(pow((x2 -x1),2) + pow((y2 -y1),2));
            d02 = sqrt(pow((x2 -x0),2) + pow((y2 -y0),2));
            
            Z01 = z1 - z0;
            Z12 = z2 - z1;
               
      	    Cc =  (2.0 * ((x1 -x0)*(y2-y1) - (y1-y0)*(x2-x1))) / (d01 * d12 * d02);	
             
            PHI1C = 2.0 * asin((d01 * Cc) / 2.0);
            PHI2C = 2.0 * asin((d12 * Cc) / 2.0);
            
            C3D1C = PHI1C / sqrt(Z01 * Z01 + d01 * d01 * (PHI1C * PHI1C)/(4.0 * pow(sin(PHI1C/2.0),2))); 
            C3D2C = PHI2C / sqrt(Z12 * Z12 + d12 * d12 * (PHI2C * PHI2C)/(4.0 * pow(sin(PHI2C/2.0),2)));

            Theta1C = acos((Z01 * C3D1C)/PHI1C);
            Theta2C = acos((Z12 * C3D2C)/PHI2C); 
            
            n1C = 1.0/((PHI1C/2.0) * (cos(PHI1C/2.0)/sin(PHI1C/2.0)) * pow(sin(Theta1C),2) + pow(cos(Theta1C),2));
            n2C = 1.0/((PHI2C/2.0) * (cos(PHI2C/2.0)/sin(PHI2C/2.0)) * pow(sin(Theta2C),2) + pow(cos(Theta2C),2));      
            
            PHI0_tilde = 0.5* (PHI1C * n1C + PHI2C * n2C);
            
            Theta0_tilde = Theta2C - Theta1C + ((1.0 - n2C) * (cos(Theta2C)/sin(Theta2C)) - (1.0 - n1C) * (cos(Theta1C)/sin(Theta1C)));
            
            RhoPHI = -0.5 * ((PHI1C * n1C / C3D1C) + (PHI2C * n2C / C3D2C));
            rhoPhi_vec.push_back(RhoPHI);
            
            RhoTheta = (1.0 - n1C) * (cos(Theta1C)/sin(Theta1C)) / C3D1C  -  (1.0 - n2C) * (cos(Theta2C)/sin(Theta2C)) / C3D2C; 
            rhoTheta_vec.push_back(RhoTheta);
  
            h_transPlane_t = calculate_h_transPlane_t( d01, d12, d02);                              
            h_longPlane_z  = calculate_h_longPlane_z( d01, d12, d02, Z01, Z12);
            
            gamma_trans2 = pow((h_transPlane_t[0]*sigma_t),2) + pow((h_transPlane_t[1]*sigma_t),2) + pow((h_transPlane_t[2]*sigma_t),2);
            //gamma_trans2 = gamma_trans2 * pow(sin((Theta1C + Theta2C) / 2.0),2);
            gamma_trans_vec.push_back(gamma_trans2);
            
            gamma_long2  = pow((h_longPlane_z[0]*sigma_z),2) + pow((h_longPlane_z[1]*sigma_z),2) + pow((h_longPlane_z[2]*sigma_z),2);
            gamma_long_vec.push_back(gamma_long2);
            
            
            C3D_hit_min = - (Theta0_tilde * RhoTheta * gamma_trans2 + PHI0_tilde * RhoPHI * gamma_long2) / (RhoTheta * RhoTheta * gamma_trans2 + RhoPHI * RhoPHI * gamma_long2);
            
            chi2_3D_min_hit = pow((Theta0_tilde * RhoPHI - PHI0_tilde * RhoTheta), 2) / (RhoTheta * RhoTheta * gamma_trans2 + RhoPHI * RhoPHI * gamma_long2);
            
            sigma_C3D_hit = sqrt( (gamma_trans2 * gamma_long2) / (RhoTheta * RhoTheta * gamma_trans2 + RhoPHI * RhoPHI * gamma_long2) );
            sigmaC3D_vec.push_back(sigma_C3D_hit);
            
            
            P3D = 0.3 * B_avg / C3D_hit_min;
            
            sum_chi2 += chi2_3D_min_hit / num_events;
              
            sigma_P3D = P3D * sigma_C3D_hit / C3D_hit_min;
            //cout<<sigma_P3D<<endl;       
       
            sum_sigma += sigma_P3D/ num_events;
               
            r = P3D / P_gen;
            
            h1 -> Fill(r);
            h2 -> Fill(sigma_P3D);
            h3 -> Fill(chi2_3D_min_hit);
            pull_parameter = (P3D - P_gen) / sigma_P3D;
            h4 -> Fill(pull_parameter);              		  
        
        }                     
                    
     avg_chi2 = sum_chi2; 
      
     //cout<<sum_sigma<<endl;
           
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
   
   TCanvas *c5 = new TCanvas();
   c5->SetCanvasSize(1800, 1200);
   c5->SetWindowSize(1800, 1000);
   c5->Divide(2,2);
   
   TCanvas *c6 = new TCanvas();
   c6->SetCanvasSize(1800, 1200);
   c6->SetWindowSize(1800, 1000);
   c6->Divide(2,2);
   
   TCanvas *c7 = new TCanvas();
   c7->SetCanvasSize(1800, 1200);
   c7->SetWindowSize(1800, 1000);
   //c6->Divide(2,2);
   
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
   
   c5->cd(1);
   h1->Draw("hist PLC");
   c5->Modified(); c5->Update();
   
   c5->cd(2);
   h2->Draw("hist PLC");
   c5->Modified(); c5->Update();
   
   c5->cd(3);
   h3->Draw("hist PLC");
   c5->Modified(); c5->Update();
   
   c5->cd(4);
   h4->Draw("hist PLC");
   c5->Modified(); c5->Update();
  
   c6->cd(1);   
   TGraph *gr2 = new TGraph(num_events, &theta_vec[0], &rhoPhi_vec[0]);
   gr2->SetTitle("Control Plot (#rho_{#phi})");
   gr2->GetYaxis()->SetTitle("#rho_{#phi}");
   gr2->GetXaxis()->SetTitle("#theta [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr2->SetLineColor(kRed);
   gr2->Draw("AP*");
   
   c6->cd(2);   
   TGraph *gr3 = new TGraph(num_events, &theta_vec[0], &rhoTheta_vec[0]);
   gr3->SetTitle("Control Plot (#rho_{#theta})");
   gr3->GetYaxis()->SetTitle("#rho_{#theta}");
   gr3->GetXaxis()->SetTitle("#theta [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr3->SetLineColor(kRed);
   gr3->Draw("AP*");
   
   c6->cd(3);   
   TGraph *gr4 = new TGraph(num_events, &theta_vec[0], &gamma_trans_vec[0]);
   gr4->SetTitle("Control Plot (#Gamma^{2}_{#perp})");
   gr4->GetYaxis()->SetTitle("#Gamma^{2}_{#perp}");
   gr4->GetXaxis()->SetTitle("#theta [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr4->SetLineColor(kRed);
   gr4->Draw("AP*");
   
   c6->cd(4);   
   TGraph *gr5 = new TGraph(num_events, &theta_vec[0], &gamma_long_vec[0]);
   gr5->SetTitle("Control Plot (#Gamma^{2}_{#parallel})");
   gr5->GetYaxis()->SetTitle("#Gamma^{2}_{#parallel}");
   gr5->GetXaxis()->SetTitle("#theta [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr5->SetLineColor(kRed);
   gr5->Draw("AP*");
   
   c7->cd();   
   TGraph *gr1 = new TGraph(num_events, &theta_vec[0], &sigmaC3D_vec[0]);
   gr1->SetTitle("Control Plot (#sigma_{C3D})");
   gr1->GetYaxis()->SetTitle("#sigma_{C3D}");
   gr1->GetXaxis()->SetTitle("#theta [degrees]");
   //gr1->GetYaxis()->SetLimits(0.0,150.0);
   //gr1->GetXaxis()->SetLimits(-1.0, 1.0);
   gr1->SetLineColor(kRed);
   gr1->Draw("AP*");
   
    
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
   
   c6->Write();
   c6->Clear();
  
   c7->Write();
   c7->Clear();
  
   fileout.Close();
       
   return 0;
   
  }










