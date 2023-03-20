#include "Libraries.h"
#include <bits/stdc++.h>

using namespace std;

int main()
{  
   gROOT->SetStyle("Plain");
                                           
   TLatex Tl;                          
   Tl.SetTextSize(0.03);        
   
   TMultiGraph  *mg1  = new TMultiGraph();
   TMultiGraph  *mg2  = new TMultiGraph();
   TMultiGraph  *mg3  = new TMultiGraph();
   TMultiGraph  *mg4  = new TMultiGraph();
   
   vector<double> Xeff;
   vector<double> chi2;
   
   vector<double> B_avg;
   vector<double> chi2_2;
   int n1, n2;
   
   vector<double> Xeff_con;
   vector<double> chi2_con;
   
   vector<double> B_avg_con;
   vector<double> chi2_2_con;
   int n1con, n2con;
     
 /* --------------------------------------------------------------chi2 dist ---------------------------------------------------------*/   
   Xeff.push_back(0.01); chi2.push_back(1.004);
   Xeff.push_back(0.02); chi2.push_back(1.070);
   Xeff.push_back(0.05); chi2.push_back(1.151);
   Xeff.push_back(0.08); chi2.push_back(1.204);
   Xeff.push_back(0.10); chi2.push_back(1.261);
   Xeff.push_back(0.20); chi2.push_back(1.301);
   Xeff.push_back(0.50); chi2.push_back(1.328);
   Xeff.push_back(2.0); chi2.push_back(0.850);
   Xeff.push_back(5.0); chi2.push_back(0.427);
   Xeff.push_back(10.0); chi2.push_back(0.233);
   
   B_avg.push_back(2.0); chi2_2.push_back(1.000);
   B_avg.push_back(1.0); chi2_2.push_back(1.043);
   B_avg.push_back(0.5); chi2_2.push_back(1.032);
   B_avg.push_back(0.1); chi2_2.push_back(0.466);
   B_avg.push_back(0.05); chi2_2.push_back(0.142);
   B_avg.push_back(0.01); chi2_2.push_back(0.006);
   
   n1= Xeff.size();
   n2= B_avg.size();
   
   Xeff_con.push_back(0.01); chi2_con.push_back(0.9935);
   Xeff_con.push_back(0.02); chi2_con.push_back(1.036);
   Xeff_con.push_back(0.05); chi2_con.push_back(1.012);
   Xeff_con.push_back(0.08); chi2_con.push_back(1.066);
   Xeff_con.push_back(0.10); chi2_con.push_back(1.087);
   Xeff_con.push_back(0.20); chi2_con.push_back(1.24);
   Xeff_con.push_back(0.30); chi2_con.push_back(5.02);
   Xeff_con.push_back(0.50); chi2_con.push_back(45.45);
  
   B_avg_con.push_back(2.0); chi2_2_con.push_back(1.007);
   B_avg_con.push_back(1.0); chi2_2_con.push_back(1.044);
   B_avg_con.push_back(0.5); chi2_2_con.push_back(1.130);
   B_avg_con.push_back(0.4); chi2_2_con.push_back(1.276);
   B_avg_con.push_back(0.2); chi2_2_con.push_back(129.74);
   
   n1con = Xeff_con.size();
   n2con = B_avg_con.size();

//----------------------------------------------------------------------------------------------------------------------------------------------//
  
   vector<double> Xeff_un;
   vector<double> sigma_p_un;
   
   vector<double> B_avg_un;
   vector<double> sigma2_p_un;
   int n1_un, n2_un;
   
   vector<double> Xeff_co;
   vector<double> sigma_p_co;
   
   vector<double> B_avg_co;
   vector<double> sigma2_p_co;
   int n1_co, n2_co;

   Xeff_un.push_back(0.01); sigma_p_un.push_back(0.90089);
   Xeff_un.push_back(0.05); sigma_p_un.push_back(2.0234);
   Xeff_un.push_back(0.08); sigma_p_un.push_back(2.56868);
   Xeff_un.push_back(0.1); sigma_p_un.push_back(2.88979);
   Xeff_un.push_back(0.2); sigma_p_un.push_back(4.15653);
   Xeff_un.push_back(0.5); sigma_p_un.push_back(7.04862);
   Xeff_un.push_back(2.0); sigma_p_un.push_back(14.2599);
   Xeff_un.push_back(5.0); sigma_p_un.push_back(16.6219);
   Xeff_un.push_back(10.0); sigma_p_un.push_back(17.6365);
   
   B_avg_un.push_back(2.0); sigma2_p_un.push_back(0.901108);
   B_avg_un.push_back(1.0); sigma2_p_un.push_back(1.80798);
   B_avg_un.push_back(0.5); sigma2_p_un.push_back(3.66223);
   B_avg_un.push_back(0.1); sigma2_p_un.push_back(18.56300);
   B_avg_un.push_back(0.05); sigma2_p_un.push_back(23.0353);
   B_avg_un.push_back(0.01); sigma2_p_un.push_back(24.4765);
   
   n1_un= Xeff_un.size();
   n2_un= B_avg_un.size();
   
   Xeff_co.push_back(0.01); sigma_p_co.push_back(0.90518);
   Xeff_co.push_back(0.05); sigma_p_co.push_back(2.19399);
   Xeff_co.push_back(0.08); sigma_p_co.push_back(2.8582);
   Xeff_co.push_back(0.1); sigma_p_co.push_back(3.24197);
   Xeff_co.push_back(0.20); sigma_p_co.push_back(4.90706);
   Xeff_co.push_back(0.50); sigma_p_co.push_back(9.70885);
   Xeff_co.push_back(2.0); sigma_p_co.push_back(18.9973);
   Xeff_co.push_back(5.0); sigma_p_co.push_back(17.5137);
   Xeff_co.push_back(10.0); sigma_p_co.push_back(27.5393);
   
   B_avg_co.push_back(2.0); sigma2_p_co.push_back(0.906035);
   B_avg_co.push_back(1.0); sigma2_p_co.push_back(1.82232);
   B_avg_co.push_back(0.5); sigma2_p_co.push_back(3.75283);
   B_avg_co.push_back(0.1); sigma2_p_co.push_back(19.7226);
   B_avg_co.push_back(0.05); sigma2_p_co.push_back(25.4188);
   B_avg_co.push_back(0.01); sigma2_p_co.push_back(59.5686);
      
   n1_co = Xeff_co.size();
   n2_co = B_avg_co.size();

//-----------------------------------------------------------------------------------------------------------------------------------------//
  
   TGraph *gr1 = new TGraph(n1, &Xeff[0], &chi2[0]);
   TGraph *gr2 = new TGraph(n1con, &Xeff_con[0], &chi2_con[0]);
   gr1->SetLineColor(kBlue);                                                        //constrained - red,   unconstrained - blue
   gr2->SetLineColor(kRed);
   mg1->Add(gr1);   
   mg1->Add(gr2);
      
   TGraph *gr3 = new TGraph(n2, &B_avg[0], &chi2_2[0]);
   TGraph *gr4 = new TGraph(n2con, &B_avg_con[0], &chi2_2_con[0]);
   gr3->SetLineColor(kBlue);
   gr4->SetLineColor(kRed);
   mg2->Add(gr3);   
   mg2->Add(gr4);

//-------------------------------------------------------------------------
   
   TGraph *gr5 = new TGraph(n1_un, &Xeff_un[0], &sigma_p_un[0]);
   TGraph *gr6 = new TGraph(n1_co, &Xeff_co[0], &sigma_p_co[0]);
   gr5->SetLineColor(kBlue);
   gr6->SetLineColor(kRed);
   mg3->Add(gr5);   
   mg3->Add(gr6);
      
   TGraph *gr7 = new TGraph(n2_un, &B_avg_un[0], &sigma2_p_un[0]);
   TGraph *gr8 = new TGraph(n2_co, &B_avg_co[0], &sigma2_p_co[0]);
   gr7->SetLineColor(kBlue);     
   gr8->SetLineColor(kRed);
   mg4->Add(gr7);   
   mg4->Add(gr8);

     
   TFile fileout("triplet_fit_output.root","recreate");
   fileout.cd();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   //c1->Divide(1,2);
   
   TCanvas *c2 = new TCanvas();
   c2->SetCanvasSize(1800, 1200);
   c2->SetWindowSize(1800, 1000);
   
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(1800, 1200);
   c3->SetWindowSize(1800, 1000);
   
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(1800, 1200);
   c4->SetWindowSize(1800, 1000);
   
   c1->cd();   
   mg1->SetTitle("#chi^{2} values for varying material thickness");
   mg1->GetXaxis()->SetTitle("X_{eff}");
   mg1->GetYaxis()->SetTitle("<#chi^{2}>");
   mg1->Draw("AL*");
   
   c2->cd();   
   mg2->SetTitle("#chi^{2} values for varying magnetic field");
   mg2->GetXaxis()->SetTitle("B_{avg} [T]");
   mg2->GetYaxis()->SetTitle("<#chi^{2}>");
   mg2->Draw("AL*");
   
   c3->cd();   
   mg3->SetTitle("<#sigma_{p}> values for varying material thickness (P3D = 20GeV)");
   mg3->GetXaxis()->SetTitle("X_{eff}");
   mg3->GetYaxis()->SetTitle("<#sigma_{p}> [GeV]");
   mg3->Draw("AL*");
   
   c4->cd();   
   mg4->SetTitle("<#sigma_{p}> values for varying magnetic field (P3D = 20GeV)");
   mg4->GetXaxis()->SetTitle("B_{avg} [T]");
   mg4->GetYaxis()->SetTitle("<#sigma_{p}> [GeV]");
   mg4->Draw("AL*");
   
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
