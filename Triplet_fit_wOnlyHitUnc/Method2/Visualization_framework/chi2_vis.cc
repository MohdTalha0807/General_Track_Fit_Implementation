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
   
   TGraph *gr1 = new TGraph(n1, &Xeff[0], &chi2[0]);
   TGraph *gr2 = new TGraph(n1con, &Xeff_con[0], &chi2_con[0]);
   gr1->SetLineColor(kBlue);
   gr2->SetLineColor(kRed);
   mg1->Add(gr1);   
   mg1->Add(gr2);
      
   TGraph *gr3 = new TGraph(n2, &B_avg[0], &chi2_2[0]);
   TGraph *gr4 = new TGraph(n2con, &B_avg_con[0], &chi2_2_con[0]);
   gr3->SetLineColor(kBlue);
   gr4->SetLineColor(kRed);
   mg2->Add(gr3);   
   mg2->Add(gr4);
     
   TFile fileout("triplet_fit_output.root","recreate");
   fileout.cd();
   
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(1800, 1200);
   c1->SetWindowSize(1800, 1000);
   //c1->Divide(1,2);
   
   TCanvas *c2 = new TCanvas();
   c2->SetCanvasSize(1800, 1200);
   c2->SetWindowSize(1800, 1000);
   
   c1->cd();   
   mg1->SetTitle("#chi^{2} values for varying material thickness");
   mg1->GetXaxis()->SetTitle("X_{eff}");
   mg1->GetYaxis()->SetTitle("<#chi^{2}>");
   mg1->Draw("AL*");
   //mg1->GetXaxis()->SetRangeUser(-0.1, 11.0);
   //mg1->GetYaxis()->SetRangeUser(-0.1, 1.8);
   
   c2->cd();   
   mg2->SetTitle("#chi^{2} values for varying magnetic field");
   mg2->GetXaxis()->SetTitle("B_{avg} [T]");
   mg2->GetYaxis()->SetTitle("<#chi^{2}>");
   mg2->Draw("AL*");
   //mg2->GetXaxis()->SetRangeUser(-0.1, 2.5);
   //mg2->GetYaxis()->SetRangeUser(-0.1, 1.8);
   
   c1->Write();
   c1->Clear();
   
   c2->Write();
   c2->Clear();
    
   fileout.Close();
       
   return 0;
   
  }
