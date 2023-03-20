
using namespace std;

double calculate_x_eff(vector<double> &p_dir, vector<double> &det_dir)                   // function to calculate the effective path length of the particle inside the material
{
   double x_eff;
   x_eff = 0.01/ (p_dir[0]*det_dir[0] + p_dir[1]*det_dir[1] + p_dir[2]*det_dir[2]);
   
   //cout<<x_eff<<"\t"<<p_dir[0]<<"\t"<<p_dir[1]<<"\t"<<p_dir[2]<<endl;
     
   return x_eff;   
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







