#include <stdlib.h>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/EulerAngles>
#include <iostream>
#include <chrono>
#include <unistd.h>


class Denny : public Eigen::MatrixXd
{
    public:
        Denny(void):Eigen::MatrixXd() {} 

        template<typename OtherDerived>
        Denny(const Eigen::MatrixBase<OtherDerived>& other) : Eigen::MatrixXd(other) {}

        template<typename OtherDerived>
        Denny& operator =(const Eigen::MatrixBase <OtherDerived>& other) 
        {
            this->Eigen::MatrixXd::operator = (other);
            return *this;
        }
        void SetCart(double X, double Y, double Z);
        void SetDH(double t,double d,double a,double i);
        Eigen::Vector3d GetZ(void);    
        Eigen::Vector3d GetC(void);
        Eigen::Matrix4d H;
    private:

};

void Denny::SetCart(double X, double Y, double Z)
{
    this->block(0,3,3,1) << X, Y, Z;   
}

void Denny::SetDH(double t, double d, double a, double i)
{
   Eigen::MatrixXd M(4,4);
 
   M(0,0) =  std::cos(t);
   M(0,1) =  (-1*std::sin(t)*std::cos(i));
   M(0,2) =  (std::sin(t)*std::sin(i));
   M(0,3) =   a*std::cos(t);
   M(1,0) =   std::sin(t);
   M(1,1) =  (std::cos(t)*std::cos(i));
   M(1,2) =  (-1*std::cos(t)*std::sin(i));
   M(1,3) =   a*std::sin(t);
   M(2,0) =   0;
   M(2,1) =   std::sin(i);
   M(2,2) =   std::cos(i);
   M(2,3) =   d;
   M(3,0) = 0;
   M(3,1) = 0;
   M(3,2) = 0;
   M(3,3) = 1;
    
   *this = M; 
}

Eigen::Vector3d Denny::GetZ(void)
{
    Eigen::Vector3d Z;
    Z << 0, 0, 1;
    return this->block(0,0,3,3)*Z;
}

Eigen::Vector3d Denny::GetC(void)
{
    return this->block(0,3,3,1);
}

Eigen::Vector3d makeUnit(Eigen::Vector3d Input)
{
   double Mag = sqrt(std::abs(pow(Input(0),2) + pow(Input(1),2) + pow(Input(2),2)));
   return Input/Mag;
}

int main(int argc, char** argv)
{
    bool percentCartesianSetPoint = false, percentEulerSetPoint = false;
    char* pEnd;
    char icon[4][5] = {"*","**","***","****"};
    auto start  = std::chrono::high_resolution_clock::now();
    double Criteria = 0, CriteriaLast = 0, linDelta, angDelta, initDelta;   
    long Count = 0;
    Eigen::AngleAxisd angleAxis;
    Eigen::MatrixXd Jv(3,6), Jw(3,6), J(6,6), q(6,1), qdot(6,1) ,Xi(6,1);
    Eigen::Matrix3d Roti, Rotk, Rotd;
    Eigen::Vector3d Angi, Angk, Angd, DeltaNorm, Xiv, Xiw;
   
    Denny O0, O1, O2, O3, O4, O5, O6;
    Denny H1, H2, H3, H4, H5, H6;
    Denny Delta, Goal, End, Tool, THD; 

    // Set final position for each joint
    Goal =  Eigen::MatrixXd::Identity(4,4);
    O0 = Eigen::MatrixXd::Identity(4,4);
    Goal(0,0) = 1;
    Goal(1,1) = 1;
    Goal(2,2) = 1;

    // Set Cartesian origin
    Goal(0,3) =  strtod(argv[1],&pEnd);
    Goal(1,3) =  strtod(argv[2],&pEnd);
    Goal(2,3) =  strtod(argv[3],&pEnd);

    // Set inital conditions for joints
    q(0) = atan2(Goal(0,3),Goal(1,3));    
    q(1) = 45 * M_PI/180;       
    q(2) = 45 * M_PI/180;       
    q(3) = 45 * M_PI/180;       
    q(4) = 45 * M_PI/180;       
    q(5) = 45 * M_PI/180;        

    do 
    { 

	// Update the D-H matrixes
	H1.SetDH(q(0)       ,    0.1273,      0,   M_PI/2);
	H2.SetDH(q(1)-M_PI/2,         0, -0.612,           0);    
	H3.SetDH(q(2)       ,         0, -0.5723,          0);    
	H4.SetDH(q(3)-M_PI/2,  0.163941,       0,  M_PI/2);    
	H5.SetDH(q(4)       ,    0.1157,       0, -M_PI/2);    
	H6.SetDH(q(5)       ,    0.0922,       0,          0); 
	THD.SetDH(0         ,       0.1,       0,          0 );

	// Calculate the link positions         
	O1 = H1;
	O2 = O1 * H2;
	O3 = O2 * H3;
	O4 = O3 * H4;
	O5 = O4 * H5;
	O6 = O5 * H6;
	Tool = O6 * THD;

	// Update body velocity vector
	Delta = (Goal - Tool);

	Xiv(0) = Delta(0,3);
	Xiv(1) = Delta(1,3); 
	Xiv(2) = Delta(2,3); 

	Roti << Goal.block(0,0,3,3);
	// Angi = Roti.eulerAngles(2,1,2);

	Rotk << Tool.block(0,0,3,3);
	Angk = Rotk.eulerAngles(2,1,2);
        
        Xiw = Angi - Angk;
        
        Rotd = Roti.transpose() * Rotk;
        angleAxis = Rotd;
        Angd = Rotd.eulerAngles(2,1,2);

        if(percentCartesianSetPoint != true)
            Xiw = Xiw; 
	Xi << Xiv, Xiw;

	// Calc Displacement
	linDelta = sqrt(std::abs(
		    pow(Delta(0,3),2) +
		    pow(Delta(1,3),2) +
		    pow(Delta(2,3),2)));
/*
	angDelta = sqrt(std::abs(
		    pow(Angi(0)-Angk(0),2) +
		    pow(Angi(1)-Angk(1),2) +
		    pow(Angi(2)-Angk(2),2)));
*/
	angDelta = std::abs(angleAxis.angle());

// Remember initial deltas
	if(Count == 0)
        {
	    initDelta = angDelta + linDelta;
        }
	// Update the jacobian both Jv, Jw
	Jv.col(0) = O0.GetZ().cross(Tool.GetC() - O0.GetC());
	Jv.col(1) = O1.GetZ().cross(Tool.GetC() - O1.GetC());
	Jv.col(2) = O2.GetZ().cross(Tool.GetC() - O2.GetC());
	Jv.col(3) = O3.GetZ().cross(Tool.GetC() - O3.GetC());
	Jv.col(4) = O4.GetZ().cross(Tool.GetC() - O4.GetC());
	Jv.col(5) = O5.GetZ().cross(Tool.GetC() - O5.GetC());

	Jw.col(0) = O0.GetZ(); 
	Jw.col(1) = O1.GetZ();
	Jw.col(2) = O2.GetZ();
	Jw.col(3) = O3.GetZ();
	Jw.col(4) = O4.GetZ();
	Jw.col(5) = O5.GetZ();

	// The jacobian... 
	J.block(0,0,3,6) = Jv;
	J.block(3,0,3,6) = Jw;

        qdot = J.colPivHouseholderQr().solve(Xi);
        //qdot = J.transpose() * (Xi);
	q  += qdot*((linDelta + angDelta*percentCartesianSetPoint)/initDelta)*.001; 
        /* 
        for(int i = 0; i < 6; i++)
        {
          q(i) =  ((q(i) < -M_PI)? -M_PI : q(i));
          q(i) =  ((q(i) >  M_PI)?  M_PI : q(i));
        }
        */

        ///* 
	if(Count % 1000 == 1 ) 
	{   
	     std::cout<<"linDelta: "<<linDelta<<"\t angDelta: "<<angDelta<<std::endl;
	     std::cout<<" q: "<<q.transpose()*180/M_PI<<std::endl;
             std::cout<<"%: "<<((linDelta + angDelta) / (initDelta))*100 <<std::endl;
             std::cout<<"End Effector:\n"<< Tool <<std::endl;
	}       
       /*
	   std::cout<<"J(x): "<<q.transpose()*180/M_PI<<std::endl;
	   std::cout<<"qdot: "<<qdot.transpose()*180/M_PI<<std::endl;
	   std::cout<<"Xi  : "<<Xi.transpose()<<std::endl;
	*/
        if(linDelta < .25)
            percentCartesianSetPoint = true;
	Count++;
	// usleep(500000);
	//    system("clear");
	//    std::cout<<icon[Count%4];
    }
    while((linDelta > .001) || (angDelta > .05)); 
    
    std::cout<<"Criteria: "<<Criteria<<std::endl;
    std::cout<<"END of Solver..."<<std::endl; 
    q = q *180/M_PI;
    for(int i = 0; i < 6; i++)
    {
        q(i) =  (q(i) >  180.0)?    (fmod((float)q(i),180.0)) : q(i);
        q(i) =  (q(i) < -180.0)? -1*(fmod((float)q(i),180.0)) : q(i);
    } 
    std::cout<<"Threshold; "<<percentCartesianSetPoint<<std::endl;  
    std::cout<<"Joint angles (deg.):\n"<< q <<std::endl;
    std::cout<<"Current Orientation (.deg):\n"<< Angk.transpose() * 180/M_PI <<std::endl;
    std::cout<<"Goal Orientation (rad):\n"<< Angi.transpose() * 180/M_PI <<std::endl;
    std::cout<<"End Effector:\n"<< Tool <<std::endl;
    std::cout<<"Goal Position:\n"<< Goal <<std::endl;
    std::cout<<"Iteration count: "<<Count<<std::endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = finish - start;
    std::cout<<"Time in (s): "<<duration.count()<<std::endl;
    
    return 0;
}





