/* -------------------------------------------------------------------------
EduLB - The educational lattice Boltzmann solver

Copyright (C) 2013 Andreas Hantsch (edulb@gmx-topmail.de)
Version 0.4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------*/

#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

/* -------------------------------------------------------------------------
0. Definition of global variables and functions
--------------------------------------------------------------------------*/
int lx,ly;	// domain size in lattice nodes
int t_0=0,t_max;	// maximum time steps
double intdensity0=0.;
int thickness_pml = 0;
double sigma_max = 0.;
double delta_t;
int time_step;
int time_push;
string version="0.4";

// Definition of D2Q9 lattice
/*
6   2   5
  \ | /
3 - 0 - 1
  / | \
7   4   8
*/
const int ex[]={0, 1, 0,-1, 0, 1,-1,-1, 1}; // D2Q9
const int ey[]={0, 0, 1, 0,-1, 1, 1,-1,-1}; // D2Q9
const int D=2, Q=9; // D2Q9
const double rt[]={4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.}; // D2Q9
const double cs2=1./3.; // D2Q9
#include "mrt.hh"

string filename(int&);
bool read_parameter(double&, double&, double&,double&, int&, int&, int&, int&, string&);
bool init_vectors(vector<bool>&, vector<double>&, vector<double>&);
bool write_results(vector<double>&, vector<double>&, vector<double>&, vector<bool>&, int);
bool read_geometry(vector<bool>&, string);
bool initialisation(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&,vector<double>&, double, double&, double);
bool check_density(vector<double>&, vector<double>&, vector<double>&, int, double, int);
bool propagation(vector<double>&, vector<double>&);
bool pml_advection(vector<double>&, const vector<double>&);
bool boundary(double, vector<double>&, vector<double>&, vector<bool>&);
void boundary_stable_fluid(double, const vector<double>&, const vector<double>& , const vector<double>&, const vector<bool>&, vector<double>&, vector<double>& );
bool calc_macr_quantities(vector<double>&, vector<double>&, vector<double>&, vector<double>&,const vector<double>&, vector<bool>&);
bool collision_srt(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<bool>&, double);
bool update_pml(vector<double>&, bool);
void collision_pml(const vector<double>&,const vector<double>&, const vector<double>&,const vector<double>&, const vector<bool>&, vector<double>&, vector<double>& );


int main(void)
{
	/* -------------------------------------------------------------------------
	1. Definition of all other variables
	--------------------------------------------------------------------------*/
	double density0, ux0; // initial values
	double omegaf;
	int dt_write, dt_density;
	int t;
	string geomfile;

	/* -------------------------------------------------------------------------
	2. Preparations
	--------------------------------------------------------------------------*/
	cout  << "/*------------------------------------------------------------" << endl;
	cout  << "EduLB 2D - The educational 2D lattice Boltzmann solver" << endl;
	cout  << "Copyright (C) 2013 Andreas Hantsch (edulb@gmx-topmail.de)" << endl;
	cout  << "Version " << version << endl << endl;
	cout  << "This program is distributed in the hope that it will be useful," << endl;
	cout  << "but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
	cout  << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
	cout  << "GNU General Public License for more details." << endl;
	cout  << "--------------------------------------------------------------*/" << endl << endl;
	read_parameter(density0, ux0, omegaf,sigma_max, dt_write, dt_density,thickness_pml,time_push, geomfile);
	
	vector<bool> obst(lx*ly); // obstacle
	vector<double> f(Q*lx*ly); // node
	vector<double> ftemp(Q*lx*ly); // temp. node
	vector<double> density(lx*ly); // mass density
	vector<double> ux(lx*ly); // x-velocity
	vector<double> uy(lx*ly); // y-velocity

	vector<double> pml(lx*ly); //pml damping value. usually called sigma
	vector<double> f_eq_old(Q*lx*ly);//save the old equilibrium state.
	
	init_vectors(obst, f, ftemp);
	read_geometry(obst,geomfile);
	initialisation(density, ux, uy, f, ftemp,f_eq_old, density0, intdensity0, ux0);
	update_pml(pml, false);
	write_results(density, ux, uy, obst, t_0);

	
	
	/* -------------------------------------------------------------------------
	3. Numerics
	--------------------------------------------------------------------------*/
	cout << "starting numerics..." << endl;

	for (t=t_0; t<t_max+1; ++t)
	{
		time_step +=1;
		//set the pml also to the left side.
		if(time_step==950){
			//update_pml(pml,true);
		}

		
		propagation(f, ftemp);  // LHS of Boltzmann equation
		//boundary(ux0, f, ftemp, obst);	// boundary conditions
		boundary_stable_fluid(ux0,density,ux,uy,obst,f,ftemp);
		//pml_advection(ftemp,pml);
		calc_macr_quantities(density, ux, uy, ftemp,pml, obst); // calculation of macroscopic quantities
		collision_srt(density, ux, uy, f, ftemp, obst, omegaf); // RHS of Boltzmann equation: SRT
		//collision_mrt(f, ftemp, obst, omegaf); // RHS of Boltzmann equation: SRT
		//collision_pml(density, ux, uy, pml, obst, f_eq_old, f);


		// io stuff
		if (t%dt_density==0)
		{
			check_density(density, ux, uy, t, intdensity0, dt_density);
		}
		if ((t>0) && (t%dt_write==0))
		{
			write_results(density, ux, uy, obst, t);
		}
	}

	/* -------------------------------------------------------------------------
	4. Some i/o
	--------------------------------------------------------------------------*/
	if ((t-1)%dt_write!=0){write_results(density, ux, uy, obst, t-1);}
	cout << "done.\n" << endl;
	return(0);
}


/* -------------------------------------------------------------------------
5. All the functions
--------------------------------------------------------------------------*/

/* -------------------------------------------------------------------------
function that creates file names without extension
--------------------------------------------------------------------------*/
string filename(int& t)
{
	int j=1;
	string space = "";
	for (j=1; j<t_max+1; j=j*10)
	{
		if (t<j)
		{
			space+="0";
		}
	}
	if (t!=t_max && t!=0 && t==10*j){space+="0";}
	if (t==0){space.erase(space.length()-1);}
	stringstream timestrings;
	timestrings << space << t;
	string timestring=timestrings.str();
	return (timestring);
}

/* -------------------------------------------------------------------------
function that reads a parameter file
--------------------------------------------------------------------------*/
bool read_parameter(double& density0, double& ux0, double& omegaf,double& sigma_max, int& dt_write, int& dt_density,int& thickness_pml,int& time_push, string& geomfile)
{
	string ab;
	double reynolds;
	double kin_visc_lb;
	char geomfilec[200];
	ifstream controlDict;
	controlDict.open("./system/controlDict", ios::in);
	cout<<"reading controlDict...";
	while(!controlDict.eof())
	{
		controlDict >> ab;
		if (ab == "TMAX")			{controlDict >> t_max;}
		if (ab == "DENSIT0")	{controlDict >> density0;}
		if (ab == "KINVISC")	{controlDict >> kin_visc_lb;}
		if (ab == "RE")				{controlDict >> reynolds;}
		if (ab == "DT_WRIT")	{controlDict >> dt_write;}
		if (ab == "DT_DENS")	{controlDict >> dt_density;}
		if (ab == "TH_PML")	{controlDict >> thickness_pml;}
		if (ab == "SIGMA_MAX")	{controlDict >> sigma_max;}
		if (ab == "T_PUSH")	{controlDict >> time_push;}
		if (ab == "GEOFILE")	{controlDict >> geomfile;}
	}
	controlDict.close();
	
	strcpy(geomfilec,geomfile.c_str()); // convert string to array of characters
	
	ifstream geometry;
	geometry.open(geomfilec, ios::in);
	geometry >> ab;
	if (ab == "LX")	{geometry >> lx	;}
	geometry >> ab;
	if (ab == "LY")	{geometry >> ly	;}
	geometry.close();
	
	cout << "calculation of further quantities...";
	//ux0=reynolds*kin_visc_lb/double(ly-1);
	ux0=reynolds*kin_visc_lb/double(600-1);
	omegaf=1./(3.*kin_visc_lb+0.5);

	double delta_x=1./double(ly-1);
	delta_t=ux0*delta_x;
	double mach=ux0/sqrt(cs2); // Mach number
	cout << "done." << endl;

	cout  << "Output of all read and calculated quantities:\n"
				<< "LX      " << lx 					<< endl
				<< "LY      " << ly 					<< endl
				<< "TMAX    " << t_max 				<< endl
				<< "DENSIT0 " << density0 		<< endl
				<< "KINVISC " << kin_visc_lb	<< endl
				<< "RE      " << reynolds			<< endl
				<< "DT_WRIT " << dt_write 		<< endl
				<< "DT_DENS " << dt_density 	<< endl
				<< "GEOFILE " << geomfile			<< endl
				<< endl
				<< "domain  " << lx << " x " << ly << " = " << lx*ly << endl
				<< "ux0     " << ux0					<< endl
				<< "omegaf  " << omegaf				<< endl
				<< "delta_x " << delta_x			<< endl
				<< "delta_t " << delta_t			<< endl
				<< "mach    " << mach					<< endl;
	if (mach>0.1){cout<<"\nWarning: Mach number = " << mach << ". Beware of compressibility error: O(Mach^2) = " << mach*mach << ".\nYou propably like to increase the resolution. \n\n";}
	cout 	<< endl;
	return(0);
}

/* -------------------------------------------------------------------------
function that writes results into a file
--------------------------------------------------------------------------*/
bool write_results(vector<double>& density, vector<double>& ux, vector<double>&uy, vector<bool>& obst, int t)
{
	int x,y;
	int pos;
	char fileresults[100];
	FILE *fp;
	string file_m="../result_folder/"+filename(t)+"_m.ssv";
	strcpy(fileresults,file_m.c_str());
	fp=fopen(fileresults, "w+");
	fprintf(fp, "x\ty\tux\tuy\tpress\trho\tobsval\n");
	bool save_compressed_for_reference = false;
	

	if(save_compressed_for_reference){
		for (y=0; y<ly; ++y)
		{
			for (x=0; x<lx; ++x)
			{
				pos=x+lx*y;
				if((x<600)&&((y>699)&&(y<1300))){
					if(x==33 && y==1000){
						std::cout<<" The values are: "<< double(x)<<" "<< double(y)<<" "<< ux[pos]<<" "<< uy[pos]<<std::endl;
					}
					
					fprintf(fp, "%.6e\t%.6e\t%.3e\t%.3e\t%.3e\t%.3e\t%d\n", double(x), double(y-700), ux[pos], uy[pos], density[pos]*cs2, density[pos], int(obst[pos]));
				}	
			}
			if((y>699)&&(y<1300)){
				fprintf(fp, "\n");
			}
		}
	}else{
		for (y=0; y<ly; ++y)
		{
			for (x=0; x<lx; ++x)
			{
				pos=x+lx*y;	
				fprintf(fp, "%.6e\t%.6e\t%.3e\t%.3e\t%.3e\t%.3e\t%d\n", double(x), double(y), ux[pos], uy[pos], density[pos]*cs2, density[pos], int(obst[pos]));
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	return(0);
}

/* -------------------------------------------------------------------------
function for initialising the particle density distribution functions
--------------------------------------------------------------------------*/
bool init_vectors(vector<bool>& obst, vector<double>& f, vector<double>& ftemp)
{
	cout << "initialise vectors...";
	int x,y,i;
	int pos;
	for (y=0; y<ly; ++y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;
			obst[pos]=false;
			for (i=0; i<Q; ++i)
			{
				f[Q*pos+i] = 0.;
				ftemp[Q*pos+i] = 0.;
			}
		}
	}
	cout << "done.\n";
	return(0);
}

/* -------------------------------------------------------------------------
function that reads geometry file neglecting domain boundaries
--------------------------------------------------------------------------*/
bool read_geometry(vector<bool>& obst, string geomfile)
{

	cout << "read geometry file...";
	int x,y;
	int pos;
	int count=0;
	double data;
	char geomfilec[200];
	string datadump;
	
	strcpy(geomfilec,geomfile.c_str()); // convert string to array of characters
	ifstream geometry;
	geometry.open(geomfilec, ios::in);
	geometry >> datadump >> datadump;
	geometry >> datadump >> datadump;
	geometry >> x;
	geometry >> y;
	geometry >> data;
	while(!geometry.eof())
	{
		if ((x<lx-1) && (y<ly-1) && (x>0) && (y>0) && data>0.5)
		{
			pos=x+lx*y;
			obst[pos]=true;
			count++;
		}
	geometry >> x;
	geometry >> y;
	geometry >> data;
	}
	geometry.close();
	cout << "done.";
	if (count>0){cout << " "	<<	count << " solid nodes read."	<< endl;}else{cout<<endl;}
	return(0);
}

/* -------------------------------------------------------------------------
function that initialises the density and velocity field
--------------------------------------------------------------------------*/
bool initialisation(vector<double>& density, vector<double>& ux, vector<double>& uy, vector<double>& f, vector<double>& ftemp,vector<double>& f_eq_old, double density0, double& intdensity0, double ux0)
{
	cout << "initialise density and velocity field, particle density distribution functions...";
	int x,y,i;
	int pos;
	double u2;
	for (y=0; y<ly; ++y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;
			density[pos]	=	density0;
 			ux[pos]				=	0.;
			//ux[pos]				=	ux0;
// 		ux[pos]				=	ux0*1.5*(4.0*(y-0.5)/(ly-2) - (2.0*(y-0.5)/(ly-2))*(2.0*(y-0.5)/(ly-2))); // parabolic profile
// 		ux[pos]				=	ux0*1.5*(4.0*y/(ly-1) - (2.0*y/(ly-1))*(2.0*y/(ly-1))); // parabolic profile
// 		ux[pos]				=	ux0*3./2.*(2.*y/(ly-1)-(1.*y/(ly-1))*(1.*y/(ly-1))); // Nusselt's velocity profile
			uy[pos]				=	0.;
			u2=ux[pos]*ux[pos] + uy[pos]*uy[pos];

			for (i=0;i<Q;++i)
			{
				f[Q*pos+i]    =density[pos]*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));;
				ftemp[Q*pos+i]=density[pos]*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));;
				f_eq_old[Q*pos+i]=density[pos]*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			}
			intdensity0+=density[pos];
		}
	}
	cout << "done.\n";
	return(0);
}

/* -------------------------------------------------------------------------
function that performs the LB propagation step (lhs of Boltzmann equation)
--------------------------------------------------------------------------*/
bool propagation(vector<double>& f, vector<double>& ftemp)
{
	int x,y;
	int x_e,x_w,y_n,y_s;
	int pos;
	for (y=0; y<ly; ++y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;

			y==ly-1 ? y_n=-1 : y_n=y+1; // avoiding periodic bc
			x==lx-1 ? x_e=-1 : x_e=x+1; // avoiding periodic bc
			y==0    ? y_s=-1 : y_s=y-1; // avoiding periodic bc
			x==0    ? x_w=-1 : x_w=x-1; // avoiding periodic bc

			ftemp[Q*(x    + y   *lx) + 0] = f[Q*pos];
			if(x_e!=-1){ftemp[Q*(x_e  + y   *lx) + 1] = f[Q*pos+1];}
			if(y_n!=-1){ftemp[Q*(x    + y_n *lx) + 2] = f[Q*pos+2];}
			if(x_w!=-1){ftemp[Q*(x_w  + y   *lx) + 3] = f[Q*pos+3];}
			if(y_s!=-1){ftemp[Q*(x    + y_s *lx) + 4] = f[Q*pos+4];}
			if(x_e!=-1 && y_n!=-1)	{ftemp[Q*(x_e  + y_n *lx) + 5] = f[Q*pos+5];}
			if(x_w!=-1 && y_n!=-1)	{ftemp[Q*(x_w  + y_n *lx) + 6] = f[Q*pos+6];}
			if(x_w!=-1 && y_s!=-1)	{ftemp[Q*(x_w  + y_s *lx) + 7] = f[Q*pos+7];}
			if(x_e!=-1 && y_s!=-1)	{ftemp[Q*(x_e  + y_s *lx) + 8] = f[Q*pos+8];}
		}
	}
	return(0);
}

/* -------------------------------------------------------------------------
function that performs boundary interactions
--------------------------------------------------------------------------*/
bool boundary(double ux0, vector<double>& f, vector<double>& ftemp, vector<bool>& obst)
{
	int x,y,i;
	int pos;
	double vf, ru, density_loc;
	//double aa=1., ab=0.; // only first order accuracy
	//double aa=2., ab=-1.; // second order accuracy

// 	Momentum at left end of domain (inlets)
	x=0; // velocity inlet at top of the domain, Zou and He
	for (y=1; y<ly-1; ++y)
	{
		pos=x+lx*y;
		//vf=ux0;
		vf = 0;
// 		vf=ux0*1.5*(4.0*(y-0.5)/(ly-2) - (2.0*(y-0.5)/(ly-2))*(2.0*(y-0.5)/(ly-2))); // parabolic profile
// 		vf=ux0*1.5*(4.0*y/(ly-1) - (2.0*y/(ly-1))*(2.0*y/(ly-1))); // parabolic profile
// 		vf=ux0*3./2.*(2.*y/(ly-1)-(1.*y/(ly-1))*(1.*y/(ly-1))); // Nusselt's velocity profile
		if(time_step<time_push){
			//here specify the size of the outlet, resp the blast.
			//if((y>5.5*thickness_pml) && (y<lx-5.5*thickness_pml)){ //good for pml 50
			if((y>975) && (y<1025)){ //inlet for the simulation on the reference solution
				
				vf = ux0;
				ru=(ftemp[Q*pos+0]+ftemp[Q*pos+2]+ftemp[Q*pos+4]+2.*(ftemp[Q*pos+3]+ftemp[Q*pos+6]+ftemp[Q*pos+7]))/(1.-vf)*vf;
				ftemp[Q*pos+1]=ftemp[Q*pos+3]+2./3.*ru;
				ftemp[Q*pos+5]=ftemp[Q*pos+7]+1./6.*ru+0.5*(ftemp[Q*pos+4]-ftemp[Q*pos+2]);
				ftemp[Q*pos+8]=ftemp[Q*pos+6]+1./6.*ru+0.5*(ftemp[Q*pos+2]-ftemp[Q*pos+4]);
			}else{
				ftemp[Q*pos + 5] = ftemp[Q*pos + 6];
				ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
				ftemp[Q*pos + 8] = ftemp[Q*pos + 7];
			}
		}else{
			ftemp[Q*pos + 5] = ftemp[Q*pos + 6];
			ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
			ftemp[Q*pos + 8] = ftemp[Q*pos + 7];
		}
		
	}

// 	Momentum at right end of the domain (outlet)
	x=lx-1;
	for (y=0; y<ly; ++y)
	{
		pos=x+lx*y;
		/*
		ftemp[Q*pos+0]=aa*ftemp[Q*pos-Q+0]+ab*ftemp[Q*pos-2*Q+0]; // takes values from x=lx-2
		ftemp[Q*pos+1]=aa*ftemp[Q*pos-Q+1]+ab*ftemp[Q*pos-2*Q+1]; // takes values from x=lx-2
		ftemp[Q*pos+2]=aa*ftemp[Q*pos-Q+2]+ab*ftemp[Q*pos-2*Q+2]; // takes values from x=lx-2
		ftemp[Q*pos+3]=aa*ftemp[Q*pos-Q+3]+ab*ftemp[Q*pos-2*Q+3]; // takes values from x=lx-2
		ftemp[Q*pos+4]=aa*ftemp[Q*pos-Q+4]+ab*ftemp[Q*pos-2*Q+4]; // takes values from x=lx-2
		ftemp[Q*pos+5]=aa*ftemp[Q*pos-Q+5]+ab*ftemp[Q*pos-2*Q+5]; // takes values from x=lx-2
		ftemp[Q*pos+6]=aa*ftemp[Q*pos-Q+6]+ab*ftemp[Q*pos-2*Q+6]; // takes values from x=lx-2
		ftemp[Q*pos+7]=aa*ftemp[Q*pos-Q+7]+ab*ftemp[Q*pos-2*Q+7]; // takes values from x=lx-2
		ftemp[Q*pos+8]=aa*ftemp[Q*pos-Q+8]+ab*ftemp[Q*pos-2*Q+8]; // takes values from x=lx-2
		*/
		ftemp[Q*pos + 6] = ftemp[Q*pos + 8];
		ftemp[Q*pos + 3] = ftemp[Q*pos + 1];
		ftemp[Q*pos + 7] = ftemp[Q*pos + 5];
	
	}

	// Momentum bounce back at top wall
	y=ly-1;
	for(x=1; x<lx-1; ++x)
	{
		pos=x+y*lx;
		ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
		ftemp[Q*pos + 7] = ftemp[Q*pos + 5];
		ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
	}

	// Momentum bounce back at bottom wall
	y=0;
	for(x=1;x<lx-1;++x)
	{
		pos=x+y*lx;
		ftemp[Q*pos + 2] = ftemp[Q*pos + 4];
		ftemp[Q*pos + 5] = ftemp[Q*pos + 7];
		ftemp[Q*pos + 6] = ftemp[Q*pos + 8];
	}

	// south-west corner of the inlet has to be defined
	pos = lx; // = 0+lx*1
	density_loc = 0.0;
	for(i=0; i<Q; ++i)
	{
		density_loc += ftemp[Q*pos + i];
	}
	pos = 0;
	ftemp[Q*pos + 2] = ftemp[Q*pos + 4];
	ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
	ftemp[Q*pos + 5] = ftemp[Q*pos + 7];
	ftemp[Q*pos + 6] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 1] + ftemp[Q*pos + 5]));
	ftemp[Q*pos + 8] = ftemp[Q*pos + 6];

	// north-west corner of the inlet has to be defined
	pos = lx*(ly-2); // = 0+lx*(ly-2)
	density_loc = 0.0;
	for(i=0; i<Q; ++i)
	{
		density_loc += ftemp[Q*pos + i];
	}
	pos = lx*(ly-1);
	ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
	ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
	ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
	ftemp[Q*pos + 7] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 3] + ftemp[Q*pos + 6]));
	ftemp[Q*pos + 5] = ftemp[Q*pos + 7];

// 	south-east corner of the outlet has to be defined
	pos = (lx-1)+lx; // = (lx-1)+lx*1
	density_loc = 0.0;
	for(i=0; i<Q; ++i)
	{
		density_loc += ftemp[Q*pos + i];
	}
	pos = lx-1;// = (lx-1)+lx*0
	ftemp[Q*pos + 2] = ftemp[Q*pos + 4];
	ftemp[Q*pos + 3] = ftemp[Q*pos + 1];
	ftemp[Q*pos + 6] = ftemp[Q*pos + 8];
	ftemp[Q*pos + 5] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 3] + ftemp[Q*pos + 6]));
	ftemp[Q*pos + 7] = ftemp[Q*pos + 5];

	// north-east corner of the inlet has to be defined
	pos = (lx-1)+lx*(ly-2); // = (lx-1)+lx*(ly-2)
	density_loc = 0.0;
	for(i=0; i<Q; ++i)
	{
		density_loc += ftemp[Q*pos + i];
	}
	pos = (lx-1)+lx*(ly-1);
	ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
	ftemp[Q*pos + 3] = ftemp[Q*pos + 1];
	ftemp[Q*pos + 7] = ftemp[Q*pos + 5];
	ftemp[Q*pos + 8] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 4] + ftemp[Q*pos + 3] + ftemp[Q*pos + 7]));
	ftemp[Q*pos + 6] = ftemp[Q*pos + 8];


	// default bounce back at all inner wall nodes
	for (y=1; y<ly-1; ++y)
	{
		for (x=1; x<lx-1; ++x)
		{
			pos=x+lx*y;
			if (obst[pos]) // bounce-back at all inner obstacle nodes
			{
				f[Q*pos+1] = ftemp[Q*pos+3];
				f[Q*pos+2] = ftemp[Q*pos+4];
				f[Q*pos+3] = ftemp[Q*pos+1];
				f[Q*pos+4] = ftemp[Q*pos+2];
				f[Q*pos+5] = ftemp[Q*pos+7];
				f[Q*pos+6] = ftemp[Q*pos+8];
				f[Q*pos+7] = ftemp[Q*pos+5];
				f[Q*pos+8] = ftemp[Q*pos+6];
			}
		}
	}
	return(0);
}

/* -------------------------------------------------------------------------
function that calculates the macroscopic quantities
--------------------------------------------------------------------------*/
bool calc_macr_quantities(vector<double>& density, vector<double>& ux, vector<double>& uy, vector<double>& ftemp, const vector<double>& pml, vector<bool>& obst)
{
	int x,y,pos,i;
	double density_loc,ux_loc,uy_loc;
	for (y=0; y<ly; ++y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;
			if (!obst[pos])
			{
				density_loc=0.;
				ux_loc=0.;
				uy_loc=0.;
				for (i=0; i<Q; ++i)
				{
					density_loc+=ftemp[Q*pos+i];
					if(i>0){
						ux_loc+=ex[i]*ftemp[Q*pos+i];
						uy_loc+=ey[i]*ftemp[Q*pos+i];
					}
				}
				density[pos]=density_loc;
				ux[pos]=ux_loc/density_loc;
				uy[pos]=uy_loc/density_loc;
			}
			else
			{
				density[pos]=0.;
				ux[pos]=0.;
				uy[pos]=0.;
			}
		}
	}
	return(0);
}

/* -------------------------------------------------------------------------
function that performs the LB collision step (rhs of Boltzmann equation)
--------------------------------------------------------------------------*/
bool collision_srt(vector<double>& density, vector<double>& ux, vector<double>& uy, vector<double>& f, vector<double>& ftemp, vector<bool>& obst, double omegaf)
{
	int x,y,pos,i;
	double density_loc, ux_loc, uy_loc, u2;
	double check_f, check_ftemp;
	double feq[Q];

	for (y=0; y<ly; ++y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;
			if (!obst[pos])
			{
				density_loc=density[pos];
				ux_loc=ux[pos];
				uy_loc=uy[pos];
				u2=ux_loc*ux_loc+uy_loc*uy_loc; // square of velocity

				check_f=0.;
				check_ftemp=0.;
				for (i=0; i<Q; ++i)
				{
					// calculating equilibrium distribution, e. (2)
					feq[i]=density_loc*rt[i]*(1. + (ex[i]*ux_loc+ey[i]*uy_loc)/cs2 + (ex[i]*ux_loc+ey[i]*uy_loc)*(ex[i]*ux_loc+ey[i]*uy_loc)/(2.*cs2*cs2) - u2/(2.*cs2));

					// solving rhs of Boltzmann equation
					f[Q*pos+i]=ftemp[Q*pos+i]+omegaf*(feq[i]-ftemp[Q*pos+i]);

					// summ up density distribution functions to check for negative densities
					check_f+=f[Q*pos+i];
					check_ftemp+=ftemp[Q*pos+i];
				}
				if (check_f < 0 || check_ftemp < 0){cout << "fatal error: density < 0 at (" << x << "," << y << ")"<< endl;exit(1);}
			}
		}
	}
	return(0);
}


/*-------------------------------------------------------------------------
function that performes the LBM PML collision term.
-------------------------------------------------------------------------*/
void collision_pml(const vector<double>& density,const vector<double>& ux,const vector<double>& uy,const vector<double>& pml, const vector<bool>& obst, vector<double>& f_eq_old, vector<double>& f ){

	int pos, x, y;
	double density_loc, ux_loc, uy_loc, u2;
	double omega_pml;
	
	vector<double> f_eq(Q*lx*ly);
	//compute the equilibrium distribution
	for(x= 0; x<lx; x++){
		for(y = 0; y<lx; y++){
			pos = lx*y + x;
			//only when we are away from obstacle
			if (!obst[pos])
			{
				density_loc=density[pos];
				ux_loc=ux[pos];
				uy_loc=uy[pos];
				u2=ux_loc*ux_loc+uy_loc*uy_loc; // square of velocity
				for (int i=0; i<Q; ++i)
				{
					// calculating equilibrium distribution, e. (2)
					f_eq[Q*pos + i]=density_loc*rt[i]*(1. + (ex[i]*ux_loc+ey[i]*uy_loc)/cs2 + (ex[i]*ux_loc+ey[i]*uy_loc)*(ex[i]*ux_loc+ey[i]*uy_loc)/(2.*cs2*cs2) - u2/(2.*cs2));
				}
			}
		}
	}

	for( x = 0; x<lx; x++){
		for( y = 0; y<ly; y++){
			pos = lx*y + x;
			if(!obst[pos]){
				//check if we are in the pml area
				if(pml[pos]>0){
					omega_pml = 0.;
					for(int i=0; i<Q;i++){
						
						//first handle the non boundary case. 
						if(((x>0)&&(y>0))&&((x<lx-1)&&(y<ly-1))){
							//add the initial equilibrium distribution. Check whether we need to subtract the starting mean
							omega_pml+=2*f_eq[Q*pos + i];

							//add the quantity for the  Q_i, i.e the integrated
							omega_pml += 0.5*(f_eq[Q*pos + i] + f_eq_old[Q*pos + i]);

							//add the quantity for the gradient of the Q_i
							//i==0 not interesing, as c0 = (0,0). Compute now $c_i\cdot \nabla Q_i $
							//note that in our case dx = 1. Hence for the gradient
							
							//compute the auxiliary quantity for all neighbors. Thus define center, north, south, east, west.
							double psi_n, psi_w, psi_e,psi_s;
							//psi_c = 0.5*(f_eq[Q*pos + i] + f_eq_old[Q*pos + i]);
							psi_e = 0.5*(f_eq[Q*(lx*y+x+1) + i] + f_eq_old[Q*(lx*y+x+1) + i]);
							psi_w = 0.5*(f_eq[Q*(lx*y+x-1) + i] + f_eq_old[Q*(lx*y+x-1) + i]);
							psi_n = 0.5*(f_eq[Q*(lx*(y+1)+x) + i] + f_eq_old[Q*(lx*(y+1)+x) + i]);
							psi_s = 0.5*(f_eq[Q*(lx*(y-1)+x) + i] + f_eq_old[Q*(lx*(y-1)+x) + i]);
							//compute the partial derivatives with a central difference scheme
							double dpsi_x, dpsi_y;
							dpsi_x = 0.5*(psi_e - psi_w);
							dpsi_y = 0.5*(psi_n - psi_s);
							//perform the central difference approximation of the gradient of psi
							{
								if(i==1){
									//c1 = (1,0)
									omega_pml+= dpsi_x;
								}
								if(i==2){
									//c2 = (0,1)								
									omega_pml+=dpsi_y;
								}
								if(i==3){
									//c3 = (-1,0)
									omega_pml -=dpsi_x;
								}
								if(i==4){
									//c4 = (0,-1)
									omega_pml-=dpsi_y;
								}
								if(i==5){
									//c5 = (1,1)
									omega_pml+=dpsi_x + dpsi_y;
								}
								if(i==6){
									//c6 =(-1,1)
									omega_pml +=dpsi_y - dpsi_x;
								}
								if(i==7){
									//c7 = (-1,-1)
									omega_pml += -dpsi_y -dpsi_x;
								}
								if(i==8){
									//c8 = (1,-1)
									omega_pml += dpsi_x - dpsi_y;
								}
							}
							//multiplicatin of the damping parameter
							omega_pml*=-pml[pos];

							f[Q*pos+i]+=omega_pml;
							
						}
					}
				}
			}
		}
	}
	//overwrite the old equilibrium state with the newer one
	f_eq_old = f_eq;
	
}

/*-------------------------------------------------------------------------
function that performes the stable fluid boundary treatment
-------------------------------------------------------------------------*/
void boundary_stable_fluid(double ux0, const vector<double>& density, const vector<double>& ux, const vector<double>& uy, const vector<bool>& obst, vector<double>& f, vector<double>& ftemp ){
	int x,y,i;
	int pos;
	double vf, ru, density_loc;

	//first consider left inlet
	x=0; // velocity inlet at top of the domain, Zou and He
	for (y=1; y<ly-1; ++y)
	{
		pos=x+lx*y;
		//vf=ux0;
		vf = 0;
		// vf=ux0*1.5*(4.0*(y-0.5)/(ly-2) - (2.0*(y-0.5)/(ly-2))*(2.0*(y-0.5)/(ly-2))); // parabolic profile
		// vf=ux0*1.5*(4.0*y/(ly-1) - (2.0*y/(ly-1))*(2.0*y/(ly-1))); // parabolic profile
		// vf=ux0*3./2.*(2.*y/(ly-1)-(1.*y/(ly-1))*(1.*y/(ly-1))); // Nusselt's velocity profile
		if(time_step<time_push){
			//here specify the size of the outlet, resp the blast.
			//if((y>5.5*thickness_pml) && (y<lx-5.5*thickness_pml)){ //good for pml 50 and size 600
			//if((y>975) && (y<1025)){
			if(true){
				
				vf = ux0;
				ru=(ftemp[Q*pos+0]+ftemp[Q*pos+2]+ftemp[Q*pos+4]+2.*(ftemp[Q*pos+3]+ftemp[Q*pos+6]+ftemp[Q*pos+7]))/(1.-vf)*vf;
				ftemp[Q*pos+1]=ftemp[Q*pos+3]+2./3.*ru;
				ftemp[Q*pos+5]=ftemp[Q*pos+7]+1./6.*ru+0.5*(ftemp[Q*pos+4]-ftemp[Q*pos+2]);
				ftemp[Q*pos+8]=ftemp[Q*pos+6]+1./6.*ru+0.5*(ftemp[Q*pos+2]-ftemp[Q*pos+4]);
			}else{
				ftemp[Q*pos + 5] = ftemp[Q*pos + 7];
				ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
				ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
			}
		}else{
			//here use the stable fluid
			ftemp[Q*pos + 5] = ftemp[Q*pos + 7];
			ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
			ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
		}
		
	}

	// Stable fluid at top wall
	y=ly-1;
	for(x=1; x<lx-1; ++x)
	{
		pos=x+y*lx;
		double ux_loc, uy_loc;
		double x_back, y_back, ux_back, uy_back, rho_back, u2_back;
		int x_left, x_right, y_up, y_down;
		double ux_down, ux_up, uy_down, uy_up, rho_down, rho_up;
		ux_loc = ux[pos];
		uy_loc = uy[pos];
		//Euler Step for backtracking
		x_back = (double(x) - delta_t*ux_loc);
		y_back = double(y) - delta_t*uy_loc;
		//check if we are in the physical domain
		if(x==90){
			
			/*std::cout<< "delta t "<<delta_t<<std::endl;
			std::cout<< "ux_loc "<<ux_loc<<std::endl;
			std::cout<< " x-back: "<<x_back<<std::endl;
			std::cout<< " y-back: "<<y_back<<std::endl;
			*/
		}

		if(((x_back >= 0.)&&(x_back<double(lx-1)))&&((y_back>=0.)&&(y_back<double(ly-1)))){
			//here hidden error
			x_left = trunc(x_back);
			x_right = x_left+1;
			y_down = trunc(y_back);
			y_up = y_down+1;
			
			//bilinear interpolation of the quantities
			ux_down = (double(x_right)-x_back)*ux[x_right + y_down*ly] + (x_back - double(x_left))*ux[x_left + y_down*ly ];
			ux_up = (double(x_right)-x_back)*ux[x_right + y_up*ly] + (x_back - double(x_left))*ux[x_left + y_up*ly ];
			ux_back = (double(y_up) - y_back)*ux_up + (y_back - double(y_down))*ux_down;
			
			uy_down = (double(x_right)-x_back)*uy[x_right + y_down*ly] + (x_back - double(x_left))*uy[x_left + y_down*ly ];
			uy_up = (double(x_right)-x_back)*uy[x_right + y_up*ly] + (x_back - double(x_left))*uy[x_left + y_up*ly ];
			uy_back = (double(y_up) - y_back)*uy_up + (y_back - double(y_down))*uy_down;

			rho_down = (double(x_right)-x_back)*density[x_right + y_down*ly] + (x_back - double(x_left))*density[x_left + y_down*ly ];
			rho_up = (double(x_right)-x_back)*density[x_right + y_up*ly] + (x_back - double(x_left))*density[x_left + y_up*ly ];
			rho_back = (double(y_up) - y_back)*rho_up + (y_back - double(y_down))*rho_down;

			u2_back = ux_back*ux_back + uy_back*uy_back;
			//use these quantities for the equilibrium distribution
			//rho_back*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 8] = rho_back*rt[8]*(1. + (ex[8]*ux_back+ey[8]*uy_back)/cs2 + (ex[8]*ux_back+ey[8]*uy_back)*(ex[8]*ux_back+ey[8]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 4] = rho_back*rt[4]*(1. + (ex[4]*ux_back+ey[4]*uy_back)/cs2 + (ex[4]*ux_back+ey[4]*uy_back)*(ex[4]*ux_back+ey[4]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 7] = rho_back*rt[7]*(1. + (ex[7]*ux_back+ey[7]*uy_back)/cs2 + (ex[7]*ux_back+ey[7]*uy_back)*(ex[7]*ux_back+ey[7]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));

			//std::cout<<"time: "<<time_step<<" We are inside the domain"<<std::endl;
		}else{//case when we leave the domain. In this case just take the former value.
			/*
			std::cout<< "delta t "<<delta_t<<std::endl;
			std::cout<< "ux_loc "<<ux_loc<<std::endl;
			std::cout<< " x-back: "<<x_back<<std::endl;
			std::cout<< " y-back: "<<y_back<<std::endl;
			std::cout<<"time: "<<time_step<<" We are outside the domain"<<std::endl;
			*/
			double u2=ux[pos]*ux[pos]+uy[pos]*uy[pos];
			ftemp[Q*pos + 7] = density[pos]*rt[7]*(1. + (ex[7]*ux[pos]+ey[7]*uy[pos])/cs2 + (ex[7]*ux[pos]+ey[7]*uy[pos])*(ex[7]*ux[pos]+ey[7]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 4] = density[pos]*rt[4]*(1. + (ex[4]*ux[pos]+ey[4]*uy[pos])/cs2 + (ex[4]*ux[pos]+ey[4]*uy[pos])*(ex[4]*ux[pos]+ey[4]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 8] = density[pos]*rt[8]*(1. + (ex[8]*ux[pos]+ey[8]*uy[pos])/cs2 + (ex[8]*ux[pos]+ey[8]*uy[pos])*(ex[8]*ux[pos]+ey[8]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
		}
	
		//ftemp[Q*pos + 7] = ftemp[Q*pos + 5];
		//ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
		//ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
	}

	// Stable fluid at bottom wall
	y=0;
	for(x=1;x<lx-1;++x)
	{
		pos=x+y*lx;
		double ux_loc, uy_loc;
		double x_back, y_back, ux_back, uy_back, rho_back, u2_back;
		int x_left, x_right, y_up, y_down;
		double ux_down, ux_up, uy_down, uy_up, rho_down, rho_up;
		ux_loc = ux[pos];
		uy_loc = uy[pos];
		//Euler Step for backtracking
		x_back = double(x) - delta_t*ux_loc;
		y_back = double(y) - delta_t*uy_loc;
		//check if we are in the physical domain
		/*
		std::cout<< "delta t "<<delta_t<<std::endl;
		std::cout<< "ux_loc "<<ux_loc<<std::endl;
		std::cout<< " x-back: "<<x_back<<std::endl;
		std::cout<< " y-back: "<<y_back<<std::endl;
		*/

		if(((x_back >= 0.)&&(x_back<=double(lx-1)))&&((y_back>=0.)&&(y_back<=double(ly-1)))){
			x_left = trunc(x_back);
			x_right = x_left+1;
			y_down = trunc(y_back);
			y_up = y_down+1;
			
			//bilinear interpolation of the quantities
			ux_down = (double(x_right)-x_back)*ux[x_right + y_down*ly] + (x_back - double(x_left))*ux[x_left + y_down*ly ];
			ux_up = (double(x_right)-x_back)*ux[x_right + y_up*ly] + (x_back - double(x_left))*ux[x_left + y_up*ly ];
			ux_back = (double(y_up) - y_back)*ux_up + (y_back - double(y_down))*ux_down;
			
			uy_down = (double(x_right)-x_back)*uy[x_right + y_down*ly] + (x_back - double(x_left))*uy[x_left + y_down*ly ];
			uy_up = (double(x_right)-x_back)*uy[x_right + y_up*ly] + (x_back - double(x_left))*uy[x_left + y_up*ly ];
			uy_back = (double(y_up) - y_back)*uy_up + (y_back - double(y_down))*uy_down;

			rho_down = (double(x_right)-x_back)*density[x_right + y_down*ly] + (x_back - double(x_left))*density[x_left + y_down*ly ];
			rho_up = (double(x_right)-x_back)*density[x_right + y_up*ly] + (x_back - double(x_left))*density[x_left + y_up*ly ];
			rho_back = (double(y_up) - y_back)*rho_up + (y_back - double(y_down))*rho_down;

			u2_back = ux_back*ux_back + uy_back*uy_back;
			//use these quantities for the equilibrium distribution
			//rho_back*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 6] = rho_back*rt[6]*(1. + (ex[6]*ux_back+ey[6]*uy_back)/cs2 + (ex[6]*ux_back+ey[6]*uy_back)*(ex[6]*ux_back+ey[6]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 2] = rho_back*rt[2]*(1. + (ex[2]*ux_back+ey[2]*uy_back)/cs2 + (ex[2]*ux_back+ey[2]*uy_back)*(ex[2]*ux_back+ey[2]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 5] = rho_back*rt[5]*(1. + (ex[5]*ux_back+ey[5]*uy_back)/cs2 + (ex[5]*ux_back+ey[5]*uy_back)*(ex[5]*ux_back+ey[5]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));

			//std::cout<<"time: "<<time_step<<" We are inside the domain"<<std::endl;
		}else{//case when we leave the domain. In this case just take the former value.
			/*
			std::cout<< "delta t "<<delta_t<<std::endl;
			std::cout<< "ux_loc "<<ux_loc<<std::endl;
			std::cout<< " x-back: "<<x_back<<std::endl;
			std::cout<< " y-back: "<<y_back<<std::endl;
			std::cout<<"time: "<<time_step<<" We are outside the domain"<<std::endl;
			*/
			double u2=ux[pos]*ux[pos]+uy[pos]*uy[pos];
			ftemp[Q*pos + 6] = density[pos]*rt[6]*(1. + (ex[6]*ux[pos]+ey[6]*uy[pos])/cs2 + (ex[6]*ux[pos]+ey[6]*uy[pos])*(ex[6]*ux[pos]+ey[6]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 2] = density[pos]*rt[2]*(1. + (ex[2]*ux[pos]+ey[2]*uy[pos])/cs2 + (ex[2]*ux[pos]+ey[2]*uy[pos])*(ex[2]*ux[pos]+ey[2]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 5] = density[pos]*rt[5]*(1. + (ex[5]*ux[pos]+ey[5]*uy[pos])/cs2 + (ex[5]*ux[pos]+ey[5]*uy[pos])*(ex[5]*ux[pos]+ey[5]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
		}
	}
	
	//Stable fluid at the back wall
	x = lx-1;
	for(y = 1; y<ly-1;y++){

		pos=x+lx*y;
		//define local and interpolated quantities.
		double ux_loc, uy_loc;
		double x_back, y_back, ux_back, uy_back, rho_back, u2_back;
		int x_left, x_right, y_up, y_down;
		double ux_down, ux_up, uy_down, uy_up, rho_down, rho_up;
		ux_loc = ux[pos];
		uy_loc = uy[pos];
		//Euler Step for backtracking
		x_back = double(x) - delta_t*ux_loc;
		y_back = double(y) - delta_t*uy_loc;
		//check if we are in the physical domain
		/*
		std::cout<< "delta t "<<delta_t<<std::endl;
		std::cout<< "ux_loc "<<ux_loc<<std::endl;
		std::cout<< " x-back: "<<x_back<<std::endl;
		std::cout<< " y-back: "<<y_back<<std::endl;
		*/

		if(((x_back >= 0.)&&(x_back<double(lx-1)))&&((y_back>=0.)&&(y_back<=double(ly-1)))){
			x_left = trunc(x_back);
			x_right = x_left+1;
			y_down = trunc(y_back);
			y_up = y_down+1;
			
			//bilinear interpolation of the quantities
			ux_down = (double(x_right)-x_back)*ux[x_right + y_down*ly] + (x_back - double(x_left))*ux[x_left + y_down*ly ];
			ux_up = (double(x_right)-x_back)*ux[x_right + y_up*ly] + (x_back - double(x_left))*ux[x_left + y_up*ly ];
			ux_back = (double(y_up) - y_back)*ux_up + (y_back - double(y_down))*ux_down;
			
			uy_down = (double(x_right)-x_back)*uy[x_right + y_down*ly] + (x_back - double(x_left))*uy[x_left + y_down*ly ];
			uy_up = (double(x_right)-x_back)*uy[x_right + y_up*ly] + (x_back - double(x_left))*uy[x_left + y_up*ly ];
			uy_back = (double(y_up) - y_back)*uy_up + (y_back - double(y_down))*uy_down;

			rho_down = (double(x_right)-x_back)*density[x_right + y_down*ly] + (x_back - double(x_left))*density[x_left + y_down*ly ];
			rho_up = (double(x_right)-x_back)*density[x_right + y_up*ly] + (x_back - double(x_left))*density[x_left + y_up*ly ];
			rho_back = (double(y_up) - y_back)*rho_up + (y_back - double(y_down))*rho_down;

			u2_back = ux_back*ux_back + uy_back*uy_back;
			//use these quantities for the equilibrium distribution
			//rho_back*rt[i]*(1. + (ex[i]*ux[pos]+ey[i]*uy[pos])/cs2 + (ex[i]*ux[pos]+ey[i]*uy[pos])*(ex[i]*ux[pos]+ey[i]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 6] = rho_back*rt[6]*(1. + (ex[6]*ux_back+ey[6]*uy_back)/cs2 + (ex[6]*ux_back+ey[6]*uy_back)*(ex[6]*ux_back+ey[6]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 3] = rho_back*rt[3]*(1. + (ex[3]*ux_back+ey[3]*uy_back)/cs2 + (ex[3]*ux_back+ey[3]*uy_back)*(ex[3]*ux_back+ey[3]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));
			ftemp[Q*pos + 7] = rho_back*rt[7]*(1. + (ex[7]*ux_back+ey[7]*uy_back)/cs2 + (ex[7]*ux_back+ey[7]*uy_back)*(ex[7]*ux_back+ey[7]*uy_back)/(2.*cs2*cs2) - u2_back/(2.*cs2));

			//std::cout<<"time: "<<time_step<<" We are inside the domain"<<std::endl;
		}else{//case when we leave the domain. In this case just take the former value.
			/*
			std::cout<< "delta t "<<delta_t<<std::endl;
			std::cout<< "ux_loc "<<ux_loc<<std::endl;
			std::cout<< " x-back: "<<x_back<<std::endl;
			std::cout<< " y-back: "<<y_back<<std::endl;
			std::cout<<"time: "<<time_step<<" We are outside the domain"<<std::endl;
			*/
			double u2=ux[pos]*ux[pos]+uy[pos]*uy[pos];
			ftemp[Q*pos + 6] = density[pos]*rt[6]*(1. + (ex[6]*ux[pos]+ey[6]*uy[pos])/cs2 + (ex[6]*ux[pos]+ey[6]*uy[pos])*(ex[6]*ux[pos]+ey[6]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 3] = density[pos]*rt[3]*(1. + (ex[3]*ux[pos]+ey[3]*uy[pos])/cs2 + (ex[3]*ux[pos]+ey[3]*uy[pos])*(ex[3]*ux[pos]+ey[3]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
			ftemp[Q*pos + 7] = density[pos]*rt[7]*(1. + (ex[7]*ux[pos]+ey[7]*uy[pos])/cs2 + (ex[7]*ux[pos]+ey[7]*uy[pos])*(ex[7]*ux[pos]+ey[7]*uy[pos])/(2.*cs2*cs2) - u2/(2.*cs2));
		}


	}

	//corner treatment
	{
		/*
		// south-west corner of the inlet has to be defined
		pos = lx; // = 0+lx*1
		density_loc = 0.0;
		for(i=0; i<Q; ++i)
		{
			density_loc += ftemp[Q*pos + i];
		}
		pos = 0;
		ftemp[Q*pos + 2] = ftemp[Q*pos + 4];
		ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
		ftemp[Q*pos + 5] = ftemp[Q*pos + 7];
		ftemp[Q*pos + 6] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 1] + ftemp[Q*pos + 5]));
		ftemp[Q*pos + 8] = ftemp[Q*pos + 6];

		// north-west corner of the inlet has to be defined
		pos = lx*(ly-2); // = 0+lx*(ly-2)
		density_loc = 0.0;
		for(i=0; i<Q; ++i)
		{
			density_loc += ftemp[Q*pos + i];
		}
		pos = lx*(ly-1);
		ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
		ftemp[Q*pos + 1] = ftemp[Q*pos + 3];
		ftemp[Q*pos + 8] = ftemp[Q*pos + 6];
		ftemp[Q*pos + 7] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 3] + ftemp[Q*pos + 6]));
		ftemp[Q*pos + 5] = ftemp[Q*pos + 7];

		//south-east corner of the outlet has to be defined
		pos = (lx-1)+lx; // = (lx-1)+lx*1
		density_loc = 0.0;
		for(i=0; i<Q; ++i)
		{
			density_loc += ftemp[Q*pos + i];
		}
		pos = lx-1;// = (lx-1)+lx*0
		ftemp[Q*pos + 2] = ftemp[Q*pos + 4];
		ftemp[Q*pos + 3] = ftemp[Q*pos + 1];
		ftemp[Q*pos + 6] = ftemp[Q*pos + 8];
		ftemp[Q*pos + 5] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 2] + ftemp[Q*pos + 3] + ftemp[Q*pos + 6]));
		ftemp[Q*pos + 7] = ftemp[Q*pos + 5];

		// north-east corner of the inlet has to be defined
		pos = (lx-1)+lx*(ly-2); // = (lx-1)+lx*(ly-2)
		density_loc = 0.0;
		for(i=0; i<Q; ++i)
		{
			density_loc += ftemp[Q*pos + i];
		}
		pos = (lx-1)+lx*(ly-1);
		ftemp[Q*pos + 4] = ftemp[Q*pos + 2];
		ftemp[Q*pos + 3] = ftemp[Q*pos + 1];
		ftemp[Q*pos + 7] = ftemp[Q*pos + 5];
		ftemp[Q*pos + 8] = 0.5*(density_loc - ftemp[Q*pos] - 2.*(ftemp[Q*pos + 4] + ftemp[Q*pos + 3] + ftemp[Q*pos + 7]));
		ftemp[Q*pos + 6] = ftemp[Q*pos + 8];
		*/
	}
	//obstacle treatment
	for (y=1; y<ly-1; ++y)
	{
		for (x=1; x<lx-1; ++x)
		{
			pos=x+lx*y;
			if (obst[pos]) // bounce-back at all inner obstacle nodes
			{
				f[Q*pos+1] = ftemp[Q*pos+3];
				f[Q*pos+2] = ftemp[Q*pos+4];
				f[Q*pos+3] = ftemp[Q*pos+1];
				f[Q*pos+4] = ftemp[Q*pos+2];
				f[Q*pos+5] = ftemp[Q*pos+7];
				f[Q*pos+6] = ftemp[Q*pos+8];
				f[Q*pos+7] = ftemp[Q*pos+5];
				f[Q*pos+8] = ftemp[Q*pos+6];
			}
		}
	}

}

/* -------------------------------------------------------------------------
function that checks the density
--------------------------------------------------------------------------*/
bool check_density(vector<double>& density, vector<double>& ux, vector<double>&uy, int t, double intdensity0, int dt_density) {
	double intdensity=0.;
	int x,y;
	long pos;
	double ma=0.;
	for (y=0; y<ly; ++y) {
		for (x=0; x<lx; ++x) {
			pos=x+lx*y;
			intdensity+=density[pos];
			ma=max(ma,sqrt(3.*(ux[pos]*ux[pos] + uy[pos]*uy[pos])));
		}
	}
	if (t%(dt_density*20)==0) {
		printf("# time step\t mass dev\tmax(Ma)\n");
	}
	printf("%11d\t% .5f %%\t%.8f", t, (intdensity-intdensity0)/intdensity0*100., ma);
	cout << endl;
	return(0);
}

//function that initializes the pml

bool update_pml(vector<double>& pml, bool with_x){

	int pos;

	for(int x = 0; x<lx; x++){
		for(int y = 0; y<ly;y++){
			
			pos = x + y*ly;
			
			pml[pos] = 0.;
			if(thickness_pml>0){
				if(y<thickness_pml){
					pml[pos] = sigma_max*double(std::pow(thickness_pml-y,2))/double(std::pow(thickness_pml,2));
				}
				if(y>ly-thickness_pml){
					pml[pos] = sigma_max*double(std::pow(thickness_pml- (ly-y),2))/double(std::pow(thickness_pml,2)); 
				}
				if(x> lx-thickness_pml){
					pml[pos] = sigma_max*double(std::pow(thickness_pml - (lx-x),2))/double(std::pow(thickness_pml,2));
				}
				//check the corners
				if((y<thickness_pml)&&(x> lx-thickness_pml)){
					pml[pos] = 0.95;
					pml[pos] = sigma_max*0.5*(double(std::pow(thickness_pml-y,2))/double(std::pow(thickness_pml,2)) + double(std::pow(thickness_pml - (lx-x),2))/double(std::pow(thickness_pml,2)));
				}
				if((y>ly-thickness_pml)&&(x> lx-thickness_pml)){
					pml[pos] = 0.5*sigma_max*(double(std::pow(thickness_pml- (ly-y),2))/double(std::pow(thickness_pml,2)) + double(std::pow(thickness_pml - (lx-x),2))/double(std::pow(thickness_pml,2)));
				}
				if((y<thickness_pml)&&(x<thickness_pml)){
					pml[pos] = 0.95;
					pml[pos] = 0.5*sigma_max*(double(std::pow(thickness_pml-y,2))/double(std::pow(thickness_pml,2)) + double(std::pow(thickness_pml-x,2))/double(std::pow(thickness_pml,2)));
				}
				if((y>ly-thickness_pml)&&(x<thickness_pml)){
					pml[pos] = 0.95;
					pml[pos] = 0.5*sigma_max*(double(std::pow(thickness_pml-x,2))/double(std::pow(thickness_pml,2)) + double(std::pow(thickness_pml- (ly-y),2))/double(std::pow(thickness_pml,2)));
				}
				if(with_x){
					if(x<thickness_pml ){
						pml[pos] = sigma_max*double(std::pow(thickness_pml-x,2))/double(std::pow(thickness_pml,2));
					}
					/*if((y>ly-thickness_pml)&&(x<thickness_pml)){
						pml[pos] = 0.95;
						pml[pos] = 0.5*(0.9/thickness_pml * x + 0.1 -0.9/thickness_pml * y + (0.1 + 0.9/thickness_pml * ly) );
					}
					if((y<thickness_pml)&&(x<thickness_pml)){
						pml[pos] = 0.95;
						pml[pos] = 0.5*(0.9/thickness_pml * y + 0.1 + 0.9/thickness_pml * x + 0.1);
					}*/
				}
			}
		}
	}
	return(0);
}

bool pml_advection(vector<double>& ftemp, const vector<double>& pml){
	int pos;
	double damping;
	for(int x = 0; x<lx; x++){
		for( int y = 0; y<ly; y++){
			pos = x + lx*y;
			if(pml[pos]<1){
				damping = 0.;
				
				{
					damping +=(1-pml[pos])*ftemp[Q*pos + 1];
					ftemp[Q*pos + 1]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 2];
					ftemp[Q*pos + 2]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 3];
					ftemp[Q*pos + 3]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 4];
					ftemp[Q*pos + 4]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 5];
					ftemp[Q*pos + 5]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 6];
					ftemp[Q*pos + 6]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 7];
					ftemp[Q*pos + 7]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 8];
					ftemp[Q*pos + 8]*=pml[pos];
				}
				/*
				else if((y>ly-thickness_pml) &&((x>thickness_pml)&&(x<lx-thickness_pml))){
					damping +=(1-pml[pos])*ftemp[Q*pos + 1];
					ftemp[Q*pos + 1]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 2];
					ftemp[Q*pos + 2]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 3];
					ftemp[Q*pos + 4]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 5];
					ftemp[Q*pos + 5]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 6];
					ftemp[Q*pos + 6]*=pml[pos];
					
				
				}
				else if((x<thickness_pml)&&((y>thickness_pml)&&(y<ly-thickness_pml))){
					
					damping +=(1-pml[pos])*ftemp[Q*pos + 2];
					ftemp[Q*pos + 2]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 3];
					ftemp[Q*pos + 3]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 4];
					ftemp[Q*pos + 4]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 6];
					ftemp[Q*pos + 6]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 7];
					ftemp[Q*pos + 7]*=pml[pos];
				
				}
				else if(x>lx-thickness_pml &&((y>thickness_pml)&&(y<ly-thickness_pml))){
					damping +=(1-pml[pos])*ftemp[Q*pos + 1];
					ftemp[Q*pos + 1]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 2];
					ftemp[Q*pos + 2]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 4];
					ftemp[Q*pos + 4]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 5];
					ftemp[Q*pos + 5]*=pml[pos];
					damping +=(1-pml[pos])*ftemp[Q*pos + 8];
					ftemp[Q*pos + 8]*=pml[pos];
				
				}
				*/
				ftemp[Q*pos + 0] +=damping;
			}
		}
	}

	return(0);
}