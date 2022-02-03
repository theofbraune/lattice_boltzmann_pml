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

const int M[Q][Q]={
 { 1, 1, 1, 1, 1, 1, 1, 1, 1},
 {-4,-1,-1,-1,-1, 2, 2, 2, 2},
 { 4,-2,-2,-2,-2, 1, 1, 1, 1},
 { 0, 1, 0,-1, 0, 1,-1,-1, 1},
 { 0,-2, 0, 2, 0, 1,-1,-1, 1},
 { 0, 0, 1, 0,-1, 1, 1,-1,-1},
 { 0, 0,-2, 0, 2, 1, 1,-1,-1},
 { 0, 1,-1, 1,-1, 0, 0, 0, 0},
 { 0, 0, 0, 0, 0, 1,-1, 1,-1}
};

const int N[Q][Q]={
 { 4,-4, 4, 0, 0, 0, 0, 0, 0},
 { 4,-1,-2, 6,-6, 0, 0, 9, 0},
 { 4,-1,-2, 0, 0, 6,-6,-9, 0},
 { 4,-1,-2,-6, 6, 0, 0, 9, 0},
 { 4,-1,-2, 0, 0,-6, 6,-9, 0},
 { 4, 2, 1, 6, 3, 6, 3, 0, 9},
 { 4, 2, 1,-6,-3, 6, 3, 0,-9},
 { 4, 2, 1,-6,-3,-6,-3, 0, 9},
 { 4, 2, 1, 6, 3,-6,-3, 0,-9}
};



/* -------------------------------------------------------------------------
function that performs the LB collision step (rhs of Boltzmann equation)
--------------------------------------------------------------------------*/
bool	collision_mrt(vector<double>& f, vector<double>& ftemp, vector<bool>& obst, double omegaf)// RHS of Boltzmann equation
{
	int x,y,pos,i,j;
	double pt1=1.63, pt2=1.14, pt4=1.92, pt6=1.92, pt7=omegaf, pt8=omegaf; // following lallemand_2000_pre_61_6546
// 	double pt1=1.64, pt2=1.54, pt4=1.7, pt6=1.7, pt7=omegaf, pt8=omegaf; // following mukherjee_2007_cf_36_1149
	double usq;
	double check_f;
	double mf[Q];

	for (y=ly-1; y>=0; --y)
	{
		for (x=0; x<lx; ++x)
		{
			pos=x+lx*y;
			if (!obst[pos]) // bounce-back at inner wall nodes; no slip
			{
				for (j=0; j<Q; ++j)
				{
					mf[j]=0.;
					for (i=0; i<Q; ++i)
					{
						mf[j]+=M[j][i]*ftemp[Q*pos+i];
					}
				}
				usq=mf[3]*mf[3]+mf[5]*mf[5];
				mf[1]=mf[1]-pt1*(mf[1]-(-2.*mf[0]+3.*usq));
				mf[2]=mf[2]-pt2*(mf[2]-(    mf[0]-3.*usq));
				mf[4]=mf[4]-pt4*(mf[4]+mf[3]);
				mf[6]=mf[6]-pt6*(mf[6]+mf[5]);
				mf[7]=mf[7]-pt7*(mf[7]-(mf[3]*mf[3]-mf[5]*mf[5]));
				mf[8]=mf[8]-pt8*(mf[8]-mf[3]*mf[5]);

				check_f=0.;
				for (i=0; i<Q; ++i)
				{
					f[Q*pos+i]=0.;
					for (j=0; j<Q; ++j)
					{
						f[Q*pos+i]+=N[i][j]*mf[j]*1./36.;
					}
					check_f+=f[Q*pos+i];
				}
				if (check_f < 0){cout << "fatal error: density < 0" << endl;exit(1);}
			}
		}
	}
	return (0);
}

