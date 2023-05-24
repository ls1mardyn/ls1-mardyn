#include "Graphit.h"
#include "Random.h"

#include <math.h>
#include <iostream>

int Graphit::getNumberOfAtoms()
{
   return numberOfAtoms;
}

void Graphit::calculateCoordinatesOfAtoms(
   int numberOfLayers, double xLength, double zLength, double A, double wo_wall
) {
	double B = 2.0 * A;
	double C = 0.86602540378 * A;  // sin(pi/3)

	double xCoor = 0.05*A;
	double yCoor = 0.0;
	double zCoor = 0.05*C;
	numberOfAtoms=0;
	for(int Index=1; Index<=numberOfLayers; Index++)
	{
		/*******************
		* Creates the odd layer (1,3,5,..)
		*******************/
		if(Index%2==1)
		{
			xCoor = 0.1*A;
			zCoor = 0.1*C;
			int i = 0;
			int j = 0;
			int temp = 0;
                        if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
                        {
                        	coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
				coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
				coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
                                componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
				numberOfAtoms++;
			}
			while(zCoor < zLength)
			{
				while(xCoor < xLength)
				{
					if (i%2==0)
					{
						xCoor += B;
						i=1;
					}
					else
					{
						xCoor += A;
						i=0;
						temp=0;
					}
					if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
					{
						coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
						coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
						coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
                             			componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
						numberOfAtoms++;
					}
				}
				zCoor += C;
                                j++;
                                if(j == 14) j = 0;
				if(j%2==1)
				{
					i = 1;
					xCoor = 0.6*A;
				}
				else
				{
					i = 0;
					xCoor = 0.1*A;
				}
				if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
				{
					coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
					coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
					coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
                             		componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
					numberOfAtoms++;
				}
				temp = j;
			}
		}
		/********************
		* Creates the even layer
		*********************/
		else
		{
			// Odd Layer has a different starting point
			xCoor = B-0.4*A;
			zCoor = 0.1*C;
			int i = 0; //Makes it possible to change the different lengths in one direction
			int j = 7; //Necessary for changing the different startpoints of the xCoor
			
			// Code for writing the Coordinates together
                        if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
                        {
				coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
				coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
				coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
              			componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
				numberOfAtoms++;
			}
			while(zCoor < zLength)
			{
				while(xCoor < xLength)
				{
					if (i%2==0)
					{
						xCoor += A;
						i=0;
					}
					else
					{
						xCoor += B;
						i=1;
					}
					if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
					{
						coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
						coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
						coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
                             			componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
						numberOfAtoms++;
					}
			
				}
				zCoor += C;
                                j++;
                                if(j == 14) j = 0;
				if(j%2==1)
				{
					i = 1;
					xCoor = 0.1*A;
				}
				else
				{
					i = 0;
					xCoor = B-0.4*A;
				}	
				if((xCoor < xLength && zCoor < zLength) && ((zCoor < (0.5 - 0.5*wo_wall)*zLength) || (zCoor > (0.5 + 0.5*wo_wall)*zLength)))
				{
					coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
					coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
					coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
                      			componentsOfAtoms[numberOfAtoms] = this->comp(i, j);
					numberOfAtoms++;
				}
			}
		}
		yCoor = yCoor + Z;
	}
}

double Graphit::getX(int number)
{
	return coordinatesOfAtoms[0][number];
}

double Graphit::getY(int number)
{
	return coordinatesOfAtoms[1][number];
}

double Graphit::getZ(int number)
{
	return coordinatesOfAtoms[2][number];
}

void Graphit::calculateVelocities(double T, double U)
{	
   double absoluteVelocity = sqrt(3.0*T / ATOMIC_MASS_C);
   double phi, omega;

   Random* r = new Random();
   r->init((int)(10000000.0*T) - (int)(3162300.0*U));
	
   for(int i=0; i < this->numberOfAtoms; i++)
   {
      phi = 6.283185 * r->rnd();
      omega = 6.283185 * r->rnd();
      velocitiesOfAtoms[0][i] = absoluteVelocity*cos(phi)*cos(omega);
      velocitiesOfAtoms[1][i] = absoluteVelocity*cos(phi)*sin(omega);
      velocitiesOfAtoms[2][i] = U + absoluteVelocity*sin(phi);
   }
}

double  Graphit::getVelocityX(int number)
{
   return velocitiesOfAtoms[0][number];
}

double Graphit::getVelocityY(int number)
{
   return velocitiesOfAtoms[1][number];
}

double Graphit::getVelocityZ(int number)
{
   return velocitiesOfAtoms[2][number];
}

void Graphit::reset()
{
   this->numberOfAtoms = 0;
   for(int d=0; d<4; d++)
   {
      this->coordinatesOfAtoms[d].clear();
      this->velocitiesOfAtoms[d].clear();
   }
}

unsigned Graphit::comp(int ti, int tj)
{
   int i = ti % 2;
   int j = tj % 14;

   // std::cout << "\t\t(i, j) \t=\t (" << i << ", " << j << ")\n";
   // std::cout.flush();

   if((j == 2) && (i == 1)) return CID_I;
   else if((j ==  3) && (i == 0)) return CID_I;
   else if((j ==  3) && (i == 1)) return CID_ZI;
   else if((j ==  7) && (i == 0)) return CID_I;
   else if((j ==  7) && (i == 1)) return CID_Z;
   else if((j == 11) && (i == 0)) return CID_I;
   else if((j == 11) && (i == 1)) return CID_ZI;
   else if((j == 12) && (i == 1)) return CID_I;
   else return CID_C;
}

