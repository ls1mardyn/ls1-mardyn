/*
 * GNU GPL version 2, original version by D. Dylus
 */

#include "Graphit.h"
#include "Random.h"

#include <math.h>
#include <iostream>
using namespace std;

int Graphit::getNumberOfAtoms()
{
   return numberOfAtoms;
}

void Graphit::calculateCoordinatesOfAtoms(
   int numberOfLayers, double xLength, double zLength, double A
) {
	double B = 2.0 * A;
	double C = 0.86602540378 * A;  // sin(pi/3)

	double xCoor = 0.1*A;
	double yCoor = 0.0;
	double zCoor = 0.1*C;
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
			coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
			coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
			coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
			numberOfAtoms++;
			while(zCoor < zLength)
			{
				while(xCoor < xLength)
				{
					if (i%2==0 && temp%2==0)
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
					if (xCoor < xLength && zCoor < zLength)
					{
						coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
						coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
						coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
						numberOfAtoms++;
					}
				}
				i=0;
				zCoor += C;
				if(j%2==0)
				{
					xCoor = 0.6*A;
					j=1;
				}
				else
				{
					xCoor = 0.1*A;
					j=0;
				}
				if (xCoor < xLength && zCoor < zLength)
				{
					coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
					coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
					coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
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
			int j = 0; //Necessary for changing the different startpoints of the xCoor
			
			// Code for writing the Coordinates together
			coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
			coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
			coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
			numberOfAtoms++;
			while(zCoor < zLength)
			{
				while(xCoor < xLength)
				{
					if (i%2==0)
					{
						xCoor += A;
						i=1;
					}
					else
					{
						xCoor += B;
						i=0;
					}
					if (xCoor < xLength && zCoor < zLength)
					{
						coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
						coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
						coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
						numberOfAtoms++;
					}
			
				}
				i=0;
				zCoor=zCoor+C;
				if(j%2==0)
				{
					xCoor = 0.1*A;
					j=1;
				}
				else
				{
					xCoor = B-0.4*A;
					j=0;
				}	
				if (xCoor < xLength && zCoor < zLength)
				{
					coordinatesOfAtoms[0][numberOfAtoms] = xCoor;
					coordinatesOfAtoms[1][numberOfAtoms] = yCoor;
					coordinatesOfAtoms[2][numberOfAtoms] = zCoor;
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
