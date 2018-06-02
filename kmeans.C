//--------------------------------------------
// Implementation of k-means clustering
// using the Hartigan-Wong algorithm.
// For now, only in two dimensions...
//
// JDOK
// 06-01-18
//--------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdlib>

//----------------------
// Structs
//----------------------

struct point {
	float x;
	float y;
};

//----------------------
// Variables
//----------------------

//Number of clusters
const int K = 3;

//Points to be clustered
std::vector<point> points;

//Array to contain all clusters
std::vector<std::map<int, point>> S;

//----------------------
// Functions
//----------------------


void calculateClusterMean(std::map<int, point> s, float mean_x, float mean_y)
{
	int clusterSize = s.size();
	float x_cm = 0.0;
	float y_cm = 0.0;

	map <int, point> :: iterator itr;
	for (itr = s.begin(); itr != s.end(); ++itr)
	{
		int index = itr->first;
		point p   = itr->second;
		x_cm += p.x;
		y_cm += p.y;
	}

	x_cm = (float) x_cm / clusterSize;
	y_cm = (float) y_cm / clusterSize;

	mean_x = x_cm;
	mean_y = y_cm;
}


/*
 * Calculate the "cost" or spread of a given cluster
 * phi = \sum(x-\mu)^2
 */
float calculateCost(std::map<int, point> s)
{
	//Compute centroid of cluster
	float x_cm, y_cm;
	calculateClusterMean(s, x_cm, y_cm);

	//Loop over all points in cluster and find distance to the centroid
	float cost = 0.0;
	map <int, point> :: iterator itr;
	for (itr = s.begin(); itr != s.end(); ++itr)
	{
		int index = itr->first;
		point p   = itr->second;

		float distsq = TMath::Power(p.x - x_cm, 2.0) + TMath::Power(p.y - y_cm, 2.0);
		cost += distsq;
	}

	return cost;
}


/*
 *
 */


/*
 * Load a set of 2D points corresponding to the initial nucleon
 * coordinates from AMTP event record (npart-xy.dat)
 */
void loadPoints()
{
	string fname = "npart-xy.dat";
	ifstream infile(fname.c_str());

	string line;

	//Header information
	int eventnumber = 0;
	float iterflag = 0;
	int aproj = 0;
	int atarg = 0;
	float impactparameter = 0;

	//Nucleon information
	float x = 0;
	float y = 0;
	int parentnucleus = 0;
	int collcategory = 0;
	float nucleonimpactparameter = 0;
	int initialflavor = 0;
	int finalflavor = 0;

	while (getline(infile, line))
	{
		istringstream iss(line);

		//If a new event is found...
		if (iss >> eventnumber >> iterflag >> aproj >> atarg >> impactparameter)
		{
			//This file lists all nucleons in the collision (spectator + wounded)
			for (int i = 1; i <= aproj + atarg; i++)
			{
				point nucleon;

				getline(infile, line);
				istringstream iss(line);
				iss >> x >> y >> parentnucleus >> collcategory >> nucleonimpactparameter >> initialflavor >> finalflavor;

				nucleon.x = x;
				nucleon.y = y;

				//Was the nucleon a participant?
				if (collcategory > 0)
				{
					points.push_back(nucleon);
				}
			}
		}
	}
}


/*
 * Take the points and randomly assign them to clusters
 * There must be at least as many points as there are clusters
 */
void initializeClustersRandomly()
{
	//Require at least as many points as clusters
	if (points.size() < K) return;

	//Assign the first K points to the first K clusters
	for (int i = 0; i < K; i++)
	{
		std::map<int, point> cluster;
		cluster[0] = points[i];
		S.push_back(cluster);
	}

	//For the rest of the points, assign them randomly to an existing cluster
	for (int i = K; i < points.size(); i++)
	{
		int clus = rand() % K;
		S[clus][i] = points[i];
	}
}


void kmeans()
{
	//Initialize points and store them in a vector
	loadPoints();

	//Randomly assign points to clusters
	initializeClustersRandomly();

}