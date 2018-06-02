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

/*
 * Calculates the centroid location of a given cluster
 */
void calculateClusterMean(std::map<int, point> s, float &mean_x, float &mean_y)
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
 * Evaluates the function Delta which should be minimized iteratively
 */
float evaluateDelta(point p, std::map<int, point> sn, std::map<int, point> sm)
{
	float x_cm_n, x_cm_m, y_cm_n, y_cm_m;
	calculateClusterMean(sn, x_cm_n, y_cm_n);
	calculateClusterMean(sm, x_cm_m, y_cm_m);

	float distsq_n = TMath::Power(p.x - x_cm_n, 2.0) + TMath::Power(p.y - y_cm_n, 2.0);
	float distsq_m = TMath::Power(p.x - x_cm_m, 2.0) + TMath::Power(p.y - y_cm_m, 2.0);

	float size_n = sn.size();
	float size_m = sm.size();

	return (size_n / (size_n + 1)) * distsq_n - (size_m / (size_m + 1)) * distsq_m;
}


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
 * Transfer point in cluster 1 to cluster 2
 */
void transferPoint(int point_index, int clus_index1, int clus_index2)
{
	//Add point to cluster 2
	S[clus_index2][point_index] = points[point_index];

	//Remove from cluster 1
	std::map<int, point>::iterator it;
	it = S[clus_index1].find(point_index);
	S[clus_index1].erase (it);
}


/*
 * Carry out Hartigan-Wong algorithm by reassigning points to
 * clusters until the function Delta is minimized
 */
void findClusters()
{
	//Point index and clusters that minimize Delta
	int minPointIndex = -999;
	int min_clus_index_n;
	int min_clus_index_m;
	float minDeltaVal = 1E10;

	//Loop over all existing clusters
	for (int iclus = 0; iclus < S.size(); iclus++)
	{
		std::map<int, point> cluster_n = S[iclus];

		//Loop over all points in the cluster at hand
		map <int, point> :: iterator itr;
		for (itr = cluster_n.begin(); itr != cluster_n.end(); ++itr)
		{
			int index = itr->first;
			point p   = itr->second;

			//Loop over all of the other clusters, trying to find a cluster 'm' that
			//minimizes Delta along with the current cluster 'n' and point 'p'
			//Clusters 'n' and 'm' have to be different!
			for (int jclus = 0; jclus < S.size(); jclus++)
			{
				if (iclus == jclus) continue;

				std::map<int, point> cluster_m = S[jclus];

				//Evaluate function Delta with cluster 'n', 'm', and point 'p'
				float deltaVal =  evaluateDelta(p, cluster_n, cluster_m);

				//Does this minimize Delta?
				if (deltaVal < minDeltaVal)
				{
					minDeltaVal      = deltaVal;
					minPointIndex    = index;
					min_clus_index_n = iclus;
					min_clus_index_m = jclus;
				}
			}
		}
	}

	//After the iteration, transfer point in cluster N that minimizes Delta to set M
	transferPoint(minPointIndex, min_clus_index_n, min_clus_index_m);
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

	//Carry out actual k-means algorithm
	findClusters();

}