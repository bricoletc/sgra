#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <map>
#include <vector>

using namespace std;

map<double, char> load_monoisotopic_masses(string fname){
	/*
	 *Loads monoisotopic masses from file containing it
	 */

	map<double, char> m;

	string line;
	ifstream file(fname);

	while(getline(file,line)){
	
		double weight;
		char residue;

		stringstream linestream(line);
		linestream >> residue >> weight;		

		m[round(weight*1000.0)/1000.0]=residue;
	}

	file.close();

	return m;
}

vector <pair<double,char> > vectorial_load_monoisotopic_masses(string fname){
	/*
	 *Loads monoisotopic masses from file containing it
	 */

	vector < pair <double,char> > vec;

	string line;
	ifstream file(fname);

	while(getline(file,line)){
	
		double weight;
		char residue;

		stringstream linestream(line);
		linestream >> residue >> weight;		

		vec.push_back(make_pair(weight,residue));
	}

	file.close();

	return vec;
}
