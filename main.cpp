#include <map>
#include "load_mono.c"
#include <vector>
#include <cmath> //For exponentiation
#include <stack>
using namespace std;

#define MAX_MASSES 100
#define MONOISOTOPIC_PRECISION 3 //The number of decimal places the monoisotopic masses are stored at

/* Disclaimer: Much code reuse from Geeks for Geeks:
 * -longest path between any pair of vertices
 * -https://www.geeksforgeeks.org/find-longest-path-directed-acyclic-graph/ 
 *  /

/* Part I:
 *
 * Flesh the structure of a Graph class. 
 * The graph by definition is directed (residues are path from one mass to the next) and acyclic (edges are necessarily between increasing masses). 
 * Whereas computing the longest path in a general graph is NP hard, a DAG has an efficient longest path computation.
 */

class Graph{
	
	int V; // Number of edges in graph

	/* Graph represented using adjacency array */
	vector <pair <int,char> >* adj; 

	stack <int> Stack; //Stores the topological sorting of the graph (Topmost elements of stack are parents).

	void TopologicalSortUtil(int v, bool visited[]);

public:
	Graph(int V); //Constructor

	void addEdge(int u, int v, char residue);

	void TopologicalSort();

	void LongestPath();
};


/* Part II:
 *
 * Implement graph populating, topological sorting, and longest path determination.
 */

Graph::Graph(int V){
	this->V=V;
	adj=new vector <pair <int,char> >[V];
	this->Stack=Stack;
}

void Graph::addEdge(int u, int v, char residue){
	adj[u].push_back(make_pair(v,residue));
}

void Graph::TopologicalSortUtil(int v, bool visited[]){
	/*
	 * Recursive algorithm which descends to graph leaves before adding any node to the Stack.
	 * Is a DFS with stack insertion of a node only once all its children have been explored.
	 */

	visited[v]=true;

	vector <pair <int,char> >::iterator it;

	for (it=adj[v].begin(); it!=adj[v].end();++it){
		if (!visited[it->first]) TopologicalSortUtil(it->first,visited);

	}	

	Stack.push(v);


}

void Graph::TopologicalSort(){
	/*
	 * Calls topological sorting for each node of the graph. This ensures a disconnected graph still ends up with all its nodes in the stack.
	 */
	bool* visited=new bool[V];
	for (int i=0;i<V;++i){
		visited[i]=false;
	}

	for (int i=0;i<V;++i){
		if (!visited[i]) TopologicalSortUtil(i,visited);
	}
	
}

void Graph::LongestPath(){
	continue;
}

/* Part III:
 *
 * Load up data, and process it.
 */

int main(int argc, char* argv[]){

	/* Load monoisotopic amino acid masses */
	string fname="monoisotopic_table.txt";

	map<double, char> load_monoisotopic_masses(string);

	map<double, char> residue_masses=load_monoisotopic_masses(fname);



	/* Load masses */
	fname=argv[1];
	double masses_array[MAX_MASSES];
	double mass;
	int num_masses=0;

	ifstream input(fname);

	while(input >> mass){
		masses_array[num_masses++]=mass;	
	}

	/*Populate the graph with valid residues
	 * Residues are edges between two mass peaks.
	 * */
	int i,j;
	double first_mass,second_mass,tested_mass;
	map<double,char>::const_iterator residue;

	Graph graph(num_masses);

	for (i=0;i<num_masses-1;++i){
		first_mass=masses_array[i];
		for (j=i+1;j<num_masses;++j){
			second_mass=masses_array[j];
			tested_mass=round(abs(second_mass-first_mass)*pow(10.0,MONOISOTOPIC_PRECISION))/pow(10.0,MONOISOTOPIC_PRECISION);

			residue=residue_masses.find(tested_mass);

			if (residue!=residue_masses.end()){
				if (second_mass>first_mass) graph.addEdge(i,j,residue->second);
				
				else graph.addEdge(j,i,residue->second); 
			}
		}
	
	}

	return 0;
}