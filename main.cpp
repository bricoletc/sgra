#include <map>
#include "load_mono.c"
#include <vector>
#include <cmath> //For exponentiation
#include <stack>
#include <limits.h> //For INT_MIN 
using namespace std;

#define MAX_MASSES 100
#define MONOISOTOPIC_PRECISION 3 //The number of decimal places the monoisotopic masses are stored at
#define NINF INT_MIN
#define EPSILON 0.0001 //The 'error' about the monoisotopic mass tolerated for finding a residue.

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

	void addEdge(int u, int v, char residue);

	void TopologicalSortUtil(int v, bool visited[]);

	void LongestPathUtil(int v, int* max_dist, string* longest_peptide);

public:
	Graph(int V); //Constructor

	void Populate(double masses_array[],vector <pair <double,char> > residue_masses);

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

void Graph::Populate(double masses_array[],vector <pair <double,char> > residue_masses){
	/*Populate the graph with valid residues
	 * Residues are edges between two mass peaks.
	 * */

	int i,j;
	double first_mass,second_mass,tested_mass;

	vector <pair<double,char> >::const_iterator residue;

	for (i=0;i<V-1;++i){
		first_mass=masses_array[i];
		for (j=i+1;j<V;++j){
			second_mass=masses_array[j];
			tested_mass=abs(second_mass-first_mass);

			for (residue=residue_masses.begin();residue<residue_masses.end(); ++residue){
				if (abs(tested_mass-residue->first)<EPSILON) {
				
					//cout << "Residue found: " << residue->second << endl;
					if (second_mass>first_mass) this->addEdge(i,j,residue->second);
					
					else this->addEdge(j,i,residue->second); 

					break; //Shouldn't have any residues below EPSILON apart in mass.
				}
			}
		}
	
	}
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

void Graph::LongestPathUtil(int v, int* max_dist, string* longest_peptide){
	int dist[V];
	string peptides[V];

	for (int i=0;i<V;++i){
	       	dist[i]=NINF;
		peptides[i]="";
	}
	dist[v]=0;


	vector <pair <int,char> >::iterator it;
	stack <int> stack_copy=Stack; 

	while(stack_copy.empty()==false){
		int node=stack_copy.top();
		if (dist[node]>NINF){
			for (it=adj[node].begin();it<adj[node].end();++it){
				if (dist[it->first]<dist[node]+1){ //1 because edges are unweighted here; all edge weights are 1.
					dist[it->first]=dist[node]+1;
					peptides[it->first]=peptides[node]+it->second;	
				}
			}	
		}
		stack_copy.pop();
	}
	for (int i=0;i<V;++i){
		if (dist[i]>*max_dist){
			*max_dist=dist[i];
			*longest_peptide=peptides[i];
			//cout << peptides[i] << endl;
		}
	}
}

void Graph::LongestPath(){
	int max_dist=0;
	string longest_peptide="";

	for (int v=0;v<V;++v){
		LongestPathUtil(v,&max_dist,&longest_peptide);
	}

	cout << "Longest path: " <<max_dist << endl;
	cout << "Longest peptide: " << longest_peptide << endl;
}

/* Part III:
 *
 * Load up data, and process it.
 */

int main(int argc, char* argv[]){

	/* Load monoisotopic amino acid masses */
	string fname="monoisotopic_table.txt";

	vector <pair<double,char> > vectorial_load_monoisotopic_masses(string fname);

	vector <pair<double,char> > residue_masses=vectorial_load_monoisotopic_masses(fname);



	/* Load masses */
	fname=argv[1];
	double masses_array[MAX_MASSES];
	double mass;
	int num_masses=0;

	ifstream input(fname);

	while(input >> mass){
		masses_array[num_masses++]=mass;	
	}


	Graph graph(num_masses);

	graph.Populate(masses_array,residue_masses);


	graph.TopologicalSort();

	graph.LongestPath();

	return 0;
}
