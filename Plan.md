### Work plan

1. Building the graph

	For each possible pair of masses, check if their difference (largest-smallest) corresponds to a residue. If it does, draw a directed edge between them.

2. Graph structure

	A list where each element represents a node. Each element is a vector of pair<int,string>. Each element of the vector is a neighbour: the int represents the index of that node, and the string represents the residue produced.


3. Graph traversal

	Two alternatives:

	* DFS from each node
	* Topological sort
