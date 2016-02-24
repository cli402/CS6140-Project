#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include <stack>
#include <unordered_set>
#include <time.h>
#include <cstddef>         // std::size_t
using namespace std;

///////////////////////////////////////////////////////////////////////
//
// Graph struct
//
///////////////////////////////////////////////////////////////////////

struct Node {
	int src;            // source
	vector<int> dest;   // destination
	Node(int a) {
		src = a;
	};
};

struct Graph {
	int V, E;                  // V, # of node; E, # of edges
	vector<struct Node> nodes; // nodes
	Graph(int a, int b) {
		V = a;
		E = b;
	}
};
// split the string by space
void splitString(string &s, vector<int> &result);
// generate the graph
struct Graph parseEdges(string graph_file);
// parse the filename
string SplitFilename (string str);




//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for BnB begin  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
unordered_set<int> branchAndBound(string graph_file, struct Graph G, int cutoff);
void BnB(string graph_file, int cutoff);
void printNode(vector<Node> &nodes);
int vc_approx_bnb(vector<Node> nodes);
bool verify(vector<Node> &nodes);
vector<Node> generate_graph(vector<Node> nodes, unordered_set<int> &path);
//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for BnB end  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for Approx begin  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void approx_ED(string graph_file, int cutoff);
//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for Approx  end  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for LS1 begin  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
vector<int> vc_approx(struct Graph G);
int cost(struct Graph gra, vector<int> candidate);
vector<int> construct(struct Graph gra, vector<int> candidate);
void ESA(string graph_file, int cutoff, int seed);
//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for LS1  end  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for LS2 begin  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// greedy
unordered_map<int, int> VC_Greedy(struct Graph gra);
void local_search(string graph,int cutoff,int seed);
struct Node compute_lowest_score(struct Graph gra, unordered_map<int, int> cover, vector<unordered_map<int, int>> edgeweight);
void add_vertex(struct Graph *uncovered_edge, unordered_map<int, int> *cover, struct Graph gra);
unordered_map<int, int> ls2_vc_approx(struct Graph G);
bool check_VC(struct Graph gra, unordered_map<int, int> cover);
//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// function for LS2  end  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


////////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/////////////////////////




///////////////////////////////////////////////////////////////////////
//
// Main function
//
///////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

	string graph_file;
	string alg;
	int cutoff;
	int seed;

	if (argc < 9) {
		cout << "Usage: " << argv[0] << "-inst <input_file> −alg [BnB|Approx|LS1|LS2] -time <cutoff_in_seconds> -seed <random_seed>" << endl;
		return 1;
	}


  for(int i=1; i<8; ++i) {
    if( strcmp(argv[i],"-inst") ==0) 
      graph_file = argv[i+1];
    else if(strcmp(argv[i],"-alg")==0) 
      alg = argv[i+1];
    else if(strcmp(argv[i], "-time")==0) 
      cutoff = atoi(argv[i+1]);
    else if(strcmp(argv[i], "-seed")==0) 
      seed = atoi(argv[i+1]);
  }


	// run algorithm
	if (alg == "BnB") {
    	BnB(graph_file, cutoff);
	} else if (alg == "Approx") {
      approx_ED(graph_file, cutoff);
	} else if (alg == "LS1") {
		ESA(graph_file, cutoff, seed);
	} else if (alg == "LS2") {
		local_search(graph_file,cutoff,seed);
	}


	return 0;
}

/////////////////////////////////////////////////////
//
// Implement your function here
//
/////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
//////function for approx begin ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void approx_ED(string graph_file, int cutoff) {

	  struct Graph G = parseEdges(graph_file);
    clock_t start = clock();
	  vector<int> cover = vc_approx(G);
    clock_t end = clock();

	  string graph_name = SplitFilename(graph_file);
    string output_file = "../Output/Solution/Approx/" + graph_name + "_Approx_" + to_string(cutoff) + ".sol";
    ofstream output;
    output.open(output_file);
    output<<cover.size()<<endl;
	  for (auto each : cover) {
		     output << each << ",";
	  }

	  ofstream trace;
	  string trace_file = "../Output/Trace/Approx/" + graph_name + "_Approx_" + to_string(cutoff) + ".trace";
	  trace.open(trace_file);
    trace<<(end-start)/(float)CLOCKS_PER_SEC<<", "<<cover.size()<<endl;

}
////////////////////////////////////////////////////////////////////////////
//////function for approx end ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////
//////function for BnB begin ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void BnB(string graph_file, int cutoff){
    
    struct Graph G = parseEdges(graph_file);
	  string graph_name = SplitFilename(graph_file);
    
    string output_file1 = "../Output/Solution/BnB/" + graph_name + "_BnB_" + to_string(cutoff) + ".sol";
    
    ofstream output1;
    
    auto result = branchAndBound(graph_name, G, cutoff);
    
    //output the result into the file
    output1.open(output_file1);
    output1<<result.size()<<endl;
    for (auto i : result) {
        output1<<i<<",";
    }
    output1<<endl;
    output1.close();
}

unordered_set<int> branchAndBound(string graph_file, struct Graph G, int cutoff) {
    int m = G.V;
    stack< unordered_set<int> > stk;
    unordered_set<int> path;
    int record = vc_approx_bnb(G.nodes);
    vector<Node> nodes = G.nodes;
    unordered_set<int> result;
    stk.push(path);
    string output_file = "../Output/Trace/BnB/" + graph_file + "_BnB_" + to_string(cutoff) + ".trace";
    ofstream output;
    output.open(output_file);
    clock_t start = clock();
    clock_t end = clock();
    
    while (!stk.empty()) {
        unordered_set<int> new_path = stk.top();
        stk.pop();
        
        if (new_path.size() > record) continue;
        // time to cutoff
        end = clock();
        if ((end-start)/(float)CLOCKS_PER_SEC >= cutoff) {
            output.close();
            return result;
        }
        
        vector<Node> new_graph = generate_graph(nodes, new_path);
        if (verify(new_graph)) {
            if (new_path.size() < record) {
                record = new_path.size();
                result.clear();
                for (auto &m : new_path) {
                    result.insert(m);
                }
                output<<to_string((end-start)/(float)CLOCKS_PER_SEC)<<","<<to_string(result.size())<<endl;
            }
            continue;
        }
        // select a node in this graph
        int next = 1;
        int index = 0;
        for (auto mm : new_graph) {
            if (mm.dest.size() > index){
                index = mm.dest.size();
                next = mm.src;
            }
        }
        new_graph.clear();
        //case 1: choose this node
        new_path.insert(next);
        new_graph = generate_graph(nodes, path);
        int lower_bound = path.size() +  vc_approx_bnb(new_graph)/2;
        new_graph.clear();
        if (lower_bound < record) {
            stk.push(new_path);
        }
        new_path.erase(next);
        //case 2: do not choose this node
        for (auto k : nodes[next-1].dest) {
            new_path.insert(k);
        }
        new_graph = generate_graph(nodes, path);
        lower_bound = path.size() +  vc_approx_bnb(new_graph)/2;
        new_graph.clear();
        if (lower_bound < record) {
            stk.push(new_path);
        }
        for (auto k : nodes[next-1].dest) {
            new_path.erase(k);
        }
    }
    output.close();
    return result;
}

bool verify(vector<Node> &nodes) {
    for (auto &i : nodes){
        //cout<<i.src<<" "<<i.dest.size()<<endl;
        if (i.dest.size() > 0)
            return false;
    }
    return true;
}

vector<Node> generate_graph(vector<Node> nodes, unordered_set<int> &path) {
    if (path.size() == 0) return nodes;
    for (auto &i : path) {
        for (int m = 0; m<nodes.size(); m++) {
            if (i == nodes[m].src){
                nodes[m].dest.clear();
            }else {
                for (int n = 0; n<nodes[m].dest.size(); n++){
                    if (i == nodes[m].dest[n]) {
                        nodes[m].dest.erase(nodes[m].dest.begin() + n);
                        break;
                    }
                }
            }
        }
    }
    return nodes;
}

int vc_approx_bnb(vector<Node> nodes) {
    vector<int> VC;
    for(auto node_v: nodes) {
        int v = node_v.src;
        for(int u: node_v.dest) {
            if(u < v || nodes[u-1].dest.size()==0)
                continue;
            node_v.dest.clear();
            nodes[u-1].dest.clear();
            VC.push_back(v);
            VC.push_back(u);
            break;
        }
    }
    return VC.size();
}

////////////////////////////////////////////////////////////////////////////
//////function for BnB begin ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////
//////function for LS1 begin ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//
// cost function
//
//////////////////////////////////////////////////////

vector<int> vc_approx(struct Graph G) {

	vector<int> VC;

	for(auto node_v: G.nodes) {
	   	int v = node_v.src;
		for(int u: node_v.dest) {
			if(u < v || G.nodes[u - 1].dest.size()==0)
			continue;
			node_v.dest.clear();
			G.nodes[u - 1].dest.clear();
			VC.push_back(v);
			VC.push_back(u);
			break;
		}
	}

  	return VC;
}


//////////////////////////////////////////////////////
//
// cost function
//
//////////////////////////////////////////////////////

int cost(struct Graph gra, vector<int> candidate) 
{
	int count = 0;

	for (int i = 0; i < gra.V; i++) {
		if (find(candidate.begin(), candidate.end(), gra.nodes[i].src) == candidate.end()) {	
			for (auto each : gra.nodes[i].dest) {
				if (each > i && find(candidate.begin(), candidate.end(), each) == candidate.end())
					count++;
			}	
		}
	}

	return candidate.size() + count;
}

//////////////////////////////////////////////////////
//
// construct vertex cover
//
//////////////////////////////////////////////////////
vector<int> construct(struct Graph gra, vector<int> candidate) {
	for (auto each : candidate) {
		gra.nodes[each - 1].dest.clear();
	}

	for (int i = 0; i < gra.V; i++) {
		if (gra.nodes[i].dest.size() > 0) {
			for (auto each : gra.nodes[i].dest){
				if (gra.nodes[each - 1].dest.size() > 0) {
					candidate.push_back(i + 1);
					break;
				}
			}
			gra.nodes[i].dest.clear();
		}
	}

	return candidate;
}

////////////////////////////////////////////////////////////////////
//
// Efficient Simulate Annealing 
//
////////////////////////////////////////////////////////////////////

void ESA(string graph_file, int cutoff, int seed)
{
	struct Graph gra = parseEdges(graph_file);
	vector<int> candidate = vc_approx(gra);

	string graph_name = SplitFilename(graph_file);

	ofstream trace;
	string trace_file = "../Output/Trace/LS1/" + graph_name + "_LS1_" + to_string(cutoff) + "_" + to_string(seed) + ".trace";
	trace.open(trace_file);

	srand(seed);

	int step = 0;
	int iter = 20;
	double T = 1;

	int bestSofar = cost(gra, candidate);
	vector<int> best = candidate;

	clock_t start = clock(); 
	clock_t end = clock();
	float rt= (end - start) / (float) CLOCKS_PER_SEC;

	while (rt < cutoff && T > 0)  {

		for (int i = 0; i < iter; i++) {

			vector<int> neighbor(candidate);
			int selected = rand() % gra.V;

			auto iter = find(neighbor.begin(), neighbor.end(), selected);
			if (iter != neighbor.end()){
				neighbor.erase(iter);
				int delta = cost(gra, neighbor) - cost(gra, candidate);
				if (delta <= 0)
					candidate = neighbor;
				else {
					double P = exp(-delta * (1 + gra.nodes[selected].dest.size() / float(gra.V)) / T);
					double R = float(rand()) / RAND_MAX;

					if (R < P) {
						candidate = neighbor;
					}
				}
			}
			else {
				neighbor.push_back(selected);
				int delta = cost(gra, neighbor) - cost(gra, candidate);
				if (delta <= 0)
					candidate = neighbor;
				else {
					double P = exp(-delta * (1 - gra.nodes[selected].dest.size() / float(gra.V)) / T);
					double R = float(rand()) / RAND_MAX;

					if (R < P) {
						candidate = neighbor;
					}
				}
			}

			end = clock();
			rt= (end - start) / (float) CLOCKS_PER_SEC;

			int temp = cost(gra, candidate);
			if (temp < bestSofar){
				bestSofar = temp;
				best = candidate;
				trace << rt<< ", " << temp << endl;
			}

			if (rt > cutoff)
				break;
			
		}

		T = 0.95 * T;
	}

	trace.close();
	best = construct(gra, best);

	ofstream sol;
	string sol_file = "../Output/Solution/LS1/" + graph_name + "_LS1_" + to_string(cutoff) + "_" + to_string(seed) + ".sol";
	sol.open(sol_file);

	sol << best.size() << endl;
	for (auto each : best) {
		// cout << each << ",";
		sol << each << ",";
	}
	
	sol.close();

}


////////////////////////////////////////////////////////////////////////////
//////function for LS1 end   ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
//////function for LS2 begin ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
unordered_map<int, int> VC_Greedy(struct Graph gra)
{
	unordered_map<int, int> cover;

	while (gra.E > 0) {
		int next = 1;
		for (auto each : gra.nodes) {
			if (each.dest.size() > gra.nodes[next - 1].dest.size())
				next = each.src;
		}
		//cout<< "next "<<next<< " size of "<< gra.nodes[next-1].dest.size()<<endl;
		//cout<< "the cover size " << cover.size()<<endl;
		cover.insert({ next, next });

		for (auto each : gra.nodes[next - 1].dest) {
			for (int i = 0; i < gra.nodes[each - 1].dest.size(); i++) {
				if (gra.nodes[each - 1].dest[i] == next) {
					gra.nodes[each - 1].dest.erase(gra.nodes[each - 1].dest.begin() + i);
					break;
				}
			}
		}
		gra.E = gra.E - gra.nodes[next - 1].dest.size();
		gra.nodes[next - 1].dest.clear();

	}
	return cover;
}

void add_vertex(struct Graph &uncovered_edge, unordered_map<int, int> &cover, struct Graph gra){
	//choose a entering vertex randomly, then add it to the cover, remove edges from the uncoverd edge
	//cout << uncovered_edge.nodes.size() << "uncover" << endl;
	int v = rand() % uncovered_edge.nodes.size();
	int rand_next = rand() % uncovered_edge.nodes[v].dest.size();
	int next = uncovered_edge.nodes[v].dest[rand_next];
	//if there are multiple nodes in it, decide if the added would help reduce them
	uncovered_edge.nodes[v].dest.erase(uncovered_edge.nodes[v].dest.begin()+rand_next);
	if (uncovered_edge.nodes.size()>=1){
		int x = 1;
		int y = 0;
		for (int i = 0; i < uncovered_edge.nodes.size(); i++){
			if (uncovered_edge.nodes[i].src == next){
				x += uncovered_edge.nodes[i].dest.size();
				y += 1;
				uncovered_edge.nodes.erase(uncovered_edge.nodes.begin() + i);
			}
			else{
			for (int j = 0; j < uncovered_edge.nodes[i].dest.size(); j++){
				if (next == uncovered_edge.nodes[i].dest[j]){
					uncovered_edge.nodes[i].dest.erase(uncovered_edge.nodes[i].dest.begin() + j);
					x += 1;
					break;
				}
			}
		}
		}
		//cout << "x is " << x<<endl;
		uncovered_edge.V -= y;
		uncovered_edge.E -= x;
		//cout << "Edge number is "<< uncovered_edge.E<<endl;
	}
	// if no node left, just erase the node it self
	int x = 0;
	while(x<uncovered_edge.nodes.size()){
		if (uncovered_edge.nodes[x].dest.size() == 0){
			uncovered_edge.nodes.erase(uncovered_edge.nodes.begin()+x);
			uncovered_edge.V -= 1;
			//uncovered_edge.E -= 1;
		}
		else{ x++; }
	}

	//add the next to the cover
	int i = cover.size();
	cover.insert({ next, next });
	if (cover.size() == i){
		//cout << "wrong!!" << endl;
		//cout << "the node number" << next << endl;
	}
}
void initialize_edge_weight(struct Graph gra, vector< unordered_map<int, int> > *edgeweight){
	for (int i = 0; i< gra.nodes.size(); i++){
		unordered_map<int, int> node;
		for (auto j : gra.nodes[i].dest){
			//initialize the edge weight to 1, key to be the edge. the i-1,j ->> 1
			node.insert({ j, 1 });
		}
		edgeweight->push_back(node);
	}
}

void add_weight_uncoverd_edge(vector< unordered_map<int, int> > *edgeweight, struct Graph uncovered_edge, int a){
	vector< unordered_map<int, int> >& weight = *edgeweight;

	for (int i = 0; i<uncovered_edge.nodes.size(); i++){
		for (int j = 0; j<uncovered_edge.nodes[i].dest.size(); j++){
			weight[uncovered_edge.nodes[i].src - 1][uncovered_edge.nodes[i].dest[j]] += 1;
		}
	}
	//change the total weight if the total weight exceed the total
	int sum = 0;
	for (int i = 0; i < weight.size(); i++){
		for (auto j : weight[i]){
			sum += j.second;
		}
	}
	//if the sum greater 
	if (sum / a >= 100){
		for (int i = 0; i < weight.size(); i++){
			for (auto j : weight[i]){
				j.second *= 0.3;
			}
		}
	}
}

bool check_VC(struct Graph gra, unordered_map<int, int> cover){
	
	for (auto i : cover){
		for (int x = 0; x < gra.nodes[i.first - 1].dest.size();x++){
			for (int j = 0; j < gra.nodes[gra.nodes[i.first -1].dest[x]-1].dest.size();j++){
				if (gra.nodes[gra.nodes[i.first -1].dest[x] -1].dest[j] == i.first){
					gra.nodes[gra.nodes[i.first - 1].dest[x] - 1].dest.erase(gra.nodes[gra.nodes[i.first - 1].dest[x] - 1].dest.begin() + j);
				}
			}
		}
		gra.nodes[i.first - 1].dest.clear();
	}
	int result = 0;
	for (auto i : gra.nodes){
		result += i.dest.size();
	}
	return 0 == result;
}

void local_search(string graph_file,int cutoff,int seed){
	//initialize the VC and graph
	struct Graph gra = parseEdges(graph_file);
	string graph_name = SplitFilename(graph_file);
	//unordered_map<int, int> cover = VC_Greedy(gra);
	unordered_map<int, int> cover = ls2_vc_approx(gra);
	unordered_map<int,int> best_solution;
	struct Graph uncovered_edge(0, 0);
	vector< unordered_map<int, int> > edge_weight;
	initialize_edge_weight(gra, &edge_weight);
	srand(seed);
	clock_t start = clock();
	int iter = 0;

	//open io to output the trace file
	ofstream output_trace;
	string trace = "../Output/Trace/LS2/" + graph_name + "_LS2_" + to_string(cutoff)+"_" + to_string(seed) + ".trace";
	output_trace.open(trace, ios::app);
	while (float(clock() - start) / (float)CLOCKS_PER_SEC < cutoff){
	//while (iter<200){
		//if the result is vetex cover
		if (uncovered_edge.E == 0){
			//choose the vertex to remove
			//cout << check_VC(gra, cover)<<endl;
			output_trace<<float(clock() - start) / (float)CLOCKS_PER_SEC <<','<< cover.size()<<endl;
			best_solution = cover;
			struct Node a = compute_lowest_score(gra, cover, edge_weight);
			cover.erase(a.src);
			if (a.dest.size() != 0){
				uncovered_edge.nodes.push_back(a);
				uncovered_edge.V += 1;
				//cout << "uncovered_edge added "<<a.dest.size()<<endl;
				uncovered_edge.E += a.dest.size();
				//cout << "erased "<<a.src<<endl;
			}
		}
		//choose a node from it

		struct Node u = compute_lowest_score(gra, cover, edge_weight);
		//cout << "erased "<< u.src<< endl;
		cover.erase(u.src);

		if (u.dest.size() != 0){
			uncovered_edge.nodes.push_back(u);
			uncovered_edge.V += 1;
			//cout << " the size of added to the uncovered_edge" << u.dest.size()<<endl;
			uncovered_edge.E += u.dest.size();
		}
		else{
			//cout << "diminished by 1" << endl;
			continue;
		}
		//choose a entering vertex randomly, then add it to the cover, remove edges from the uncoverd edge
		add_vertex(uncovered_edge, cover, gra);
		// configure the edge weight, update the edge weight
		add_weight_uncoverd_edge(&edge_weight, uncovered_edge, gra.V);
		iter++;
	}
	output_trace.close();
	//output to the solution
	ofstream output;
	string output_file = "../Output/Solution/LS2/"+graph_name + "_LS2_"+to_string(cutoff)+"_"+to_string(seed)+".sol";
	output.open(output_file);
	// output the result into the file
	output<<best_solution.size()<<endl;
	for (auto i : best_solution){
		output << i.first << ",";
	}
	output.close();
}

unordered_map<int,int> ls2_vc_approx(struct Graph G) {

	unordered_map<int,int> VC;
	for (auto node_v : G.nodes) {
		int v = node_v.src;
		for (int u : node_v.dest) {
			if (u < v || G.nodes[u - 1].dest.size() == 0)
				continue;
			node_v.dest.clear();
			G.nodes[u - 1].dest.clear();
			VC.insert({v,v});
			VC.insert({u,u});
			break;
		}
	}
	return VC;
}

struct Node compute_lowest_score(struct Graph gra, unordered_map<int, int> cover, vector< unordered_map<int, int> > edgeweight){
	int result = 1000;
	int index = 0;
	struct Node a(0);
	for (auto i : cover) {
		int score = 0;
		vector<int> dest;
		for (int j = 0; j<gra.nodes[i.first - 1].dest.size(); j++){
			if (cover.find(gra.nodes[i.first - 1].dest[j]) == cover.end()){
				score += edgeweight[i.first - 1][gra.nodes[i.first - 1].dest[j]];
				dest.push_back(gra.nodes[i.first - 1].dest[j]);
			}
		}
		if (score < result) {
			a.dest = dest;
			a.src = i.first;
			result = score;
		}
		//cout<< dest.size()<<" edges added"<<endl;
	}
	return a;
}

////////////////////////////////////////////////////////////////////////////
//////function for LS2 end ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////
//
// Help functions to construct the graph
//
//////////////////////////////////////////////////////

void splitString(string &s, vector<int> &result){
	string ss;
	int index = 0;
	for (int i = 0; i<s.size(); i++) {
		if (s[i] == ' ') {
			ss = s.substr(index, i);
			index = i + 1;
			result.push_back(std::stoi(ss));
		}
	}
}

struct Graph parseEdges(string graph_file) {
	string path = graph_file;
	std::ifstream infile(path);
	int V, E, W;
	vector<int> input;
	string s;
	getline(infile, s);
	splitString(s, input);
	V = (int)input[0];
	E = (int)input[1];
	//W = (int)input[2];
	struct Graph gra(V,E);
	int index = 1;
	input.clear();
	while (!infile.eof()) {
		if (index > V) break;
		getline(infile, s);
		splitString(s, input);
		struct Node newNode = Node(index++);
		for (auto i : input) {
			newNode.dest.push_back(i);
		}
		gra.nodes.push_back(newNode);
		input.clear();
	}
	return gra;
}

string SplitFilename (string str)
{
  std::size_t first = str.find_last_of("/");
  std::size_t last = str.find_last_of(".");
  return str.substr(first + 1, last - 1 - first);
}
