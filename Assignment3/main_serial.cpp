// write the same code as in main.cpp but this time in a serial manner 
#include <bits/stdc++.h>
#include <unistd.h>
using namespace std;


class EdgeData{
public:
    int u;
    int v;
    // we are going to assume u is the smaller one
    int sup;
    int truss_number;
    int g;
    set<pair<int,int>> triangles;
    vector<int> histogram;
};


class Node{
public:
    int id;
    int degree;
    int rank;
    set<int> adjlist;

    Node(){
        this->id = -1;
        this->degree = -1;
        this->rank = -1;
        adjlist = set<int>();
    }
    Node(int id, int degree){
        this->id = id;
        this->degree = degree;
        this->rank = -1;
        adjlist = set<int>();
    }

    // make the default comparator for the priority queue, which uses rank 
    bool operator<(const Node& other) const{
        return this->rank < other.rank;
    }

};


void insertTriangle(int u, int v, int w, map<pair<int, int>, EdgeData>& edgeList){
    if (edgeList.find(make_pair(u, v)) == edgeList.end()){
        edgeList[make_pair(u, v)] = EdgeData();
        edgeList[make_pair(u, v)].u = u;
        edgeList[make_pair(u, v)].v = v;
        edgeList[make_pair(u, v)].sup = 1;
        edgeList[make_pair(u, v)].triangles.insert({w, INT_MAX});
    }
    else{
        edgeList[make_pair(u, v)].sup += 1;
        edgeList[make_pair(u, v)].triangles.insert({w, INT_MAX});
    }
}


string result_parse(string& args, string finding){
    string result = "";
    int start = args.find(finding);
    if (start == string::npos){
        if (finding == "--p") return "0";
        else return result;
    }
    start += finding.size() + 1;
    for (int i=start;i<args.size();i++){
        if (args[i] == '-' && args[i+1] == '-'){
            break;
        }
        result += args[i];
    }
    return result;
}


int main( int argc, char** argv ){

    string args = "";
    for (int i = 1; i < argc; i++){
        args += argv[i];
    }

    int taskid = stoi(result_parse(args, "--taskid"));
    string inputpath = result_parse(args, "--inputpath");
    string headerpath = result_parse(args, "--headerpath");
    string outputpath = result_parse(args, "--outputpath");
    int k_minimum = stoi(result_parse(args, "--startk"));
    int k_maximum = stoi(result_parse(args, "--endk"));
    int verbose = stoi(result_parse(args, "--verbose"));
    int p = stoi(result_parse(args, "--p"));

    ifstream input(inputpath, ios::binary);
    ifstream header(headerpath, ios::binary);
    ofstream output;
    output.open(outputpath);

    int n, m;
    input.read((char*)&n, sizeof(int));
    input.read((char*)&m, sizeof(int));

    int* offsets = new int[n];
    header.read((char*)offsets, sizeof(int)*n);

    vector<Node> nodes(n);

    int sum=0;
    for(int i=0; i<n; i++){
        int degree;
        if (i==n-1){
            degree = 2*m-sum;
        }
        else{
            degree = (offsets[i+1]-offsets[i]-8)/4;
        }
        sum += degree;
        nodes[i].id = i;
        nodes[i].degree = degree;
    }

    // cout <<"i1\n";

    // sort the nodes in order of their degree
    sort(nodes.begin(), nodes.end(), [](Node a, Node b){
        if (a.degree == b.degree){
            return a.id < b.id;
        }
        return a.degree < b.degree;
    });

    // assign ranks to the nodes
    for(int i=0; i<n; i++){
        nodes[i].rank = i;
    }

    // sort the nodes in order of their id 
    sort(nodes.begin(), nodes.end(), [](Node a, Node b){
        return a.id < b.id;
    });

    // cout << "i2\n";
    for (int u = 0; u<n; u++ ){
        // read the edges of the nodes assigned to this processor using the offsets
        input.seekg(offsets[u]+8);
        int degree = nodes[u].degree;
        set<int> adjlist;
        for(int j=0; j<degree; j++){
            int neighbour;
            input.read((char*)&neighbour, sizeof(int));
            if (nodes[neighbour].rank > nodes[u].rank){
                adjlist.insert(neighbour);
            }
        }
        nodes[u].adjlist = adjlist;
    }

    map<pair<int, int>, EdgeData> edgeList;
    // cout << "i3\n";
    // cout << "n = " << n << "\n";
    for (int u = 0; u<n; u++ ){
        // if (u%500 == 0)cout << u << endl;
        for (auto it1 = nodes[u].adjlist.begin(); it1!= nodes[u].adjlist.end(); it1++){
            for (auto it2 = next(it1); it2!= nodes[u].adjlist.end(); it2++){
                int v = *it1;
                int w = *it2;   
                if (nodes[w] < nodes[v]){
                    swap(v, w);
                }

                if (nodes[v].adjlist.find(w) != nodes[v].adjlist.end()){
                    // cout << u << " " << v << " " << w << endl;
                    insertTriangle(u, v, w, edgeList);
                    insertTriangle(v, w, u, edgeList);
                    insertTriangle(u, w, v, edgeList);
                }
            }
        }
    }
    // cout << "i4\n";
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        it->second.truss_number = it->second.sup + 2;
    }    

    cout << "Triangle Enumeration done ...\n";

    set<pair<int,int>> settled;

    while(settled.size() < edgeList.size()){
        // find the edge with minimum truss number
        int min_truss = INT_MAX;
        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
            if (settled.find(it->first) == settled.end()){
                if (it->second.truss_number < min_truss){
                    min_truss = it->second.truss_number;
                }
            }
        }
        set<pair<int,int>> to_settle;
        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
            if ( settled.count(it->first)==0 && edgeList[it->first].truss_number == min_truss){
                to_settle.insert(it->first);
            }
        }

        for (auto edge: to_settle){
            settled.insert(edge);
	    int u = edge.first;
            int v = edge.second;
            for (auto &triangle: edgeList[edge].triangles){
                int w = triangle.first;
                int uu, ww;
                int vvv, www;
                if (nodes[w] < nodes[u]){
                    uu = w;
                    ww = u;
                }
                else{
                    uu = u;
                    ww = w;
                }
                if (nodes[w] < nodes[v]){
                    vvv = w;
                    www = v;
                }
                else{
                    vvv = v;
                    www = w;
                }
                bool settled_uu_ww = settled.count(make_pair(uu,ww));   
                bool settled_vvv_www = settled.count(make_pair(vvv,www));

                if (!settled_uu_ww && !settled_vvv_www){
                    if (edgeList[make_pair(uu,ww)].truss_number > min_truss) edgeList[make_pair(uu,ww)].truss_number-=1;
                    if (edgeList[make_pair(vvv,www)].truss_number > min_truss)edgeList[make_pair(vvv,www)].truss_number-=1;
                }
            }
        }
    }

    // lets print the truss number of each edge in the edgeList
    // for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
    //     cout << it->first.first << " " << it->first.second << " " << it->second.truss_number << endl;
    // }

    
    // now make a trussList using the edgeList 
    vector<tuple<int,int,int>> trussList;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        trussList.push_back(make_tuple(it->first.first, it->first.second, it->second.truss_number));
    }


    cout << "Truss Decomposition done ..." << endl;


    for (int k = k_minimum; k <= k_maximum; k++){
        bool exists = false;
        for (auto it = trussList.begin(); it!= trussList.end(); it++){
            if (get<2>(*it) >= k+2){
                exists = true;
                break;
            }
        }
        if (exists){
            // find all the edges with truss number >= k+2
            vector<tuple<int,int,int>> trussList_k;
            for (auto it = trussList.begin(); it!= trussList.end(); it++){
                if (get<2>(*it) >= k+2){
                    trussList_k.push_back(*it);
                }
            }
            // cout << "k = " << k << " size = " << trussList_k.size() << endl;
            // now using this trussList_k, find the connected components using bfs 
            // and store the connected components in a vector of vectors
            vector<vector<int>> connectedComponents;
            vector<bool> visited(n, false);
            // make adjacency list of the trussList_k
            vector<vector<int>> adjList(n);
            for (auto it = trussList_k.begin(); it!= trussList_k.end(); it++){
                adjList[get<0>(*it)].push_back(get<1>(*it));
                adjList[get<1>(*it)].push_back(get<0>(*it));
            }
            // now do bfs on the adjList
            vector<int> color(n, -1);
            int cur_color = 0;
            for (int i=0;i<n;i++){
                if (!visited[i]){
                    vector<int> component;
                    queue<int> q;
                    q.push(i);
                    color[i] = cur_color;
                    visited[i] = true;
                    while(!q.empty()){
                        int u = q.front();
                        q.pop();
                        component.push_back(u);
                        color[u] = cur_color;
                        for (auto v: adjList[u]){
                            if (!visited[v]){
                                q.push(v);
                                visited[v] = true;
                            }
                        }
                    }
                    cur_color++;
                    if (component.size()>1) connectedComponents.push_back(component);
                    else{
                        color[component[0]] = -1;
                    }
                }
            }

            // for (auto i=0;i<n;i++){
            //     cout << color[i] << " ";
            // }cout << endl;

            vector<int> influencer_vertices;
            set<int> components[n];
            for(int i = 0; i < n; i ++)
            {   
                for(auto u:nodes[i].adjlist)
                {
                    if ( color[u]!=-1 ) components[i].insert(color[u]);
                    if ( color[i]!=-1 ) components[u].insert(color[i]);
                }
            }
            for(int i = 0; i < n; i ++)
            {
                if (components[i].size() >= p) influencer_vertices.push_back(i);
            }
            
            output << influencer_vertices.size() << endl;
            for (auto it = influencer_vertices.begin(); it!= influencer_vertices.end(); it++){
                output << *it << " ";
            } output << endl;
            
            // output <<"1\n";
            // output << connectedComponents.size() << endl;
            // for (auto it = connectedComponents.begin(); it!= connectedComponents.end(); it++){
            //     sort(it->begin(), it->end());
            //     for (auto it2 = it->begin(); it2!= it->end(); it2++){
            //         output << *it2 << " ";
            //     }
            //     output << endl;
            // }

        }
        else {
            // cout << k << endl;
            // cout << "Number of connected components = 0" << endl;
            cout <<"0\n";
        }
    }

}
