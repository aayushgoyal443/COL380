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


// void do_window_expansion(int k_max, pair<int,int> & window, int &gamma_active, int gamma_max, vector<pair<int, int>> & active_set, map<pair<int, int>, EdgeData>& edgeList, vector<vector<pair<int,int>>> &buckets){
//     int k = window.second;
//     double delta = 0.1;
//     while (gamma_active <= delta * gamma_max && k != k_max){
//         for (auto u: buckets[k+1]){
//             active_set.push_back(u);
//             gamma_active += edgeList[u].triangles.size();
//         }
//         k++;
//     }
//     window.second = k;
// }


// bool updateTriangle(pair<int, int> edge, int val_old, int val_new, map<pair<int, int>, EdgeData>& edgeList){
//     int t_num = edgeList[edge].truss_number;
//     if (val_old >= t_num && val_new >= t_num ){
//     }
//     else if (val_old >= t_num && val_new < t_num) {
//         edgeList[edge].g -= 1;
//         edgeList[edge].histogram[val_new]+=1;
//     }
//     else if (val_old < t_num && val_new < t_num) {
//         edgeList[edge].histogram[val_old]-=1;
//         edgeList[edge].histogram[val_new]+=1;
//     }
//     if (edgeList[edge].g < t_num-2){
//         edgeList[edge].truss_number-=1;
//         edgeList[edge].g += edgeList[edge].histogram[edgeList[edge].truss_number];   
//         return true;
//     }
//     else return false;
// }


// void update_triangles_and_histograms(int u, int v, int w, int val_old, int val_new, map<pair<int, int>, EdgeData>& edgeList, set<pair<int,int> > &changed_edges, vector<Node> &nodes){
    
//     if (nodes[u].rank > nodes[w].rank)
//     {
//         bool updated = updateTriangle(make_pair(w, u), val_old, val_new, edgeList);  
//         if (updated) changed_edges.insert(make_pair(w, u));
//     }
//     else{
//         bool updated = updateTriangle(make_pair(u, w), val_old, val_new, edgeList);
//         if (updated) changed_edges.insert(make_pair(u, w));
//     }
//     if (nodes[v].rank > nodes[w].rank)
//     {
//         bool updated = updateTriangle(make_pair(w, v), val_old, val_new, edgeList);   
//         if (updated) changed_edges.insert(make_pair(w, v));
//     }
//     else{ 
//         bool updated = updateTriangle(make_pair(v, w), val_old, val_new, edgeList);  
//         if (updated) changed_edges.insert(make_pair(v, w)); 
//     }

// }


int main( int argc, char** argv ){

    string inputpath = argv[1];
    string headerpath = argv[2];
    string outputpath = argv[3];
    int k_minimum = atoi(argv[4]);
    int k_maximum = atoi(argv[5]);

    ifstream input(inputpath, ios::binary);
    ifstream header(headerpath, ios::binary);

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
    
    for (int u = 0; u<n; u++ ){
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
    
    // // lets print the support for each edge in the edgeList
    // for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
    //     cout << it->first.first << " " << it->first.second << " " << it->second.sup << endl;
    // }

    // print all the triangles 
    // for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
    //     for (auto it1 = it->second.triangles.begin(); it1!= it->second.triangles.end(); it1++){
    //         cout << it->first.first << " " << it->first.second << " " << it1->first << endl;
    //     }
    // }

    // remember to deallocate the memory

    // initialize the truss values in the edgeList

    // int k_min = INT_MAX;
    // int k_max = INT_MIN;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        it->second.truss_number = it->second.sup + 2;
        // it->second.g = it->second.sup;
        // it->second.histogram = vector<int>(it->second.truss_number , 0);
        // k_max = max(k_max, it->second.truss_number);
        // k_min = min(k_min, it->second.truss_number);
    }    

    /*
    vector<vector<pair<int,int>>> buckets(k_max+1);
    pair<int,int> window = make_pair(k_min,k_min);

    vector<pair<int,int>> active_set;

    int gamma_active = 0;
    int gamma_max = 0;
    // let add elements to the buckets and also update the active set and gamma_active
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        buckets[it->second.truss_number].push_back(it->first);
        if (it->second.truss_number == k_min){
            active_set.push_back(it->first);
            gamma_active += it->second.sup;
        }
    }
    // // lets print the size of each bucket, size of the active set and gamma_active
    // for (int i=0; i<=k_max; i++){
    //     cout << "Bucket " << i << " size: " << buckets[i].size() << endl;
    // }
    // cout << "Active set size: " << active_set.size() << endl;
    // cout << "Gamma active: " << gamma_active << endl;
    // cout << "Gamma max: " << gamma_max << endl;
    // cout << "Window: " << window.first << " " << window.second << endl;
    // cout << "--------------------------------------------" << endl;

    while (true){
        
        if (active_set.size() == 0 and window.second == k_max){
            break;
        }
        do_window_expansion(k_max, window, gamma_active, gamma_max, active_set, edgeList, buckets);
        set<pair<int, int>> changed_edges;
        for(auto edge: active_set){
            for (auto &triangle: edgeList[edge].triangles){
                if (edgeList[edge].truss_number >= triangle.second){
                    continue;
                }
                int val_old = triangle.second;
                triangle.second = edgeList[edge].truss_number;
                int val_new = triangle.second;
                int u = edge.first;
                int v = edge.second;
                int w = triangle.first;
                update_triangles_and_histograms(u, v, w, val_old, val_new, edgeList, changed_edges, nodes);
            }
        }

        // update active set with the edges which were changed in the current iteration
        active_set.clear();
        gamma_active = 0;
        for (auto &edge: changed_edges){
            active_set.push_back(edge);  
            gamma_active += edgeList[edge].sup; 
        }
        gamma_max = max(gamma_max, gamma_active);
        // TODO: DEallocate memory to avoid segmentation fault
    }   

    // lets print the truss number of each edge in the edgeList
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        cout << it->first.first << " " << it->first.second << " " << it->second.truss_number << endl;
    }
    */

    // implement the minTruss algorithm
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

    // now rank 0 will have the trussList of all the edges
        // cout << "rank = " << rank << " size = " << recvList.size() << endl;
        // for (auto it = recvList.begin(); it!= recvList.end(); it++){
        //     cout << get<0>(*it) << " " << get<1>(*it) << " " << get<2>(*it) << endl;
        // }

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
            for (int i=0;i<n;i++){
                if (!visited[i]){
                    vector<int> component;
                    queue<int> q;
                    q.push(i);
                    visited[i] = true;
                    while(!q.empty()){
                        int u = q.front();
                        q.pop();
                        component.push_back(u);
                        for (auto v: adjList[u]){
                            if (!visited[v]){
                                q.push(v);
                                visited[v] = true;
                            }
                        }
                    }
                    if (component.size()>1) connectedComponents.push_back(component);
                }
            }
            // now connectedComponents has all the connected components
            // print this connectedComponents
            // cout << "k = " << k << endl;
            cout <<"1\n";
            // cout << "Number of connected components = " << connectedComponents.size() << endl;
            cout << connectedComponents.size() << endl;
            for (auto it = connectedComponents.begin(); it!= connectedComponents.end(); it++){
                sort(it->begin(), it->end());
                for (auto it2 = it->begin(); it2!= it->end(); it2++){
                    cout << *it2 << " ";
                }
                cout << endl;
            }

        }
        else {
            // cout << k << endl;
            // cout << "Number of connected components = 0" << endl;
            cout <<"0\n";
        }
    }

}
