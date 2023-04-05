#include <bits/stdc++.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>
#include <chrono>
//[]#include <ctime>
using namespace std;


// Implement a DSU class and write a function which takes an edge list and prints the connected components of the graph

class DSU{
public:
    vector<int> parent;
    vector<int> rank;
    vector<int> size;
    int num_components;

    DSU(int n){
        this->parent = vector<int>(n);
        this->size = vector<int>(n, 1);
        this->num_components = n;
        for (int i=0; i<n; i++){
            this->parent[i] = i;
        }
    }

    int find(int u){
        if (u == this->parent[u]){
            return u;
        }
        return this->parent[u] = find(this->parent[u]);
    }

    void unite(int u, int v){
        u = find(u);
        v = find(v);
        if (u != v){
            if (this->size[u] < this->size[v]){
                swap(u, v);
            }
            this->parent[v] = u;
            this->size[u] += this->size[v];
            this->num_components--;
        }
    }

    // take edgeList as input and complete the parent array
    void complete_parent(pair<int, int>* edgeList, int arrSize){
        for (int i=0; i<arrSize; i++){
            int u = edgeList[i].first;
            int v = edgeList[i].second;
            unite(u, v);
        }
    }

    void print_connected_components(ofstream& out){
        // cout << this->num_components << endl;
        map<int,vector<int>> components;
        for (int i=0; i<this->parent.size(); i++){
            if (size[find(i)]> 1) components[find(i)].push_back(i);
        }
        out << components.size() << endl;
        for (auto component: components){
            for (auto node: component.second){
                out << node << " ";
            }
            out << endl;
        }
    }


};



class EdgeData{
public:
    // we are going to assume u is the smaller one
    int sup;
    int truss_number;
    vector<pair<int,int>> triangles;
};


class Node{
public:
    int id;
    int degree;
    int rank;
    int owner;
    set<int> adjlist;

    Node(){
        this->id = -1;
        this->degree = -1;
        this->rank = -1;
        this->owner = -1;
        adjlist = set<int>();
    }
    Node(int id, int degree){
        this->id = id;
        this->degree = degree;
        this->rank = -1;
        this->owner = -1;
        adjlist = set<int>();
    }

    // make the default comparator for the priority queue, which uses rank 
    bool operator<(const Node& other) const{
        return this->rank < other.rank;
    }

};


struct Query{
    int u, v, w;
    //default constructor
    Query(){
        this->u = -1;
        this->v = -1;
        this->w = -1;
    }
    Query(int a, int b, int c){
        this->u = a;
        this->v = b;
        this->w = c;
    }
};


struct EdgeQuery{
    pair<int, int> edge_to_change;
    pair<int, int> other_edge;
    int settled;
    int tag;

    //add default constructor
    EdgeQuery(){
        this->edge_to_change = make_pair(-1, -1);
        this->other_edge = make_pair(-1, -1);
        this->settled = -1;
        this->tag = -1;
    }
    EdgeQuery(pair<int, int> a, pair<int, int> b, int c, int d){
        this->edge_to_change = a;
        this->other_edge = b;
        this->settled = c;
        this->tag = d;
    }
};


void insertTriangle(int u, int v, int w, map<pair<int, int>, EdgeData>& edgeList){
    if (edgeList.find(make_pair(u, v)) == edgeList.end()){
        edgeList[make_pair(u, v)] = EdgeData();
        edgeList[make_pair(u, v)].sup = 1;
        edgeList[make_pair(u, v)].triangles.push_back({w, INT_MAX});
    }
    else{
        edgeList[make_pair(u, v)].sup += 1;
        edgeList[make_pair(u, v)].triangles.push_back({w, INT_MAX});
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

void print_time(string s, chrono::high_resolution_clock::time_point& start, chrono::high_resolution_clock::time_point& end, int rank){
    // end = chrono::high_resolution_clock::now();
    // // print the time duration in seconds
    // cout <<"[" << rank << "]" << s << " " << chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-6 << " mseconds" << endl;
    // start = chrono::high_resolution_clock::now();
}

int main( int argc, char** argv ){

    // initialize begin time high resolution clock
    auto start_init = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto end = chrono::high_resolution_clock::now();
    
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
    int rank, size;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    vector<Node> nodes(n);
    print_time("init", start, end, rank );
    if (rank == 0){
        int sum=0;
        // [OMP_PARALLEL_FOR]
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
        print_time("read", start, end, rank );
        // sort the nodes in order of their degree
        sort(nodes.begin(), nodes.end(), [](Node a, Node b){
            if (a.degree == b.degree){
                return a.id < b.id;
            }
            return a.degree < b.degree;
        });

        // assign ranks to the nodes
        // [OMP_PARALLEL_FOR]
        for(int i=0; i<n; i++){
            nodes[i].rank = i;
            nodes[i].owner = (i%size);
        }

        // sort the nodes in order of their id 
        sort(nodes.begin(), nodes.end(), [](Node a, Node b){
            return a.id < b.id;
        });
        print_time("sort", start, end, rank );
        // now they are back in the oringinal order

    }

    // broadcast the nodes to all the processors
    MPI_Bcast(nodes.data(), n * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    // broadcast the offsets to all the processors
    MPI_Bcast(offsets, n, MPI_INT, 0, MPI_COMM_WORLD);
    print_time("broadcast", start, end, rank );     
    for (int u = 0; u<n; u++ ){
        if (nodes[u].owner!= rank) continue;
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
    print_time("readlist", start, end, rank );
    map<pair<int, int>, EdgeData> edgeList;
    
    vector<vector<Query>> sendQueries(size);
    for (int u = 0; u<n; u++ ){ // OMP_PARALLEL_FOR
        if (nodes[u].owner!= rank) continue;
        for (auto it1 = nodes[u].adjlist.begin(); it1!= nodes[u].adjlist.end(); it1++){
            for (auto it2 = next(it1); it2!= nodes[u].adjlist.end(); it2++){
                int v = *it1;
                int w = *it2;   
                if (nodes[w] < nodes[v]){
                    swap(v, w);
                }
                int owner = nodes[v].owner;
                if (owner == rank){
                    // if the owner is the same as the current processor, then we can answer the query
                    // and we don't need to send it to anyone
                    if (nodes[v].adjlist.find(w) != nodes[v].adjlist.end()){
                        // cout << u << " " << v << " " << w << endl;
                        insertTriangle(u, v, w, edgeList); // [OMP_PARALLEL_TASK]
                        insertTriangle(v, w, u, edgeList);
                        insertTriangle(u, w, v, edgeList);
                    }
                }
                else{
                    sendQueries[owner].push_back(Query(u, v, w));
                }
            }
        }
    }
    // cout << "owner of 4492 is " << nodes[4492].owner << endl;
    // cout << "owner of 4490 is " << nodes[4492].owner << endl;
    print_time("query_buffer", start, end, rank );
    // u we need to find out how many queries we are going to send to each processor
    int* sendCounts = new int[size];
    int* sendOffsets = new int[size];
    int* recvCounts = new int[size];
    int* recvOffsets = new int[size];
    for(int i=0; i<size; i++){
        sendCounts[i] = sendQueries[i].size();
        sendOffsets[i] = 0;
        recvCounts[i] = 0;
        recvOffsets[i] = 0;
    }

    // now we need to find out how many queries we are going to receive from each processor
    MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);


    // now we need to find out the offsets for the send and receive buffers
    for(int i=1; i<size; i++){
        sendOffsets[i] = sendOffsets[i-1] + sendCounts[i-1];
        recvOffsets[i] = recvOffsets[i-1] + recvCounts[i-1];
    }

    // now we need to create the send and receive buffers
    vector<Query> sendBuffer = vector<Query>(sendOffsets[size-1] + sendCounts[size-1]);
    vector<Query> recvBuffer = vector<Query>(recvOffsets[size-1] + recvCounts[size-1]);
    
    // now we need to copy the queries into the send buffer
    for(int i=0; i<size; i++){ // [OMP_PARALLEL_FOR]
        for(int j=0; j<sendCounts[i]; j++){
            sendBuffer[sendOffsets[i] + j] = sendQueries[i][j];
        }
    }

    // multiply them by sizeof(Query) to get the number of bytes
    for (int i=0;i<size;i++){
        sendCounts[i] *= sizeof(Query);
        sendOffsets[i] *= sizeof(Query);
        recvCounts[i] *= sizeof(Query);
        recvOffsets[i] *= sizeof(Query);
    }


    // now we need to send the queries to the appropriate processors
    MPI_Alltoallv(sendBuffer.data(), sendCounts, sendOffsets, MPI_BYTE, recvBuffer.data(), recvCounts, recvOffsets, MPI_BYTE, MPI_COMM_WORLD);

    // divide the values by sizeof(Query)
    for (int i=0;i<size;i++){
        sendCounts[i] /= sizeof(Query);
        sendOffsets[i] /= sizeof(Query);
        recvCounts[i] /= sizeof(Query);
        recvOffsets[i] /= sizeof(Query);
    }

    // now we need to create the send and receive buffers
    char* sendBuffer2 = new char[recvOffsets[size-1] + recvCounts[size-1]];
    char* recvBuffer2 = new char[sendOffsets[size-1] + sendCounts[size-1]];

    // now we need to answer the queries and copy the answers into the send buffer
    for(int i=0; i<size; i++){
        for(int j=0; j<recvCounts[i]; j++){ // [OMP_PARALLEL_FOR] Maybe
            Query q = recvBuffer[recvOffsets[i] + j];
            if (nodes[q.v].adjlist.find(q.w) != nodes[q.v].adjlist.end()){
                sendBuffer2[recvOffsets[i] + j] = '1';
                insertTriangle(q.v, q.w, q.u, edgeList);   
            }
            else{
                sendBuffer2[recvOffsets[i] + j] = '0';
            }
        }
    }

    print_time("query_answer", start, end, rank );
    // // now we need to send the answers to the appropriate processors
    MPI_Alltoallv(sendBuffer2, recvCounts, recvOffsets, MPI_CHAR, recvBuffer2, sendCounts, sendOffsets, MPI_CHAR, MPI_COMM_WORLD);

    //now we have received the answers from the other processors
    //we need to insert the triangles into the edgeList
    for(int i=0; i<size; i++){
        for(int j=0; j<sendCounts[i]; j++){ // [OMP_PARALLEL_FOR] Maybe
            if (recvBuffer2[sendOffsets[i] + j]=='1'){
                Query q = sendBuffer[sendOffsets[i] + j];
                insertTriangle(q.u, q.v, q.w, edgeList);
                insertTriangle(q.u, q.w, q.v, edgeList);
            }
        }
    }
    // test 1,2,3 -> 8 sec
    print_time("query_insert", start, end, rank );
   // OMP_PARALLEL_FOR
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        it->second.truss_number = it->second.sup+2;
    }
    print_time("truss_init", start, end, rank );
    int edge_list_count = edgeList.size();
    MPI_Allreduce(MPI_IN_PLACE, &edge_list_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    print_time("edge_list_count", start, end, rank );
    set<pair<int, int>> settled;
    int iter = 0;
    while(true){
        // if (iter == 1) break;
        iter++;

        int count = settled.size();
        // print_time("while loop iter = " + to_string(iter), start, end, rank);
        MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (count == edge_list_count){
            break;
        }

        int min_truss = INT_MAX; // OMP_PARALLEL_FOR        
        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){            
            if (settled.find(it->first) == settled.end()){
                if (it->second.truss_number < min_truss){
                    min_truss = it->second.truss_number;
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &min_truss, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        set<pair<int,int>> to_settle;
        /// OMP_PARALLEL_FOR
        // #pragma omp parallel for
        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
            if ( settled.count(it->first)==0 && edgeList[it->first].truss_number == min_truss){
                to_settle.insert(it->first);
            }
        }
        int to_settle_count = to_settle.size();
        int global_to_settle_count=INT_MIN;
        MPI_Allreduce(&to_settle_count, &global_to_settle_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        auto it = to_settle.begin();
        // print_time("before follop to_settle", start, end, rank);

        for (int i=0;i<global_to_settle_count;i++){
            pair<int,int> edge = {-1,-1};
            int u=-1,v=-1;
            vector<vector<EdgeQuery>> sendEdgeQuery(size);
            
            if (it != to_settle.end()){
                edge = *it;
                it++;
                settled.insert(edge);
                u = edge.first;
                v = edge.second;
                for (auto& triangle: edgeList[edge].triangles){ // OMP_PARALLEL_FOR MAybe

                    int w = triangle.first;
                    int uu=u, ww=w;
                    int vvv=v, www=w;
                    if (nodes[w] < nodes[u]){
                        uu = w;
                        ww = u;
                    }
                    if (nodes[w] < nodes[v]){
                        vvv = w;
                        www = v;
                    }
                    int owner_uu = nodes[uu].owner;
                    int owner_vvv = nodes[vvv].owner;

                    if (owner_uu == rank && owner_vvv == rank ){
                        string uu_ww = to_string(uu) + "_" + to_string(ww);
                        bool settled_uu_ww = settled.count(make_pair(uu,ww));   
                        bool settled_vvv_www = settled.count(make_pair(vvv,www));
                        // #pragma omp critical (uu_ww)
                        if (!settled_uu_ww && !settled_vvv_www){
                            if (edgeList[make_pair(uu,ww)].truss_number > min_truss){
                                edgeList[make_pair(uu,ww)].truss_number-=1;
                            }
                            if (edgeList[make_pair(vvv,www)].truss_number > min_truss){
                                edgeList[make_pair(vvv,www)].truss_number-=1;
                            }
                        }                   
                    }
                    else if (owner_uu != rank && owner_vvv != rank){
                        string owner_uu_name = to_string(owner_uu);
                        // #pragma omp critical (owner_uu_name)
                        sendEdgeQuery[owner_uu].push_back(EdgeQuery({uu,ww},{vvv,www}, 0, true ));
                    }
                    else if (owner_uu != rank){
                        string owner_uu_name = to_string(owner_uu);
                        // #pragma omp critical (owner_uu_name)
                        sendEdgeQuery[owner_uu].push_back(EdgeQuery({uu,ww},{vvv,www}, (owner_vvv== rank)? settled.count(make_pair(vvv,www)) : 1 ,  false ));
                    }
                    else if (owner_vvv != rank){
                        string owner_vvv_name = to_string(owner_vvv);
                        // #pragma omp critical (owner_vvv_name)
                        sendEdgeQuery[owner_vvv].push_back(EdgeQuery({vvv,www},{uu,ww}, (owner_uu== rank)? settled.count(make_pair(uu,ww)) : 1, false ));
                    }
                }
            }

            vector<pair<int,int>> send_main_edge(size, {u,v});
            vector<pair<int,int>> recv_main_edge(size);
            MPI_Alltoall( send_main_edge.data() , sizeof(pair<int,int>), MPI_BYTE, recv_main_edge.data(), sizeof(pair<int,int>), MPI_BYTE, MPI_COMM_WORLD);

            // now send the sendEdgeQuery to the appropriate processors
            int* sendCounts = new int[size];
            int* sendOffsets = new int[size];
            int* recvCounts = new int[size];
            int* recvOffsets = new int[size];
            for(int i = 0; i < size; i++){
                sendCounts[i] = sendEdgeQuery[i].size();
                sendOffsets[i] = 0;
                recvCounts[i] = 0;
                recvOffsets[i] = 0;
            }
            MPI_Alltoall(sendCounts, 1, MPI_INT, recvCounts, 1, MPI_INT, MPI_COMM_WORLD);
            for(int i = 1; i < size; i++){
                sendOffsets[i] = sendOffsets[i-1] + sendCounts[i-1];
                recvOffsets[i] = recvOffsets[i-1] + recvCounts[i-1];
            }
            int totalSend = sendOffsets[size-1] + sendCounts[size-1];
            int totalRecv = recvOffsets[size-1] + recvCounts[size-1];

            EdgeQuery* sendBuffer = new EdgeQuery[totalSend];
            EdgeQuery* recvBuffer = new EdgeQuery[totalRecv];

            for(int i = 0; i < size; i++){
                for(int j = 0; j < sendCounts[i]; j++){
                    sendBuffer[sendOffsets[i] + j] = sendEdgeQuery[i][j];
                }
            }

            // multiply by counts and offsets by sizeof
            for (int i = 0; i < size; i++){
                sendCounts[i] *= sizeof(EdgeQuery);
                sendOffsets[i] *= sizeof(EdgeQuery);
                recvCounts[i] *= sizeof(EdgeQuery);
                recvOffsets[i] *= sizeof(EdgeQuery);
            }

            MPI_Alltoallv(sendBuffer, sendCounts, sendOffsets, MPI_BYTE, recvBuffer, recvCounts, recvOffsets, MPI_BYTE, MPI_COMM_WORLD);

            // now divide by the size 
            for (int i = 0; i < size; i++){
                sendCounts[i] /= sizeof(EdgeQuery);
                sendOffsets[i] /= sizeof(EdgeQuery);
                recvCounts[i] /= sizeof(EdgeQuery);
                recvOffsets[i] /= sizeof(EdgeQuery);
            }

            int* sendResponseBuffer = new int[totalRecv];
            int* recvResponseBuffer = new int[totalSend];
            for(int i = 0; i < size; i++){
                for(int j = 0; j < recvCounts[i]; j++){
                    EdgeQuery q = recvBuffer[recvOffsets[i] + j];
                    int uu = q.edge_to_change.first;
                    int w = q.edge_to_change.second;
                    int vv = q.other_edge.first;
                    int ww = q.other_edge.second;
                    int other_owner = nodes[vv].owner;
                    int other_settled = q.settled;
                    int tag  = q.tag;
                    if (other_owner == rank){
                        other_settled = settled.count(make_pair(vv,ww));
                    }
                    int settled_me  = settled.count(make_pair(uu,w));
                    
                    if (!other_settled && !settled_me){
                        if (edgeList[{uu,w}].truss_number > min_truss) edgeList[{uu,w}].truss_number-=1;
                        if (tag && edgeList[{vv,ww}].truss_number > min_truss) edgeList[{vv,ww}].truss_number-=1;
                    }
                    sendResponseBuffer[recvOffsets[i] + j] = settled_me;
                    // if (u == uu && v == w && !lied ) 
                    // {      
                    //     lied = true;                
                    //     sendResponseBuffer[recvOffsets[i] + j] = false;
                    // }
                }
            }

            MPI_Alltoallv(sendResponseBuffer, recvCounts, recvOffsets, MPI_INT, recvResponseBuffer, sendCounts, sendOffsets, MPI_INT, MPI_COMM_WORLD);
            // now use thne recvResponseBuffer to update the other_edge in the sendEdgeQuery
            // print_time("[594]before updating the other edge", start, end, rank);
            for(int i = 0; i < size; i++){
                for(int j = 0; j < sendCounts[i]; j++){
                    EdgeQuery q = sendEdgeQuery[i][j];
                    int u = q.other_edge.first;
                    if (nodes[u].owner != rank) continue;
                    
                    int w = q.other_edge.second;
                    int vv  = q.edge_to_change.first;
                    int ww = q.edge_to_change.second;
                    if (nodes[vv].owner == rank) continue;
                    int other_owner = nodes[vv].owner;
                    int other_settled = recvResponseBuffer[sendOffsets[i] + j];
                    if (recv_main_edge[other_owner].first == vv && recv_main_edge[other_owner].second == ww){
                        other_settled = false;
                    }
                    int settled_me  = settled.count(make_pair(u,w));
                    if (!other_settled && !settled_me){
                        if (edgeList[{u,w}].truss_number > min_truss) edgeList[{u,w}].truss_number-=1;
                    }
                }
            }

        }
        // print_time("[594]after updating the other edge", start, end, rank);
    }
    print_time("query_send_m", start, end, rank );
            
    // now make a trussList using the edgeList 
    vector<tuple<int,int,int>> trussList;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        trussList.push_back(make_tuple(it->first.first, it->first.second, it->second.truss_number));
    }
    print_time("query_trussList_m", start, end, rank );
    for (int k = k_minimum; k <= k_maximum; k++){
        int exists = 0;
        print_time("loop start", start, end, rank);
        for (auto it = trussList.begin(); it!= trussList.end(); it++){
            if (get<2>(*it) >= k+2){
                exists = true;
                break;
            }
        }
        print_time("before allreduce", start, end, rank);
        MPI_Allreduce(MPI_IN_PLACE, &exists, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        print_time("after allreduce", start, end, rank);
        if (verbose ==0 && taskid==1 ){
            if (rank == 0){
                if (exists) output << "1 ";
                else output << "0 ";
            }
            continue;
        }
        if (exists == 0){
            if (rank == 0) output << "0\n";
            continue;
        }

        // find all the edges with truss number >= k+2
        vector<pair<int,int>> trussList_k;
        for (auto it = trussList.begin(); it!= trussList.end(); it++){
            if (get<2>(*it) >= k+2){
                trussList_k.push_back({get<0>(*it), get<1>(*it)});
            }
        }

        // receive the trussList_k from all other processes in the rank 0

        // gather the number of edges in each processor
        int* numEdges = new int[size];
        for(int i = 0; i < size; i ++)
            numEdges[i] = 0;
        int my_edges = trussList_k.size();
        MPI_Gather(&my_edges, 1, MPI_INT, numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // gather the edges from all the processors
        int* recvEdgeOffsets = new int[size];
        recvEdgeOffsets[0] = 0;
        for (int i=1;i<size;i++){
            recvEdgeOffsets[i] = recvEdgeOffsets[i-1] + numEdges[i-1];
        }
        int total_edges = recvEdgeOffsets[size-1] + numEdges[size-1];
        pair<int,int>* recvEdges = new pair<int,int>[total_edges];

        int pair_size = sizeof(pair<int,int>);
        // now multiply by pair_size
        for (int i=0;i<size;i++){
            numEdges[i] *= pair_size;
            recvEdgeOffsets[i] *= pair_size;
        }

        MPI_Gatherv(trussList_k.data(), pair_size*my_edges, MPI_BYTE, recvEdges, numEdges, recvEdgeOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);

        // divide 
        for (int i=0;i<size;i++){
            numEdges[i] /= pair_size;
            recvEdgeOffsets[i] /= pair_size;
        }

        // now rank 0 has all the edges in the recvEdges
        // now do DSU on the edges and find the connected components
        if (rank ==0){
            
            DSU dsu(n);
            dsu.complete_parent(recvEdges, total_edges);

            if (taskid ==1 && verbose==1){
                output << "1\n";
                dsu.print_connected_components(output);
                continue;
            }

            vector<int> influencer_vertices;
            set<int> components[n];
            for (int u = 0; u<n; u++ ){
                // read the edges of the nodes assigned to this processor using the offsets
                input.seekg(offsets[u]+8);
                int degree = nodes[u].degree;
                for(int j=0; j<degree; j++){
                    int neighbour;
                    input.read((char*)&neighbour, sizeof(int));
                    if (dsu.size[dsu.find(neighbour)]>1) components[u].insert(dsu.find(neighbour));
                }
                if ( components[u].size() >= p ) {
                    influencer_vertices.push_back(u);
                }
            }
            output << influencer_vertices.size() << endl;
            if (verbose ==0){
                for (auto it = influencer_vertices.begin(); it!= influencer_vertices.end(); it++){
                    output << *it << " ";
                }
                output << endl;
            }
            else{
                for (auto it = influencer_vertices.begin(); it!= influencer_vertices.end(); it++){
                    output << *it << endl;
                    for (int i=0;i<n;i++){
                        if (components[*it].count(dsu.find(i))){
                            output << i << " ";
                        }
                    }
                    output <<"\n";
                }
            }          
            print_time("query_influencer_m", start, end, rank );    
        }
    }

    MPI_Finalize();
    print_time("query_total_m", start_init, end, rank );
}