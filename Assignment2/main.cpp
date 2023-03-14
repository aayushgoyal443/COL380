#include <bits/stdc++.h>
#include <unistd.h>
#include <mpi.h>
using namespace std;


class EdgeData{
public:
    int u;
    int v;
    // we are going to assume u is the smaller one
    int sup;
    int truss_number;
    int g;
    vector<pair<int,int>> triangles;
    vector<int> histogram;
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


struct TriangleQuery{
    int u, v, val_old, val_new;
    //default constructor
    TriangleQuery()
    {
        this->u = -1;
        this->v = -1;
        this->val_old = -1;
        this->val_new = -1;
    }
    TriangleQuery(int a, int b, int c, int d){
        this->u = a;
        this->v = b;
        this->val_old = c;
        this->val_new = d;
    }
};


void insertTriangle(int u, int v, int w, map<pair<int, int>, EdgeData>& edgeList){
    if (edgeList.find(make_pair(u, v)) == edgeList.end()){
        edgeList[make_pair(u, v)] = EdgeData();
        edgeList[make_pair(u, v)].u = u;
        edgeList[make_pair(u, v)].v = v;
        edgeList[make_pair(u, v)].sup = 1;
        edgeList[make_pair(u, v)].triangles.push_back({w, INT_MAX});
    }
    else{
        edgeList[make_pair(u, v)].sup += 1;
        edgeList[make_pair(u, v)].triangles.push_back({w, INT_MAX});
    }
}


void do_window_expansion(int k_max, pair<int,int> & window, int &gamma_active, int gamma_max, vector<pair<int, int>> & active_set, map<pair<int, int>, EdgeData>& edgeList, vector<vector<pair<int,int>>> &buckets){
    int k = window.second;
    double delta = 0.1;
    while (gamma_active <= delta * gamma_max && k != k_max){
        for (auto u: buckets[k+1]){
            active_set.push_back(u);
            gamma_active += edgeList[u].triangles.size();
        }
        k++;
    }
    window.second = k;
}


bool updateTriangle(pair<int, int> edge, int val_old, int val_new, map<pair<int, int>, EdgeData>& edgeList){
    int t_num = edgeList[edge].truss_number;
    if (val_old >= t_num && val_new >= t_num ){
    }
    else if (val_old >= t_num && val_new < t_num) {
        edgeList[edge].g -= 1;
        edgeList[edge].histogram[val_new]+=1;
    }
    else if (val_old < t_num && val_new < t_num) {
        edgeList[edge].histogram[val_old]-=1;
        edgeList[edge].histogram[val_new]+=1;
    }
    if (edgeList[edge].g < t_num-2){
        edgeList[edge].truss_number-=1;
        edgeList[edge].g += edgeList[edge].histogram[edgeList[edge].truss_number];   
        return true;
    }
    else return false;
}


void update_triangles_and_histograms(int u, int v, int w, int val_old, int val_new, int rank,  int size, map<pair<int, int>, EdgeData>& edgeList, set<pair<int,int> > &changed_edges, vector<Node> &nodes, vector<vector<TriangleQuery>>& sendQueries){
    
    if (nodes[u].rank > nodes[w].rank)
    {
        int owner = nodes[w].owner;
        if (owner == rank)
        {
            bool updated = updateTriangle(make_pair(w, u), val_old, val_new, edgeList);  
            if (updated) changed_edges.insert(make_pair(w, u));
        }
        else{
            sendQueries[owner].push_back(TriangleQuery(w, u, val_old, val_new));
        }
    }
    else{
        bool updated = updateTriangle(make_pair(u, w), val_old, val_new, edgeList);
        if (updated) changed_edges.insert(make_pair(u, w));
    }
    if (nodes[v].rank > nodes[w].rank)
    {
        int owner = nodes[w].owner;
        if (owner == rank)
        {
            bool updated = updateTriangle(make_pair(w, v), val_old, val_new, edgeList);   
            if (updated) changed_edges.insert(make_pair(w, v));
        }
        else{
            sendQueries[owner].push_back(TriangleQuery(w, v, val_old, val_new));
        }
    }
    else{ 
        int owner = nodes[v].owner;
        if (owner == rank)
        {
            bool updated = updateTriangle(make_pair(v, w), val_old, val_new, edgeList);  
            if (updated) changed_edges.insert(make_pair(v, w)); 
        }
        else{
            sendQueries[owner].push_back(TriangleQuery(v, w, val_old, val_new));
        }
    }

}


// void print_connected_components(int source, int rank, int size, vector<Node> &nodes, map<pair<int, int>, EdgeData>& edgeList, vector<vector<TriangleQuery>>& sendQueries, vector<vector<TriangleQuery>>& recvQueries, vector<vector<pair<int,int>>> &buckets, vector<pair<int, int>> & active_set, int k_max, int gamma_max, int &gamma_active, pair<int,int> & window, set<pair<int,int> > &changed_edges, vector<int> &trussness, vector<int> &trussness_new, vector<int> &trussness_old, vector<int> &trussness_old2, vector<int> &trussness_old3, vector<int> &trussness_old4, vector<int> &trussness_old5, vector<int> &trussness_old6, vector<int> &trussness_old7, vector<int> &trussness_old8, vector<int> &trussness_old9, vector<int> &trussness_old10, vector<int> &trussness_old11, vector<int> &trussness_old12, vector<int> &trussness_old13, vector<int> &trussness_old14, vector<int> &trussness_old15, vector<int> &trussness_old16, vector<int> &trussness_old17, vector<int> &trussness_old18, vector<int> &trussness_old19, vector<int> &trussness_old20, vector<int> &trussness_old21, vector<int> &trussness_old22, vector<int> &trussness_old23, vector<int> &trussness_old24, vector<int> &trussness_old25, vector<int> &trussness_old26, vector<int> &trussness_old27, vector<int> &trussness_old28, vector<int> &trussness_old29, vector<int> &trussness_old30, vector<int> &trussness_old31, vector<int> &trussness_old32, vector<int> &trussness_old33, vector<int> &trussness_old34, vector<int> &trussness_old35, vector<int> &trussness_old36, vector<int> &trussness_old37, vector<int> &trussness_old38, vector<int> &trussness_old39, vector<int>)
// {   
//     vector<int> frontier, next_frontier;
//     vector<int> visited(nodes.size(), 0);
//     int owner = nodes[source].owner;
//     if (owner == rank)
//     {
//         frontier.push_back(source);
//         visited[source] = 1;
//     }
//     while (true)
//     {
//         next_frontier.clear();
//         vector<vector<int>> sendUpdates(size);
//         vector<vector<int>> recvUpdates(size);
//         vector<int> sendCounts(size, 0);
//         vector<int> recvCounts(size, 0);
//         vector<int> sendOffsets(size, 0);
//         vector<int> recvOffsets(size, 0);
        
//         for (int i=0; i<frontier.size(); i++)
//         {
//             int cur = frontier[i];
//             for(auto u: nodes[cur].adj)
//             {
//                 if (visited[u] == 0)
//                 {
//                     int owner = nodes[u].owner;
//                     if (owner == rank)
//                     {
//                         next_frontier.push_back(u);
//                         visited[u] = 1;
//                     }
//                     else{
//                         sendUpdates[owner].push_back(u);
//                         sendCounts[owner]++;
//                     }
//                 }
//             }
            
//         }
//         MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
//         for (int i=0; i<size; i++)
//         {
//             sendOffsets[i] = (i==0)?0:(sendOffsets[i-1] + sendCounts[i-1]);
//             recvOffsets[i] = (i==0)?0:(recvOffsets[i-1] + recvCounts[i-1]);
//         }

//         for (int i=0; i<size; i++)
//         {
//             sendCounts[i] *= sizeof(int);
//             recvCounts[i] *= sizeof(int);
//             sendOffsets[i] *= sizeof(int);
//             recvOffsets[i] *= sizeof(int);
//         }
//         vector<int> sendBuffer(sendOffsets[size-1] + sendCounts[size-1]);
//         vector<int> recvBuffer(recvOffsets[size-1] + recvCounts[size-1]);
//         for (int i=0; i<size; i++)
//         {
//             memcpy(sendBuffer.data() + sendOffsets[i], sendUpdates[i].data(), sendCounts[i]);
//         }

//         MPI_Alltoallv(sendBuffer.data(), sendCounts.data(), sendOffsets.data(), MPI_BYTE, recvBuffer.data(), recvCounts.data(), recvOffsets.data(), MPI_INT, MPI_COMM_WORLD);
//         for (int i=0; i<size; i++)
//         {
//             for (int j=0; j<recvCounts[i]/sizeof(int); j++)
//             {
//                 int u = recvBuffer[recvOffsets[i]/sizeof(int) + j];
//                 if (visited[u] == 0)
//                 {
//                     next_frontier.push_back(u);
//                     visited[u] = 1;
//                 }
//             }
//         }
//         frontier = next_frontier;
//         int fsize = frontier.size();
//         int global_fsize = 0;
//         MPI_Allreduce(&global_fsize, &fsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//         if (global_fsize == 0) break;
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
    int rank, size;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    vector<Node> nodes(n);

    if (rank == 0){
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
            nodes[i].owner = (i%size);
        }

        // sort the nodes in order of their id 
        sort(nodes.begin(), nodes.end(), [](Node a, Node b){
            return a.id < b.id;
        });

        // now they are back in the oringinal order

    }

    // broadcast the nodes to all the processors
    MPI_Bcast(nodes.data(), n * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    // broadcast the offsets to all the processors
    MPI_Bcast(offsets, n, MPI_INT, 0, MPI_COMM_WORLD);
            
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

    map<pair<int, int>, EdgeData> edgeList;
    
    vector<vector<Query>> sendQueries(size);
    for (int u = 0; u<n; u++ ){
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
                        insertTriangle(u, v, w, edgeList);
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
    for(int i=0; i<size; i++){
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
        for(int j=0; j<recvCounts[i]; j++){
            Query q = recvBuffer[recvOffsets[i] + j];
            // cout << "Node num: "<< q.v << " | Node id: "  << nodes[q.v].id  <<" | rank: " << rank << " | owner: " <<  nodes[q.v].owner<<"\n";
            if (nodes[q.v].adjlist.find(q.w) != nodes[q.v].adjlist.end()){
                sendBuffer2[recvOffsets[i] + j] = '1';
                insertTriangle(q.v, q.w, q.u, edgeList);   
            }
            else{
                sendBuffer2[recvOffsets[i] + j] = '0';
            }
        }
    }


    // // now we need to send the answers to the appropriate processors
    MPI_Alltoallv(sendBuffer2, recvCounts, recvOffsets, MPI_CHAR, recvBuffer2, sendCounts, sendOffsets, MPI_CHAR, MPI_COMM_WORLD);

    // // now we have received the answers from the other processors
    // // we need to insert the triangles into the edgeList
    for(int i=0; i<size; i++){
        for(int j=0; j<sendCounts[i]; j++){
            if (recvBuffer2[sendOffsets[i] + j]=='1'){
                Query q = sendBuffer[sendOffsets[i] + j];
                insertTriangle(q.u, q.v, q.w, edgeList);
                insertTriangle(q.u, q.w, q.v, edgeList);
            }
        }
    }
    
    // lets print the support for each edge in the edgeList
    // for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
    //     cout << it->first.first << " " << it->first.second << " " << it->second.sup << endl;
    // }
    // remember to deallocate the memory

    // initialize the truss values in the edgeList
    int k_min = INT_MAX;
    int k_max = INT_MIN;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        it->second.truss_number = it->second.sup + 2;
        it->second.g = it->second.sup;
        it->second.histogram = vector<int>(it->second.truss_number , 0);
        k_max = max(k_max, it->second.truss_number);
        k_min = min(k_min, it->second.truss_number);
    }    
    
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
    int iter = 0;
    int completed_count = 0;
    while (completed_count < size){
        iter++;
        int my_status = 0;
        if (active_set.size() == 0 and window.second == k_max){
            my_status = 1;
        }
        do_window_expansion(k_max, window, gamma_active, gamma_max, active_set, edgeList, buckets);
        set<pair<int, int>> changed_edges;
        vector<vector<TriangleQuery>> sendQueries(size);
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
                update_triangles_and_histograms(u, v, w, val_old, val_new,rank, size, edgeList, changed_edges, nodes, sendQueries);
            }
        }

        int* numSendQueries = new int[size];
        int* numRecvQueries = new int[size];
        int* sendOffsets = new int[size];
        int* recvOffsets = new int[size];
        for (int i=0;i<size;i++){
            numSendQueries[i] = sendQueries[i].size();
            numRecvQueries[i] = 0;
            sendOffsets[i] = 0;
            recvOffsets[i] = 0;
        }

        MPI_Alltoall(numSendQueries, 1, MPI_INT, numRecvQueries, 1, MPI_INT, MPI_COMM_WORLD);

        for (int i=1;i<size;i++){
            sendOffsets[i] = sendOffsets[i-1] + numSendQueries[i-1];
            recvOffsets[i] = recvOffsets[i-1] + numRecvQueries[i-1];
        }

        vector<TriangleQuery> sendBuffer(sendOffsets[size-1] + numSendQueries[size-1]);
        vector<TriangleQuery> recvBuffer(recvOffsets[size-1] + numRecvQueries[size-1]);

        for (int i=0;i<size;i++){
            for (int j=0;j<numSendQueries[i];j++){
                sendBuffer[sendOffsets[i] + j] = sendQueries[i][j];
            }
        }

        // multiply the values by size of datatype 
        for (int i=0;i<size;i++){
            numSendQueries[i] *= sizeof(TriangleQuery);
            numRecvQueries[i] *= sizeof(TriangleQuery);
            sendOffsets[i] *= sizeof(TriangleQuery);
            recvOffsets[i] *= sizeof(TriangleQuery);
        }

        MPI_Alltoallv(sendBuffer.data(), numSendQueries, sendOffsets, MPI_BYTE, recvBuffer.data(), numRecvQueries, recvOffsets, MPI_BYTE, MPI_COMM_WORLD);

        // divide the values by size of datatype
        for (int i=0;i<size;i++){
            numSendQueries[i] /= sizeof(TriangleQuery);
            numRecvQueries[i] /= sizeof(TriangleQuery);
            sendOffsets[i] /= sizeof(TriangleQuery);
            recvOffsets[i] /= sizeof(TriangleQuery);
        }

        for (int i=0;i<size;i++){
            for (int j=0;j<numRecvQueries[i];j++){
                TriangleQuery q = recvBuffer[recvOffsets[i] + j];
                updateTriangle({q.u, q.v}, q.val_old, q.val_new, edgeList);
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
        MPI_Allreduce(&my_status, &completed_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        // TODO: DEallocate memory to avoid segmentation fault
    }   
    
    
    // now make a trussList using the edgeList 
    vector<tuple<int,int,int>> trussList;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        trussList.push_back(make_tuple(it->first.first, it->first.second, it->second.truss_number));
    }

    // now gather this trussList from all the processors and gather in rank 0 using MPI_Gatherv
    // cout << "rank = " << rank << " size = " << trussList.size() << endl;
    // gather the number of edges in each processor
    int* numEdges = new int[size];
    for(int i = 0; i < size; i ++)
        numEdges[i] = 0;
    int my_edges = trussList.size();
    MPI_Gather(&my_edges, 1, MPI_INT, numEdges, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // gather the trussList from all the processors
    int* recvTrussOffsets = new int[size];
    recvTrussOffsets[0] = 0;
    for (int i=1;i<size;i++){
        recvTrussOffsets[i] = recvTrussOffsets[i-1] + numEdges[i-1];
    }

    int tuple_size = sizeof(tuple<int,int,int>);
    vector<tuple<int,int,int>> recvList(recvTrussOffsets[size-1] + numEdges[size-1]);

    // multiply the values by size of datatype
    for (int i=0;i<size;i++){
        numEdges[i] *= tuple_size;
        recvTrussOffsets[i] *= tuple_size;
    }

    MPI_Gatherv(trussList.data(), my_edges*tuple_size, MPI_BYTE, recvList.data(), numEdges, recvTrussOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);

    // now rank 0 will have the trussList of all the edges
    if (rank == 0){
        // cout << "rank = " << rank << " size = " << recvList.size() << endl;
        // for (auto it = recvList.begin(); it!= recvList.end(); it++){
        //     cout << get<0>(*it) << " " << get<1>(*it) << " " << get<2>(*it) << endl;
        // }

        for (int k = k_minimum; k <= k_maximum; k++){
            bool exists = false;
            for (auto it = recvList.begin(); it!= recvList.end(); it++){
                if (get<2>(*it) >= k+2){
                    exists = true;
                    break;
                }
            }
            if (exists){
                // find all the edges with truss number >= k+2
                vector<tuple<int,int,int>> trussList_k;
                for (auto it = recvList.begin(); it!= recvList.end(); it++){
                    if (get<2>(*it) >= k+2){
                        trussList_k.push_back(*it);
                    }
                }
                cout << "k = " << k << " size = " << trussList_k.size() << endl;
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
                        connectedComponents.push_back(component);
                    }
                }
                // now connectedComponents has all the connected components
                // print this connectedComponents
                cout << "k = " << k << endl;
                cout << "Number of connected components = " << connectedComponents.size() << endl;
                for (auto it = connectedComponents.begin(); it!= connectedComponents.end(); it++){
                    sort(it->begin(), it->end());
                    for (auto it2 = it->begin(); it2!= it->end(); it2++){
                        cout << *it2 << " ";
                    }
                    cout << endl;
                }

            }
            else {
                cout << "k = " << k << endl;
                cout << "Number of connected components = 0" << endl;
            }
        }

    }
    

    MPI_Finalize();
}