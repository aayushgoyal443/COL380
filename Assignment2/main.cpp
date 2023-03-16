#include <bits/stdc++.h>
#include <unistd.h>
#include <mpi.h>
using namespace std;


struct hashFunction
{
  size_t operator()(const pair<int , 
                    int> &x) const
  {
    return x.first ^ x.second;
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


vector<int> print_connected_components(int source, int rank, int size, vector<Node> &nodes, vector<int>& visited, vector<vector<int>>& adjList)
{   
    vector<int> component;
    vector<int> frontier, next_frontier;
    int owner = nodes[source].owner;
    if (owner == rank)
    {
        frontier.push_back(source);
        visited[source] = 1;
        component.push_back(source);
    }
    while (true)
    {
        next_frontier.clear();
        vector<vector<int>> sendUpdates(size);
        vector<vector<int>> recvUpdates(size);
        vector<int> sendCounts(size, 0);
        vector<int> recvCounts(size, 0);
        vector<int> sendOffsets(size, 0);
        vector<int> recvOffsets(size, 0);

        for (int i=0; i<frontier.size(); i++)
        {
            int cur = frontier[i];
            for(auto u: adjList[cur])
            {
                if (visited[u] == 0)
                {
                    int owner = nodes[u].owner;
                    if (owner == rank)
                    {
                        next_frontier.push_back(u);
                        visited[u] = 1;
                        component.push_back(u);
                    }
                    else{
                        sendUpdates[owner].push_back(u);
                        sendCounts[owner]++;
                    }
                }
            }
        }

        MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i=0; i<size; i++)
        {
            sendOffsets[i] = (sendOffsets[i-1] + sendCounts[i-1]);
            recvOffsets[i] = (recvOffsets[i-1] + recvCounts[i-1]);
        }
        vector<int> sendBuffer(sendOffsets[size-1] + sendCounts[size-1]);
        vector<int> recvBuffer(recvOffsets[size-1] + recvCounts[size-1]);
        for (int i=0; i<size; i++)
        {
            for (int j=0; j<sendCounts[i]; j++)
            {
                sendBuffer[sendOffsets[i] + j] = sendUpdates[i][j];
            }
        }
        MPI_Alltoallv(sendBuffer.data(), sendCounts.data(), sendOffsets.data(), MPI_INT, recvBuffer.data(), recvCounts.data(), recvOffsets.data(), MPI_INT, MPI_COMM_WORLD);
        for (int i=0; i<size; i++)
        {
            for (int j=0; j<recvCounts[i]; j++)
            {
                int u = recvBuffer[recvOffsets[i] + j];
                if (visited[u] == 0)
                {
                    next_frontier.push_back(u);
                    visited[u] = 1;
                    component.push_back(u);
                }
            }
        }
        frontier = next_frontier;
        int fsize = frontier.size();
        int global_fsize = 0;
        MPI_Allreduce(&fsize, &global_fsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (global_fsize == 0) break;
    }
    return component;
}


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
    // cout << "owner of 4492 is " << nodes[4492].owner << endl;
    // cout << "owner of 4490 is " << nodes[4492].owner << endl;
    
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

    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        it->second.truss_number = it->second.sup+2;
    }

    int edge_list_count = edgeList.size();
    MPI_Allreduce(MPI_IN_PLACE, &edge_list_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    unordered_set<pair<int, int>, hashFunction> settled;
    int iter = 0;
    while(true){
        iter++;
        int count = settled.size();
        MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (count == edge_list_count){
            break;
        }

        int min_truss = INT_MAX;
        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
            if (settled.find(it->first) == settled.end()){
                if (it->second.truss_number < min_truss){
                    min_truss = it->second.truss_number;
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &min_truss, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
        set<pair<int,int>> to_settle;

        for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
            if ( settled.count(it->first)==0 && edgeList[it->first].truss_number == min_truss){
                to_settle.insert(it->first);
            }
        }
        int to_settle_count = to_settle.size();
        int global_to_settle_count=INT_MIN;
        MPI_Allreduce(&to_settle_count, &global_to_settle_count, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        auto it = to_settle.begin();
        for (int i=0;i<global_to_settle_count;i++){
            pair<int,int> edge = {-1,-1};
            int u,v;
            vector<vector<EdgeQuery>> sendEdgeQuery(size);
            
            if (it != to_settle.end()){
                edge = *it;
                it++;
                settled.insert(edge);
                u = edge.first;
                v = edge.second;

                for (auto& triangle: edgeList[edge].triangles){

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
                        bool settled_uu_ww = settled.count(make_pair(uu,ww));   
                        bool settled_vvv_www = settled.count(make_pair(vvv,www));
                        if (!settled_uu_ww && !settled_vvv_www){
                            if (edgeList[make_pair(uu,ww)].truss_number > min_truss) edgeList[make_pair(uu,ww)].truss_number-=1;
                            if (edgeList[make_pair(vvv,www)].truss_number > min_truss) edgeList[make_pair(vvv,www)].truss_number-=1;
                        }                   
                    }
                    else if (owner_uu != rank && owner_vvv != rank){
                        sendEdgeQuery[owner_uu].push_back(EdgeQuery({uu,ww},{vvv,www}, 0, true ));
                    }
                    else if (owner_uu != rank) sendEdgeQuery[owner_uu].push_back(EdgeQuery({uu,ww},{vvv,www}, (owner_vvv== rank)? settled.count(make_pair(vvv,www)) : 1 ,  false ));
                    else if (owner_vvv != rank) sendEdgeQuery[owner_vvv].push_back(EdgeQuery({vvv,www},{uu,ww}, (owner_uu== rank)? settled.count(make_pair(uu,ww)) : 1, false ));

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

    }

    // // print the truss number of each edge
    // for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
    //     cout << it->first.first << " " << it->first.second << " " << it->second.truss_number << endl;
    // }
    
    // TODO: Now we just need to print the components directly
    // now make a trussList using the edgeList 
    vector<tuple<int,int,int>> trussList;
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        trussList.push_back(make_tuple(it->first.first, it->first.second, it->second.truss_number));
    }
    
    for (int k = k_minimum; k <= k_maximum; k++){
        int exists = 0;
        for (auto it = trussList.begin(); it!= trussList.end(); it++){
            if (get<2>(*it) >= k+2){
                exists = true;
                break;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &exists, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        vector<vector<int>> connectedComponents;
        if (exists){
            // find all the edges with truss number >= k+2
            vector<pair<int,int>> trussList_k;
            for (auto it = trussList.begin(); it!= trussList.end(); it++){
                if (get<2>(*it) >= k+2){
                    trussList_k.push_back({get<0>(*it), get<1>(*it)});
                }
            }

            vector<vector<int>> adjList(nodes.size());
            vector<vector<pair<int, int>>> sendTrussList(size);
            for (auto p : trussList_k){
                adjList[p.first].push_back(p.second);
                int other_owner = nodes[p.second].owner;
                if (other_owner != rank){
                    sendTrussList[other_owner].push_back({p.second, p.first});
                }
                else{
                    adjList[p.second].push_back(p.first);
                }
            }
            vector<int> sendCounts(size, 0);
            vector<int> recvCounts(size, 0);
            vector<int> sendOffsets(size, 0);
            vector<int> recvOffsets(size, 0);
            for (int i=0; i<size; i++)
            {
                sendCounts[i] = sendTrussList[i].size();
            }
            MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
            for (int i=0; i<size; i++)
            {
                sendOffsets[i] = (sendOffsets[i-1] + sendCounts[i-1]);
                recvOffsets[i] = (recvOffsets[i-1] + recvCounts[i-1]);
            }
            vector<pair<int, int>> sendTrussListBuffer(sendOffsets[size-1] + sendCounts[size-1]);
            vector<pair<int, int>> recvTrussListBuffer(recvOffsets[size-1] + recvCounts[size-1]);

            for (int i=0; i<size; i++)
            {
                for (int j=0; j<sendCounts[i]; j++)
                {
                    sendTrussListBuffer[sendOffsets[i] + j] = sendTrussList[i][j];
                }
            }

            for (int i = 0; i < size; i++)
            {
                sendCounts[i] *= 2;
                recvCounts[i] *= 2;
                sendOffsets[i] *= 2;
                recvOffsets[i] *= 2;
            }


            MPI_Alltoallv(sendTrussListBuffer.data(), sendCounts.data(), sendOffsets.data(), MPI_INT, recvTrussListBuffer.data(), recvCounts.data(), recvOffsets.data(), MPI_INT, MPI_COMM_WORLD);

            for (int i=0; i<recvTrussListBuffer.size(); i++)
            {
                adjList[recvTrussListBuffer[i].first].push_back(recvTrussListBuffer[i].second);
            }

            // now using this trussList_k, find the connected components using bfs 
            // and store the connected components in a vector of vectors
            vector<int> visited(n,0);
            for (int i=0;i<n;i++){
                if (!visited[i]){
                    vector<int> component = print_connected_components(i, rank, size, nodes, visited, adjList);

                    // now gather all the components from all the processors
                    // gather the number of components in each processor
                    int* numComponents = new int[size];
                    for(int i = 0; i < size; i ++)
                        numComponents[i] = 0;
                    int my_components = component.size();
                    MPI_Gather(&my_components, 1, MPI_INT, numComponents, 1, MPI_INT, 0, MPI_COMM_WORLD);

                    // gather the components from all the processors
                    int* recvComponentOffsets = new int[size];
                    recvComponentOffsets[0] = 0;
                    for (int i=1;i<size;i++){
                        recvComponentOffsets[i] = recvComponentOffsets[i-1] + numComponents[i-1];
                    }

                    vector<int> recvComponents(recvComponentOffsets[size-1] + numComponents[size-1]);

                    MPI_Gatherv(component.data(), my_components, MPI_INT, recvComponents.data(), numComponents, recvComponentOffsets, MPI_INT, 0, MPI_COMM_WORLD);
                    if (recvComponents.size()>1) connectedComponents.push_back(recvComponents);
                    for (auto it = recvComponents.begin(); it!= recvComponents.end(); it++){
                        visited[*it] = 1;
                    }
                    // now broadcast this visited to all processors
                    MPI_Bcast(visited.data(), n, MPI_INT, 0, MPI_COMM_WORLD);
                }
            }
        }
        if (rank ==0){
            if (exists){
                cout <<"1\n";
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
                cout << "0\n";
            }
        }
    }
    

    MPI_Finalize();
}