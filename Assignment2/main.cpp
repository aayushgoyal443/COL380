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
    int truss;
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



int main( int argc, char** argv ){

    string inputpath = argv[1];
    string headerpath = argv[2];
    string outputpath = argv[3];

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
    for (auto it = edgeList.begin(); it!= edgeList.end(); it++){
        cout << it->first.first << " " << it->first.second << " " << it->second.sup << endl;
    }
    // remember to deallocate the memory
    
    MPI_Finalize();

}