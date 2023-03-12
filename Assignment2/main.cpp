#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

class Node{
public:
    int id;
    int degree;
    int rank;
    int owner;
    Node(){
        this->id = -1;
        this->degree = -1;
        this->rank = -1;
        this->owner = -1;
    }
    Node(int id, int degree){
        this->id = id;
        this->degree = degree;
        this->rank = -1;
        this->owner = -1;
    }

};

struct Query{
    int first, second;
    bool exists;
    Query(int a, int b){
        this->first = a;
        this->second = b;
        this->exists = false;
    }
};

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

    vector<Node> nodes(n);

    int rank, size;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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
            return a.degree < b.degree;
        });

        // assign ranks to the nodes
        for(int i=0; i<n; i++){
            nodes[i].rank = i;
            nodes[i].owner = i%size;
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

    // for(int i  = 0; i < n ; i ++){
    //     cout << rank << " " << nodes[i].id << " " << nodes[i].degree << " " << nodes[i].rank << " " << offsets[i] << endl;
    // }
    // for (int i=0; i<size; i++){
    //     // read the edges of the nodes assigned to this processor using the offsets
    //     for(auto it=processors[i]->adjlist.begin(); it!=processors[i]->adjlist.end(); it++){
    //         int node = it->first;
    //         int offset = offsets[node];
    //         input.seekg(offset+8);
    //         int degree = nodes[node]->degree;
    //         for(int j=0; j<degree; j++){
    //             int neighbour;
    //             input.read((char*)&neighbour, sizeof(int));
    //             if (nodes[neighbour]->rank > nodes[node]->rank){
    //                 processors[i]->adjlist[node].insert(neighbour);
    //             }
    //         }
    //     }
    // }

    

    // MPI_Barrier(MPI_COMM_WORLD);
    
    // vector<vector<Query>> queries(size);
    // for (auto u:processors[rank]->adjlist){
    //     int node = u.first;
    //     for (auto v=u.second.begin(); v!=u.second.end(); v++){
    //         for (auto w= next(v); w!=u.second.end(); w++){
    //             if (nodes[*v] -> rank < nodes[*w] -> rank){
    //                 // send query to processor with rank = nodes[*v]->owner
    //                 queries[nodes[*v]->owner].push_back(Query(*v, *w));
    //             }
    //             else{
    //                 // send query to processor with rank = nodes[*w]->owner
    //                 queries[nodes[*w]->owner].push_back(Query(*w, *v));
    //             }
    //         }
    //     }
    // }
    
    // // send queries to the processors
    // for (int i=0; i<size; i++){
    //     if (i == rank) continue;
    //     int size = queries[i].size();
    //     MPI_Send(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //     for (auto q:queries[i]){
    //         MPI_Send(&q.first, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //         MPI_Send(&q.second, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    //     }
    // }

    // // receive queries from the processors
    // for (int i=0; i<size; i++){
    //     if (i == rank) continue;
    //     int size;
    //     MPI_Recv(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     vector<bool> answers(size);
    //     for (int j=0; j<size; j++){
    //         int a, b;
    //         MPI_Recv(&a, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         MPI_Recv(&b, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         // check if the query exists
    //         if (processors[rank]->adjlist[a].find(b) != processors[rank]->adjlist[a].end()){
    //             // send the answer to the processor with rank = nodes[a]->owner
    //             answers[j] = true;
    //         }
    //         else{
    //             // send the answer to the processor with rank = nodes[b]->owner
    //             answers[j] = false;
    //         }
    //     }
    //     // send the answers to the processors
    //     MPI_Send(&answers[0], size, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);

    // }

    // // receive answers from the processors
    // for (int i=0; i<size; i++){
    //     if (i == rank) continue;
    //     int size;
    //     MPI_Recv(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     vector<bool> answers(size);
    //     MPI_Recv(&answers[0], size, MPI_C_BOOL, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     // update the queries
    //     for (int j=0; j<size; j++){
    //         if (queries[i][j].first == nodes[queries[i][j].first]->owner){
    //             queries[i][j].exists = answers[j];
    //         }
    //         else{
    //             queries[i][j].exists = answers[j];
    //         }
    //     }
    // }

    MPI_Finalize();

}