#include <bits/stdc++.h>
#include <omp.h>
#include "library.hpp"
using namespace std;

const int MAX_VAL = pow(2,16)-1;


class Block{
public:
    int i;
    int j;
    int** data;
    int m;

    Block(int i, int j, int** data, int m){
        this->i = i;
        this->j = j;
        this->data = data;
        this->m = m;
    }

    Block(int i, int j, int m){
        this->i = i;
        this->j = j;
        this->m = m;
        this->data = new int*[m];
        for (int a=0;a<m;a++){
            this->data[a] = new int[m];
        }
    }

    void show(){
        cout << this->i << " " << this->j << endl;
        for (int a=0;a<this->m;a++){
            for (int b=0;b<this->m;b++){
                cout << this->data[a][b] << " ";
            }
            cout << "\n";
        }
    }

    Block* inner(Block* b){
        Block* result = new Block(this->i, b->j, this->m);
        for (int p=0;p<this->m;p++){
            for (int q=0;q<this->m;q++){
                result->data[p][q] = Inner(this->data[p][0], b->data[0][q]);
                for (int k=1;k<this->m;k++){
                    int temp = Inner(this->data[p][k], b->data[k][q]);
                    if (temp > MAX_VAL){
                        temp = MAX_VAL;
                    }
                    result->data[p][q] = Outer( result->data[p][q] , temp);
                    if (result->data[p][q] > MAX_VAL){
                        result->data[p][q] = MAX_VAL;
                    }
                }
            }
        }
        return result;
    }

    void outer(Block* b){
        for (int p=0;p<this->m;p++){
            for (int q=0;q<this->m;q++){
                this->data[p][q] = Outer(this->data[p][q], b->data[p][q]);
                if (this->data[p][q] > MAX_VAL){
                    this->data[p][q] = MAX_VAL;
                }
            }
        }
    }

};


class SparseMatrixInput{

public:
    int n;
    int m;
    int k;
    vector<Block*> blocks;

    SparseMatrixInput(int n, int m, int k){
        this->n = n;
        this->m = m;
        this->k = k;
    }

    void insert(Block* block){
        this->blocks.push_back(block);
    }

    void show(){
        for (int i=0;i<this->blocks.size();i++){
            this->blocks[i]->show();
            cout <<"\n";
        }
    }

};


class SparseMatrixOutput{

public:
    int n;
    int m;
    int k;
    map<pair<int,int>, Block*> blocks;


    SparseMatrixOutput(int n, int m){
        this->n = n;
        this->m = m;
    }

    void insert(Block* block){
        pair<int,int> key = {block->i, block->j};
        string name = to_string(key.first) + "_" + to_string(key.second);
        #pragma omp critical (name)
        if (this->blocks.find(key) == this->blocks.end()){
            this->blocks[key] = block;
        }
        else{
            this->blocks[key]->outer(block);
        }
    }

    void show(){
        cout << this->n << " " << this->m << " " << this->k << endl;
        for (auto it = this->blocks.begin(); it != this->blocks.end(); ++it){
            it->second->show();
            cout <<"\n";
        }
    }


};


// read 4 bytes of data from the input file and return the integer 
int readInt(ifstream &input){
    int result;
    input.read((char*)&result, sizeof(int));
    return result;
}

int** transpose(int** data, int m){
    int** result = new int*[m];
    for (int a=0;a<m;a++){
        result[a] = new int[m];
    }
    for (int a=0;a<m;a++){
        for (int b=0;b<m;b++){
            result[a][b] = data[b][a];
        }
    }
    return result;
}

void read_input(SparseMatrixInput* matrix, ifstream &input){
    int k = matrix->k;
    int m = matrix->m;
    for (int z=0;z<k;z++){
        int i = readInt(input);
        int j = readInt(input);
        int** data = new int*[m];
        for (int a=0;a<m;a++){
            data[a] = new int[m];
            for (int b=0;b<m;b++){
                unsigned char val;
                input.read(( char *)&val, sizeof(unsigned char));
                data[a][b] = (int)val;
            }
        }
        Block* block = new Block(i,j,data,m);
        matrix->insert(block);
        if (i!=j){
            Block* block2 = new Block(j,i,transpose(data,m),m);
            matrix->insert(block2);
        }
    }
}

void write_output(SparseMatrixOutput* result, ofstream &output){
    output.write((char*)&result->n, sizeof(int));
    output.write((char*)&result->m, sizeof(int));
    int k = result->blocks.size();
    output.write((char*)&k, sizeof(int));
    cout << result->n << " " << result->m << " " << k << endl;
    for (auto it = result->blocks.begin(); it != result->blocks.end(); ++it){
        output.write((char*)&it->second->i, sizeof(int));
        output.write((char*)&it->second->j, sizeof(int));
        for (int a=0;a<it->second->m;a++){
            for (int b=0;b<it->second->m;b++){
                unsigned short val = (unsigned short) it->second->data[a][b];
                output.write((char*)&val, sizeof(unsigned short));
            }
        }
    }
}



int main( int argc, char** argv ){
    if (argc !=3){
        cout <<"Wrong input format";
    }

    string inputFile = argv[1];
    string outputFile = argv[2];
    
    // open a byte stream file for reading and writing

    ifstream input(inputFile, ios::binary);
    ofstream output(outputFile, ios::binary);

    // print the content of the input file 
    int n = readInt(input);
    int m = readInt(input);
    int k= readInt(input);

    cout << n << " " << m << " " << k << endl;

    SparseMatrixInput* matrix  = new SparseMatrixInput(n,m,k);
    read_input(matrix, input);
    
    // matrix->show();

    SparseMatrixOutput* result = new SparseMatrixOutput(n,m);

    auto start3 = chrono::system_clock::now();

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (Block* block1: matrix->blocks){
                #pragma omp task
                {
                    for (Block* block2: matrix->blocks){
                        if (block1->j == block2->i && block1->i <= block2->j){
                            Block* block3 = block1->inner(block2);
                            result->insert(block3);
                        }
                    }
                }
            }
        }
    }

    auto stop3 = chrono::system_clock::now();
    auto duration3 = chrono::duration_cast<chrono::nanoseconds>(stop3 - start3);
    cout << "Multiply time: " << (1e-6)*duration3.count() << " ms" << endl;
    // result->show();

    // write the result to the output file
    write_output(result, output);

}