#include <bits/stdc++.h>
using namespace std;

#define int unsigned int

class SparseMatrixInput{

public:
    int n;
    int m;
    int k;
    int new_n;
    int* block_data;
    int* block_index;
    int* cuda_block_data;
    int* cuda_block_index;

    SparseMatrixInput(){
        this->n = 0;
        this->m = 0;
        this->k = 0;
        this->new_n = 0;
        this->block_data = NULL;
        this->block_index = NULL;
    }

    void initialize_arrays(){
        this->new_n = this->n/this->m;
        this->block_data = new int[this->k*this->m*this->m];
        this->block_index = new int[(this->new_n)*(this->new_n)];
        // intialize the block index to -1
        for (int i=0;i<(this->new_n)*(this->new_n);i++){
            this->block_index[i] = this->n*this->n;
        }
        // copy the data to the gpu
        cudaMalloc((void**)&this->cuda_block_data, this->k*this->m*this->m*sizeof(int));
        cudaMalloc((void**)&this->cuda_block_index, (this->new_n)*(this->new_n)*sizeof(int));
    }

    void copy_to_cuda(){
        cudaMemcpy(this->cuda_block_data, this->block_data, this->k*this->m*this->m*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(this->cuda_block_index, this->block_index, (this->new_n)*(this->new_n)*sizeof(int), cudaMemcpyHostToDevice);
    }


    // print the entire matrix
    void print_matrix(){
        cout << "n: " << this->n << " m: " << this->m << " k: " << this->k << endl;
        for (int i=0;i<this->new_n*this->new_n;i++){
            if (this->block_index[i] == this->n*this->n)
                continue;
            cout << "Block index: " << i/this->new_n << " " << i%this->new_n << endl;
            for (int j=0;j<this->m*this->m;j++){
                cout << this->block_data[this->block_index[i]+j] << " ";
            }
            cout << endl;
        }
    }


    // make a destructor which deallocate memory from cuda as well
    ~SparseMatrixInput(){
        delete[] this->block_data;
        delete[] this->block_index;
        cudaFree(this->cuda_block_data);
        cudaFree(this->cuda_block_index);
    }


};


class SparseMatrixOutput{

public:
    int n;
    int m;
    int k;
    int new_n;
    int* block_data;
    int* block_index;
    int* cuda_block_data;
    int* cuda_block_index;

    SparseMatrixOutput(int n, int m){
        this->n = n;
        this->m = m;
        this->k = 0;
        this->new_n = n/m;
        this->block_data = new int[this->new_n*this->new_n*this->m*this->m];
        this->block_index = new int[this->new_n*this->new_n];
        // intialize the block index to -1
        for (int i=0;i<(this->new_n)*(this->new_n);i++){
            this->block_index[i] = i*m*m;
            // data values to zero 
            for (int j=0;j<this->m*this->m;j++){
                this->block_data[i*m*m+j] = 0;
            }
        }
        // copy the data to the gpu
        cudaMalloc((void**)&this->cuda_block_data, this->new_n*this->new_n*this->m*this->m*sizeof(int));
        cudaMalloc((void**)&this->cuda_block_index, (this->new_n)*(this->new_n)*sizeof(int));
        cudaMemcpy(this->cuda_block_data, this->block_data, this->new_n*this->new_n*this->m*this->m*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(this->cuda_block_index, this->block_index, (this->new_n)*(this->new_n)*sizeof(int), cudaMemcpyHostToDevice);
    }

    void check_zeros(){
        for (int i=0;i<this->new_n*this->new_n;i++){
            int flag = 0;
            for (int j=0;j<this->m*this->m;j++){
                if (this->block_data[this->block_index[i]+j] != 0){
                    flag = 1;
                    break;
                }
            }
            if (flag == 0){
                this->block_index[i] = this->n*this->n;
            }
            else{
                this->k++;
            }
        }
    }

    // print the entire matrix
    void print_matrix(){
        cout << "n: " << this->n << " m: " << this->m << " k: " << this->k << endl;
        for (int i=0;i<this->new_n*this->new_n;i++){
            if (this->block_index[i] == this->n*this->n)
                continue;
            cout << "Block index: " << i/this->new_n << " " << i%this->new_n << endl;
            for (int j=0;j<this->m*this->m;j++){
                cout << this->block_data[this->block_index[i]+j] << " ";
            }
            cout << endl;
        }
    }

    // make a destructor which deallocate memory from cuda as well
    ~SparseMatrixOutput(){
        delete[] this->block_data;
        delete[] this->block_index;
        cudaFree(this->cuda_block_data);
        cudaFree(this->cuda_block_index);
    }

};


// read 4 bytes of data from the input file and return the integer 
int readInt(ifstream &input){
    int result;
    input.read((char*)&result, sizeof(int));
    return result;
}


void read_input(SparseMatrixInput* matrix, ifstream &input){

    int n = matrix->n = readInt(input);
    int m = matrix->m = readInt(input);
    int k = matrix->k= readInt(input);
    matrix->initialize_arrays();
    int nn = n/m;

    int iter =0;
    for (int z=0;z<k;z++){
        int i = readInt(input);
        int j = readInt(input);
        // cout << i << " " << j << endl;
        matrix->block_index[i*nn + j] = iter;
        for (int a=0;a<m*m;a++){
            unsigned short val;
            input.read(( char *)&val, sizeof(unsigned short));
            matrix->block_data[iter++] = (int)val;
        }
    }
    // copy the data to the gpu
    matrix->copy_to_cuda();
}

void write_output(SparseMatrixOutput* result, ofstream &output){
    output.write((char*)&result->n, sizeof(int));
    output.write((char*)&result->m, sizeof(int));
    output.write((char*)&result->k, sizeof(int));

    for (int i=0; i<result->new_n*result->new_n; i++){
        if (result->block_index[i] != result->n*result->n){
            int ii = i/result->new_n;
            int jj = i%result->new_n;
            output.write((char*)&ii, sizeof(int));
            output.write((char*)&jj, sizeof(int));
            output.write((char*)&result->block_data[result->block_index[i]], sizeof(int)*result->m*result->m);
        }
    }
}

// multiply the two matrices and store the result in the result matrix using cuda 
__global__ void multiply(int* index1, int* index2, int* data1, int* data2, int* result_data, int n, int m){
    
    // printf("Hello\n");
    
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int nn = n/m;
    int mm = m*m;

    int i = index/nn;
    int j = index%nn;


    // cout << index << " " << i << " " << j << endl;

    if (index < nn*nn ){
        for (int idx = 0; idx < nn; idx++){
            int block1_index = index1[i*nn + idx];
            int block2_index = index2[idx*nn + j];
            // printf("%d %d %d %d\n",i ,j, block1_index, block2_index);
            if (block1_index != n*n && block2_index != n*n){
                // printf("hello\n");
                for (int a=0;a<m;a++){
                    for (int b=0;b<m;b++){
                        for (int c=0;c<m;c++){
                            unsigned long val = result_data[i*nn*mm + j*mm + a*m + b] + (unsigned long) data1[block1_index + a*m + c]*data2[block2_index + c*m + b];
                            result_data[i*nn*mm + j*mm + a*m + b] = min(val,(unsigned long) 0xffffffff);
                        }
                    }
                }
            }
        }
    }
}

void print_time(string s, chrono::high_resolution_clock::time_point& start, chrono::high_resolution_clock::time_point& end){
    if (s == "Total time"){
    end = chrono::high_resolution_clock::now();
    // print the time duration in mili seconds
    cout << s << " " << chrono::duration_cast<chrono::nanoseconds>(end - start).count()*1e-6 << " mseconds" << endl;
    start = chrono::high_resolution_clock::now();
    }
}


int32_t main( int32_t argc, char** argv ){
    if (argc !=4){
        cout <<"Wrong input format";
    }
    auto start_total = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto end = chrono::high_resolution_clock::now();

    string inputFile1 = argv[1];
    string inputFile2 = argv[2];
    string outputFile = argv[3];

    ifstream input1(inputFile1, ios::binary);
    ifstream input2(inputFile2, ios::binary);
    ofstream output(outputFile, ios::binary);

    print_time("Opening the input files", start, end);

    SparseMatrixInput* matrix1  = new SparseMatrixInput();
    read_input(matrix1, input1);

    SparseMatrixInput* matrix2  = new SparseMatrixInput();
    read_input(matrix2, input2);

    int n = matrix1->n;
    int m = matrix1->m;
    SparseMatrixOutput* result = new SparseMatrixOutput(n,m);

    print_time("Reading the input files", start, end);

    multiply<<<(n * n + 1023) / 1024, 1024>>>(matrix1->cuda_block_index, matrix2->cuda_block_index, matrix1->cuda_block_data, matrix2->cuda_block_data, result->cuda_block_data, n, m);

    // synchronize the threads
    cudaDeviceSynchronize();

    print_time("multiplication", start, end);

    // copy the result back to the host
    cudaMemcpy(result->block_data, result->cuda_block_data, result->new_n*result->new_n*result->m*result->m*sizeof(int), cudaMemcpyDeviceToHost);
    result->check_zeros();
    write_output(result, output);

    print_time("writing output", start, end);

    print_time("Total time", start_total, end);

}