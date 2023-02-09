#include <bits/stdc++.h>
// #include <library.hpp>
using namespace std;


// read 4 bytes of data from the input file and return the integer 
int readInt(ifstream &input){
    int result;
    input.read((char*)&result, sizeof(int));
    return result;
}



int main( int argc, char** argv ){
    if (argc !=3){
        cout <<"Wrong input format";
    }

    string inputFile = argv[1];
    string outputFile = argv[2];
    
    // open a byte stream file for reading and writing

    ifstream input(inputFile, ios::in | ios::binary);
    // ofstream output(outputFile, ios::out | ios::binary);

    // print the content of the input file 
    int n = readInt(input);
    int m = readInt(input);
    int k= readInt(input);
    // cout <<  n<< " " << m  <<" " << k << endl;

    


}