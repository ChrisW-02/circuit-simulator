#include <iostream>
#include <complex>

typedef std::complex<float> element;

struct Matrix{
    
    int dim;
    element** mat;
    element** inv;

    void delete(element** a){
        for(int i = 0; i < dim; ++i)
        {
            delete []a[i];
        }
        delete []a;
    }

    void set_dim(){
        mat = new element* [dim];
        inv = new element* [dim];
        for(int i = 0; i < dim; i++)
        {
               mat[i] = new element [dim];
               inv[i] = new element [dim];
        }
    } 


    void inputmat(){
        set_dim();
        for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            std::cin >> mat[i][j]; 
        }
    }
    }

    void print(element** a){
        for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            std::cout << a[i][j]<< "\t";
        }
        std::cout << "\n";
    }
    }

    
    //get inverse


    element* matsolv(element* b, element* x){
        /*if(dim != size){
            std::cout<<"Error: dimensions mismatch!!!" <<std::endl;
        }*/
        
        for(int i=0; i<dim; i++){
            for (int j=0; j<dim; j++){
                b[i] += inv[i][j] * x[j];
            }
        }
        delete_inv();
        return b;
    }



};

int main(){
    Matrix A;
    std::cout<<"Please enter the dimension of the matrix"<<std::endl;
    std::cin >> A.dim;

    std::cout << "please enter the entries of the matrix" << std::endl;
    A.inputmat();
    A.print(A.mat);

    element b[2] = {1,2};

    
    element* x;
    x = new element;

    x = A.matsolv(x,b);


    std::cout<<std::endl;


    std::cout << x[0] <<x[1]<< std::endl;


    delete[] x;
    A.delete_mat();

}


//possible enhancement:
//combine vector and matrix as one class, ie. share the same input/print function
//use 1D array instead of 2D to increase time efficiency
//check the dimension of the vector and the matrix