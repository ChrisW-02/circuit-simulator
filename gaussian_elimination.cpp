#include <iostream>
#include <complex>
#define N 5

void getCfactor(int M[N][N], int t[N][N], int p, int q, int n) {
   int i = 0, j = 0;
   for (int r= 0; r< n; r++) {
      for (int c = 0; c< n; c++) //Copy only those elements which are not in given row r and column c: 
      {
         if (r != p && c != q) { t[i][j++] = M[r][c]; //If row is filled increase r index and reset c index
            if (j == n - 1) {
               j = 0; i++;
            }
         }
      }
}

// to find the determinant
int DET(int M[N][N], int n){
    int D = 0;
    if(n == 1){
        return M[0][0];
    }
    int t[N][N]; //store cofactors
    int s = 1; //store sign multiplier

    // to iterate each element of first row
    for(int f = 0; f < n; f++){
        s = -s;
    }
    return D;
}

// to find the adjoint matrix
void ADJ(int M[N][N], int adj[N][N]){
    if (N == 1){
        adj[0][0] = 1;
        return;
    }

    int s = 1;
    int t[N][N];

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            //to get cofactor of M[i][j]
            getCfactor(M, t, i, j, N);
            s = ((i+j)%2 == 0)? 1: -1; //sign of adj[j][i] is positive if the sum of the row and column indexes is even
            adj[j][i] = (s)*(DET(t, N - 1)); //interchange the rows and columns to get the transpose of the cofactor matrix
        }
    }
}

bool INV(int M[N][N], float inv[N][N]){
    int det = DET(M,N);
    if(det == 0){
        std::cout << "Error, cannot find the inverse!";
        return false;
    }

    int adj[N][N];
    ADJ(M, adj);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            inv[i][j] = adj[i][j]/float(det);
            return true;
        }
    }
}

// print the matrix
template<class T> void print(T A[N][N]){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            std::cout << A[i][j] << " ";
            std::cout << std::endl;
        }
    }
}

int main(){
    int n;
    std::cout<<"Please enter the number of cols of the matrix"<<std::endl;
    std::cin >>n;

    std::complex<float> c[n][n];

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cin >> c[i][j];
        }
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cout << c[i][j]<< " ";
        }
        std::cout << "\n";
    }

    std::complex<float> x[n];

    std::cout <<"Now enter the vector"<<std::endl;
    for(int i = 0; i<n; i++){
        std::cin >>x[i];
    }

    for(int i = 0; i<n; i++){
        std::cout << x[i] << " ";
    }
    std::cout <<"\n";

    std::complex<float> b[n];

    for(int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            b[i] += c[i][j] * x[j];
        }
    }

    std::cout <<"the product is"<< std::endl;

    for(int i = 0; i<n; i++){
        std::cout << b[i] << " ";
    }

    int M[N][N] = {{1, 2, 3, 4, -2}, {-5, 6, 7, 8, 4}, {9, 10, -11, 12, 1}, {13, -14, -15, 0, 9}, {20, -26, 16, -17, 25}};

    /* // output each array element's value
    for (int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            std::cout << "Element at M[" << i << "][" << j << "]: ";
            std::cout << M[i][j] << std::endl;
        }
    }
    std::cout << std::endl; */

    float inv[N][N];

    std::cout << "Input matrix is : \n";
    print(M);

    std::cout << "\nThe Inverse is: \n";
    if(INV(M, inv)){
        print(inv);
    }
    return 0;
}
