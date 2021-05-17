#include <iostream>
#include <complex>

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


}
