#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <Eigen/Dense>

#include "component.h"

typedef std::complex<float> element;
typedef std::vector<element> netvec;


int main(){
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Welcome to the AC Circuit Simulation!" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl << std::endl;

    std::string inputline;

    std::vector<std::string> firstlet;
    std::vector<std::string> line;
    std::string filename;

    std::cout << "Please enter name of the input file: ";
    std::cin >> filename;

    std::ifstream infile;
    infile.open(filename);

    if(!infile.is_open()){
        std::cout << "Error, can not open input file." << std::endl;
        return EXIT_FAILURE;
    }

    while(getline(infile, inputline)){
        if(inputline.empty()){
        }
        else{
            line.push_back(inputline);
        }
    }

    for(int i = 0; i < line.size(); i++){
           std::cout << line[i] << std::endl;
    }

    std::cout << std::endl;

    std::string modeline = line[line.size()-2];
    std::vector<std::string> info2 = Getinfo(modeline);

    float n_steps, f_n;

    n_steps = 1;
    f_n = 0;

    std::cout << "simulation mode: " << info2[0] << std::endl;

    if(info2[0] == ".ac"){
        n_steps = log10(string_to_float(info2[4])/(string_to_float(info2[3])))*(string_to_float(info2[2]));
    }
    else{
        throw std::invalid_argument("File not compatible because simulation does not use AC analysis.");
    }

    std::vector<netvec> results;

    //std::cout << "No. of steps: " << n_steps << std::endl;

    for(int i = 0; i < line.size(); i++){
        line = newline(line);
    }

    // vector 'line' only contains component now
    for(int i = 0; i < line.size(); i++){
           std::cout << line[i] << std::endl;
    }

    //std::cout << "Get Nodes" << std::endl;

    //std::cout << "Line Size: " << line.size() << std::endl;

   // get the number of nodes in the circuit
   int tmp = 0; // record highest node

   for(int i = 0; i < line.size(); i++){
        std::vector<std::string> info;
        info = Getinfo(line[i]);
        int a, b;
        a = label_to_idx(info[1]);
        b = label_to_idx(info[2]);

        if(a > tmp){
            tmp = a;
        }
        if(b > tmp){
            tmp = b;
        }
        //std::cout << "tmp: "<< tmp << std::endl;
    }

    int n = tmp; // number of nodes
    //std::cout << "n: " << n << std::endl;

    // if file is doing ac analysis, then continue
    for(float i = 0; i <= n_steps; i++){
        // find the angular frequency conditions
        float f_n = (string_to_float(info2[3])) * pow(10, i/(string_to_float(info2[2])));
        float w = 2* M_PI * f_n;

        std::cout << "frequency: " << f_n << std::endl;
        std::cout << "w (angular frequency): " << w << std::endl;

        // building class and insert values to the matrix for each frequency

        // creation of a (complex number) matrix using number of nodes
        Eigen::MatrixXcf A_mat(n,n);
        A_mat.setZero(); // set all matrix values to zero
        //std::cout << A_mat << std::endl;

        const std::complex<float> Re(1.0f, 0.0f); // set real part
        const std::complex<float> Im(0.0f, 1.0f); // set imaginary part

        Eigen::MatrixXcf b_vec(n,1); // creation of b matrix
        Eigen::MatrixXcf x_vec(n,1); // creation of x matrix
        b_vec.setZero(); // set values of b matrix to zero

        for(int i = 0; i < line.size(); i++){
            // determine the first letter
            char firstLetter = line[i].at(0);
            //std::cout << "currently taking in: "<< firstLetter << std::endl;

            std::vector<std::string> info;
            info = Getinfo(line[i]);

            // classify depending on the first letter
            // order of code: R, C, L, D, Q, M

            if(firstLetter == 'R'){
                Resistor R(info);
                int i = R.getNodeA();
                int j = R.getNodeB();

                if (i == 0){
                    j--;
                    A_mat(j,j) += R.G;
                }
                else if(j == 0){
                    i--;
                    A_mat(i,i) += R.G;
                }
                else{
                    i--;
                    j--;
                    A_mat(i,i) += R.G;
                    A_mat(j,j) += R.G;
                    A_mat(i,j) -= R.G;
                    A_mat(j,i) -= R.G;
                }
            }

            if(firstLetter == 'C'){
                Capacitor C(info,w);
                int i = C.getNodeA();
                int j = C.getNodeB();

                if (i == 0){
                    j--;
                    A_mat(j,j) += C.G;
                }
                else if(j == 0){
                    i--;
                    A_mat(i,i) += C.G;
                }
                else{
                    i--;
                    j--;
                    A_mat(i,i) += C.G;
                    A_mat(j,j) += C.G;
                    A_mat(i,j) -= C.G;
                    A_mat(j,i) -= C.G;
                }
            }

            if(firstLetter == 'L'){
                Inductor L(info, w);
                int i = L.getNodeA();
                int j = L.getNodeB();

                if (i == 0){
                    j--;
                    A_mat(j,j) += L.G;
                }

                else if(j == 0){
                    i--;
                    A_mat(i,i) += L.G;
                }
                else{
                    i--;
                    j--;
                    A_mat(i,i) += L.G;
                    A_mat(j,j) += L.G;
                    A_mat(i,j) -= L.G;
                    A_mat(j,i) -= L.G;
                }
            }

            if(firstLetter == 'D'){
                Diode D(info);
                line.push_back(D.resistor_rd);
                line.push_back(D.Id);
            }

            if(firstLetter == 'Q'){
                BJT Q(info);
                line.push_back(Q.resistor_rbe);
                line.push_back(Q.resistor_r0);
                line.push_back(Q.Ic);
            }

            if(firstLetter == 'M'){
                MOSFET M(info);
                line.push_back(M.resistor_rgs);
                line.push_back(M.resistor_r0);
                line.push_back(M.Id);
            }
            

        }

        // iterate again for current source(s)
        for(int i = 0; i < line.size(); i++){
            //std::cout << "entering second loop" << std::endl;

            // this part is repeated so it needs to be simplified
            char firstLetter = line[i].at(0);
            //std::cout << "currently taking in: " << firstLetter << std::endl;

            std::vector<std::string> info;
            info = Getinfo(line[i]);

            if(firstLetter == 'I'){
                Current_source I1(info);
                int i = I1.getNodeA(); // in (tail of arrow)
                int j = I1.getNodeB(); // out (head of arrow)

                if(i ==0){
                    j--;
                    b_vec(j,0) += I1.I;
                }
                else if(j==0){
                    i--;
                    b_vec(i,0) -= I1.I;
                }
                else{
                    i--;
                    j--;
                    b_vec(i,0) -= I1.I;
                    b_vec(j,0) += I1.I;
                }
            }

            if(firstLetter == 'G'){
                V_Current_source G1(info);
                int i = G1.getNodeA()-1; // in
                int j = G1.getNodeB()-1; // out
                int k = G1.getNodeC()-1; // control +
                int l = G1.getNodeD()-1; // control -


                //for none ground nodes
                //for a, b, c, d nodes
                //A(a,c)+=G
                //A(a,d)-=G
                //A(b,c)-=G
                //A(b,d)+=G
                if(i!= -1){              //!=-1 means not ground
                    if(k!=-1){
                        A_mat(i,k) +=G1.G;
                          }
                    if(l!=-1){
                        A_mat(i,l) -=G1.G;
                          }
                }

                if(j!= -1){
                    if(k!=-1){
                        A_mat(j,k) -=G1.G;
                    }
                    if(l!=-1){
                        A_mat(j,l) +=G1.G;
                    }
                }
            }
        }

        // iterate again for voltage source(s)
        for(int i = 0; i < line.size(); i++){
            //std::cout << "entering second loop" << std::endl;

            // repeated again
            char firstLetter = line[i].at(0);
            //std::cout << "currently taking in: "<< firstLetter << std::endl;

            std::vector<std::string> info;
            info = Getinfo(line[i]);

            if(firstLetter == 'V'){
                Voltage_source V1(info);
                int i = V1.getNodeA(); // anode
                int j = V1.getNodeB(); // cathode

                if(i ==0){
                    j--;
                    b_vec(j,0) = V1.V;

                    for(int k = 0; k < n; k++){
                        //changing row j to represent vj = -vsrc
                        A_mat(j,k) = 0;
                    }

                    A_mat(j,j) = -1;
                }
                else if(j ==0){
                    i--;
                    b_vec(i,0) = V1.V;

                    for(int k = 0; k < n; k++){
                        //changing row i to represent vi = vsrc
                        A_mat(i,k) = 0;
                    }

                    A_mat(i,i) = 1;
                }

                else{
                    i--;
                    j--;
                    b_vec(j,0) += b_vec(i,0); //i_supernode
                    b_vec(i,0) = V1.V; //v_src
                    A_mat(j,j) += A_mat(i,i);


                    for(int k = 0; k < n; k++){
                        //changing row i to represent vi = vj + vsrc
                        A_mat(i,k) = 0;
                    }

                    A_mat(i,i) = 1;
                    A_mat(i,j) = -1;
                }
            }
        }

        std::cout << "----------------------------------------" << std::endl;

        // print the A matrix and b vector (visualization purposes)
        std::cout << std::endl << "the A matrix " << std::endl;
        std::cout << A_mat << std::endl;

        std::cout << std::endl << "the b vector" << std::endl;
        std::cout << b_vec << std::endl;

        Eigen::MatrixXcf A_inv(n,n); // creation of matrix to compute the inverse
        A_inv = A_mat.inverse(); // compute inverse
        std::cout << std::endl << "the inverse A matrix" << std::endl;
        std::cout << A_inv << std::endl;

        std::cout << std::endl << "test the inverse matrix (should produce identity matrix)" << std::endl;
        std::cout << A_mat * A_inv << std::endl;

        std::cout << std::endl << "the x vector (solution)" << std::endl;
        x_vec = A_inv * b_vec;
        std::cout << x_vec << std::endl;

        std::cout << std::endl;

        // store values of x vector (matrix) into results vector to save to a file
        netvec voltage_at_nodes;

        for(int i = 0; i < n; i++){
            voltage_at_nodes.push_back(x_vec(i,0));
        }

        results.push_back(voltage_at_nodes);
    }
    
    

    std::ofstream outfile;
    outfile.open("voltage.txt");

    if(!outfile.is_open()){
        std::cout << "error opening file" << std::endl;
        return EXIT_FAILURE;
    }

    for(int i =0; i<results.size(); i++){
        for(int j =0; j<results[0].size(); j++){
            outfile << results[i][j] <<" ";
        }
        outfile<< "\n";
    }

    outfile.close();
    
    std::string in;
    std::cout<<"Generating Transfer Function"<<std::endl;
    
    while(in!="END"){
        std::vector<std::string> bode_results;

    // user input values for the reference node and output node
        int n_ref, n_out;
        std::cout << "Enter the reference node: ";
        std::cin >> n_ref;
        while(n_ref >= (n+1)){
            std::cout << "Error, invalid reference node value. Enter a valid value for the reference node: ";
            std::cin >> n_ref;
        }
            std::cout << "Enter the output node: ";
        std::cin >> n_out;
        
        while((n_out >= (n+1)) || (n_out == n_ref)){
            std::cout << "Error, invalid output node value. Enter a valid value for the output node: ";
            std::cin >> n_out;
            
        }
    
        std::string heading;
        heading.append("Freq.   ");
        heading.append("V(n");
        heading.append(std::to_string(n_out));
        heading.append(")/V(n");
        heading.append(std::to_string(n_ref));
        heading.append(")");
        bode_results.push_back(heading);
    
        n_ref = n_ref - 1;
        n_out = n_out - 1;

        std::cout << "Reference Node Values are: " << std::endl;
        if(n_ref == -1){
            for(int i = 0; i < results.size(); i++){
                std::cout << "(0,0)" << std::endl;
            }
        }
        else{
            for(int i = 0; i < results.size(); i++){
                std::cout << results[i][n_ref] << std::endl;
            }
        }

        std::cout << std::endl << "Output Node Values are: " << std::endl;
        for(int i = 0; i < results.size(); i++){
            std::cout << results[i][n_out] << std::endl;
        }
    
    

        std::cout << std::endl << "Calculation of Bode Plot Entries: " << std::endl;
        for(int i = 0; i < results.size(); i++){
            element a = results[i][n_ref];
            element b = results[i][n_out];
            element c; // phasor

            if(n_ref == -1){
            // if n_ref = -1, then value is referring to ground
                c = b;
            }
            else{
                c = b / a; // transfer function = output node / reference node
            }
        
        //-----------debug using cartesian form
        
        /*
            float re = c.real();
            float im = c.imag();
            float f_n = (string_to_float(info2[3])) * pow(10, i/(string_to_float(info2[2])));
            std::string values;
            values.append(std::to_string(f_n));
            values.append("    ");
            values.append(std::to_string(re));
            values.append(", ");
            values.append(std::to_string(im));
        */
        
        

        
        float phase = get_phase(c);
            float magnitude = 20*log10(get_magnitude(c));

        float f_n = (string_to_float(info2[3])) * pow(10, i/(string_to_float(info2[2])));
        //float w = 2* M_PI * f_n;

        // store bode plot results into bode_results in CSV format
        std::string values;
        values.append(std::to_string(f_n));
        values.append(", ");
        values.append(std::to_string(magnitude));
        values.append("dB, ");
        values.append(std::to_string(phase));
       
         

            bode_results.push_back(values);
        }

    // print bode plot results (visualization purposes)
    for(int i = 0; i < bode_results.size(); i++){
        std::cout << bode_results[i] << std::endl;
    }
    std::cout << std::endl;
        
    std::string title;
        title.append("bode_plot_");
        title.append("n");
        title.append(std::to_string(n_out+1));
        title.append("onn");
        title.append(std::to_string(n_ref+1));
        title.append(".txt");
        
        std::cout<<title<<std::endl;

    std::ofstream outfile2;
    outfile2.open(title);

    if(!outfile2.is_open()){
        std::cout << "error opening file" << std::endl;
        return EXIT_FAILURE;
    }

    for(int i = 0; i < bode_results.size(); i++){
        outfile2 << bode_results[i] << std::endl;
    }

    outfile2.close();
        
    std::cout<<"Enter 'END' to end the program or enter anything to generate another bode plot"<<std::endl;
    std::cin >> in;
    }
    
    
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "The AC Circuit Simulation is completed." << std::endl;
    std::cout << "Please view the results for the bode plot in the text file." << std::endl;
    std::cout << std::endl;
    std::cout << "Thank you for using the AC Circuit Simulator software package." << std::endl;
    std::cout << "Credits to Aman Narain, Christina Wang, and Yuhe Zhang" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl << std::endl;
}

