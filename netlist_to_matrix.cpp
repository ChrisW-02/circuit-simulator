#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sstream>

///////////////////////general component class and functions

float string_to_float(const std::string& label){
 
    if (label.find("p") != std::string::npos) {
        return std::stof(label) * pow(10,-12);
    }

    if (label.find("n") != std::string::npos) {
        return std::stof(label) * pow(10,-9);
    }

    if (label.find("u") != std::string::npos || label.find("μ") != std::string::npos) {
        return std::stof(label) * pow(10,-6);
    }

    if (label.find("m") != std::string::npos) {
        return std::stof(label) * pow(10,-3);
    }

    if (label.find("k") != std::string::npos) {
        return std::stof(label) * pow(10,3);
    }

    if (label.find("Meg") != std::string::npos) {
        return std::stof(label) * pow(10,6);
    }

    if (label.find("G") != std::string::npos) {
        return std::stof(label) * pow(10,9);
    }
    
    return std::stof(label);
}

//extracted from the bt assignment

int label_to_idx(const std::string& label){
    
    std::string out;
    
    for(int i = 1; i < label.size(); i++){
        out.push_back(label[i]);
    }
    
    return std::stoi(out) - 1;
}
 
class Component{
    

public:

    std::complex<float> G;
    std::complex<float> I;
    std::complex<float> V;

    float w;

    virtual ~Component(){
        std::cout << "Base class component destructor." << std::endl;
    }

    //virtual float getFreq();

    virtual std::complex<float> getI(){
        return I;
    }

    std::string getCompName(){
        return CompName;
    };

    virtual int getNodeA(){
        return label_to_idx(NodeA);
    };

    virtual int getNodeB(){
        return label_to_idx(NodeB);
    };

    virtual int getNodeC(){
        return label_to_idx(NodeC);
    };
        

    virtual int getNodeD(){
        return label_to_idx(NodeD);
    };
       

    virtual bool AC(std::string value){
        char mode = value[0];
        if(mode == 'A'){
            return true;
        }
        return false;
    }

protected:

    std::string CompName;
    std::string NodeA;
    std::string NodeB;
    std::string NodeC;
    std::string NodeD;
    std::string value;
    

};

std::vector<std::string> Getinfo(std::string aline)
{
    // Used to split string around spaces.
    std::stringstream ss(aline);
    std::vector<std::string> info;
  
    std::string word; // for storing each word
  
    // Traverse through all words
    // while loop till we get 
    // strings to store in string word
    while (ss >> word) 
    {
        // print the read word
        info.push_back(word);
    }
    return info;
}

std::complex<float> get_AC(std::string str1, std::string str2){
    str1.erase(0,3);
    int idx = str2.size()-1;
    str2.erase(idx,1);
    std::istringstream is('(' + str1 + "," + str2 + ')');
    std::complex<float> c;
    is >> c;
    return c;
}

std::vector<std::string> newline(std::vector<std::string> line){
    
    std::string ast = "*";
    std::string c = "c";
    std::string end = "end";

    std::string::size_type idx;

    for(int i = 0; i < line.size(); i++){
        idx = line[i].find(ast);
        if(idx == std::string::npos){}
        else{
            line.erase(line.begin() + (i));
        }
    }

    for(int i = 0; i < line.size(); i++){
        idx = line[i].find(c);
        if(idx == std::string::npos){}
        else{
            line.erase(line.begin() + i);
        }
    }

    for(int i = 0; i < line.size(); i++){
        idx = line[i].find(end);
        if(idx == std::string::npos){}
        else{
            line.erase(line.begin() + i);
        }
    }

    return line;
}

/////////////////////////class for each type of components
class Voltage_source: public Component{
public: 

    Voltage_source(){
        V = 0;
    }

    Voltage_source(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        G = 0;

        if(AC(info[3])){
            //get the value after AC
            V = get_AC(info[3], info[4]);
        }
        else{
            V = 0; 
        }
    }
};

class Current_source: public Component{
public: 
    
    Current_source(){
        I = 0;
    }

    Current_source(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        G = 0;

        if(AC(info[3])){
            //get the value after AC
            I = get_AC(info[3], info[4]);
        }
        else{
            I = 0; 
        }
    }
};

class Resistor: public Component{
public:

    Resistor(){
        //default constructor
        G = 0;
    }

    Resistor(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        if(info[3]=="INFINITE"){
            G = 0;
        }
        else{
            G.real(1/string_to_float(info[3]));
        }
    }

};

class Capacitor: public Component{
public:

    Capacitor(){
        G = 0;
    }

    Capacitor(std::vector<std::string> info ){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        G.imag(w * string_to_float(info[3]));
        
    }
    
};

class Inductor: public Component{
public:


    Inductor(){
        //default constructor
        G = 0;
    }

    Inductor(std::vector<std::string> info){

        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        G.imag(-1/(w * string_to_float(info[3])));
    }
    
};

class Diode: public Component{
public:

    Diode(){
        //default constructor
        G = 0;
    }

    Diode(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        V = 0.7; //VD
        G = is_const * (exp(V / VT)) /VT ; //dynamic conductance of diode = Id/VT, in small signal model a diode equals to a resistance
    }

    float is_const = 20e-9; // Reverse saturation current, 20 nA
    float VT = 0.025; // can we just use VT

};


class BJT: public Component{
public:

//Q1 N003 N001 0 NPN
    
    std::string resistor_rbe;
    std::string resistor_r0;
    std::string Ic;

    BJT(std::vector<std::string> info){
        resistor_rbe.append("Rbe_");
        resistor_rbe.append(info[0]);
        resistor_rbe.push_back(' ');
        resistor_rbe = resistor_rbe + (info[2]);
        resistor_rbe.push_back(' ');
        resistor_rbe = resistor_rbe + (info[3]);
        resistor_rbe.push_back(' ');

        resistor_r0.append("R0_");
        resistor_r0.append(info[0]);
        resistor_r0.push_back(' ');
        resistor_r0 = resistor_r0 + (info[1]);
        resistor_r0.push_back(' ');
        resistor_r0 = resistor_r0 + (info[3]);
        resistor_r0.push_back(' ');

        Ic.append("Ic_");
        Ic.append(info[0]);
        Ic.push_back(' ');
        Ic = Ic + (info[1]);
        Ic.push_back(' ');
        Ic = Ic + (info[3]);
        Ic.push_back(' ');

        if(info[4] == "NPN"){//"2N2222"

            //gm = Ic/VT = 800mA/25mV = 32;
            //Rbe = beta/gm = 200/32 = 6.25;
            //r0 = VA/Ic = 100/0.8 = 125;

            resistor_rbe.append("6.25");

            resistor_r0.append("125");

            Ic.append("AC(0.8 0)");

        }

        if(info[4] == "PNP"){//"2N2907"
            //gm = Ic/VT = 600mA/25mV = 24;
            //Rbe = beta/gm = 250/24 = 10.4167;
            //r0 = VA/Ic = 120/0.6 = 200;

            resistor_rbe.append("10.4167");

            resistor_r0.append("200");

            Ic.append("AC(0.6 0)");

        }
        std::cout << "resistor_rbe: " << resistor_rbe <<std::endl;
        std::cout << "resistor_r0: " << resistor_r0 <<std::endl;
        std::cout << "Ic: " << Ic << std::endl;

    }

};

class MOSFET: public Component{
public:
    
    std::string resistor_rgs;
    std::string resistor_r0;
    std::string Id;

    MOSFET(std::vector<std::string> info){

        //betweem gate and source is open circuit with potential difference of Vgs
        //represent as a resistor with infinite resistance;

        resistor_rgs.append("Rgs_");
        resistor_rgs.append(info[0]);
        resistor_rgs.push_back(' ');
        resistor_rgs = resistor_rgs + (info[2]);
        resistor_rgs.push_back(' ');
        resistor_rgs = resistor_rgs + (info[3]);
        resistor_rgs.push_back(' ');
        resistor_rgs.append("infinite");

        resistor_r0.append("R0_");
        resistor_r0.append(info[0]);
        resistor_r0.push_back(' ');
        resistor_r0 = resistor_r0 + (info[1]);
        resistor_r0.push_back(' ');
        resistor_r0 = resistor_r0 + (info[3]);
        resistor_r0.push_back(' ');

        Id.append("Id_");
        Id.append(info[0]);
        Id.push_back(' ');
        Id = Id + (info[1]);
        Id.push_back(' ');
        Id = Id + (info[3]);
        Id.push_back(' ');

        if(info[4] == "NMOS"){//"FDS5680"
           
            resistor_r0.append("0.025");

            Id.append("AC(0.8 0)");

            //gm = 2K(Vgs-Vt)
            //r0 = VA/Id = ???;
            //Id = gm * Vgs;


        }

        if(info[4] == "PMOS"){//"A06407"
            
            resistor_r0.append("200");

            Id.append("AC(0.6 0)");

            //r0 = VA/Id = 120/0.6 = 200;
            //Id = gm * Vgs

        }
        std::cout << "resistor_rgs: " << resistor_rgs <<std::endl;
        std::cout << "resistor_r0: " << resistor_r0 <<std::endl;
        std::cout << "Id: " << Id << std::endl;

    }
    
};

class V_Current_source: public Component{
public: 

    V_Current_source(){
        I = 0;
    }

    V_Current_source(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        NodeC = info[3];
        NodeD = info[4];
        
        G = 0;

        //write a program to obtain the voltage between node c and d
        //need to finish this
        //Vcd = ??;

        I = Vcd * string_to_float(info[5]);

    }
};


////////////////////////matrix

typedef std::complex<float> element;
typedef std::vector<std::string> netvec;

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
        
        for(int i=0; i<dim; i++){
            for (int j=0; j<dim; j++){
                b[i] += inv[i][j] * x[j];
            }
        }
        delete_inv();
        return b;
    }

};





/////////////////////////////////main

int main(){
    std::string inputline;
    
    std::vector<std::string> firstlet;
    std::vector<std::string> line;
    std::string filename;

    std::cout << "please enter name of the input file" << std::endl;
    std::cin >> filename;
 
    std::ifstream infile;
    infile.open(filename);
 
    if(!infile.is_open()){
        std::cout << "error, can't open input file" << std::endl;
        return EXIT_FAILURE;
    }

    std::string tmp;

    while(getline(infile, inputline)){
        line.push_back(inputline);
    }
    
    for(int i = 0; i < line.size(); i++){
           std::cout << line[i] << std::endl;
    }
    
    std::cout << std::endl;

    for(int i = 0; i < line.size(); i++){
        line = newline(line);
    }

    //write a program to get the number of nodes in the circuit

    int n = 4; //number of nodes

    //creating a matrix
    Matrix A;
    A.dim = n;

    //creating the right hand side vector
    element b[n] = {0};

    
    //vector 'line' only contanis component now
    for(int i = 0; i < line.size(); i++){
           std::cout << line[i] << std::endl;
    }

    ////////////////////////////////building class and insert values to the matrix
    
    for(int i = 0; i < line.size(); i++){

        //determine the first letter
        char firstLetter = line[i].at(0);
        std::cout << firstLetter << std::endl;

        std::vector<std::string> info;
        info = Getinfo(line[i]);

        //to classify...
        if(firstLetter == 'R'){

            Resistor R(info);

            int i = R.getNodeA;
            int j = R.getNodeB;
            if(j<i){
                int tmp = j;
                j = i;
                i = tmp;
            }

            A[i][i] += R.G; //eg. G11
            A[j][j] += R.G; //eg. G22

            for(int k =i+1; k<n; k++){   //eg. G1n and Gn1
                A[i][k] -= R.G;
                A[k][i] -= R.G;
            }
            for(int k =0; k<j; k++){   //eg. G2n and Gn2
                if(k!=i){
                    A[j][k] -= R.G;
                    A[k][j] -= R.G;
                }
            }
            
        }

        if(firstLetter == 'C'){
            
            Capacitor C(info);
            int i = C.getNodeA;
            int j = C.getNodeA;
            
            if(j<i){
                int tmp = j;
                j = i;
                i = tmp;
            }

            A[i][i] += C.G; //eg. G11
            A[j][j] += C.G; //eg. G22

            for(int k =i+1; k<n; k++){   //eg. G1n and Gn1
                A[i][k] -= C.G;
                A[k][i] -= C.G;
            }
            for(int k =0; k<j; k++){   //eg. G2n and Gn2
                if(k!=i){
                    A[j][k] -= C.G;
                    A[k][j] -= C.G;
                }
            }
            
        }

        if(firstLetter == 'L'){
            
            Inductor L(info);
            int i = L.getNodeA;
            int j = L.getNodeA;
            
            if(j<i){
                int tmp = j;
                j = i;
                i = tmp;
            }

            A[i][i] += L.G; //eg. G11
            A[j][j] += L.G; //eg. G22

            for(int k =i+1; k<n; k++){   //eg. G1n and Gn1
                A[i][k] -= L.G;
                A[k][i] -= L.G;
            }
            for(int k =0; k<j; k++){   //eg. G2n and Gn2
                if(k!=i){
                    A[j][k] -= L.G;
                    A[k][j] -= L.G;
                }
            }
            
        }

        if(firstLetter == 'D'){
        
            Diode D(info);
            int i = D.getNodeA;
            int j = D.getNodeA;

            //??????a diode is a one-way resistance, in the other direction it has infinite resistance
            
           if(j<i){
                int tmp = j;
                j = i;
                i = tmp;
            }

            A[i][i] += L.G; //eg. G11
            A[j][j] += L.G; //eg. G22

            for(int k =i+1; k<n; k++){   //eg. G1n and Gn1
                A[i][k] -= L.G;
                A[k][i] -= L.G;
            }
            for(int k =0; k<j; k++){   //eg. G2n and Gn2
                if(k!=i){
                    A[j][k] -= L.G;
                    A[k][j] -= L.G;
                }
            }
            
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

    //iterate again

    for(int i = 0; i < line.size(); i++){
        if(firstLetter == 'I'){
            Current_source I1(info);
            int i = I1.getNodeA; //in
            int j = I1.getNodeB; //out
            b[i] += I1.I;
            b[j] -= I1.I;
            
        }

        if(firstLetter == 'G'){
            V_Current_source G1(info);
            int i = G1.getNodeA; //in
            int j = G1.getNodeB; //out
            b[i] += G1.I;
            b[j] -= G1.I;
        }
    }

    //iterate again
    for(int i = 0; i < line.size(); i++){

        if(firstLetter == 'V'){

            Voltage_source V1(info);

            int i = V1.getNodeA; //anode
            int j = V1.getNodeB; //cathode
            b[j]+=b[i]; //i_supernode
            b[i] = V1.V; //v_src

            A[j][j] += A[i][i];

            for(int k=0; k<n; k++){ //changing row i to represent vi = vj + vsrc
                A[i][k] = 0;
            }
            A[i][i] = 1;
            A[i][j] = -1;
            
        }

    }
}