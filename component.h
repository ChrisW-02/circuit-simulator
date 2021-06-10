#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iomanip>

// ********************************************************************************
// Component header to provide capability to all component classes and functions
// ********************************************************************************

// function to convert the unit multipliers
float string_to_float(const std::string& label){

    if(label.find("p") != std::string::npos){
        return std::stof(label) * pow(10,-12);
    }

    if(label.find("n") != std::string::npos){
        return std::stof(label) * pow(10,-9);
    }

    if((label.find("u") != std::string::npos) || (label.find("μ") != std::string::npos)){
        return std::stof(label) * pow(10,-6);
    }

    if(label.find("m") != std::string::npos){
        return std::stof(label) * pow(10,-3);
    }

    if(label.find("k") != std::string::npos){
        return std::stof(label) * pow(10,3);
    }

    if(label.find("Meg") != std::string::npos){
        return std::stof(label) * pow(10,6);
    }

    if(label.find("G") != std::string::npos){
        return std::stof(label) * pow(10,9);
    }

    return std::stof(label);
}

int label_to_idx(const std::string& label){
    if(label == "0"){
        int zero = 0;
        return zero;
    }
    std::string out;

    for(int i = 1; i < label.size(); i++){
        out.push_back(label[i]);
    }

    return std::stoi(out);
}

// ********************************************************************************
// COMPONENT CLASS (general)
// ********************************************************************************
class Component{
public:
    std::complex<float> G;
    std::complex<float> I;
    std::complex<float> V;
    float w;

    virtual ~Component(){};

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

std::vector<std::string> Getinfo(std::string aline){
    std::stringstream ss(aline); // used to split string around spaces
    std::vector<std::string> info;
    std::string word; // for storing each word

    // while loop to traverse through all words till we get strings to store in string word
    while (ss >> word){
        // print the read word
        info.push_back(word);
    }
    return info;
}

// AC function is in form <amplitude> <phase>
std::complex<float> get_AC(std::string str1, std::string str2){
    str1.erase(0,3);
    int idx = str2.size()-1;
    str2.erase(idx,1);
    float A = std::stof(str1);
    float theta = std::stof(str2);
    theta = theta / 180 * M_PI;

    std::complex<float> c;

    // A*exp(j*theta) or A*(cos(theta)+isin(theta))
    // store sine and cosine in real and imaginary parts
    c.real(A*cos(theta));
    c.imag(A*sin(theta));
    return c;
}

// function to solve the magnitude of the phasor
float get_magnitude(std::complex<float> val){
    float A = std::abs(val); // absolute value of complex
    return A;
}

// function to solve the phase of the phasor (theta)
float get_phase(std::complex<float> val){
    float angle = std::arg(val); // phase angle of complex in radians
    float theta = angle * (180.0/M_PI); // convert to radians
    return theta;
}

// function to ignore excess information given in the netlist file
std::vector<std::string> newline(std::vector<std::string> line){
    std::string ast = "*";
    std::string c = "c";
    std::string endFile = "end";

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
        idx = line[i].find(endFile);
        if(idx == std::string::npos){}
        else{
            line.erase(line.begin() + i);
        }
    }

    return line;
}

// ********************************************************************************
// CLASSES FOR EACH TYPE OF COMPONENT (specific)
// ********************************************************************************
// class for voltage source
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
            // get the value after AC
            V = get_AC(info[3], info[4]);
        }
        else{
            V = 0;
        }
    }
};

// class for current source
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

// class for resistor
class Resistor: public Component{
public:
    Resistor(){
        G = 0;
    }

    Resistor(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        if(info[3] == "INFINITE"){
            G = 0;
        }
        else{
            G.real(1/string_to_float(info[3]));
        }
    }

};

// class for capacitor
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

// class for inductor
class Inductor: public Component{
public:
    Inductor(){
        G = 0;
    }

    Inductor(std::vector<std::string> info){
        CompName = info[0];
        NodeA = info[1];
        NodeB = info[2];
        G.imag(-1/(w * string_to_float(info[3])));
    }
};

// class for diode
class Diode: public Component{
public:
    std::string resistor_rd;
    std::string Id;

    Diode(std::vector<std::string> info){
        CompName = info[0];

        resistor_rd.append("Rd_");
        resistor_rd.append(info[0]);
        resistor_rd.push_back(' ');
        resistor_rd = resistor_rd + (info[1]);
        resistor_rd.push_back(' ');
        resistor_rd = resistor_rd + (info[2]);
        resistor_rd.push_back(' ');
        resistor_rd.append("1.6186e9"); // rd = VT/Id

        Id.append("Id_");
        Id.append(info[0]);
        Id.push_back(' ');
        Id = Id + (info[1]);
        Id.push_back(' ');
        Id = Id + (info[2]);
        Id.push_back(' ');
        Id.append("AC(1.5445e-11 0)"); // Id = 10^-12*(e^(0.7/0.025)-1)
    }
};

// class for BJT
class BJT: public Component{
public:
    //Q1 N003 N001 0 NPN

    std::string resistor_rbe;
    std::string resistor_r0;
    std::string Ic;

    BJT(std::vector<std::string> info){
        CompName = info[0];

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


        // details for NPN (2N2222) BJT
        if(info[4] == "NPN"){
            // gm = Ic/VT = 800mA/25mV = 32;
            // Rbe = beta/gm = 200/32 = 6.25;
            // r0 = VA/Ic = 100/0.8 = 125;
            Ic.append("Ic_");
            Ic.append(info[0]);
            Ic.push_back(' ');
            Ic = Ic + (info[1]);
            Ic.push_back(' ');
            Ic = Ic + (info[3]);
            Ic.push_back(' ');

            resistor_rbe.append("6.25");
            resistor_r0.append("125");
            Ic.append("AC(0.8 0)");
        }

        // details for PNP (2N2907) BJT
        if(info[4] == "PNP"){
            // gm = Ic/VT = 600mA/25mV = 24;
            // Rbe = beta/gm = 250/24 = 10.4167;
            // r0 = VA/Ic = 120/0.6 = 200;
            
            Ic.append("Ic_");
            Ic.append(info[0]);
            Ic.push_back(' ');
            Ic = Ic + (info[3]);
            Ic.push_back(' ');
            Ic = Ic + (info[1]);
            Ic.push_back(' ');

            resistor_rbe.append("10.4167");
            resistor_r0.append("200");
            Ic.append("AC(0.6 0)");
        }
    }
};

// class for MOSFET
class MOSFET: public Component{
public:
    std::string resistor_rgs;
    std::string resistor_r0;
    std::string Id;

    MOSFET(std::vector<std::string> info){
        CompName = info[0];

        // between gate and source is open circuit with potential difference of Vgs
        // represent as a resistor with infinite resistance;

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

        // details for NMOS (BSB012N03LX3) MOSFET
        if(info[4] == "NMOS"){
            // gm = 2K(Vgs-Vt)
            // r0 = VA/Id = ???;
            // Id = gm * Vgs;

            // Id = 0.2A
            // Vt = 2V
            // k = u*Cou*(W/L) = 0.713m
            // gm = 2k(Vgs - Vt) = 2 ( k * Id)^0.5
            // Id = k * (Vgs - Vt)^2 = 1.2m
            // r0 = 1/(lamtha * Id) = 1/0.09 *1.2m = 9260
            Id.append("Id_");
            Id.append(info[0]);
            Id.push_back(' ');
            Id = Id + (info[1]);
            Id.push_back(' ');
            Id = Id + (info[3]);
            Id.push_back(' ');

            resistor_r0.append("9260");
            Id.append("AC(0.001 0)");
        }

        // details for PMOS (FDS4435A) MOSFET
        if(info[4] == "PMOS"){
            // Id = -0.2A (continuous drain current)
            // rds(on) = 14 ohms when Vgs = -10V
            
            Id.append("Id_");
            Id.append(info[0]);
            Id.push_back(' ');
            Id = Id + (info[3]);
            Id.push_back(' ');
            Id = Id + (info[1]);
            Id.push_back(' ');

            resistor_r0.append("16670");
            Id.append("AC(0.0012 0)");
        }
    }
};

// class for voltage controlled current source
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

        G = string_to_float(info[5]);
    }
};
