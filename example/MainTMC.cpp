#include "TMC.h"
#include <vector>

using namespace std;


class F123M0_dummy : public F123M0 {
    public : 
        F123M0_dummy (){};
        ~ F123M0_dummy(){};
        
        double getMasslessF1(double x, double Q) override {
            return pow(x, 0.6)*pow(1-x, 3)*log(Q);
        }
        double getMasslessF2(double x, double Q) override {
            return pow(x, 0.6)*pow(1-x, 3)*log(Q)/(2*x);
        };
        double getMasslessF3(double x, double Q) override{
            return pow(x, 0.6)*pow(1-x, 3)*log(Q);
        };
};


int main() {
    F123M0_dummy FDummy;
    TMCs dummyTMCs (&FDummy); 

    std::vector<double> xl = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.575, 0.60, 0.625, 
                0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95};
     std::vector<double> Qlist = {1.3, 1.5, 2., 3., 4., 6.};
     std::string pref = "TMC_Dummy_";
     ////////////////////////////////////////////////////////////////////////////////////////////////
    
    for (auto& Q : Qlist){
        std::cout << pref << "_Q_"<< Q << " = {";
        for (auto& x : xl){
           // double Ftmc = dummyTMCs.getOPE_TMCF2(x, Q);
           // double Ftmc = dummyTMCs.getOPELeading_TMCF2(x, Q);
            double Ftmc = dummyTMCs.getPARAM_TMCF2 (x, Q);
            std::cout<< Ftmc<<", "; 
        }; 
        std::cout << "}\n";
    };
}
