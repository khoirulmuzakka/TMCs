/**
 * @file TMC.h
 *
 * Created by Faiq on 10.03.2022
 * 
 */


#ifndef TMCS_H
#define TMCS_H
#include <iostream>
#include "cubature.h"
#include <map>
#include <cmath>

using namespace std;

//Forward declaration for the integrands, so that cubature integrator know these functions in advance
int integrandForH1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int integrandForH2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int integrandForH3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
int integrandForG2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

double SF_2F1(double a, double b, double c, double x);

/**
 * @brief This represents a purely virtual class which provide common interface to the TMC class. 
 * In order to use the TMCs class, one must create a class which inherit F123M0 and implements all the virtula functions.
 * 
 */
class F123M0 {
    public : 
        F123M0(){}; 
        ~F123M0(){};

        /**
         * @brief Virtual function to quary massless (Mp=0) F1 structure function
         */
        virtual double getMasslessF1(double x, double Q)=0;

        /**
         * @brief Virtual function to quary massless (Mp=0) F2 structure function
         */
        virtual double getMasslessF2(double x, double Q) =0;

        /**
         * @brief Virtual function to quary massless (Mp=0) F3 structure function
         */
        virtual double getMasslessF3(double x, double Q)=0;
};


/**
 * @brief Class that implements the 2F1-based TMC parameterizations
 */
class TMCParameterization{
    private : 
        int jmax= 4; //maximal of j values. jmax=4 is enough, increasing the value of jmax>4 WILL NOT improve the parameterizations
        map<std::string, double> pars = { //the value of the parameters from the TMC paper
            {"lambda1" , 2.352}, {"delta1" , -0.122}, 
            {"lambda2", 2.264} , {"delta2", -0.074}, 
            {"lambda3", 2.090} , {"delta3", 0.035}
        }; 
        double Mp = 0.938; //proton mass
        
    public : 
        TMCParameterization() {}; 
        ~TMCParameterization(){};

        /**
         * @brief Query the ratio of TMC F1 over leaidng TMC F1 obtained from the 2F1-parameterizations
         */
        double getF1TMCoF1Leading (double x, double Q){
            double r = pow(1.0+4.*pow(x*Mp/Q, 2), 0.5); 
            double xi = 2.*x/(1.+r);
            double ft = pow(Mp/Q, 2) * (x*xi/r) * H2oF10(xi); 
            double st = pow(Mp/Q, 4) * (2*xi* pow(x/r, 2)) * G2oF10(xi); 
            return 1.0+ft+st;
        }

        /**
         * @brief Query the ratio of TMC F2 over leading TMC F2 obtained from the 2F1-parameterizations
         */
        double getF2TMCoF2Leading (double x, double Q){
            double r = pow(1.0+4.*pow(x*Mp/Q, 2), 0.5); 
            double xi = 2.*x/(1.+r);
            double ft = pow(Mp/Q, 2) * (pow(xi, 2)*6*x/r) * H2oF20(xi); 
            double st = pow(Mp/Q, 4) * (12*pow(xi*x/r, 2)) * G2oF20(xi); 
            return 1.0+ft+st;
        } 

        /**
         * @brief Query the ratio of TMC F3 over leading TMC F3 obtained from the 2F1-parameterizations
         */
        double getF3TMCoF3Leading (double x, double Q){
            double r = pow(1.0+4.*pow(x*Mp/Q, 2), 0.5); 
            double xi = 2.*x/(1.+r);
            double ft = pow(Mp/Q, 2) * (2.0*xi*x/r) * H3oF30(xi); 
            return 1.0+ft;
        }

    protected :
        //The full implementation of 2F1-parameterizations 
        double H2oF10 (double xi);
        double G2oF10 (double xi);
        double H2oF20 (double xi);
        double G2oF20 (double xi);
        double H3oF30 (double xi);
};   

/**
 * @brief This class provides all the functionality to query : leading TMC, the full TMC, and the parameterized TMC
 */
class TMCs {
    private : 
        F123M0* F123;
        double M = 0.938;
        TMCParameterization tmcParams;

    public : 
        F123M0* getF123engine() { return F123;};
        double xbj; 
        double Qv;

    public : 
        TMCs (F123M0* f123) : F123(f123) {}; 
        ~TMCs(){}; 

        /**
         * @brief Function to query the leading TMC F1 structure function
         */
        double getOPELeading_TMCF1 (double x, double Q);

        /**
         * @brief Function to query the leading TMC F2 structure function
         */
        double getOPELeading_TMCF2 (double x, double Q);

        /**
         * @brief Function to query the leading TMC F3 structure function
         */
        double getOPELeading_TMCF3 (double x, double Q);

        /**
         * @brief Function to query the full TMC F1 structure function
         */
        double getOPE_TMCF1 (double x, double Q);
    
        /**
         * @brief Function to query the full TMC F2 structure function
         */
        double getOPE_TMCF2 (double x, double Q);

        /**
         * @brief Function to query the full TMC F3 structure function
         */
        double getOPE_TMCF3 (double x, double Q);


        /**
         * @brief Function to query the parameterized TMC F1 structure function
         */
        double getPARAM_TMCF1 (double x, double Q);

        /**
         * @brief Function to query the parameterized TMC F2 structure function
         */
        double getPARAM_TMCF2 (double x, double Q);

        /**
         * @brief Function to query the parameterized TMC F3 structure function
         */
        double getPARAM_TMCF3 (double x, double Q);



    private : 
        double integrateSF (int mode, double x, double Q){
            const double xMin[]={x};
            const double xMax[] = {1.0};
            int maxEval = 50000;
            double error = 0.0;
            double relError = 1e-4;
            double res=0.0;
            double resErr=0.0;
            if (mode==1){
                hcubature(1, integrandForH1, this, 1, xMin, xMax,
                                maxEval,//< no limit number of calls
                                error, relError, ERROR_INDIVIDUAL,//< accuracy
                                &res,  &resErr);

                return res;
            } else if (mode==2){
                hcubature(1, integrandForH2, this, 1, xMin, xMax,
                                maxEval,//< no limit number of calls
                                error, relError, ERROR_INDIVIDUAL,//< accuracy
                                &res,  &resErr);

                return res;
            } else if (mode==3){
                hcubature(1, integrandForH3, this, 1, xMin, xMax,
                                maxEval,//< no limit number of calls
                                error, relError, ERROR_INDIVIDUAL,//< accuracy
                                &res,  &resErr);

                return res;
            } else if (mode==4){
                hcubature(1, integrandForG2, this, 1, xMin, xMax,
                                maxEval,//< no limit number of calls
                                error, relError, ERROR_INDIVIDUAL,//< accuracy
                                &res,  &resErr);

                return res;
            } else {
                std::cout << "Unknown mode!"; 
                exit(1);
            };
        };


        double getH1 ( double x, double Q){
            return integrateSF(1, x, Q);
        }
        double getH2 ( double x, double Q){
            return integrateSF(2,x, Q);
        }
        double getH3 ( double x, double Q){
            return integrateSF(3,x, Q);
        }
        double getG2 ( double x, double Q){
            return integrateSF(4,x, Q);
        }

}; 





#endif