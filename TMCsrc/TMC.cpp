#include "TMC.h"
#include "hyp_2F1.h"

double SF_2F1(double a, double b, double c, double x){
    // std::cout << a<< " " << b << " " << c<< " "<< x<< " ";
    complex<double> compa (a, 0), compb(b, 0), compc(c, 0), compx(x, 0), compres;

    compres = hyp_2F1(a, b, c, x); 
    //std::cout << compres.real () << " "<<compres.imag()<< "\n";
    return compres.real();
};

double TMCParameterization::H2oF10 (double xi){
    double ft = (1.0-xi)/xi; 
    double cst = (pars.at("lambda1") + pars.at("delta1")/(log(1.0+xi)) ) * (1.0-xi)/(xi*xi) ;
    double ser = 0.0; 
    for (int j=1; j<=jmax; j++ ){
        ser = ser+ (pow(-1, j)/(tgamma(j+1)*(j+1)) ) * SF_2F1(2, j+1, j+2, 1.0-1.0/xi); //note n! = gamma(n+1)
    }; 
    return 2.0*xi*(ft+ cst * ser);
}

double TMCParameterization::G2oF10 (double xi) {
    double ft = -log(xi)-(1.0-xi);
    double cst = (pars.at("lambda1") + pars.at("delta1")/(log(1.0+xi)) ) * pow((1.0-xi)/xi, 2) ;
    double ser = 0.0; 
    for (int j=1; j<=jmax; j++ ){
        ser = ser+ (pow(-1, j)/(tgamma(j+1)*(j+2)) ) * SF_2F1(2, j+2, j+3, 1.0-1.0/xi); //note n! = gamma(n+1)
    }; 
    return 2.0*xi*(ft+ cst * ser);
}

double TMCParameterization::H2oF20 (double xi) {
    double ft = (1.0-xi)/xi; 
    double cst = (pars.at("lambda2") + pars.at("delta2")/(log(1.0+xi)) ) * (1.0-xi)/(xi*xi) ;
    double ser = 0.0; 
    for (int j=1; j<=jmax; j++ ){
        ser = ser+ (pow(-1, j)/(tgamma(j+1)*(j+1)) ) * SF_2F1(2, j+1, j+2, 1.0-1.0/xi); //note n! = gamma(n+1)
    }; 
    return ft+ cst * ser;
}

double TMCParameterization::G2oF20 (double xi) {
    double ft = -log(xi)-(1.0-xi);
    double cst = (pars.at("lambda2") + pars.at("delta2")/(log(1.0+xi)) ) * pow((1.0-xi)/xi, 2) ;
    double ser = 0.0; 
    for (int j=1; j<=jmax; j++ ){
        ser = ser+ (pow(-1, j)/(tgamma(j+1)*(j+2)) ) * SF_2F1(2, j+2, j+3, 1.0-1.0/xi); //note n! = gamma(n+1)
    }; 
    return ft+ cst * ser;

}

double TMCParameterization::H3oF30 (double xi) {
    double ft = -log(xi);
    double cst = (pars.at("lambda3") + pars.at("delta3")/(log(1.0+xi)) ) * (1.0-xi)/xi ;
    double ser = 0.0; 
    for (int j=1; j<=jmax; j++ ){
        ser = ser+ (pow(-1, j)/(tgamma(j+1)*(j+1)) ) * SF_2F1(1, j+1, j+2, 1.0-1.0/xi); //note n! = gamma(n+1)
    }; 
    return ft+ cst * ser;

}


// OPE Implementations
double TMCs::getOPELeading_TMCF1 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F1= F123->getMasslessF1(xi, Q);
    return (x/(xi*r))*F1;
};

double TMCs::getOPELeading_TMCF2 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F2= F123->getMasslessF2(xi, Q);
    return (pow(x/xi, 2)/pow(r, 3))*F2;
};

double TMCs::getOPELeading_TMCF3 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F3 = F123->getMasslessF3(xi, Q);
    return (x/(xi*r*r))*F3;
};

double TMCs::getPARAM_TMCF1(double x, double Q) {
    double F1l = TMCs::getOPELeading_TMCF1(x, Q); 
    return F1l * tmcParams.getF1TMCoF1Leading(x, Q);
}

double TMCs::getPARAM_TMCF2(double x, double Q) {
    double F2l = TMCs::getOPELeading_TMCF2(x, Q); 
    return F2l * tmcParams.getF2TMCoF2Leading(x, Q);
}

double TMCs::getPARAM_TMCF3(double x, double Q) {
    double F3l = TMCs::getOPELeading_TMCF3(x, Q); 
    return F3l * tmcParams.getF3TMCoF3Leading(x, Q);
}

double TMCs::getOPE_TMCF1 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F = F123->getMasslessF1(xi, Q);
    double H2 = getH2(xi, Q); 
    double G2 = getG2 (xi, Q); 
    return (x/(xi*r))*F+ pow(x*M/(r*Q) , 2)*H2+ 2.*pow(M/Q, 4)*pow(x/r, 3)*G2;
};

double TMCs::getOPE_TMCF2 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F = F123->getMasslessF2(xi, Q);
    double H2 = getH2(xi, Q); 
    double G2 = getG2 (xi, Q); 
    return (pow(x/xi, 2)/pow(r, 3)) * F + (6. * pow(M/Q, 2)*pow(x, 3)/pow(r, 4))* H2+ (12.0*pow(M*x/Q, 4)/pow(r, 5))* G2;
};

double TMCs::getOPE_TMCF3 (double x, double Q){
    xbj = x; 
    Qv = Q;
    double r = pow(1.0+4.*pow(x*M/Q, 2), 0.5); 
    double xi = 2.*x/(1.+r);
    double F = F123->getMasslessF3(xi, Q);
    double H3 = getH3(xi, Q);
    return (x/(xi*r*r)) * F + 2.0* (pow(x*M/Q, 2)/pow(r, 3)) * H3;
};

//Implementation for the integrands
int integrandForH1 (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    TMCs *tmc = ((TMCs  *) fdata);
    fval[0] = 2.* tmc->getF123engine()->getMasslessF1(x[0], tmc->Qv)/x[0];
    return 0;
}

int integrandForH2 (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    TMCs *tmc = ((TMCs *) fdata);
    fval[0] = tmc->getF123engine()->getMasslessF2(x[0], tmc->Qv)/pow(x[0], 2);
    return 0;
}

int integrandForH3 (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    TMCs *tmc = ((TMCs  *) fdata);
    fval[0] = tmc->getF123engine()->getMasslessF3(x[0], tmc->Qv)/x[0];
    return 0;
}

int integrandForG2 (unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    TMCs *tmc = ((TMCs  *) fdata);
    fval[0] = (x[0]-tmc->xbj)*tmc->getF123engine()->getMasslessF2(x[0], tmc->Qv)/pow(x[0], 2);
    return 0;
}

