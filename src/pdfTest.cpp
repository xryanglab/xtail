#include "Rcpp.h"
#include "Rmath.h"
#include <cmath>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector probDensity(NumericVector k, NumericVector alpha,
NumericVector intercept, NumericVector sf, NumericVector betas,
NumericVector cond){
    // estimate the probability density of beta.
    double sum = 0.0;
    NumericVector z(betas.size());
    for (int beta_i = 0; beta_i < betas.size(); beta_i++){
        double temp = 1.0;
        for (int i = 0; i < k.size(); i++){
            temp *= R::dnbinom_mu(k[i], 1/alpha[i],
                 sf[i] * pow(2, intercept[0] + betas[beta_i] * cond[i]), 0);//no log
        }
        z[beta_i] = temp;
        sum += temp;
    }
    for (int i = 0; i < z.size(); i++){
        z[i] /= sum;
    }

    return z;
}

// [[Rcpp::export]]
Rcpp::NumericVector probMatrix(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector side){
    // Xtail generates a joint probability matrix by multiplying two probability density dsitributions. 
    // This return pvalue, 

    double pmax = 0;
    int pmaxi = 0; //max probability index 
    double drecord = 0;
    double s = 0; // 2*s is final pvalue
    NumericVector ret(2);

    if (side[0] >= 0){
        //lower triangle,include the diagonal
        for (int k=x.size()-1;k>-1;k--){
            int j = 0;
            int i = k;
            drecord = 0;
            while(i<x.size()){
                drecord += x[i] * y[j];
                ++j;
                ++i;
            }
            if (pmax<=drecord){ //if there many same pmax, get the morst one.
                pmax = drecord;
                pmaxi = k;
            }
        }
        //upper triangle
        for (int k=1;k<x.size();k++){
            int i = 0;
            int j = k;
            drecord = 0;
            while(j<y.size()){
                drecord += x[i] * y[j];
                ++j;
                ++i;
            }
            if (pmax<drecord){
                pmax = drecord;
                pmaxi = -k;
            }
            s += drecord;
        }
        ret[0] = pmaxi; // index of max probability
        if (s > 0.5){
            s = 1 - s;
        }
        ret[1] = 2 * s;

    }else{
        // upper triangle, include the diagonal
        for (int k=y.size()-1;k>-1;k--){
            int j=k;
            int i=0;
            drecord = 0;
            while(j<y.size()){
                drecord += x[i]*y[j];
                ++j;
                ++i;
            }
            if (pmax<=drecord){
                pmax = drecord;
                pmaxi = -k;
            }
        }
        // lower triangle
        for (int k=1;k<x.size();k++){
            int i=k;
            int j=0;
            drecord = 0;
            while(i<x.size()){
                drecord += x[i]*y[j];
                ++j;
                ++i;
            }
            if (pmax<drecord){
                pmax = drecord;
                pmaxi = k;
            }
            s += drecord;
        }
        ret[0] = pmaxi;
        if (s>0.5){
            s = 1 - s;
        }
        ret[1] = 2 * s;
    }
    return ret;  // first is index of betas, thecond is pvalue.

}

// [[Rcpp::export]]
NumericVector probMatrixCI(NumericVector x, NumericVector y, NumericVector side, NumericVector ci){
    // Xtail generates a joint probability matrix by multiplying two probability density dsitributions. 
    // This return pvalue and CI

    double pmax = 0;
    int pmaxi = 0; //max probability index 
    double drecord = 0;
    double s = 0; // 2*s is final pvalue
    double drecords = 0; 

    NumericVector ret(4); // return index of betas, pvalue,lowCI index, highCI index.
    int lowCI = -1;
    int highCI = -1;

    if (side[0] >= 0){
        //lower triangle,include the diagonal
        for (int k=x.size()-1;k>-1;k--){
            int j = 0;
            int i = k;
            drecord = 0;
            while(i<x.size()){
                drecord += x[i] * y[j];
                ++j;
                ++i;
            }
            if (pmax<=drecord){ //if there many same pmax, get the morst one.
                pmax = drecord;
                pmaxi = k;
            }
            drecords += drecord;
            if (lowCI == -1){
                if (drecords >= (1-ci[0])/2){
                    lowCI = k;
                }
            }
            if (highCI == -1){
                if (drecords >= (ci[0]+(1-ci[0])/2)){
                    highCI = k;
                }
            }

        }
        //upper triangle
        for (int k=1;k<x.size();k++){
            int i = 0;
            int j = k;
            drecord = 0;
            while(j<y.size()){
                drecord += x[i] * y[j];
                ++j;
                ++i;
            }
            if (pmax<drecord){
                pmax = drecord;
                pmaxi = -k;
            }
            s += drecord;
            drecords += drecord;
            if (lowCI == -1){
                if (drecords >= (1-ci[0])/2){
                    lowCI = -k;
                }
            }
            if (highCI == -1){
                if (drecords >= ci[0]+(1-ci[0])/2){
                    highCI = -k;
                }
            }
        }
        ret[0] = pmaxi; // index of max probability
        if (s > 0.5){
            s = 1 - s;
        }
        ret[1] = 2 * s;
        ret[2] = lowCI;
        ret[3] = highCI;

    }else{
        // upper triangle, include the diagonal
        for (int k=y.size()-1;k>-1;k--){
            int j=k;
            int i=0;
            drecord = 0;
            while(j<y.size()){
                drecord += x[i]*y[j];
                ++j;
                ++i;
            }
            if (pmax<=drecord){
                pmax = drecord;
                pmaxi = -k;
            }
            drecords += drecord;
            if(lowCI == -1){
                if (drecords >= (1-ci[0])/2){
                    lowCI = -k;
                }
            }
            if(highCI == -1){
                if (drecords >= ci[0] +(1-ci[0])/2){
                    highCI = -k;
                }
            }
        }
        // lower triangle
        for (int k=1;k<x.size();k++){
            int i=k;
            int j=0;
            drecord = 0;
            while(i<x.size()){
                drecord += x[i]*y[j];
                ++j;
                ++i;
            }
            if (pmax<drecord){
                pmax = drecord;
                pmaxi = k;
            }
            drecords += drecord;
            s += drecord;
            if(lowCI == -1){
                if(drecords >= (1-ci[0])/2){
                    lowCI = k;
                }
            }
            if(highCI == -1){
                if(drecords >= ci[0]+(1-ci[0])/2){
                    highCI = k;
                }
            }
        }
        ret[0] = pmaxi;
        if (s > 0.5){
            s = 1 - s;
        }
        ret[2] = lowCI;
        ret[3] = highCI;
        ret[1] = 2 * s;
    }
    return ret;
}

// [[Rcpp::export]]
NumericVector boundFC(NumericVector k1, NumericVector alpha1, NumericVector intercept1,
    NumericVector sf1, NumericVector FC1, NumericVector cond1, NumericVector k2,  
    NumericVector alpha2, NumericVector intercept2, NumericVector sf2, NumericVector FC2,
    NumericVector cond2, NumericVector stepFC, NumericVector toleration){
    	
    NumericVector ret(2);

    //left side
    double left_fc = min(FC1[0],FC2[0]);
    double currentFC = left_fc;
    double currentP = 0;
    double P1 = 1.0;
    double P2 = 1.0;

    while (true){
    	currentFC -= stepFC[0];
        P1 = 1.0;
        P2 = 1.0;
    	for (int i = 0; i<k1.size();i++){
    		P1 *= R::dnbinom_mu(k1[i],1/alpha1[i],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0);
    		P2 *= R::dnbinom_mu(k2[i],1/alpha2[i],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0);
        }
            if (currentP == max(P1,P2)){
                break;
            }else{
                currentP = max(P1,P2);
                if(currentP <= toleration[0]){
                    break;
                }
            }

    }
    ret[0] = currentFC;
    //right side
    double right_fc = max(FC1[0],FC2[0]);
    currentFC = right_fc;
    while (true){
        currentFC += stepFC[0];
        P1 = 1.0;
        P2 = 1.0;
        for (int i = 0; i<k1.size();i++){
            P1 *= R::dnbinom_mu(k1[i],1/alpha1[i],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0);
            P2 *= R::dnbinom_mu(k2[i],1/alpha2[i],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0);
        }
            if (currentP == max(P1,P2)){
                break;
            }else{
                currentP = max(P1,P2);
                if(currentP <= toleration[0]){
                    break;
                }
            }

    }
    ret[1] = currentFC;
    return ret;
}




// [[Rcpp::export]]
NumericVector xtail_test(NumericVector k1, NumericVector k2, NumericVector ints1,
    NumericVector ints2, NumericVector log2FC_1, NumericVector log2FC_2, NumericVector disper1,
    NumericVector disper2, NumericVector sf1, NumericVector sf2, NumericVector cond1,
    NumericVector cond2, IntegerVector bin, NumericVector ci){


    NumericVector ret(4); // deltaTE,  pval, ci_low, ci_high
    if (isnan(log2FC_1[0]) || isnan(log2FC_2[0])){
        ret[1] = NA_REAL;
        ret[0] = NA_REAL;
        ret[2] = NA_REAL;
        ret[3] = NA_REAL;
        return ret;
    }

    int bins=bin[0];
    NumericVector stepFC(1);
    stepFC[0] = 0.2;

    NumericVector tolerance(1);
    tolerance[0] = 1e-50;
    NumericVector left_right = boundFC(k1,disper1,ints1,sf1,log2FC_1,cond1,
                                        k2,disper2,ints2,sf2,log2FC_2,cond2,stepFC,tolerance);
    NumericVector betas(bins);
    for (int i = 0; i < bins; i++){
        betas[i] = left_right[0] + (left_right[1] - left_right[0])/(bins-1)*i;
    }
    NumericVector log2_density1 = probDensity(k1, disper1, ints1, sf1, betas, cond1);
    NumericVector log2_density2 = probDensity(k2, disper2, ints2, sf2, betas, cond2);

    NumericVector side(1);
    side[0] = log2FC_1[0] - log2FC_2[0];
   
    if (ci[0] == 0){
        NumericVector retValues = probMatrix(log2_density1, log2_density2, side);
        if(retValues[0] >= 0){
            ret[0] = betas[0] - betas[retValues[0]];
        }else{
            ret[0] = betas[-retValues[0]] - betas[0];
        }
        ret[1] = retValues[1];
        ret[2] = 0;
        ret[3] = 0;

    }else{

        NumericVector retValues = probMatrixCI(log2_density1, log2_density2, side,ci);
        // cout <<pvalue[1]<<endl;
        // cout<<pvalue[2]<<endl;
        double lowCI = 0;
        double highCI = 0;

        if(retValues[0] >= 0){
            ret[0] = betas[0] - betas[retValues[0]];
        }else{
            ret[0] = betas[-retValues[0]] - betas[0];
        }
        ret[1] = retValues[1];
        //lowCI = retValues[2];
        //highCI = retValues[3];
        // lowCI       
        if (retValues[2]>=0){
            highCI = betas[0] - betas[retValues[2]];
        }else{
            highCI = betas[-retValues[2]] - betas[0];
        }
        if (retValues[3]>=0){
            lowCI = betas[0] - betas[retValues[3]];
        }else{
            lowCI = betas[-retValues[3]] - betas[0];
        }

        ret[2] = min(lowCI,highCI);
        ret[3] = max(lowCI,highCI);
    }
    return ret;
}
