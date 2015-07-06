#include "Rcpp.h"
#include "Rmath.h"
#include <cmath>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericVector likelihood(NumericVector k, NumericVector alpha,
NumericVector intercept, NumericVector sf, NumericVector betas,
NumericVector cond){
    NumericVector z(betas.size());
    for (int beta_i = 0; beta_i < betas.size(); beta_i++){
        double temp = 1.0;
        for (int i = 0; i < k.size(); i++){
            temp *= R::dnbinom_mu(k[i], 1/alpha[0],
                 sf[i] * pow(2, intercept[0] + betas[beta_i] * cond[i]), 0);//no log
        }
        z[beta_i] = temp;
    }
    double sum = 0.0;
    for (int i = 0; i < betas.size(); i++){
        sum += z[i] * (betas[1] - betas[0]);
    }
    for (int i = 0; i < z.size(); i++){
        z[i] /= sum;
    }
    return z;
}

// [[Rcpp::export]]
Rcpp::NumericVector fastSum(Rcpp::NumericVector x1, Rcpp::NumericVector y1){
    Rcpp::NumericVector x(x1);
    Rcpp::NumericVector y(y1);

    double c = 0;
	// double d = 0;
	double e = 0;
	for (int i = 0; i < x.size(); i++){
		for (int j = 0; j < i; j++){
			c += x[i] * y[j];
		}
		for (int j = i+1; j < y.size(); j++){
		    e += x[i] * y[j];
		}
		//d += x[i] * y[i]; // if it is of no use, no calculation
	}
	double ret = c < e ? c : e;
    Rcpp::NumericVector ret1(1);
    ret1[0] = 2 * ret;
    return ret1;
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
    double currentP = 1.0;
    for (int i = 0; i<k1.size(); i++){
    	currentP *= min(R::dnbinom_mu(k1[i],1/alpha1[0],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0),
    				R::dnbinom_mu(k2[i],1/alpha2[0],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0));
    }
    while (true){
    	currentFC -= stepFC[0];
    	double tmp = 1.0;
    	for (int i = 0; i<k1.size();i++){
    		tmp *= min(R::dnbinom_mu(k1[i],1/alpha1[0],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0),
    			R::dnbinom_mu(k2[i],1/alpha2[0],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0));
        }
    		if( ((currentP - tmp) <= toleration[0]) & ((currentP - tmp) >=0)){
    			break;
    		}else{
    			currentP = tmp;
    		}
    }
    ret[0] = currentFC - stepFC[0];
    //right side
    double right_fc = max(FC1[0],FC2[0]);
    currentFC = right_fc;
    currentP = 1.0;
    for (int i = 0; i<k1.size(); i++){
        currentP *= min(R::dnbinom_mu(k1[i],1/alpha1[0],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0),
                    R::dnbinom_mu(k2[i],1/alpha2[0],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0));
    }
    while (true){
        currentFC += stepFC[0];
        double tmp = 1.0;
        for (int i = 0; i<k1.size();i++){
            tmp *= min(R::dnbinom_mu(k1[i],1/alpha1[0],sf1[i] * pow(2, intercept1[0] + currentFC * cond1[i]),0),
                R::dnbinom_mu(k2[i],1/alpha2[0],sf2[i] * pow(2, intercept2[0] + currentFC * cond2[i]),0));
        }
            if( ((currentP - tmp) <= toleration[0]) & ((currentP - tmp) >=0)){
                break;
            }else{
                currentP = tmp;
            }
    }
    ret[1] = currentFC + stepFC[0];
    return ret;
}

// [[Rcpp::export]]
NumericVector posterior(NumericVector k, NumericVector alpha,
NumericVector intercept, NumericVector sf, NumericVector priorSigma, NumericVector betas,
NumericVector cond){
    NumericVector z = likelihood(k, alpha, intercept, sf, betas, cond);

    double diff = betas[1] - betas[0];
    double sum = 0.0;
    for (int i = 0; i < betas.size(); i++){
        z[i] *= R::dnorm((double)(betas[i]), 0.0, (double)(priorSigma[0]), 0);
        // z[i] *= R::dnorm(1, 0.0, 1, 0);
        sum += z[i] * diff;
    }
    for (int i = 0; i < betas.size(); i++){
        z[i] = z[i] / sum * diff; //multiply diff to reduce time to consume
    }

    return z;
}


// [[Rcpp::export]]
NumericVector overlapCoef(NumericVector a, NumericVector b){
    double asum = 0.0;
    double minSum = 0.0;
    for (int i = 0; i < a.size(); i++){
        minSum += min(a[i], b[i]);
        asum += a[i];
    }
    NumericVector ret(1);
    ret[0] = minSum / asum;
    return ret;
}



// [[Rcpp::export]]
NumericVector xtail_test(NumericVector k1, NumericVector k2, NumericVector ints1,
    NumericVector ints2, NumericVector log2FC_1, NumericVector log2FC_2,
    NumericVector priorSigma1, NumericVector priorSigma2, NumericVector disper1,
    NumericVector disper2, NumericVector sf1, NumericVector sf2, NumericVector cond1,
    NumericVector cond2, IntegerVector bin){

    NumericVector ret(2);
    if (isnan(log2FC_1[0]) || isnan(log2FC_2[0])){
       ret[1] = NA_REAL;
        ret[0] = NA_REAL;
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
    NumericVector log2_density1 = posterior(k1, disper1, ints1, sf1, priorSigma1, betas, cond1);
    NumericVector log2_density2 = posterior(k2, disper2, ints2, sf2, priorSigma2, betas, cond2);
    NumericVector ovl = overlapCoef(log2_density1, log2_density2);

    NumericVector pvalue = fastSum(log2_density1, log2_density2);
    ret[0] = ovl[0];
    ret[1] = pvalue[0];

    return ret;

}
