#include "NURBSBasis.h"

void NURBSBasis::evaluateBasis(int span, double t, double* R){
    int p =bspline_->getDegree();//获得阶数
    double N[4];
    bspline_->evaluateBasis(span,t,N);//计算当前区间上的所有非零B样条基函数
    int index =span-p;
    double w_sum=0.0;
    for(int i=0;i<=p;i++){
        w_sum+=N[i]*weights_[index+i];
    }
    for (int i = 0; i <= p; ++i) {
        R[i] = (N[i] * weights_[index+ i]) / w_sum;
    }
}

void NURBSBasis::evaluateDerBasis(int span, double t, double* dR){
    int p = bspline_->getDegree();
    double N[4];
    double dN[4];
    bspline_->evaluateAllBasis(span,t,N,dN);

    double W=0.0;
    double dW=0.0;
    int index = span - p;
    for(int i=0;i<=p;i++)
    {
        double w = weights_[index+i];
        W  += N[i]  * w;
        dW += dN[i] * w;
    }
    for(int i=0;i<=p;i++)
    {
        double w = weights_[index+i];

        double A  = N[i]  * w;
        double dA = dN[i] * w;

        dR[i] = (dA * W - A * dW) / (W * W);
    }
}

void NURBSBasis::evaluateAllBasis(int span, double t, double*R,double*dR){
    int p = bspline_->getDegree();
    double N[4];
    double dN[4];
    bspline_->evaluateAllBasis(span,t,N,dN);

    double W=0.0;
    double dW=0.0;
    int index = span - p;
    for(int i=0;i<=p;i++)
    {
        double w = weights_[index+i];
        W  += N[i]  * w;
        dW += dN[i] * w;
    }
    for(int i=0;i<=p;i++)
    {
        double w = weights_[index+i];

        double A  = N[i]  * w;
        double dA = dN[i] * w;

        R[i]  = A / W;
        dR[i] = (dA * W - A * dW) / (W * W);
    }
}

double NURBSBasis::getValue(int j, double t){
    int span = bspline_->findSpan(t);
    int p = bspline_->getDegree();
    // 非零范围检查
    if (j < span - p || j > span){
        return 0.0;
    }
    double N[4];   
    evaluateBasis(span, t, N);
    int localIndex = j - (span - p);
    return N[localIndex];
}