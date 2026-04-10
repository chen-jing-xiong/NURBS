#ifndef NURBS_BASIS_H
#define NURBS_BASIS_H

#include "BSplineBasis.h"

class NURBSBasis {
private:
    BSplineBasis* bspline_; // 使用指针持有 B-Spline 对象
    double *weights_; // 存储权重
    int numWeights_;
public:
    NURBSBasis(BSplineBasis* bspline, double* weights, int numWeights)
        : bspline_(bspline),weights_(weights),numWeights_(numWeights) {
        
    }


    void evaluateBasis(int span, double t, double* R);//计算当前区域所有基函数值
    void evaluateDerBasis(int span, double t, double* dR);
    void evaluateAllBasis(int span, double t, double*R,double*dR);
    double getValue(int j, double t);
    int getDegree() const { return bspline_->getDegree(); }
    double* getKnots(){
        return bspline_->getKnots();
    }

    int getNumWeights(){ return numWeights_; }
    int findSpan(double t) const{return bspline_->findSpan(t);}
};
#endif