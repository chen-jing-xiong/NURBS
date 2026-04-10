#pragma once
#include <cmath>  
#include "BSplineBasis.h"


class BSplineCurve2D{
private:
    BSplineBasis* basis_;   // 基函数类对象
    int numPoints_;         // 控制点数量
    double* px_;            // x方向坐标数组
    double* py_;            // y方向坐标数组
    int degree_;
public:
  // 构造函数
    BSplineCurve2D(BSplineBasis* basis, int n,  double* x, double* y) 
    : basis_(basis), numPoints_(n) 
    {
        degree_= basis->getDegree();
        px_ = new double[numPoints_];
        py_ = new double[numPoints_];
        
        for (int i = 0; i < numPoints_; ++i) {
            px_[i] = x[i];
            py_[i] = y[i];
        }
    }

    
    // 析构函数
    ~BSplineCurve2D() {
        delete[] px_;
        delete[] py_;
        delete basis_;
    }
    
    void evaluateValue(double t, double& x, double& y);
    void evaluateDerivative(double t, double& dx, double& dy);
    void evaluateAll(double t, double& x, double& y, double& dx, double& dy);

    void evaluateBatch(const double* t_values, int n_samples, double* out_x, double* out_y);
    double getJacobian(double t);
    
    int getnumber(){return numPoints_;}

    BSplineBasis* getBasis() {
        return basis_;
    }

    int getDegree() const {
        return degree_;
    }
};
