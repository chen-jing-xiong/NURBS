#include "BSplineCurve2D.h"

void BSplineCurve2D:: evaluateValue(double t, double& x, double& y) {
        int span = basis_->findSpan(t);//查找范围
        double N[4]; 
        basis_->evaluateBasis(span, t, N);
        x = 0.0; y = 0.0;
        for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        x += N[i] * px_[idx];
        y += N[i] * py_[idx];
        }
    }

void BSplineCurve2D:: evaluateDerivative(double t, double& dx, double& dy) {
    int span = basis_->findSpan(t);
    double dN[4]; 
    basis_->evaluateDerBasis(span, t, dN);
    dx = 0.0;
    dy = 0.0;
    for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        dx += dN[i] * px_[idx];
        dy += dN[i] * py_[idx];
    }
}

void BSplineCurve2D::evaluateAll(double t, double& x, double& y, double& dx, double& dy) {
    int span = basis_->findSpan(t);
    double N[4] = {0.0};
    double dN[4] = {0.0};
    basis_->evaluateAllBasis(span, t, N, dN);
    x = 0.0; y = 0.0;
    dx = 0.0; dy = 0.0;
    for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        // 位置累加
        x += N[i] * px_[idx];
        y += N[i] * py_[idx];        
        // 导数累加
        dx += dN[i] * px_[idx];
        dy += dN[i] * py_[idx];
    }
}

void BSplineCurve2D::evaluateBatch(const double* t_values, int n_samples, double* out_x, double* out_y){
    double N[4]; 
    for (int i = 0; i < n_samples; ++i) {
        double t = t_values[i];
        int span = basis_->findSpan(t);
        
        basis_->evaluateBasis(span, t, N);

        double x = 0.0, y = 0.0;
        for (int j = 0; j <= degree_; ++j) {
            int idx = span - degree_ + j;
            x += N[j] * px_[idx];
            y += N[j] * py_[idx];
        }
        
        out_x[i] = x;
        out_y[i] = y;
    }
}

double BSplineCurve2D::getJacobian(double t)  {
        double dx, dy;
        evaluateDerivative(t, dx, dy);
        return std::sqrt(dx * dx + dy * dy);
    }