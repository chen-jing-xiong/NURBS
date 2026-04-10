#include "NURBSCurve2D.h"

void NURBSCurve2D:: evaluateValue(double t, double& x, double& y) {
        int span = basis_->findSpan(t);//查找范围
        double R[4]; 
        basis_->evaluateBasis(span, t, R);
        x = 0.0; y = 0.0;
        for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        x += R[i] * px_[idx];
        y += R[i] * py_[idx];
        }
    }

void NURBSCurve2D:: evaluateDerivative(double t, double& dx, double& dy) {
    int span = basis_->findSpan(t);
    double dR[4]; 
    basis_->evaluateDerBasis(span, t, dR);
    dx = 0.0;
    dy = 0.0;
    for (int i = 0; i <= degree_; ++i) {
        int idx = span - degree_ + i;
        dx += dR[i] * px_[idx];
        dy += dR[i] * py_[idx];
    }
}

void NURBSCurve2D::evaluateAll(double t, double& x, double& y, double& dx, double& dy) {
    int span = basis_->findSpan(t);
    double N[4];
    double dN[4];
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

void NURBSCurve2D::evaluateBatch(const double* t_values, int n_samples, double* out_x, double* out_y){
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

double NURBSCurve2D::getJacobian(double t)  {
        double dx, dy;
        evaluateDerivative(t, dx, dy);
        return std::sqrt(dx * dx + dy * dy);
    }