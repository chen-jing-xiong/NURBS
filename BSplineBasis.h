#ifndef BSPLINE_BASIS_H
#define BSPLINE_BASIS_H

class BSplineBasis {
private:
    int degree_;            //阶数
    int numKnots_;          //节点数
    double* knots_;         //节点向量
    bool isUniform_;        // 是否为均匀节点
    double h_;              // 节点间距（均匀节点时使用）

public:
    BSplineBasis(int degree, const double* knots, int numKnots,bool isUniform)
    :degree_(degree), numKnots_(numKnots), isUniform_(isUniform)
    {
        knots_ = new double[numKnots];
        for(int i=0; i<numKnots; ++i){
            knots_[i] = knots[i];
        }
         // 如果是均匀节点，计算步长
        if (isUniform_) {
            h_ = knots_[degree_ + 1] - knots_[degree_];
        } else {
            h_ = 0.0;
        }
    }
    ~BSplineBasis() {
    delete[] knots_; 
    }

    int findSpan(double t);//搜索值位于的区间范围
    void basisP2(int span, double t, double N[3]);//2阶样条函数
    void basisP2Der(int span, double t, double dN[3]);//2阶样条函数导数
    void basisP2WithDer(int span, double t,double N[3], double dN[3]);

    void basisP3(int span, double t, double N[4]);//3阶样条函数
    void basisP3Der(int span, double t, double dN[4]);//3阶样条函数导数
    void basisP3WithDer(int span, double t, double N[4], double dN[4]);
    void evaluateBasis(int span, double t, double* N);
    void evaluateDerBasis(int span, double t, double* dN);
    void evaluateAllBasis(int span, double t, double*N,double*dN);
    double getValue(int j, double t);
    double getDerValue(int j, double t);

    int getDegree() {return degree_;}

    int getNumKnots() {
        return numKnots_;
    }

    double* getKnots(){
        return knots_;
    }

    double getKnot(int i){
        // 实际代码中建议加越界检查
        return knots_[i];
    }
};

#endif