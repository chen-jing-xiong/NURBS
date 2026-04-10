#include "BSplineBasis.h"

int BSplineBasis::findSpan(double t) {
    if (isUniform_) {
        int span = int((t - knots_[degree_]) / h_) + degree_;
        if (span >= numKnots_ - degree_ - 1)
            span = numKnots_ - degree_ - 2;
        return span;
    }else{
        int p = degree_;             // B样条次数
        int m = numKnots_ - 1;      // 最大节点索引
        int n = m - p - 1;           // 最后一个非零基函数索引

        if (t <= knots_[p]) return p;       // 左端点
        if (t >= knots_[n + 1]) return n;   // 右端点（左连续）

        // 二分查找初始化
        int low = p;
        int high = n + 1;

        // mid计算用 safe formula 避免整数溢出
        int mid = low + (high - low) / 2;

        // 二分查找核心
        while (!(t >= knots_[mid] && t < knots_[mid + 1])) {
            if (t < knots_[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = low + (high - low) / 2;
        }
        return mid;
    }
}

void BSplineBasis::basisP2(int span, double t, double N[3]){
    const double u_im1 = knots_[span - 1]; //u_{i-1}
    const double u_i   = knots_[span];     //u_{i}
    const double u_ip1 = knots_[span + 1]; //u_{i+1}
    const double u_ip2 = knots_[span + 2]; //u_{i+2}
    // 区间长度
    const double width = (u_ip1 - u_i);

    if (width == 0.0){
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.0;
        return;
    }else{
        // 一阶基函数
        const double Ni_deg1   = ((t - u_i) / width);  //N_{i,1}
        const double Nim1_deg1 = ((u_ip1 - t) / width);//N_{i-1,1}

        // 二阶分母
        const double denom1 = (u_ip1 - u_im1);
        const double denom2 = (u_ip2 - u_i);
        // N_{span-2,2}
        if (denom1 != 0.0){
            N[0] = (((u_ip1 - t) * Nim1_deg1) / denom1);
        }else{
            N[0] = 0.0;
        }
        // N_{span-1,2}
        double term1 = 0.0;
        double term2 = 0.0;

        if (denom1 != 0.0){
            term1 = (((t - u_im1) * Nim1_deg1) / denom1);
        }

        if (denom2 != 0.0){
            term2 = (((u_ip2 - t) * Ni_deg1) / denom2);
        }
        N[1] = (term1 + term2);

        // N_{span,2}
        if (denom2 != 0.0){
            N[2] = (((t - u_i) * Ni_deg1) / denom2);
        }else{
            N[2] = 0.0;
        }
    }
}
void BSplineBasis::basisP2Der(int span, double t, double dN[3]) {
    const double u_im1 = knots_[span - 1];
    const double u_i   = knots_[span];
    const double u_ip1 = knots_[span + 1];
    const double u_ip2 = knots_[span + 2];

    const double width = (u_ip1 - u_i);

    if (width == 0.0) {
        dN[0] = 0.0;
        dN[1] = 0.0;
        dN[2] = 0.0;
        return;
    }

    // 一阶基函数
    const double Ni_deg1   = (t - u_i) / width;
    const double Nim1_deg1 = (u_ip1 - t) / width;

    const double denom1 = (u_ip1 - u_im1);
    const double denom2 = (u_ip2 - u_i);

    // N'_{span-2,2}
    if (denom1 != 0.0) {
        dN[0] = (-2.0 * Nim1_deg1) / denom1;
    } else {
        dN[0] = 0.0;
    }

    // N'_{span-1,2}
    double term1 = 0.0;
    double term2 = 0.0;

    if (denom1 != 0.0) {
        term1 = (2.0 * Nim1_deg1) / denom1;
    }

    if (denom2 != 0.0) {
        term2 = (-2.0 * Ni_deg1) / denom2;
    }

    dN[1] = term1 + term2;

    // N'_{span,2}
    if (denom2 != 0.0) {
        dN[2] = (2.0 * Ni_deg1) / denom2;
    } else {
        dN[2] = 0.0;
    }
}

void BSplineBasis::basisP2WithDer(int span, double t, double N[3], double dN[3]) {
    const double u_im1 = knots_[span - 1]; // u_{i-1}
    const double u_i   = knots_[span];     // u_{i}
    const double u_ip1 = knots_[span + 1]; // u_{i+1}
    const double u_ip2 = knots_[span + 2]; // u_{i+2}
    
    const double width = u_ip1 - u_i;

    if (width == 0.0) {
        N[0] = N[1] = N[2] = 0.0;
        dN[0] = dN[1] = dN[2] = 0.0;
        return;
    }

    // 计算一次公共的 1 阶基函数
    const double Ni_1   = (t - u_i) / width;
    const double Nim1_1 = (u_ip1 - t) / width;

    // 计算一次公共的 2 阶分母
    const double denom1 = u_ip1 - u_im1;
    const double denom2 = u_ip2 - u_i;

    // 同时计算 N[0] 和 dN[0] (N_{span-2,2})
    if (denom1 != 0.0) {
        double factor1 = Nim1_1 / denom1; // 提取公共项
        N[0]  = (u_ip1 - t) * factor1;
        dN[0] = -2.0 * factor1;
    } else {
        N[0]  = 0.0;
        dN[0] = 0.0;
    }

    // 同时计算 N[2] 和 dN[2] (N_{span,2})
    if (denom2 != 0.0) {
        double factor2 = Ni_1 / denom2;   // 提取公共项
        N[2]  = (t - u_i) * factor2;
        dN[2] = 2.0 * factor2;
    } else {
        N[2]  = 0.0;
        dN[2] = 0.0;
    }

    // 依赖于前面的 denom1 和 denom2
    double term1_N = 0.0, term2_N = 0.0;
    double term1_dN = 0.0, term2_dN = 0.0;

    if (denom1 != 0.0) {
        double factor1 = Nim1_1 / denom1;
        term1_N  = (t - u_im1) * factor1;
        term1_dN = 2.0 * factor1;
    }
    
    if (denom2 != 0.0) {
        double factor2 = Ni_1 / denom2;
        term2_N  = (u_ip2 - t) * factor2;
        term2_dN = -2.0 * factor2;
    }
    N[1]  = term1_N + term2_N;
    dN[1] = term1_dN + term2_dN;
}
void BSplineBasis::basisP3(int span, double t, double N[4]) {
    // ---- 预读取节点 ----
    const double u_im2 = knots_[span - 2];
    const double u_im1 = knots_[span - 1];
    const double u_i   = knots_[span];
    const double u_ip1 = knots_[span + 1];
    const double u_ip2 = knots_[span + 2];
    const double u_ip3 = knots_[span + 3];

    const double width = u_ip1 - u_i;
    if (width == 0.0) {
        N[0] = N[1] = N[2] = N[3] = 0.0;
        return;
    }

    // ---- 1 阶基函数 ----
    const double N1_i   = (t - u_i) / width;
    const double N1_im1 = (u_ip1 - t) / width;

    // ---- 2 阶基函数 ----
    double N_deg2[3];
    const double denom1 = u_ip1 - u_im1;
    const double denom2 = u_ip2 - u_i;

    if (denom1 != 0.0) {
        N_deg2[0] = (u_ip1 - t) * N1_im1 / denom1;
        N_deg2[1] = (t - u_im1) * N1_im1 / denom1;
    } else {
        N_deg2[0] = N_deg2[1] = 0.0;
    }

    if (denom2 != 0.0) {
        N_deg2[1] += (u_ip2 - t) * N1_i / denom2;
        N_deg2[2]  = (t - u_i) * N1_i / denom2;
    } else {
        N_deg2[2] = 0.0;
    }

    // ---- 3 阶基函数 ----
    const double denomB = u_ip1 - u_im2;
    const double denomC = u_ip2 - u_im1;
    const double denomD = u_ip3 - u_i;

    // N_{span-3,3}
    if (denomB != 0.0) {
        double factor = N_deg2[0] / denomB;
        N[0] = (u_ip1 - t) * factor;
    } else {
        N[0] = 0.0;
    }

    // N_{span-2,3}
    double term1 = 0.0, term2 = 0.0;
    if (denomB != 0.0) term1 = (t - u_im2) * N_deg2[0] / denomB;
    if (denomC != 0.0) term2 = (u_ip2 - t) * N_deg2[1] / denomC;
    N[1] = term1 + term2;

    // N_{span-1,3}
    term1 = 0.0; term2 = 0.0;
    if (denomC != 0.0) term1 = (t - u_im1) * N_deg2[1] / denomC;
    if (denomD != 0.0) term2 = (u_ip3 - t) * N_deg2[2] / denomD;
    N[2] = term1 + term2;

    // N_{span,3}
    if (denomD != 0.0) {
        N[3] = (t - u_i) * N_deg2[2] / denomD;
    } else {
        N[3] = 0.0;
    }
}

void BSplineBasis::basisP3Der(int span, double t, double dN[4]) {
    // ---- 预读取节点 ----
    const double u_im2 = knots_[span - 2];
    const double u_im1 = knots_[span - 1];
    const double u_i   = knots_[span];
    const double u_ip1 = knots_[span + 1];
    const double u_ip2 = knots_[span + 2];
    const double u_ip3 = knots_[span + 3];

    const double width = u_ip1 - u_i;
    if (width == 0.0) {
        dN[0] = dN[1] = dN[2] = dN[3] = 0.0;
        return;
    }

    // ---- 1 阶基函数 ----
    const double N1_i   = (t - u_i) / width;
    const double N1_im1 = (u_ip1 - t) / width;

    // ---- 2 阶基函数 (N'_{3} 依赖于 N_{2} 的值) ----
    double N_deg2[3];
    const double denom1 = u_ip1 - u_im1;
    const double denom2 = u_ip2 - u_i;

    if (denom1 != 0.0) {
        N_deg2[0] = (u_ip1 - t) * N1_im1 / denom1;
        N_deg2[1] = (t - u_im1) * N1_im1 / denom1;
    } else {
        N_deg2[0] = N_deg2[1] = 0.0;
    }

    if (denom2 != 0.0) {
        N_deg2[1] += (u_ip2 - t) * N1_i / denom2;
        N_deg2[2]  = (t - u_i) * N1_i / denom2;
    } else {
        N_deg2[2] = 0.0;
    }

    // ---- 3 阶基函数导数 ----
    // N'_{i,p} = p * ( N_{i,p-1}/(u_{i+p}-u_i) - N_{i+1,p-1}/(u_{i+p+1}-u_{i+1}) )
    const double denomB = u_ip1 - u_im2; 
    const double denomC = u_ip2 - u_im1; 
    const double denomD = u_ip3 - u_i;   

    // dN[0] (即 N'_{span-3, 3})
    if (denomB != 0.0) {
        dN[0] = -3.0 * N_deg2[0] / denomB;
    } else {
        dN[0] = 0.0;
    }

    // dN[1] (即 N'_{span-2, 3})
    dN[1] = 0.0;
    if (denomB != 0.0){
        dN[1] += 3.0 * N_deg2[0] / denomB;
    } 
    if (denomC != 0.0){
        dN[1] -= 3.0 * N_deg2[1] / denomC;
    } 

    // dN[2] ( N'_{span-1, 3})
    dN[2] = 0.0;
    if (denomC != 0.0) {
        dN[2] += 3.0 * N_deg2[1] / denomC;
    }
    if (denomD != 0.0){
        dN[2] -= 3.0 * N_deg2[2] / denomD;
    } 

    // dN[3] ( N'_{span, 3})
    if (denomD != 0.0) {
        dN[3] = 3.0 * N_deg2[2] / denomD;
    } else {
        dN[3] = 0.0;
    }
}

void BSplineBasis::basisP3WithDer(int span, double t, double N[4], double dN[4]) {
    // ---- 预读取节点 ----
    const double u_im2 = knots_[span - 2];
    const double u_im1 = knots_[span - 1];
    const double u_i   = knots_[span];
    const double u_ip1 = knots_[span + 1];
    const double u_ip2 = knots_[span + 2];
    const double u_ip3 = knots_[span + 3];

    const double width = u_ip1 - u_i;
    if (width == 0.0) {
        for (int i = 0; i < 4; ++i) {
            N[i] = dN[i] = 0.0;
        }
        return;
    }

    // ---- 1 阶基函数 (Linear) ----
    const double Ni_1   = (t - u_i) / width;
    const double Nim1_1 = (u_ip1 - t) / width;

    // ---- 2 阶基函数 (Quadratic) ----
    double N_deg2[3];
    const double denom1 = u_ip1 - u_im1;
    const double denom2 = u_ip2 - u_i;

    if (denom1 != 0.0) {
        N_deg2[0] = (u_ip1 - t) * Nim1_1 / denom1;
        N_deg2[1] = (t - u_im1) * Nim1_1 / denom1;
    } else {
        N_deg2[0] = N_deg2[1] = 0.0;
    }

    if (denom2 != 0.0) {
        N_deg2[1] += (u_ip2 - t) * Ni_1 / denom2;
        N_deg2[2]  = (t - u_i) * Ni_1 / denom2;
    } else {
        N_deg2[2] = 0.0;
    }

    // ---- 3 阶基函数与导数 (Cubic) ----
    const double denomB = u_ip1 - u_im2;
    const double denomC = u_ip2 - u_im1;
    const double denomD = u_ip3 - u_i;

    // 计算 N[0] 和 dN[0] (N_{span-3,3})
    if (denomB != 0.0) {
        double factorB = N_deg2[0] / denomB;
        N[0]  = (u_ip1 - t) * factorB;
        dN[0] = -3.0 * factorB;
    } else {
        N[0] = dN[0] = 0.0;
    }

    // 计算 N[3] 和 dN[3] (N_{span,3})
    if (denomD != 0.0) {
        double factorD = N_deg2[2] / denomD;
        N[3]  = (t - u_i) * factorD;
        dN[3] = 3.0 * factorD;
    } else {
        N[3] = dN[3] = 0.0;
    }

    // 计算 N[1], dN[1] 和 N[2], dN[2]
    double term1_N = 0.0, term2_N = 0.0;
    double term1_dN = 0.0, term2_dN = 0.0;
    double term3_N = 0.0, term4_N = 0.0;
    double term3_dN = 0.0, term4_dN = 0.0;

    if (denomB != 0.0) {
        double factorB = N_deg2[0] / denomB;
        term1_N  = (t - u_im2) * factorB;
        term1_dN = 3.0 * factorB;
    }

    if (denomC != 0.0) {
        double factorC = N_deg2[1] / denomC;
        // 贡献给 N[1]/dN[1] 的右半部分
        term2_N  = (u_ip2 - t) * factorC;
        term2_dN = -3.0 * factorC;
        // 贡献给 N[2]/dN[2] 的左半部分
        term3_N  = (t - u_im1) * factorC;
        term3_dN = 3.0 * factorC;
    }

    if (denomD != 0.0) {
        double factorD = N_deg2[2] / denomD;
        term4_N  = (u_ip3 - t) * factorD;
        term4_dN = -3.0 * factorD;
    }

    N[1]  = term1_N + term2_N;
    dN[1] = term1_dN + term2_dN;

    N[2]  = term3_N + term4_N;
    dN[2] = term3_dN + term4_dN;
}

void BSplineBasis::evaluateBasis(int span, double t, double* N)
{
    if (degree_ == 2){
        basisP2(span, t, N);
    }
    else if (degree_ == 3){
        basisP3(span, t, N);
    }
}

void BSplineBasis::evaluateDerBasis(int span, double t, double* N)
{
    if (degree_ == 2){
        basisP2Der(span, t, N);
    }
    else if (degree_ == 3){
        basisP3Der(span, t, N);
    }
}

void BSplineBasis::evaluateAllBasis(int span, double t, double* N, double* dN) {
    if (degree_ == 2) {
        basisP2WithDer(span, t, N, dN);
    } else if (degree_ == 3) {
        basisP3WithDer(span, t, N, dN);
    }
}

double BSplineBasis::getValue(int j, double t){
    int span = findSpan(t);
    int p = degree_;
    // 非零范围检查
    if (j < span - p || j > span){
        return 0.0;
    }
    double N[p + 1];   
    evaluateBasis(span, t, N);
    int localIndex = j - (span - p);
    return N[localIndex];
}

double BSplineBasis::getDerValue(int j, double t)
{
    int span = findSpan(t);
    int p = degree_;
    // 非零范围检查
    if (j < span - p || j > span){
        return 0.0;
    }
    double N[p + 1];   
    evaluateDerBasis(span, t, N);
    int localIndex = j - (span - p);
    return N[localIndex];
}