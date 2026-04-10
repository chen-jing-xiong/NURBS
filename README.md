

## ✨ 主要特性

- 支持 **B-Spline** 和 **NURBS** 二维曲线
- 支持 **2阶（二次）** 和 **3阶（三次）** 曲线
- 高效的基函数计算（包含组合计算 `evaluateAllBasis`）
- 节点向量支持 **均匀** 和 **非均匀** 两种模式
- 提供单点求值、导数求值、位置+导数同时求值、批量求值接口
- 实现了曲线雅可比（Jacobian）计算
- 支持多段曲线拼接（`BSplineMultiPatch` / `NURBSMultiPatch`）
- 内存管理清晰，析构函数自动释放资源

NURBS-BSpline-Curve/
├── BSplineBasis.h
├── BSplineBasis.cpp
├── BSplineCurve2D.h
├── BSplineCurve2D.cpp
├── NURBSBasis.h
├── NURBSBasis.cpp
├── NURBSCurve2D.h
├── NURBSCurve2D.cpp
├── BSplineMultiPatch.h
├── NURBSMultiPatch.h
└── README.md

### 1. 创建 B-Spline 曲线

```cpp
// 准备节点向量、控制点、权重等...
int degree = 3;
double knots[] = {0,0,0,0,1,2,3,4,4,4,4};  // 示例节点
double px[] = {...};   // x 坐标
double py[] = {...};   // y 坐标

BSplineBasis* basis = new BSplineBasis(degree, knots, sizeof(knots)/sizeof(double), false);
BSplineCurve2D curve(basis, numPoints, px, py);

// 求值
double x, y, dx, dy;
curve.evaluateAll(0.5, x, y, dx, dy);
//创建 NURBS 曲线
double weights[] = {1.0, 1.0, 1.0, 1.0, ...};

BSplineBasis* bsBasis = new BSplineBasis(...);
NURBSBasis* nurbsBasis = new NURBSBasis(bsBasis, weights, numWeights);

NURBSCurve2D nurbsCurve(nurbsBasis, numPoints, px, py);

