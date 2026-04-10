#pragma once
#include "BSplineCurve2D.h"

class BSplineMultiPatch{
private:
    BSplineCurve2D** patches_; // 指针数组，存储每一片的地址
    int numPatches_;           // 当前已添加的片数
    int totalDofs_;            // 总的控制点数，首尾相邻
public:
    BSplineMultiPatch(BSplineCurve2D** patches, int n) {
        numPatches_ = n;
        totalDofs_ = 0;
        patches_ = new BSplineCurve2D*[numPatches_];    // 分配指针数组空间
        for (int i = 0; i < numPatches_; ++i) {
            patches_[i] = patches[i];
            totalDofs_ +=patches_[i]->getnumber();
        }
        totalDofs_= totalDofs_-numPatches_;
    }
    ~BSplineMultiPatch() {
        for (int i = 0; i < numPatches_; ++i) {
            delete patches_[i];
        }
        delete[] patches_;
    }
    
    int getTotalDofs() const { return totalDofs_; }
    int getNumPatches() const { return numPatches_;}

    BSplineCurve2D* getPatch(int index) const {return patches_[index];}
};