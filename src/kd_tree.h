#ifndef KD_TREE_H_
#define KD_TREE_H_

#include <vector>
#include <stdio.h>
#include "helper_functions.h"

struct BiNode{
    LandmarkObs* pLandmark;
    BiNode* left;
    BiNode* right;
};
typedef BiNode* BiTree;

#endif