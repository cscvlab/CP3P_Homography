#include "demo_Ransac.h"
#include "demo_P3P.h"
#include "P3P_HD.h"
using namespace std;
int main() {

    //"isEPNP=true && isCeres=true" : LO-RANSAC framework which encloses the local and global optimization.
    //"isEPNP=false && isCeres=false" :  Vanilla RANSAC framework (No optimization is applied).
    bool isEPNP=false;
    bool isCeres=false;


    //Several sets of data from the Cambridge KingsCollege scene.
    string filsPath="../Cambridge/KingsCollege_1";

    demo_Ransac(hd::P3P_methods::hd_mul4, isEPNP, isCeres, filsPath);
    demo_Ransac(hd::P3P_methods::hd_mul3, isEPNP, isCeres, filsPath);
    demo_Ransac(hd::P3P_methods::hd_uni3_1v, isEPNP, isCeres, filsPath);
    demo_Ransac(hd::P3P_methods::hd_uni3_3v, isEPNP, isCeres, filsPath);

    demo_P3P(hd::P3P_methods::hd_mul3);
    demo_P3P(hd::P3P_methods::hd_mul4);
    demo_P3P(hd::P3P_methods::hd_uni3_1v);
    demo_P3P(hd::P3P_methods::hd_uni3_3v);

}
