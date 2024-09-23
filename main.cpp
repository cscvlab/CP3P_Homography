#include "demo_Ransac.h"
#include "demo_P3P.h"
#include "P3P_HD.h"
using namespace std;
int main() {

    bool isEPNP=true;
    bool isCeres=true;

    demo_Ransac(hd::P3P_methods::hd_mul4, isEPNP, isCeres);
    demo_Ransac(hd::P3P_methods::hd_mul3, isEPNP, isCeres);
    demo_Ransac(hd::P3P_methods::hd_uni3_1v, isEPNP, isCeres);
    demo_Ransac(hd::P3P_methods::hd_uni3_3v, isEPNP, isCeres);

    demo_P3P(hd::P3P_methods::hd_mul3);
    demo_P3P(hd::P3P_methods::hd_mul4);
    demo_P3P(hd::P3P_methods::hd_uni3_1v);
    demo_P3P(hd::P3P_methods::hd_uni3_3v);

}
