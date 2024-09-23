#include<iostream>
#include "P3P_HD.h"
// #include "utils/matrix.h"
#include "demo_P3P.h"
#include <Eigen/Core>
using namespace std;

template<class T>
void print_result(string method_name,int num_sols,std::vector<Eigen::Matrix<T, 3, 3>> &Rs, std::vector<Eigen::Vector3<T>> &Ts) {
    cout<<"Method: "<<method_name<<endl;
    cout<<"The num of solutions is :"<<num_sols<<endl;
    cout<<"The matrix of the P3P problems is: "<<endl;
    for(int i=0;i<num_sols;i++){
        cout<<"Solution "<<i+1<<endl;
        cout<<"R matrix: "<<endl;
        for(int t=0;t<3;t++){
            for(int j=0;j<3;j++){
                cout<<Rs[i](t*3+j)<<" ";
            }
            cout<<endl;
        }
        cout << "t matrix: " << endl;
        for (int t = 0; t < 3; t++) {
            cout << Ts[i](t) << " ";
        }
        cout << endl<<endl;
    }
}

int demo_P3P(hd::P3P_methods methods_id){

    double imagePoints[12]={650.0694006368067 ,1381.360697742374 ,1147.187644348962 ,1949.537909630508 ,
                                1187.956439904574 ,1429.31115805356 ,84.32218903568389 ,1279.093532946494  ,
                                1 ,1 ,1 ,1};

    double worldPoints[12]={7.990260726989884,0.6251698576368507,10.76873510411773,8.405694430677894 ,
                            2.444672074608359,-0.7061331846727531,-6.853648093269328,3.335931846624678,
                            3.433383852087895   ,0.02056500410048895   ,-0.7582608758809428,-7.583658175429444 };

    double Intrin[3] = {1000., 1000., 1000.};

    std::vector<Eigen::Matrix<double, 3, 3>> Rs(4);
    std::vector<Eigen::Vector3<double>> Ts(4);
    int num_sols=0;

    switch(methods_id) {
        case hd::P3P_methods::hd_mul4: {
            num_sols= hd::P3P_HD(imagePoints,worldPoints,Intrin,Rs,Ts,5,hd::P3P_methods::hd_mul4);
            print_result("Mul4",num_sols,Rs,Ts);
            break;
        }
        case hd::P3P_methods::hd_mul3: {
            num_sols= hd::P3P_HD(imagePoints,worldPoints,Intrin,Rs,Ts,5,hd::P3P_methods::hd_mul3);
            print_result("Mul3",num_sols,Rs,Ts);
            break;
        }
        case hd::P3P_methods::hd_uni3_1v: {
            num_sols= hd::P3P_HD(imagePoints,worldPoints,Intrin,Rs,Ts,5,hd::P3P_methods::hd_uni3_1v);
            print_result("Uni_3_1V",num_sols,Rs,Ts);
            break;
        }
        case hd::P3P_methods::hd_uni3_3v: {
            num_sols= hd::P3P_HD(imagePoints,worldPoints,Intrin,Rs,Ts,5,hd::P3P_methods::hd_uni3_3v);
            print_result("Mul_3_3V",num_sols,Rs,Ts);
            break;
        }
    }

    return 0;
}



