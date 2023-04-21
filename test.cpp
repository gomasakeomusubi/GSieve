#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>

NTL_CLIENT

int main(){
    ostringstream oss;
    ifstream ifs;
    mat_ZZ B;

    oss.str("");
    oss << "inpute.txt";
    ifs.open(oss.str(), ios::in);
    if(ifs.is_open()){
        ifs >> B;
        ifs.close();
    }
    else cin >> B;

    cout << B << endl;

    return 0;
}