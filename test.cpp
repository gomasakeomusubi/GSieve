#include <vector>
#include <iostream>

using namespace std;

struct TES{
    int a;
    float b;
};

void printTES(vector<TES*> a){
    for(size_t i = 0; i < a.size(); i++){
        cout << "(" << a[i]->a << ", " << a[i]->b << ") " << flush;
    }
    cout << endl;
}

void func(TES *tes){
    tes->b += 1;
}

int main(){
    vector<TES*> tes, tes2;

    // TES *c = new TES;
    for(int i = 0; i < 5; i++){
        TES *c = new TES;
        c->a = i;
        c->b = i + 0.5;
        tes.push_back(c);
    }

    cout << "original  : " << flush;
    printTES(tes);

    // tesに格納されているTES*型の変数は削除されるがvectorの要素は変わらない
    // ポインタを変数に代入すると、値をコピーしているがポインタの指すアドレスは同じなので、片方をdeleteするともう片方もdeleteされる。

    // eraseではvectorからは削除されるがポインタ自身は残る
    
    // for(int i = 0; i < tes.size(); i++){
    //     if(1 <= tes[i]->a && tes[i]->a <= 2){
    //         tes.erase(tes.begin()+i);
    //         i--;
    //     }
    // }

    TES *f = new TES;
    f = tes[2];
    for(int i = 0; i < tes.size(); i++){
        if(i == 2) func(tes[i]);
    }

    cout << "after     : " << flush;
    printTES(tes);

    cout << "(" << f->a << ", " << f->b << ")" << endl;

    for(auto &t : tes) delete t;

    return 0;
}