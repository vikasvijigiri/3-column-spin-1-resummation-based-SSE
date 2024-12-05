#include <iostream>
#include <string>
#include "Bar.h"

using namespace std;




int main()
{
    Foo theFoo;
    Bar theBar;
    theBar.PrintFoo();

    cout << theFoo.MyData << endl; // correct! accessing private members of foo from main because main is friend of Foo.

    return 0;
}
