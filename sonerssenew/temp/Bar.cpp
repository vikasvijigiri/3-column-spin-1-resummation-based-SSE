#include "Bar.h"
#include "Foo.h"
#include <iostream>
using namespace std;
void Bar::PrintFoo()const
{
    cout << fooObj.MyData << endl; // accessing private member data of class Foo from Bar
}

