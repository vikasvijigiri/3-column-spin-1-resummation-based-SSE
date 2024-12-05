#include "Foo.h"



class Foo; // forward declaration
class Bar; // forward declaration



class Bar // friend of Foo but Foo is not my Friend so it can't access mine
{
    Foo fooObj;

    public:
        Bar(){}
        void PrintFoo()const;
};
