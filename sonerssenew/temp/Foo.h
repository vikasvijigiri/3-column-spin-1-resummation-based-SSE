#pragma once
class Foo
{
    public:
        Foo() : MyData(100){}

    private:
        int MyData;
        friend class Bar; // any access scope is ok. Bar will have full access to Foo
        friend int main(); // main now can access private members
};
