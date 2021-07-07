/*class a
{
    public:
    int type;
    a(int t)
    :type(t){}
};

class b: public a
{
    public:
    int bb;
    b(int t, int bbb)
    :a(t),bb(bbb){}
};*/

#include <stdio.h>
#include <vector>
using namespace std;

int main()
{
    vector<int> a;
    a.push_back(1);
    printf("%llu",a.size());
    return 0;
}