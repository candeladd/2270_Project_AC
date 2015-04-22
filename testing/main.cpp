#include <iostream>

using namespace std;

int main()
{
    string geno = "GT";

    int genoInt = geno[0] + geno[1];

    switch(genoInt)
    {
    case 130:
        return 0;
    case 132:
        return 1;
    case 136:
        return 2;
    case 149:
        return 3;
    case 134:
        return 4;
    case 138:
        return 5;
    case 151:
        return 6;
    case 142:
        return 7;
    case 155:

        return 8;
    case 168:
        return 9;
    default:
        cout << "Invalid Genotype" << '\n';
    }

    return -1;
}
