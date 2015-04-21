#include <iostream>
#include "SNP_Fun.h"
#include <fstream>
using namespace std;

int main()
{
        string userinput;
    SNP_Fun *Analyze= new SNP_Fun();

    while (userinput != "5")// main menu loop
    {
        cout	<<	"======Main Menu====="	<<	endl;
        cout	<<	"1. Enter your SNP file name"	<<	endl;
        cout	<<	"2. Get results"	<<	endl;
         getline(cin, userinput);
        if (userinput == "1")
        {
            string infile;
            cout	<<	"Enter your SNP file name"	<<	endl;
            getline(cin, infile);
            Analyze->createVector(infile)



        }

    }
    #include <iostream>
#include "SNP_Fun.h"

using namespace std;

    cout << "Hello World!" << endl;
    return 0;
}
