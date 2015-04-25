#include <iostream>
#include "SNP_Fun.h"
#include <fstream>
using namespace std;

int main()
{
    SNP_Fun *Initial = new SNP_Fun();

    Initial->initial_data();

    string userinput;

    SNP_Fun *User= new SNP_Fun();

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
            User->createMatrix(infile);
        }

        else if(userinput == "2")
        {
            string RSID_str;
            string chromo;
            string geno;

            cout << "Enter your SNP: ";
            getline(cin, RSID_str);

            cout << "Enter your chromosome: ";
            getline(cin, chromo);

            cout << "Enter your genotype: ";
            getline(cin, geno);
            cout << '\n';

            Initial->retrieveRSID(RSID_str, stoi(chromo), geno);
        }

    }

    return 0;
}
