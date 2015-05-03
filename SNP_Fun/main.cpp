#include <iostream>
#include "SNP_Fun.h"
#include <fstream>
#include <cctype>
using namespace std;

int main()
{
    SNP_Fun *Initial = new SNP_Fun(); //Create an instance for a "default" database

    Initial->initial_data(); //Creates the "default" database

    string userinput; //For user input
    bool user_created = false; //Flag for user database creation

    SNP_Fun *User= new SNP_Fun(); //Create an instance for a "user" database

    while (userinput != "6")// main menu loop
    {
        cout	<<	"======Main Menu====="	<<	endl;
        cout	<<	"1. Create User Database"	<<	endl;
        cout	<<	"2. Search known SNP Database"	<<	endl;
        cout    <<  "3. Get all known User effects" << endl;
        cout	<<  "4. Print all User SNPs matching genotype"	<<	endl;
        cout    <<  "5. Print all User SNPs matching chromosome" << endl;
        cout    <<  "6. Quit" << endl;


        getline(cin, userinput);

        if (userinput == "1" && !user_created)
        {
            string infile;
            cout	<<	"Enter your SNP file name"	<<	endl;
            getline(cin, infile);
            User->createMatrix(infile);
            user_created = true;
        }
        else if (userinput == "1" && user_created)
        {
            cout << "User database already exists" << '\n';
        }
        else if(userinput == "2")
        {
            string RSID_str;
            string chromo;
            string geno;

            bool valid = true;
            int i = 0;

            cout << "Enter your SNP: ";
            getline(cin, RSID_str);

            cout << "Enter your chromosome: ";
            getline(cin, chromo);

            cout << "Enter your genotype: ";
            getline(cin, geno);
            cout << '\n';

            if(RSID_str == "" || chromo == "" || geno == "")
            {
                cout << "Please enter a valid SNP, chromosome, and genotype combination" << '\n';
                valid = false;
            }

            //checks to see if values for chromo are not numbers and outputs error
            while(i < chromo.length() && valid)
            {
                if(isdigit(chromo[i]) == false)
                {
                    cout << "Value entered is an invalid chromosome. Please enter valid value" << endl;

                    valid = false;
                }

                i++;
            }

            if(valid)
            {
                if(stoi(chromo) < 24)
                {
                    Initial->retrieveRSID(RSID_str, stoi(chromo), geno);
                }
                else
                {
                    cout << "Invalid chromosome value. Please choose a value from 1 - 23." << '\n';
                }
            }
        }
        else if (userinput == "3")
        {
            string infile;
            cout	<<	"Enter your SNP file name"	<<	endl;
            getline(cin, infile);
            Initial->compareData(infile);
        }
        else if (userinput == "4" && user_created)
        {
            string geno;

            cout << "Enter your genotype: ";
            getline(cin, geno);
            cout << '\n';

            User->printMatchingGeno(geno);
        }
        else if(userinput == "4" && !user_created)
        {
            cout << "Please create user database first" << '\n';
        }
        else if(userinput == "5" && user_created)
        {
            string chromo1;
            bool valid = true;
            int i = 0;

            cout << "Enter your chromosome: ";
            getline(cin, chromo1);

            while(i < chromo1.length() && valid)
            {
                if(isdigit(chromo1[i]) == false)
                {
                    cout << "Value entered is an invalid chromosome. Please enter valid value" << endl;

                    valid = false;
                }

                i++;
            }

            if(valid)
            {
                if(stoi(chromo1) < 24)
                {
                    int chromo = stoi(chromo1);
                    User->printAllForChromosome(chromo);
                }
                else
                {
                    cout << "Invalid chromosome value. Please choose a value from 1 - 23." << '\n';
                }
            }

        }
        else if(userinput == "5" && !user_created)
        {
            cout << "Please create user database first" << '\n';
        }

    }

    Initial->~SNP_Fun();
    User->~SNP_Fun();

    cout << "Goodbye!"<<endl;

    return 0;
}
