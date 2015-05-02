#include <iostream>
#include "SNP_Fun.h"
#include <fstream>
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
        cout    <<  "======Main Menu====="  <<  endl;
        cout    <<  "1. Create User Database"   <<  endl;
        cout    <<  "2. Search known SNP Database"  <<  endl;
        cout    <<  "3. Get all known User effects" << endl;
        cout    <<  "4. Print all SNPs matching genotype"   <<  endl;
        cout    <<  "5. Print all SNPs matching chromosome" << endl;
        cout    <<  "6. Quit" << endl;
        
        
        getline(cin, userinput);
        
        if (userinput == "1" && !user_created)
        {
            string infile;
            cout    <<  "Enter your SNP file name"  <<  endl;
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
            int chromo;
            string geno;
            
            cout << "Enter your SNP: ";
            getline(cin, RSID_str);
            if(cin.fail())
            {
                cout<<"input not valid"<<endl;
            }
            cout << "Enter your chromosome: ";
            cin >> chromo;
            if(cin.fail())
            {
                cout<<"input not valid"<<endl;
                break;
            }
            cout << "Enter your genotype: ";
            getline(cin, geno);
            cout << '\n';
            
            if(RSID_str == "" || geno == "")
            {
                cout << "Please enter a valid SNP, chromosome, and genotype combination" << '\n';
            }
            if (cin.fail()) {
                cout<<"Invalid input"<<endl;
            }
            else
            {
                    Initial->retrieveRSID(RSID_str, chromo, geno);
            }
        }
        else if (userinput == "3")
        {
            string infile;
            cout    <<  "Enter your SNP file name"  <<  endl;
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
            cout << "Enter your chromosome: ";
            getline(cin, chromo1);
            int chromo = stoi(chromo1);
            
            User->printAllForChromosome(chromo);
            
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
