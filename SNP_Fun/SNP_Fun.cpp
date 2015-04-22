#include <iostream>
#include <fstream>
#include "SNP_Fun.h"
using namespace std;

SNP_Fun::SNP_Fun()
{
    //ctor
}

SNP_Fun::~SNP_Fun()
{
    //dtor
}

int SNP_Fun::hashRSID(string id, int hashSize)
{
    int sum = 0;

    for(int i = 0; i < id.length(); i++)
    {
        sum = sum + id[i];
    }

    sum = sum % hashSize;

    return sum;
}

int SNP_Fun::hashGeno(std::string geno)
{
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

void SNP_Fun::add_RSID(int geno, int chromo, int id, std::string effect)
{

    if(GC_Table[geno][chromo].empty())
    {
        vector<rsid> rsid_Table (10);
        GC_Table[geno][chromo] = rsid_Table;

        rsid *new_rsid = new rsid;

        GC_Table[geno][chromo][id] = *new_rsid;
    }
    else
    {
        cout << "fatality" << '\n';
    }
}

void SNP_Fun::createMatrix(string filename)
{
    string token;
    ifstream infile;
    infile.open(filename.c_str());


    int counter = 0;
    int column = 0;

        while (!infile.eof())
    {

        string rsId;
        string Chrom;
        string pos;
        string genotype;
        int Chromosome = 0;
        int position = 0;

        getline(infile, rsId, ' ');
        getline(infile, Chrom, ' ');
        getline(infile, pos, ' ');
        getline(infile, genotype);
        Chromosome = stoi(Chrom);

        //sort_Data
    }

}

 void SNP_Fun::sort_Data(std::string, int, int, std::string)
 {

 }
