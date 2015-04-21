#include "SNP_Fun.h"

SNP_Fun::SNP_Fun()
{
    //ctor
}

SNP_Fun::~SNP_Fun()
{
    //dtor
}

void SNP_Fun::createMatrix(string filename)
{
    string token;
    ifstream infile;
    infile.open(filename.c_str());


    int counter = 0;
    int column = 0;

        while (!myfile.eof())
    {

        string rsId;
        string Chrom;
        string pos;
        string genotype
        int Chromosome = 0;
        int position = 0;

        getline(myfile, rsId, ' ');
        getline(myfile, Chrom, ' ');
        getline(myfile, pos, '   ');
        getline(myfile, genotype);
        Chromosome = stoi()

        sort_Data
    }

}

 void SNP_Fun::sort_Data(std::string, int, int, std::string)
 {

 }
