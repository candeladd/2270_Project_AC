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
        string genotype_in;
        int Chromosome_in = 0;
        int position_in = 0;

        getline(infile, rsId, ' ');
        getline(infile, Chrom, ' ');
        getline(infile, pos, ' ');
        getline(infile, genotype_in);
        Chromosome_in = stoi(Chrom);
        position_in = stoi(pos);

        rsid *temp = sort_Data(rsId, chromosome_in, position_in, genotype_in);

    }

}

 rsid* SNP_Fun::sort_Data(std::string str, int chromosome1, int position1, std::string genotype1)
 {
        rsid *rsID = new rsid;// creates a new instance of the struct
        rsID->id = str; //sets the id variable
        rsID->chromosome = chromosome1;//sets the chromosome variable for this instance
        rsID->position = position1;// sets the position variable for this instance
        rsID->genotype = genotype1;// sets the genotype variable for this instance

 }

void SNP_Fun::place_in(rsid* temp)
{
    int geno_sum = hashGeno(std::string temp->genotype)
    int
    for (int i=0;i<CHROM;i++)
    {
        if (i == temp->chromosome)
        {
            if (GC_Table[geno_sum][i] == NULL)
            {
                GC_Table[geno_sum][i] = temp;
                GC_Table[geno_sum][i]->next = NULL;
            }
            else
            {
                rsid* temp2 = GC_Table[geno_sum][i];
                if (temp2->id.compare(temp->id)>0)// if A>B A(B) you will get value greater that zero. B is what you want to insert
                {
                    temp->next = temp2;
                    GC_Table[geno_sum][i] = temp;
                    return;
                }
                while (temp2->next != NULL)
                {

                    if (temp2->next->id.compare(temp->id)>0)// if A<B A(B) alue less than 0
                    {
                       break;
                    }
                    temp2=temp2->next;

                }
                temp->next=temp2->next;
                temp2->next=temp;

            }

        }
    }
    }
