#ifndef SNP_FUN_H
#define SNP_FUN_H

#define GENE 10
#define CHROM 23

#include <vector>
#include <iostream>


struct rsid
{
    std::string id = "empty";
    std::string info = "empty";
    std::string genotype = "empty";
    int chromosome = NULL;
    int position = NULL;
    rsid *next=NULL;
};

class SNP_Fun
{
    public:
        SNP_Fun();
        ~SNP_Fun();
        void createMatrix(std::string filename);
        rsid* sort_Data(std::string, int, int, std::string);
        int hashGeno(std::string geno);
        void add_RSID(int geno, int chromo, int id, std::string effect);
        int hashRSID(std::string id, int hashSize);
        void place_in(rsid* temp);
        rsid* findRSID(std::string name);
    protected:
    private:
        std::vector<rsid*> GC_Table [GENE][CHROM]= {};

};

#endif // SNP_FUN_H
