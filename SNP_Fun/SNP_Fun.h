#ifndef SNP_FUN_H
#define SNP_FUN_H

#define GENE 10
#define CHROM 23

#include <vector>
#include <iostream>
#include <fstream>


struct RSID
{
    std::string id = "empty";
    std::string info = "empty";
    std::string genotype = "empty";
    int chromosome = NULL;
    int position = NULL;
    RSID *next=NULL;
};

class SNP_Fun
{
    public:
        SNP_Fun();
        ~SNP_Fun();
        std::ofstream outputFile;
        void createMatrix(std::string filename);
        RSID* sort_Data(std::string RSID_str, int chromosome1, int position1, std::string genotype1);
        int getGenoInt(std::string geno);
        void default_add_rsid(RSID* rsid_obj, std::string effect);
        void user_add_rsid(RSID* rsid_obj);
        void initial_data();
        int hashRSID(std::string id, int hashSize);
        void printMatchingGeno(std::string geno);
        void retrieveRSID(std::string RSID_str, int chromo, std::string geno);
        void printAllForChromosome(int chromo);
        void compareData(std::string filename);
        void retrieveData(std::string RSID_str, int chromo, int genoInt);
    protected:
    private:
        std::vector<RSID*> GC_Table [GENE][CHROM]= {};
};

#endif // SNP_FUN_H
