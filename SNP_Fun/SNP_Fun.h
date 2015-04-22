#ifndef SNP_FUN_H
#define SNP_FUN_H

#define GENE 10
#define CHROM 23

#include <vector>
#include <iostream>


struct rsid
{
    std::string info = "empty";
    std::string genotype = "empty";
    int chromosome = NULL;
};

class SNP_Fun
{
    public:
        SNP_Fun();
        ~SNP_Fun();
        void createMatrix(std::string filename);
        void sort_Data(std::string, int, int, std::string);
        int hashGeno(std::string geno);
        void add_RSID(int geno, int chromo, int id, std::string effect);
        int hashRSID(std::string id, int hashSize);
    protected:
    private:
        std::vector<rsid> GC_Table [GENE][CHROM];

};

#endif // SNP_FUN_H
