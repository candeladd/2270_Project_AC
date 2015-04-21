#ifndef SNP_FUN_H
#define SNP_FUN_H

#define GENE 10
#define CHROM 23

#include <vector>
#include <iostream>

//const int gene = 10;
//const int chrom = 23;

struct rsid
{
    std::string info;
};

class SNP_Fun
{
    public:
        SNP_Fun();
        ~SNP_Fun();
        void SNP_Fun::createMatrix(std::string filename);
        void sort_Data(std::string, int, int, std::string);
    protected:
    private:
        std::vector<rsid> foo [GENE][CHROM];

};

#endif // SNP_FUN_H
