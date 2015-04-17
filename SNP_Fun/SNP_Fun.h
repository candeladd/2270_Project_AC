#ifndef SNP_FUN_H
#define SNP_FUN_H
#include <vector>
#include <iostream>

const int gene = 10;
const int chrom = 23;

struct rsid
{
    std::string info;
};

class SNP_Fun
{
    public:
        SNP_Fun();
        ~SNP_Fun();
    protected:
    private:
        std::vector<rsid> foo [gene][chrom];

};

#endif // SNP_FUN_H
