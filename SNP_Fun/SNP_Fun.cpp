#include <iostream>
#include <fstream>
#include "SNP_Fun.h"

SNP_Fun::SNP_Fun()
{
    //ctor
}

SNP_Fun::~SNP_Fun()
{
    //dtor
}

int SNP_Fun::hashRSID(std::string id, int hashSize)
{
    int sum = 0;

    for(int i = 0; i < id.length(); i++)
    {
        sum = sum + id[i];
    }

    sum = sum % hashSize;

    return sum;
}

int SNP_Fun::getGenoInt(std::string geno)
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
        std::cout << "Invalid Genotype" << '\n';
    }

    return -1;
}

void SNP_Fun::default_add_rsid(RSID* rsid_obj, std::string effect)
{
    int genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string
    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string
    rsid_obj->info = effect; //Sets the effect of the rsid

    if(GC_Table[genoInt][rsid_obj->chromosome].empty()) //If a vector doesn't in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome] = rsid_Table;

        GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
    }
    else
    {
        //std::cout << genoInt << ":" << rsid_obj->chromosome << ":" << idInt << '\n';

        if(GC_Table[genoInt][rsid_obj->chromosome][idInt] == NULL)
        {
            GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
        }
        else
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome][idInt];

            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp;

                GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
            }
            else
            {

                while(temp->next != NULL && (stoi(rsid_obj->id.substr(2, rsid_obj->id.length())) > stoi(temp->next->id.substr(2, temp->id.length()))))
                {
                    temp = temp->next;
                }

                if(temp->next != NULL)
                {

                    rsid_obj->next = temp->next;
                    temp->next = rsid_obj;
                }

                else
                {
                    temp->next = rsid_obj;
                }
            }
        }
    }
}

void SNP_Fun::user_add_rsid(RSID* rsid_obj)
{
    int genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string
    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string

    if(GC_Table[genoInt][rsid_obj->chromosome].empty()) //If a vector doesn't in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome] = rsid_Table;

        GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
    }
    else
    {
        //std::cout << genoInt << ":" << rsid_obj->chromosome << ":" << idInt << '\n';

        if(GC_Table[genoInt][rsid_obj->chromosome][idInt] == NULL)
        {
            GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
        }
        else
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome][idInt];

            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp;

                GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
            }
            else
            {

                while(temp->next != NULL && (stoi(rsid_obj->id.substr(2, rsid_obj->id.length())) > stoi(temp->next->id.substr(2, temp->id.length()))))
                {
                    temp = temp->next;
                }

                if(temp->next != NULL)
                {

                    rsid_obj->next = temp->next;
                    temp->next = rsid_obj;
                }

                else
                {
                    temp->next = rsid_obj;
                }
            }
        }
    }
}

void SNP_Fun::initial_data()
{
    RSID *temp = sort_Data("rs53576", 3, 8762685, "AA");
    default_add_rsid(temp, "Tendency towards a lack of empathy");

    temp = sort_Data("rs53576", 3, 8762685, "AG");
    default_add_rsid(temp, "Tendency towards a lack of empathy");

    temp = sort_Data("rs53576", 3, 8762685, "AG");
    default_add_rsid(temp, "Tendency for an optimistic and empathetic attitude; handles stress well");

    temp = sort_Data("rs1815739", 11, 66560624, "CC");
    default_add_rsid(temp, "Better performing muscles; likely better sprinter when running");

    temp = sort_Data("rs1815739", 11, 66560624, "CT");
    default_add_rsid(temp, "Mixed muscle types; likely better sprinter when running");

    temp = sort_Data("rs1815739", 11, 66560624, "TT");
    default_add_rsid(temp, "Impaired muscle performance; likely endurance athlete");

    temp = sort_Data("rs7412", 19, 44908822, "CC");
    default_add_rsid(temp, "More likely to gain weight if taking olanzapine");

    temp = sort_Data("rs7412", 19, 44908822, "CT");
    default_add_rsid(temp, "More likely to gain weight if taking olanzapine");

    temp = sort_Data("rs7412", 19, 44908822, "TT");
    default_add_rsid(temp, "Normal likelihood of gaining weight if taking olanzapine");

    temp = sort_Data("rs429358", 19, 44908684, "CT");
    default_add_rsid(temp, "3x or greater increased risk for Alzheimer's");

    temp = sort_Data("rs429358", 19, 44908684, "TT");
    default_add_rsid(temp, "1.4x increased risk for heart disease");

    temp = sort_Data("rs6152", 19, 44908684, "CC");
    default_add_rsid(temp, "3x or greater increased risk for Alzheimer's");

    temp = sort_Data("rs1800497", 11, 113400106, "CC");
    default_add_rsid(temp, "Normal OCD risk; normal Tardive Diskinesia risk; lower ADHD risk; less alcohol dependence; higher risk of Postoperative Nausea; lower obesity; Bupropion is effective");

    temp = sort_Data("rs1800497", 11, 113400106, "CT");
    default_add_rsid(temp, "0.5x lower OCD risk, 0.87x lower Tardive Diskinesia risk, higher ADHD risk. More Alcohol Dependence. Lower risk of Postoperative Nausea. Increased obesity. Bupropion is not effective.");

    temp = sort_Data("rs1800497", 11, 113400106, "CT");
    default_add_rsid(temp, "0.25x lower OCD; 0.56x lower Tardive Diskinesia; higher ADHD; 1.4x Alcohol Dependence; lower Postoperative Nausea; Increased obesity; less pleasure response; Bupropion ineffective; 2.4x risk for adenoma recurrence.");

    temp = sort_Data("rs1805007", 16, 89919709, "CC");
    default_add_rsid(temp, "Normal response to anesthetics");

    temp = sort_Data("rs1805007", 16, 89919709, "CT");
    default_add_rsid(temp, "Increased response to anesthetics");

    temp = sort_Data("rs1805007", 16, 89919709, "TT");
    default_add_rsid(temp, "Increased response to anesthetics and 13-20x higher likelihood of red hair");

    temp = sort_Data("rs9939609", 16, 53786615, "AA");
    default_add_rsid(temp, "Obesity risk and 1.6x risk for Type-2 Diabetes");

    temp = sort_Data("rs9939609", 16, 53786615, "AT");
    default_add_rsid(temp, "Obesity risk and 1.3x risk for Type-2 Diabetes");

    temp = sort_Data("rs9939609", 16, 53786615, "TT");
    default_add_rsid(temp, "Lower risk of obesity and Type-2 Diabetes");

    temp = sort_Data("rs12979860", 19, 39248147, "CC");
    default_add_rsid(temp, "80% of Hepatitis C patients respond to treatment");

    temp = sort_Data("rs12979860", 19, 39248147, "CT");
    default_add_rsid(temp, "20-40% of Hepatitis C patients respond to treatment");

    temp = sort_Data("rs12979860", 19, 39248147, "TT");
    default_add_rsid(temp, "20-25% of Hepatitis C patients respond to treatment");

    temp = sort_Data("rs7903146", 10, 112998590, "TT");
    default_add_rsid(temp, "2x increased risk for Type-2 Diabetes");

    temp = sort_Data("rs7903146", 10, 112998590, "CT");
    default_add_rsid(temp, "1.4x increased risk for Type-2 Diabetes");

    temp = sort_Data("rs7903146", 10, 112998590, "CC");
    default_add_rsid(temp, "Normal/Lowered risk for Type-2 Diabetes");

    temp = sort_Data("rs4680", 22, 19963748, "AA");
    default_add_rsid(temp, "More exploratory, lower COMT enzymatic activity, therefore higher dopamine levels; lower pain threshold, enhanced vulnerability to stress, yet also more efficient at processing information under most conditions");

    temp = sort_Data("rs4680", 22, 19963748, "GG");
    default_add_rsid(temp, "Less exploratory, higher COMT enzymatic activity, therefore lower dopamine levels; higher pain threshold, better stress resiliency, albeit with a modest reduction in executive cognition performance under most conditions");

    temp = sort_Data("rs8099917", 19, 39252525, "GG");
    default_add_rsid(temp, "Lower odds of responding to peg-ifnalpha/rbv treatment");

    temp = sort_Data("rs8099917", 19, 39252525, "GT");
    default_add_rsid(temp, "Moderately lower odds of responding to peg-ifnalpha/rbv treatment");

    temp = sort_Data("rs8099917", 19, 39252525, "TT");
    default_add_rsid(temp, "Normal odds of responding to peg-ifnalpha/rbv treatment");

    temp = sort_Data("rs1800629", 6, 31575254, "AA");
    default_add_rsid(temp, "Generally higher risk for certain diseases: Allograft Rejection, Asthma, Chronic Obstructive Pulmonary Disease, Crohn's Disease, Exfoliation Glaucoma, Graves' Disease, Heart Disease, Leprosy, Liver Disease, Lymphoma, Mediterranean Spotted Fever, Migraines, Multiple Sclerosis. Nasal Polyps, Psoriasis, Rheumatoid Arthritis, Sarcoidosis, Sepsis, Systemic Lupus Erthematosus");

    temp = sort_Data("rs1800629", 6, 31575254, "AG");
    default_add_rsid(temp, "Generally higher risk for certain diseases: Allograft Rejection, Asthma, Chronic Obstructive Pulmonary Disease, Crohn's Disease, Exfoliation Glaucoma, Graves' Disease, Heart Disease, Leprosy, Liver Disease, Lymphoma, Mediterranean Spotted Fever, Migraines, Multiple Sclerosis. Nasal Polyps, Psoriasis, Rheumatoid Arthritis, Sarcoidosis, Sepsis, Systemic Lupus Erthematosus");

    temp = sort_Data("rs1800629", 6, 31575254, "GG");
    default_add_rsid(temp, "Generally normal risk for certain diseases: Allograft Rejection, Asthma, Chronic Obstructive Pulmonary Disease, Crohn's Disease, Exfoliation Glaucoma, Graves' Disease, Heart Disease, Leprosy, Liver Disease, Lymphoma, Mediterranean Spotted Fever, Migraines, Multiple Sclerosis. Nasal Polyps, Psoriasis, Rheumatoid Arthritis, Sarcoidosis, Sepsis, Systemic Lupus Erthematosus");

    temp = sort_Data("rs76265", 11, 27658369, "AA");
    default_add_rsid(temp, "Tendency towards introversion; but depression resistant");

    temp = sort_Data("rs76265", 11, 27658369, "AG");
    default_add_rsid(temp, "Impaired motor skills learning");

    temp = sort_Data("rs76265", 11, 27658369, "GG");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs6983267", 8, 127401060, "TT");
    default_add_rsid(temp, "Normal risk for cancers; aspirin reduces risk for colorectal cancer");

    temp = sort_Data("rs6983267", 8, 127401060, "GT");
    default_add_rsid(temp, "1.3x risk of prostate cancer; aspirin reduces risk for colorectal cancer");

    temp = sort_Data("rs6983267", 8, 127401060, "GT");
    default_add_rsid(temp, "1.6x risk of prostate cancer and increased risk for other cancers");

    temp = sort_Data("rs1801133", 1, 11796321, "CC");
    default_add_rsid(temp, "Common; normal homocysteine levels");

    temp = sort_Data("rs1801133", 1, 11796321, "CT");
    default_add_rsid(temp, "60% efficiency in processing folic acid, so could have elevated homocysteine and low B12 and Folate levels; Risk of poor detoxing, neurodegenerative issues, chemical sensitivities, coronary artery disease, depression, cancer");

    temp = sort_Data("rs1801133", 1, 11796321, "TT");
    default_add_rsid(temp, "60% efficiency in processing folic acid, resulting in high homocysteine and low B12 and Folate levels; Risk of poor detoxing, neurodegenerative issues, chemical sensitivities, coronary artery disease, depression, cancer");

    temp = sort_Data("rs2910164", 5, 160485411, "CC");
    default_add_rsid(temp, "Higher risk of cancer");

    temp = sort_Data("rs2910164", 5, 160485411, "CG");
    default_add_rsid(temp, "Higher/Earlier risk of cancer");

    temp = sort_Data("rs2910164", 5, 160485411, "GG");
    default_add_rsid(temp, "Normal risk of cancer");

    temp = sort_Data("rs2476601", 1, 113834946, "GG");
    default_add_rsid(temp, "Normal risk for autoimmune disorders");

    temp = sort_Data("rs2476601", 1, 113834946, "AG");
    default_add_rsid(temp, "2x Risk for Type-1 Diabetes, Rheumatoid Arthritis, and Addison's Disease");

    temp = sort_Data("rs2476601", 1, 113834946, "GG");
    default_add_rsid(temp, "2x risk of Rheumatoid Arthritis and other autoimmune disorders");

    temp = sort_Data("rs1801282", 3, 12351626, "GG");
    default_add_rsid(temp, "Avoid high fat diets");

    temp = sort_Data("rs1801282", 3, 12351626, "CG");
    default_add_rsid(temp, "Avoid high fat diets");

    temp = sort_Data("rs1801282", 3, 12351626, "CC");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs11614913", 12, 53991815, "CT");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs11614913", 12, 53991815, "TT");
    default_add_rsid(temp, "Decreased risk of colorectal and lung cancer");

    temp = sort_Data("rs11614913", 12, 53991815, "CC");
    default_add_rsid(temp, "Increased risk of various types of cancer");

    temp = sort_Data("rs13266634", 8, 117172544, "TT");
    default_add_rsid(temp, "Normal risk of Type-2 Diabetes");

    temp = sort_Data("rs13266634", 8, 117172544, "CT");
    default_add_rsid(temp, "Increased risk of Type-2 Diabetes");

    temp = sort_Data("rs13266634", 8, 117172544, "CC");
    default_add_rsid(temp, "Increased risk of Type-2 Diabetes");

    temp = sort_Data("rs17782313", 18, 60183864, "CT");
    default_add_rsid(temp, "Adults likely to be 0.22 Body Mass Index units higher");

    temp = sort_Data("rs17782313", 18, 60183864, "CC");
    default_add_rsid(temp, "Adults likely to be 0.44 Body Mass Index units higher");

    temp = sort_Data("rs17782313", 18, 60183864, "TT");
    default_add_rsid(temp, "Normal");



}

void SNP_Fun::createMatrix(std::string filename)
{
    std::string token;
    std::ifstream infile;
    infile.open(filename.c_str());


    int counter = 0;
    int column = 0;

        while (!infile.eof())
    {

        std::string rsId;
        std::string Chrom;
        std::string pos;
        std::string genotype_in;
        int Chromosome_in = 0;
        int position_in = 0;

        getline(infile, rsId, ' ');
        getline(infile, Chrom, ' ');
        getline(infile, pos, ' ');
        getline(infile, genotype_in);
        Chromosome_in = stoi(Chrom);
        position_in = stoi(pos);

        RSID *temp = sort_Data(rsId, Chromosome_in, position_in, genotype_in);

    }

}

 RSID* SNP_Fun::sort_Data(std::string rsid_str, int chromosome1, int position1, std::string genotype1)
 {
        RSID *rsID = new RSID;// creates a new instance of the struct
        rsID->id = rsid_str; //sets the id variable
        rsID->chromosome = chromosome1;//sets the chromosome variable for this instance
        rsID->position = position1;// sets the position variable for this instance
        rsID->genotype = genotype1;// sets the genotype variable for this instance

        return rsID;

 }

 void SNP_Fun::retrieveRSID(std::string RSID_str, int chromo, std::string geno)
 {
     int genoInt = getGenoInt(geno); //Gets the integer value of the genotype string
     int idInt = hashRSID(RSID_str, 10); //Gets the hashed value of the RSID string

     if(GC_Table[genoInt][chromo].empty())
     {
         std::cout << "RSID does not exist within database" << '\n';
         return;
     }
     else
     {
         if(GC_Table[genoInt][chromo][idInt]->id == RSID_str)
        {
            std::cout << "RSID: " << GC_Table[genoInt][chromo][idInt]->id << '\n';
            std::cout << "Chromosome: " << GC_Table[genoInt][chromo][idInt]->chromosome << '\n';
            std::cout << "Position: " << GC_Table[genoInt][chromo][idInt]->position << '\n';
            std::cout << "Genotype: " << GC_Table[genoInt][chromo][idInt]->genotype << '\n';
            std::cout << "Effect: " << GC_Table[genoInt][chromo][idInt]->info << '\n' << '\n';
            return;
        }
        else
        {
            RSID *temp = GC_Table[genoInt][chromo][idInt];

            while(temp->next != NULL)
            {
                temp = temp->next;

                if(temp->id == RSID_str)
                {
                    std::cout << "RSID: " << temp->id << '\n';
                    std::cout << "Chromosome: " << temp->chromosome << '\n';
                    std::cout << "Position: " << temp->position << '\n';
                    std::cout << "Genotype: " << temp->genotype << '\n';
                    std::cout << "Effect: " << temp->info << '\n' << '\n';
                    return;
                }
            }
        }

        std::cout << "RSID does not exist within database" << '\n';
     }

 }

/*void SNP_Fun::place_in(RSID* temp)
{
    int geno_sum = hashGeno(temp->genotype);
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
                RSID* temp2 = GC_Table[geno_sum][i];
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
}*/
