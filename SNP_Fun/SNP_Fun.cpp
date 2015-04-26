#include <iostream>
#include <fstream>
#include <sstream>
#include "SNP_Fun.h"

SNP_Fun::SNP_Fun()
{
    //ctor
}

SNP_Fun::~SNP_Fun()
{
    //dtor
}
/*hashRSID - A hash function for returning the
* hash key of an entered SNP rsid*/
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

/*getGenoInt - Returns a key, via a switch statement,
* for the inputed genotype*/
int SNP_Fun::getGenoInt(std::string geno)
{
    int genoInt = geno[0] + geno[1]; //Get the value of the two geno chars

    switch(genoInt)
    {
    case 130: //Genotype: AA
        return 0;
    case 132: //Genotype: AC
        return 1;
    case 136: //Genotype: AG
        return 2;
    case 149: //Genotype: AT
        return 3;
    case 134: //Genotype: CC
        return 4;
    case 138: //Genotype: CG
        return 5;
    case 151: //Genotype: CT
        return 6;
    case 142: //Genotype: GG
        return 7;
    case 155: //Genotype: GT
        return 8;
    case 168: //Genotype: TT
        return 9;
    case 90: // add to handle genotype "--"
        return 10;
    case 146: //added to handle genotype "II"
        return 11;
    case 141: //added to handle genotype "DI"
        return 12;
    default: //In case the genotype doesn't match any of the above
        std::cout << "Invalid Genotype" << '\n';
        return -1;
    }
}

/*default_add_rsid - A function for creating the initial/default database
* of SNP's and info for comparison or general use.*/
void SNP_Fun::default_add_rsid(RSID* rsid_obj, std::string effect)
{
    int genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string
    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string
    rsid_obj->info = effect; //Sets the effect of the rsid

    if(GC_Table[genoInt][rsid_obj->chromosome].empty()) //If a vector doesn't  exist in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome] = rsid_Table; //Add new vector

        GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj; //Assign the rsid object to hash table
    }
    else
    {
        //std::cout << genoInt << ":" << rsid_obj->chromosome << ":" << idInt << '\n';

        if(GC_Table[genoInt][rsid_obj->chromosome][idInt] == NULL) //If there is no collision add the rsid object
        {
            GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
        }
        else //Else there is a collision in hash table
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome][idInt]; //Temp for iterating

            //If temps rsid substr ints is greater than the one to be added, move temp
            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp; //Temps rsid obj gets assigned as the addition's next

                GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj; //The one to be added takes temp's place
            }
            else //Else the one to be added is less than temp
            {

                //While we haven't reached the end of the linked list and the one to be added is still greater than temp, iterate through
                while(temp->next != NULL && (stoi(rsid_obj->id.substr(2, rsid_obj->id.length())) > stoi(temp->next->id.substr(2, temp->id.length()))))
                {
                    temp = temp->next; //Iterating
                }

                //Out of the while, if we're not at the end of the linked list insert the rsid obj in between nodes
                if(temp->next != NULL)
                {

                    rsid_obj->next = temp->next; //The addition's next is the current pos next
                    temp->next = rsid_obj; //The current pos next becomes the addition
                }

                else //Else we're at the end of the linked list
                {
                    temp->next = rsid_obj; //The addition becomes the new end
                }
            }
        }
    }
}

/*user_add_rsid - User based database creation function that adds rsid
* objs that have been created from an inputted txt file.*/
void SNP_Fun::user_add_rsid(RSID* rsid_obj)
{
    int   genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string

    if (genoInt==10||genoInt==11 ||genoInt==12)// ignores rsid with genotype "--" "II" and "DI"
    {
        return;
    }

    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string

    if(GC_Table[genoInt][rsid_obj->chromosome].empty()) //If a vector doesn't exist in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome] = rsid_Table; //Add new vector

        GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj; //Add rsid obj to hash table
    }
    else
    {
        //std::cout << genoInt << ":" << rsid_obj->chromosome << ":" << idInt << '\n';

        if(GC_Table[genoInt][rsid_obj->chromosome][idInt] == NULL) //If there is no collision add the rsid object
        {
            GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj;
        }
        else //Else there is a collision in hash table
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome][idInt]; //Temp for iterating

            //If temps rsid substr ints is greater than the one to be added, move temp
            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp; //Temps rsid obj gets assigned as the addition's next

                GC_Table[genoInt][rsid_obj->chromosome][idInt] = rsid_obj; //The addition takes temp's place
            }
            else //Else the one to be added is less than temp
            {

                //While we haven't reached the end of the linked list and the one to be added is still greater than temp, iterate through
                while(temp->next != NULL && (stoi(rsid_obj->id.substr(2, rsid_obj->id.length())) > stoi(temp->next->id.substr(2, temp->id.length()))))
                {
                    temp = temp->next; //Iterating
                }

                //Out of the while, if we're not at the end of the linked list insert the rsid obj in between nodes
                if(temp->next != NULL)
                {

                    rsid_obj->next = temp->next; //The addition's next is the current pos next
                    temp->next = rsid_obj; //The current pos next becomes the addition
                }

                else //Else we're at the end of the linked list
                {
                    temp->next = rsid_obj; //The addition becomes the new end
                }
            }
        }
    }
}

/*initial_data - Where all of the default/initial rsid objs are
* created and stored into the default database - 31 SNPs -  */
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

    temp = sort_Data("rs1061170", 1, 196690107, "TT");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs1061170", 1, 196690107, "CT");
    default_add_rsid(temp, "2.5x risk for age related Muscular Dystrophy; higher mortality among 90 - 99yr olds");

    temp = sort_Data("rs1061170", 1, 196690107, "CT");
    default_add_rsid(temp, "5.9x risk for age related Muscular Dystrophy; higher mortality among 90 - 99yr olds");

    temp = sort_Data("rs1042522", 17, 7676154, "CC");
    default_add_rsid(temp, "Typically prolongs life by 3 years; chemotherapy is more effective");

    temp = sort_Data("rs1042522", 17, 7676154, "CG");
    default_add_rsid(temp, "Typically slightly increased lifespan");

    temp = sort_Data("rs1042522", 17, 7676154, "GG");
    default_add_rsid(temp, "Typically slightly shorter lifespan");

    temp = sort_Data("rs7574865", 2, 191099907, "TT");
    default_add_rsid(temp, "2.6x risk of Rheumatoid Arthritis; 3.10x risk of Lupus; 2.28x risk of Sjögren's syndrome; increased risk of Type-1 Diabetes; and increased risk of Primary Biliary Cirrhosis.");

    temp = sort_Data("rs7574865", 2, 191099907, "GT");
    default_add_rsid(temp, "1.3x risk of Rheumatoid Arthritis; 1.55x risk of Lupus; 1.42x risk of Sjögren's syndrome; increased risk of Type-1 Diabetes; and increased risk of Primary Biliary Cirrhosis.");

    temp = sort_Data("rs7574865", 2, 191099907, "TT");
    default_add_rsid(temp, "Normal risk of Rheumatoid Arthritis, Lupus, Sjögren's syndrome, Type-1 Diabetes, and Primary Biliary Cirrhosis.");

    temp = sort_Data("rs780094", 2, 27518370, "GG");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs780094", 2, 27518370, "AG");
    default_add_rsid(temp, "Metabolic consequences");

    temp = sort_Data("rs780094", 2, 27518370, "AA");
    default_add_rsid(temp, "Metabolic consequences");

    temp = sort_Data("rs738409", 22, 43928847, "GG");
    default_add_rsid(temp, "Most common genotype; slightly less damage from alcohol");

    temp = sort_Data("rs738409", 22, 43928847, "CG");
    default_add_rsid(temp, "Increased liver fat, increased odds of alcoholic liver disease");

    temp = sort_Data("rs738409", 22, 43928847, "CC");
    default_add_rsid(temp, "Higher odds of alcoholic liver disease, increased liver fat");

    temp = sort_Data("rs16969968", 15, 78590583, "AA");
    default_add_rsid(temp, "Higher risk for nicotine dependence, lower risk for cocaine dependence");

    temp = sort_Data("rs16969968", 15, 78590583, "AG");
    default_add_rsid(temp, "Slightly higher risk for nicotine dependence, lower risk for cocaine dependence");

    temp = sort_Data("rs16969968", 15, 78590583, "GG");
    default_add_rsid(temp, "Normal");

    temp = sort_Data("rs25531", 17, 30237328, "AA");
    default_add_rsid(temp, "Short form of 5-httlpr. lower levels of serotonin, slightly less happy, benefits from more support");

    temp = sort_Data("rs25531", 17, 30237328, "GG");
    default_add_rsid(temp, "Long form of 5-httlpr. less sensitive to pain");

    temp = sort_Data("rs11209026", 1, 67240275, "GG");
    default_add_rsid(temp, "Normal, but higher risk for certain autoimmune diseases");

    temp = sort_Data("rs11209026", 1, 67240275, "AG");
    default_add_rsid(temp, "0.26x lower risk for certain autoimmune diseases");

    temp = sort_Data("rs11209026", 1, 67240275, "AA");
    default_add_rsid(temp, "0.26x lower risk for certain autoimmune diseases");
}

void SNP_Fun::createMatrix(std::string filename)
{
    std::cout << "Processing File" << '\n';

    std::string token;
    std::ifstream infile;
    std::string pos;
    std::string genotype_in;
    std::string line;
    int Chromosome_in = 0;
    int position_in = 0;

    infile.open(filename.c_str());

    int counter = 0;
    int column = 0;

    while(getline(infile, line))
    {
        std::stringstream ss;

        if(line[0] != '#')
        {
            ss << line;

            std::string rsId;
            std::string Chrom;

            getline(ss, rsId, '\t');
            getline(ss, Chrom, '\t');

            if(Chrom != "X" && Chrom != "Y" && Chrom != "MT")
            {
                std::cout << Chrom << '\n';

                getline(ss, pos, '\t');
                getline(ss, genotype_in);

                Chromosome_in = std::stoi(Chrom);
                position_in = std::stoi(pos);

                RSID *temp = sort_Data(rsId, Chromosome_in, position_in, genotype_in);
                user_add_rsid(temp);
            }
        }
        //getline(infile, rsId,  '\t');
    }

}

/*sort_Data - Creates rsid objs for easy use in other functions*/
 RSID* SNP_Fun::sort_Data(std::string rsid_str, int chromosome1, int position1, std::string genotype1)
 {
        RSID *rsID = new RSID;// creates a new instance of the struct
        rsID->id = rsid_str; //sets the id variable
        rsID->chromosome = chromosome1;//sets the chromosome variable for this instance
        rsID->position = position1;// sets the position variable for this instance
        rsID->genotype = genotype1;// sets the genotype variable for this instance

        return rsID; //Returns the new obj

 }

/*retrieveRSID - Given user inputted information, this function checks the default database
* by genotype, chromosome, and rsid hash table position for an rsid obj and info to return to the user.*/
 void SNP_Fun::retrieveRSID(std::string RSID_str, int chromo, std::string geno)
 {
     int genoInt = getGenoInt(geno); //Gets the integer value of the genotype string
     int idInt = hashRSID(RSID_str, 10); //Gets the hashed value of the RSID string

     //If there is nothing at the genotype/chromosome pos, RSID isn't stored
     if(GC_Table[genoInt][chromo].empty())
     {
         std::cout << "RSID does not exist within database" << '\n';
         return;
     }

     //Else there exists a hash table at the genotype/chromosome pos
     else
     {
        //If there is nothing at the hash table key, RSID isn't stored
        if(GC_Table[genoInt][chromo][idInt] == NULL)
        {
            std::cout << "RSID does not exist within database" << '\n';
            return;
        }
        //Else if the hash table rsid obj is a match, return info
        else if(GC_Table[genoInt][chromo][idInt]->id == RSID_str)
        {
            std::cout << "RSID: " << GC_Table[genoInt][chromo][idInt]->id << '\n';
            std::cout << "Chromosome: " << GC_Table[genoInt][chromo][idInt]->chromosome << '\n';
            std::cout << "Position: " << GC_Table[genoInt][chromo][idInt]->position << '\n';
            std::cout << "Genotype: " << GC_Table[genoInt][chromo][idInt]->genotype << '\n';
            std::cout << "Effect: " << GC_Table[genoInt][chromo][idInt]->info << '\n' << '\n';
            return;
        }
        //Else if the hash table rsid isn't a match, check the linked list
        else
        {
            RSID *temp = GC_Table[genoInt][chromo][idInt]; //Temp for iterating

            //While there are still nodes to check, iterate
            while(temp->next != NULL)
            {
                temp = temp->next; //Iterating

                //If the current node is a match, return info
                if(temp->id == RSID_str)
                {
                    std::cout << "RSID: " << temp->id << '\n';
                    std::cout << "Chromosome: " << temp->chromosome << '\n';
                    std::cout << "Postition: " << temp->position << '\n';
                    std::cout << "Genotype: " << temp->genotype << '\n';
                    std::cout << "Effect: " << temp->info << '\n' << '\n';
                    return;
                }
            }
        }

        //Else the rsid doesn't exist in the database
        std::cout << "RSID does not exist within database" << '\n';
     }

 }

//finds all rsid associated with a particular genotype in a users DNA
void SNP_Fun::printMatchingGeno(std::string geno)
{
    int genoInt = getGenoInt(geno); //hash the genotype sent into the function
    for(int i=0; i < CHROM;i++) // loop through all rows for the particular genotype column you are wanting to print out
    {
        if(!GC_Table[genoInt][i].empty()) // check to see that a vector exists at each of these in each row for this column
        {
            for (int j=0; j<10;j++)// loop through each index in the vector at this row comlumn position
            {

                if (GC_Table[genoInt][i][j]!= NULL)  //if there is a value inside of the rsid_Table (vector) at each index print out that value and loop through
                {
                    RSID *temp = GC_Table[genoInt][i][j];

                    while (temp != NULL) // the linked list until there are no more values.
                    {

                        std::cout << "RSID: " << temp->id << " "
                         << "Chromosome: " << temp->chromosome <<" "<<  "Genotype: " << temp->genotype << '\n';

                        temp = temp->next;

                    }

                }
            }
        }
        else
        {
            std::cout << "no RSID with " <<geno<< " at Chromosome "<< i << '\n';
        }
    }
}

void SNP_Fun::printAllForChromosome(int chromo)
{

    for(int i=0; i < GENE;i++) // go down an entire column for one chromosome
    {
        if(!GC_Table[i][chromo].empty()) // check to see that a vector exists at each index of row for this column
        {
            for (int j=0; j<10;j++)// loop through each index in the vector at this row comlumn position
            {

                if (GC_Table[i][chromo][j]!= NULL)  //if there is a value inside of the rsid_Table (vector) at each index print out that value and loop through
                {
                    RSID *temp = GC_Table[i][chromo][j];

                    while (temp != NULL) // the linked list until there are no more values.
                    {

                        std::cout << "RSID: " << temp->id << " "
                         << "Chromosome: " << temp->chromosome <<" "<<  "Genotype: " << temp->genotype << '\n';

                        temp = temp->next;

                    }

                }
            }
        }
        else
        {
            std::cout << "no RSID at chromosome " <<chromo<< " with genotype "<< i << '\n';
        }
    }

}

