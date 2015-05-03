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

/*
Function Prototype:
int hashRSID(std::string, int);

Function Description:
hashRSID - A hash function for returning the
hash key of an entered SNP rsid.

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->hashRSID("rs564798", 10);

Within class - hashRSID("rs564798", 10);

Pre-Condition:
The input rsid is a string, then converted via hash function, and is
then stored into the int sum as the hash key.

Post-Condition:
Rsid has now been converted, via the hash function, and stored into
sum which has been returned for use in positioning rsid obj in
vector hash tables that have been stored in the 2D array.*/

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

/*
Function Prototype:
int getGenoInt(std::string);

Function Description:
getGenoInt - Returns a key, via a switch statement,
for the inputed genotype

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->getGenoInt("AA");

Inside class - getGenoInt("AA");

Pre-Condition:
Takes in a genotype string, two characters, and returns a value 0 - 12
depending on its ASCII value.

Post-Condition:
Returned value 0 -12, depending on the genotype ASCII value, for use
in determining a RSID obj position in the genotype/chromosome 2D array.
*/
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

/*
Function Prototype:
void default_add_rsid(RSID*, std::string)

Function Description:
default_add_rsid - A function for creating the initial/default database
of SNP's and info for comparison of user data or general use.

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->default_add_rsid(temp, "Implies a higher risk of cancer");

Inside class - default_add_rsid(temp, "Implies a higher risk of cancer");

Pre-Condition:
Default database chromosome/genotype positions are empty or have a pre-existing vector
hash table. Vector hash table is created in appropriate position if one doesn't
exist. Appropriate hash table position is filled with an rsid pointer obj. The string
effect is the effect the specific rsid gene mutation means for a person.

Post-Condition:
There is now a vector hash table in position genotype/chromosome that has a position
filled with a rsid pointer obj that has been determined by the hashRSID function.
*/
void SNP_Fun::default_add_rsid(RSID* rsid_obj, std::string effect)
{
    int genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string
    if(genoInt < 0 || genoInt > 9)
    {
        return;
    }
    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string
    rsid_obj->info = effect; //Sets the effect of the rsid

    if(GC_Table[genoInt][rsid_obj->chromosome - 1].empty()) //If a vector doesn't  exist in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome - 1] = rsid_Table; //Add new vector

        GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj; //Assign the rsid object to hash table
    }
    else
    {
        if(GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] == NULL) //If there is no collision add the rsid object
        {
            GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj;
        }
        else //Else there is a collision in hash table
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome - 1][idInt]; //Temp for iterating

            //If temps rsid substr ints is greater than the one to be added, move temp
            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp; //Temps rsid obj gets assigned as the addition's next

                GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj; //The one to be added takes temp's place
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

/*
Function Prototype:
void user_add_rsid(RSID*);

Function Description:
user_add_rsid - User database creation function that adds rsid
objs that have been sorted and created from an user inputted .txt file.

Example:
Outside of class - SNP_Fun *User;
                   User->user_add_rsid(temp);

Inside class - user_add_rsid(temp);

Pre-Condition:
User database Chromosome/genotype positions are empty or have a pre-existing vector
hash table. Vector hash table is created in appropriate position if one doesn't
exist. Appropriate hash table position is filled with an rsid pointer obj.

Post-Condition:
There is now a vector hash table in the genotype/chromosome position that has a position
filled with a rsid pointer obj that has been determined by the hashRSID function.
*/
void SNP_Fun::user_add_rsid(RSID* rsid_obj)
{
    int genoInt = getGenoInt(rsid_obj->genotype); //Gets the integer value of the genotype string

    if (genoInt < 0 || genoInt > 9)// ignores rsid with genotype "--" "II" and "DI"
    {
        return;
    }

    int idInt = hashRSID(rsid_obj->id, 10); //Gets the hashed value of the RSID string

    if(GC_Table[genoInt][rsid_obj->chromosome - 1].empty()) //If a vector doesn't exist in the gene/chrom position create one and add RSID
    {
        std::vector<RSID*> rsid_Table (10);
        GC_Table[genoInt][rsid_obj->chromosome - 1] = rsid_Table; //Add new vector

        GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj; //Add rsid obj to hash table
    }
    else
    {
        //std::cout << genoInt << ":" << rsid_obj->chromosome - 1 << ":" << idInt << '\n';

        if(GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] == NULL) //If there is no collision add the rsid object
        {
            GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj;
        }
        else //Else there is a collision in hash table
        {

            RSID *temp = GC_Table[genoInt][rsid_obj->chromosome - 1][idInt]; //Temp for iterating

            //If temps rsid substr ints is greater than the one to be added, move temp
            if(stoi(temp->id.substr(2, temp->id.length())) > stoi(rsid_obj->id.substr(2, rsid_obj->id.length())))
            {
                rsid_obj->next = temp; //Temps rsid obj gets assigned as the addition's next

                GC_Table[genoInt][rsid_obj->chromosome - 1][idInt] = rsid_obj; //The addition takes temp's place
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

/*
Function Prototype:
void initial_data();

Function Description:
initial_data - Where all of the default/initial rsid objs are
created and stored into the default database - 31 SNPs - 85 Gene Mutations

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->initial_data();

Inside class - initial_data();

Precondition:
Data is in it's "raw" or unhashed condition and then passed to the sort_Data function
to create a rsid obj of it. This object is then taken and passed to default_add_rsid
along with its mutation effect and is stored to the default database.

Post-Condition:
The default database is now full of new rsid objects that have been assigned their
relevant information and place in their relevant positions in the hash table. This
is where their stored for use in comparison of user data or general use.
*/
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

/*
Function Prototype:
void createMatrix(std::string);

Function Description:
createMatrix reads in from a file provided by the user. Parses file based on tab
spacing. Finds matching genotypes. Puts all RSID objects into a 3d matrix.

Example:
Outside of class - SNP_Fun *User;
                   User->createMatrix("test.txt");

Inside class - createMatrix("test.txt");

Pre-Condition:
Unused locations in User 3D matrix hashtable are "empty". User file is unparsed.
The input string is the filename to use. The file information is parsed per
rsid and stored with relevant info into rsid objs. Which are then stored into
3D matrix hash table for later use or analysis.

Post-Condition:
User file has been parsed an stored into relevant genotype/chromosome positions
inside of their relevant hashed position into the rsid vector hashtable. The
user database has been created.
*/
void SNP_Fun::createMatrix(std::string filename)
{

    std::string token;
    std::ifstream infile;
    std::string pos;
    std::string genotype_in;
    std::string line;
    int Chromosome_in = 0;
    int position_in = 0;
    int counter = 0;

    infile.open(filename.c_str());

    if(!infile.is_open())
    {
        std::cout << "Enter a valid filename" << '\n';
    }

    else
    {
        std::cout << "Processing File" << '\n';
        while(getline(infile, line))// reads in the file one line at a time and stores the line in a string called line
        {
            std::stringstream ss; // creates a string stream
            counter++;

            if(line[0] != '#')// checks to see that the first charecter in the line is not a #
            {
                ss << line; // sends data in line to the string stream ss

                std::string rsId;
                std::string Chrom;


                getline(ss, rsId, '\t'); // grabs all chars before the first tab sends them to rsId
                getline(ss, Chrom, '\t');//grabs all chars after first tab & before third and sends it to Chrom

                if(Chrom != "X" && Chrom != "Y" && Chrom != "MT") // checks to make sure the Chromosome is not X, Y, or mitochondria
                {
                    getline(ss, pos, '\t');//store the position of the rsid
                    getline(ss, genotype_in);// stores the genotype

                    Chromosome_in = std::stoi(Chrom); // turns the chromosome string into an int
                    position_in = std::stoi(pos);   // turns the positon string into an int

                    RSID *temp = sort_Data(rsId, Chromosome_in, position_in, genotype_in); // creates an instance of the struct rsid and add all user info to it
                    user_add_rsid(temp); // places the new rsid instance into a 3d matrix

                }
            }

            if(counter%100000 == 0)
            {
                std::cout << "..." << '\n';
            }

        }
    }

    infile.close();

}

/*
Function Prototype:
void compareData(std::string);

Function Description:
compareData reads in from a file provided by the user. Parses file based on tab
spacing. Finds matching genotypes. Then sends parsed data to retrieveData for
comparison against the default database.

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->compareData("test.txt");

Inside class - compareData("test.txt");

Pre-Condition:
Defalut database is full of rsid objs ready for comparison and use. User file is unparsed.
The input string is the filename to use. The file information is parsed per
rsid and each is compare against the default database by calling retrieveData.
Found matches are then outputted to console with relevant information displayed.

Post-Condition:
User file has been parsed and compared against the default database using the
retrieveData function. Matches have been outputted to console with relevant info
The state of the databases are the same.
*/
void SNP_Fun::compareData(std::string filename)
{
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

    if(!infile.is_open())
    {
        std::cout << "Enter a valid filename" << '\n';
    }

    else
    {
        std::cout << "Processing File" << '\n';
        while(getline(infile, line))// reads in the file one line at a time and stores the line in a string called line
        {
            std::stringstream ss; // creates a string stream

            if(line[0] != '#')// checks to see that the first charecter in the line is not a #
            {
                ss << line; // sends data in line to the string stream ss

                std::string rsId;
                std::string Chrom;


                getline(ss, rsId, '\t'); // grabs all chars before the first tab sends them to rsId
                getline(ss, Chrom, '\t');//grabs all chars after first tab & before third and sends it to Chrom

                if(Chrom != "X" && Chrom != "Y" && Chrom != "MT") // checks to make sure the Chromosome is not X, Y, or mitochondria
                {

                    getline(ss, pos, '\t');//store the position of the rsid
                    getline(ss, genotype_in);// stores the genotype

                    Chromosome_in = std::stoi(Chrom); // turns the chromosome string into an int
                    position_in = std::stoi(pos);   // turns the positon string into an int



                    int genoInt = getGenoInt(genotype_in); //find the genotype int for the current RSID being added to the matrix

                    if (genoInt != 11 && genoInt != 12 && genoInt != 10) // makes sure that genotype is not II,DI,or --
                    {
                        retrieveData(rsId, Chromosome_in, genoInt); // sends parsed data to retrieveData
                    }

                }
            }

        }
    }

    infile.close();

}



/*
Function Prototype:
RSID* sort_Data(std::string, int, int, std::string)

Function Description:
sort_Data - Creates rsid objs for easy use in other functions

Example:
Outside of class - SNP_Fun *Initial;
                   Initial->sort_Data("rs673546", 12, 76927632, "AA");

Inside class - sort_Data("rs673546", 12, 76927632, "AA");

Pre-Condition:
There doesn't exist a rsid obj. The data is taken either from entry or parsed info.
A rsid obj is created, assigned relevant info, and returned for use in other functions.

Post-Condition:
A rsid object has been created, assigned relevant info, and returned for use in other functions.
*/
 RSID* SNP_Fun::sort_Data(std::string rsid_str, int chromosome1, int position1, std::string genotype1)
 {
        RSID *rsID = new RSID;// creates a new instance of the struct
        rsID->id = rsid_str; //sets the id variable
        rsID->chromosome = chromosome1;//sets the chromosome variable for this instance
        rsID->position = position1;// sets the position variable for this instance
        rsID->genotype = genotype1;// sets the genotype variable for this instance

        return rsID; //Returns the new obj

 }
 /*
Function Prototype:
void retrieveRSID(std::string RSID_str, int chromo, std::string geno);

Function Description:

retrieveRSID - Given user inputted information, this function checks the default database
 by genotype, chromosome, and rsid hash table position for an rsid obj and info to return to the user.


Example:
outside of class-   SNP_Fun *Initial;
                    Initial->retrieveRSID(RSID_str, stoi(chromo), geno);

Pre-Condition:
The program user has a particular snp that they are interested in finding information about.  They know the rsid, the chromosome, and the genotype.  They will
select option 2 on the menu and enter the rsid, the chromosome, and the genotype

Post-Condition:
a print out of the rsid, chromosome, genotype, and effect of that particular genotype will appear in the console window.

*/
 void SNP_Fun::retrieveRSID(std::string RSID_str, int chromo, std::string geno)
 {
    int genoInt = getGenoInt(geno); //Gets the integer value of the genotype string
    if(genoInt < 0 || genoInt > 9)
    {
        std::cout << "Please enter a valid genotype" << '\n';
        return;
    }
    int idInt = hashRSID(RSID_str, 10); //Gets the hashed value of the RSID string

     //If there is nothing at the genotype/chromosome pos, RSID isn't stored
    if(GC_Table[genoInt][chromo - 1].empty())
    {
        std::cout << "RSID does not exist within database" << '\n';
        return;
    }

     //Else there exists a hash table at the genotype/chromosome pos
    else
    {
        //If there is nothing at the hash table key, RSID isn't stored
        if(GC_Table[genoInt][chromo - 1][idInt] == NULL)
        {
            std::cout << "RSID does not exist within database" << '\n';
            return;
        }
        //Else if the hash table rsid obj is a match, return info
        else if(GC_Table[genoInt][chromo - 1][idInt]->id == RSID_str)
        {
            std::cout << "RSID: " << GC_Table[genoInt][chromo - 1][idInt]->id << '\n';
            std::cout << "Chromosome: " << GC_Table[genoInt][chromo - 1][idInt]->chromosome << '\n';
            std::cout << "Position: " << GC_Table[genoInt][chromo - 1][idInt]->position << '\n';
            std::cout << "Genotype: " << GC_Table[genoInt][chromo - 1][idInt]->genotype << '\n';
            std::cout << "Effect: " << GC_Table[genoInt][chromo - 1][idInt]->info << '\n' << '\n';
            return;
        }
        //Else if the hash table rsid isn't a match, check the linked list
        else
        {
            RSID *temp = GC_Table[genoInt][chromo - 1][idInt]; //Temp for iterating

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

/*
Function Prototype:
void retrieveData(std::string RSID_str, int chromo, int genoInt);

Function Description:
- Given user inputted file it compares the users SNPs to the default database created in the Initial object.  it then
finds any RSID, gene, chromosome pairs with an effect and prints out the effect that the users data contains.

Example:
inside of class - retrieveData(rsId, Chromosome_in, genoInt);

Pre-Condition:
The program user has a SNP file and they want to know what information can be determined about their particular genetic information.  They send file name into this function.
the information is then analyzed.

Post-Condition:
a print out of all known effects that can be determined from their data is printed out to the console window.

*/
 void SNP_Fun::retrieveData(std::string RSID_str, int chromo, int genoInt)
 {
     int idInt = hashRSID(RSID_str, 10); //Gets the hashed value of the RSID string

     //If there is nothing at the genotype/chromosome pos, RSID isn't stored
     if(GC_Table[genoInt][chromo - 1].empty())
     {
         return;
     }

     //Else there exists a hash table at the genotype/chromosome pos
     else
     {
        //If there is nothing at the hash table key, RSID isn't stored
        if(GC_Table[genoInt][chromo - 1][idInt] == NULL)
        {
            return;
        }
        //Else if the hash table rsid obj is a match, return info
        else if(GC_Table[genoInt][chromo - 1][idInt]->id == RSID_str)
        {
            std::cout << "RSID: " << GC_Table[genoInt][chromo - 1][idInt]->id << '\n';
            std::cout << "Chromosome: " << GC_Table[genoInt][chromo - 1][idInt]->chromosome << '\n';
            std::cout << "Position: " << GC_Table[genoInt][chromo - 1][idInt]->position << '\n';
            std::cout << "Genotype: " << GC_Table[genoInt][chromo - 1][idInt]->genotype << '\n';
            std::cout << "Effect: " << GC_Table[genoInt][chromo - 1][idInt]->info << '\n' << '\n';
            return;
        }
        //Else if the hash table rsid isn't a match, check the linked list
        else
        {
            RSID *temp = GC_Table[genoInt][chromo - 1][idInt]; //Temp for iterating

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
     }

 }
/*
Function Prototype:
void printMatchingGeno(std::string geno);

Function Description:
- Once a user database is created, this function will take in a genotype from the user and print out all RSIDs and the chromosome they are indexed to for that
particular genotype.

Example:
outside of class -  SNP_Fun *User;
                    User->printAllForChromosome(chromo);

Pre-Condition:
The user has entered their file and created a 3d matrix holding their particular information.  They want to know all to the chromosomes and RSID associated
with a particular genotype.

Post-Condition:
every RSID and the chromosome where it is located are printed to the console window.
*/
void SNP_Fun::printMatchingGeno(std::string geno)
{
    int genoInt = getGenoInt(geno); //hash the genotype sent into the function
    if(genoInt < 0 || genoInt > 9)
    {
        std::cout << "Please enter a valid genotype" << '\n';
         return;
    }

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
/*
Function Prototype:
void printAllForChromosome(int chromo);

Function Description:
- Once a user database is created, this function will take in a chromosome from the user and print out all RSIDs and the genotype at that particular
chromosome.

Example:
outside of class -  SNP_Fun *User;
                    User->printAllForChromosome(int chromo);
Pre-Condition:
The user has entered their file and created a 3d matrix holding their particular information.  They want to know all to the genotypes and RSIDs found at
that chromosome.

Post-Condition:
every RSID and the genotype at that chromosome printed to the console window.
*/
void SNP_Fun::printAllForChromosome(int chromo)
{

    for(int i=0; i < GENE;i++) // go down an entire column for one chromosome
    {
        if(!GC_Table[i][chromo - 1].empty()) // check to see that a vector exists at each index of row for this column
        {
            for (int j=0; j<10;j++)// loop through each index in the vector at this row comlumn position
            {

                if (GC_Table[i][chromo - 1][j]!= NULL)  //if there is a value inside of the rsid_Table (vector) at each index print out that value and loop through
                {
                    RSID *temp = GC_Table[i][chromo - 1][j];

                    while (temp != NULL) // the linked list until there are no more values.
                    {

                        std::cout << "RSID: " << temp->id << " "
                         << "Chromosome: " << temp->chromosome <<" "<<  "Genotype: " << temp->genotype << '\n';//print all found information

                        temp = temp->next; //advance temp

                    }

                }
            }
        }
        else
        {
            std::cout << "no RSID at chromosome " <<chromo<< " with genotype "<< i << '\n'; // if there is no rsid at the chromosome
        }
    }

}
