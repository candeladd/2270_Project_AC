# 2270_Project_AC
Project for Data Structures

PROJECT SUMMARY:

We will be working on creating a program that allows a person to enter in their 23andMe data and the program will analyze the data.  It will give feed back on any predicatble traits or chareteristics that can be deteremined fom their genotype. We will be creating a structure for individual rsid's/SNP's, which ultimately stores the info, and storing them in hash table.  The class will contain double array, 23 X 10, that correlates to a chromosome (1 - 23) and a genotype (e.g. AA, GG, TT, CC). Each position within the double array will then contain a hash table of the rsid's to be efficiently and easily retreviable.

HOW TO RUN:

1) Clone the project file

Optional: Download raw data files (.txt) from https://opensnp.org/user_index (it tells you whether users have SNPs available) into the project file for extra fun or just use our test.txt file.

2) Open project file in Code::Blocks IDE (SNP_Fun.cbp)

3) Recompile (The Gear)

4) Run (The Green Arrow)


DEPENDENCIES:

The program depends on being compiled using C++11 or higher a standard with no 3rd party or special libraries.

The program also requires that the user enters a file name in the console window when prompted.  The file name should be formatted as follows filename.txt where file name would not actually be filename but the name of the file you want to run.  We have included a sample file named test.txt.  Users are also able to download files from opensnp.org and save them in the 2270_Project_AC directory with the .txt extension.  Then analyze these files using our program. 

SYSTEM REQUIREMENTS:

The program is meant to be ran in Linux with Code::Blocks IDE with no additional requirements.

GROUP MEMBERS:

Andrew Candelaresi, Austin Holler

CONTRIBUTORS:

Andrew Candelaresi, Austin Holler

OPEN AND KNOWN ISSUES/BUGS:

-If user inputs unexpected input (e.g. string when it expects an int or vise-versa, entering a chromosome that doesn't exist, etc.) the program doesn't like it.

-Instead of cout-ing RSIDs when using option #4 or #5 (i.e. printAllForChromosome or printMatchingGeno) it would be better to write out the info to a file.

-Also writing out to a file the info given to users, while keeping the cout-ing aspect, option #3 (i.e. retrieveData) would be nice.

-An option, and actual way of doing it, under option #1 that lets user decided to overwrite the user database with a new user .txt file if desired.

-Any additional error checking.

-A makefile for easy execution would be cool!

-Any other functionality that might be awesome...
