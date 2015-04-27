# 2270_Project_AC
Project for Data Structures

PROJECT SUMMARY:

We will be working on creating a program that allows a person to enter in their 23andMe data and the program will analyze the data.  It will give feed back on any predicatble traits or chareteristics that can be deteremined fom their genotype. We will be creating a structure for individual rsid's/SNP's, which ultimately stores the info, and storing them in hash table.  The class will contain double array, 23 X 10, that correlates to a chromosome (1 - 23) and a genotype (e.g. AA, GG, TT, CC). Each position within the double array will then contain a hash table of the rsid's to be efficiently and easily retreviable.

HOW TO RUN:



DEPENDENCIES:

The program depends on being compiled using C++11 or higher a standard with no 3rd party or special libraries.

SYSTEM REQUIREMENTS:

The program is meant to be ran in Linux with Code::Blocks IDE with no additional requirements.

GROUP MEMBERS:

Andrew Candelaresi, Austin Holler

CONTRIBUTORS:

Andrew Candelaresi, Austin Holler

OPEN ISSUES/BUGS:

-If user inputs unexpected input (e.g. string when it expects an int or vise-versa) the program doesn't like it.

-Instead of cout-ing RSIDs when using option #4 or #5 (i.e. printAllForChromosome or printMatchingGeno) it would be better to write out the info to a file.
