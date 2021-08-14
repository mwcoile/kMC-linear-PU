#pragma once
#include<iomanip>
#include<stdlib.h>
#include<stdexcept>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<algorithm>

struct chain
{
    std::vector<std::vector<int>> v;
    int frontend; // identity of species at start of vector
    int backend; // identity of species at end of vector
    double chain_mass;
};
// Need to keep track of all the polymer chains that have been created
typedef std::vector<chain> chain_pool;
void explicit_sequence_record(int whichA, int whichB, std::vector<int>& monomerA, std::vector<int>& monomerB, int& A_monomer_type, int& B_monomer_type, std::vector<int>& chainsA, std::vector<int>& chainsB, chain_pool& all_chains, chain_pool& loops, bool& isloop, bool& isnewchain, bool& ismonomerA,bool& ismonomerB, double& Mi_A, double& Mi_B, std::vector<double>& monomermassA, std::vector<double>& monomermassB);   

