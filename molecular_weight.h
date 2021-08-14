#pragma once
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include"chain_update.h"

void molecular_weight(double& Mn, double& Mw, chain_pool& all_chains, chain_pool& loops, std::vector<int>& monomerA, std::vector<int>& monomerB, std::vector<double>& monomermassA, std::vector<double>& monomermassB, bool& isloop, bool& isnewchain, bool& ismonomerA,bool& ismonomerB,double& sumNi, double& sumMiNi, double& sumMi2Ni, double& Mi_A, double& Mi_B);