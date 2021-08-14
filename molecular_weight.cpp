#include"molecular_weight.h"
#include"chain_update.h"

// This function updates the number average and weight average molecular weight after each kMC step.
// Monomers are assumed to not affect the molecular weight, but dimers, trimers, and all larger oligomers
// are included in the calculation
void molecular_weight(double& Mn, double& Mw, chain_pool& all_chains, chain_pool& loops, std::vector<int>& monomerA, std::vector<int>& monomerB, std::vector<double>& monomermassA, std::vector<double>& monomermassB, bool& isloop, bool& isnewchain, bool& ismonomerA,bool& ismonomerB,double& sumNi, double& sumMiNi, double& sumMi2Ni, double& Mi_A, double& Mi_B){
   // if a chain forms a loop, there is no change in the molecular weight distribution
   if (isloop == false) {
       // if a new chain is formed from two monomers, then:
        if (isnewchain == true) {
            sumNi++; // increase the total number of chains by 1
            sumMiNi+=Mi_A+Mi_B; // increase the total molecular weight of all chains (sumMiNi) by the masses of the two monomers forming the new chain (Mi_A+Mi_B)
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B); // (new chain weight)^2. (Note: Only one new chain with mass Mi is added, so Ni is 1)
        }
        // if an A monomer is added to the end of an existing chain, then:
        else if (ismonomerA==true){
            // sumNi does not change, because the total number of chains remains constant
            sumMiNi+=Mi_A; // add weight of additional A monomer to the total molecular weight of all chains
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_B*Mi_B; // newchainweight^2-oldchainweight^2
        }
        // if a B monomer is added to the end of an existing chain, then:
        else if (ismonomerB==true){
            // sumNi does not change
            sumMiNi+=Mi_B; // add weight of additional monomer
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_A*Mi_A; // newchainweight^2-oldchainweight^2
        } 
        // if two chains react to form one longer chain, then:
        else {
            sumNi--; // total number of chains decreases by 1
            // sumMiNi does not change
            sumMi2Ni+=(Mi_A+Mi_B)*(Mi_A+Mi_B)-Mi_A*Mi_A-Mi_B*Mi_B; // newchainweight^2-oldchainAweight^2-oldchainBweight^2
        }
    }
    // Reset all the flags each time the function is called
    isloop = false; //
    isnewchain=false;
    ismonomerA=false;
    ismonomerB=false;
    Mn=sumMiNi/sumNi; // calculate the new Mn
    Mw=sumMi2Ni/sumMiNi; // calculate the new Mw
}