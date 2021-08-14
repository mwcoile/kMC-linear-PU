#include<iomanip>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<array>
#include<limits>
#include<numeric>
#include<algorithm>
#include<sstream>
#include"chain_update.h"
#include"molecular_weight.h"

#if __cplusplus < 201103L
#error This file requires compiler and library support for the \
ISO C++ 2011 standard. This support is currently experimental, and must be \
enabled with the -std=c++11 or -std=gnu++11 compiler options.
#endif

// author: Matthew Coile

// Function to read input data from a file, courtesy of https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str)
{
    std::vector<std::string> result;
    std::string line;
    std::getline(str,line);

    std::stringstream lineStream(line);
    std::string cell;

    while(std::getline(lineStream,cell, ','))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}

// Function to select the reaction channel to execute
void reactionchannelselector(int &A_type, int &B_type, std::vector<int> monomerA,std::vector<int> monomerB, int mu) {
    int count=0;
    for (A_type=0;A_type<monomerA.size();A_type++) {
        for (B_type=0;B_type<monomerB.size();B_type++) {
            if (count>=mu) {
                return;
            }
            count++;
        }
    }
}

// This code is designed to generate polyurethane sequences via a "one-step" reaction
int main() {

    // READ IN USER INPUTS FROM FILE //
    std::ifstream file("input.txt");
    // Read in temperature (K)
    int T = std::stoi(getNextLineAndSplitIntoTokens(file)[1]);
    // Run to a specified conversion instead of reaction time?
    std::string conversionOrReactionTime = getNextLineAndSplitIntoTokens(file)[1];
    // Eliminate any whitespace
    std::string::iterator end_pos = std::remove(conversionOrReactionTime.begin(), conversionOrReactionTime.end(), ' ');
    conversionOrReactionTime.erase(end_pos, conversionOrReactionTime.end());
    bool useConversion=false;
    if (conversionOrReactionTime=="Y"){
        useConversion=true;
    }
    // Read in reaction time (seconds)
    double simulation_time = std::stod(getNextLineAndSplitIntoTokens(file)[1]);
    // Read in conversion (as a fraction, e.g., 0.95)
    double target_conversion = std::stod(getNextLineAndSplitIntoTokens(file)[1]);
    // Read in simulation size, i.e., number of starting monomers
    int num_of_molecules_in_simulation = std::stoi(getNextLineAndSplitIntoTokens(file)[1]);
    // Read in monomer A (diol) masses (g/mol)
    std::vector<std::string> monomermassA_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> monomermassA;
    for (int i=1;i<monomermassA_input.size();i++){
        monomermassA.push_back(std::stod(monomermassA_input[i]));
    }
    // Read in monomer A (diol) concentrations (M)
    std::vector<std::string> concentrationsA_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> concentrationsA;
    for (int i=1;i<concentrationsA_input.size();i++){
        concentrationsA.push_back(std::stod(concentrationsA_input[i]));
    }
    // Read in monomer B (diisocyanate) masses (g/mol) 
    std::vector<std::string> monomermassB_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> monomermassB;
    for (int i=1;i<monomermassB_input.size();i++){
        monomermassB.push_back(std::stod(monomermassB_input[i]));
    }
    // Read in monomer B (diisocyanate) concentrations (M)
    std::vector<std::string> concentrationsB_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> concentrationsB;
    for (int i=1;i<concentrationsB_input.size();i++){
        concentrationsB.push_back(std::stod(concentrationsB_input[i]));
    }
    // The below vectors of kinetic parameters should be in the order A1+B1, A1+B2, A1+B3... A1+Bn, A2+B1,A2+B2, A2+B3... etc.
    // Read in activation energies Ea (kJ/mol)
    std::vector<std::string> Ea_kf_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> Ea_kf;
    for (int i=1;i<Ea_kf_input.size();i++){
        Ea_kf.push_back(std::stod(Ea_kf_input[i]));
    }
    // Read in preexponential factors A (L/(mol s))
    std::vector<std::string> A_kf_input = getNextLineAndSplitIntoTokens(file);
    std::vector<double> A_kf;
    for (int i=1;i<A_kf_input.size();i++){
        A_kf.push_back(std::stod(A_kf_input[i]));
    }
    // Read in filename
    std::string filename = getNextLineAndSplitIntoTokens(file)[1];
    // Eliminate any whitespace in filename
    end_pos = std::remove(filename.begin(), filename.end(), ' ');
    filename.erase(end_pos, filename.end());
    filename.append(".txt");
    
    const int M = monomermassA.size()*monomermassB.size(); // number of reaction channels considered
    
    // for volume calculation
    const double total_A_concentration = std::accumulate(concentrationsA.begin(), concentrationsA.end(), decltype(concentrationsA)::value_type(0));
    const double total_B_concentration = std::accumulate(concentrationsB.begin(), concentrationsB.end(), decltype(concentrationsB)::value_type(0));
    double total_monomer_concentration=total_A_concentration+total_B_concentration;

    // NOTE: moleculeA = type 0, moleculeB = type 1 in the vectors storing sequence
    std::vector<int> monomerA; // Track number of bifunctional monomer(s) containing functional group A
    std::vector<int> monomerB; // Track number of bifunctional monomer(s) containing functional group B
    
    // initialize the numbers of A and B monomers present based on concentrations and simulation size
    for (int i=0; i<concentrationsA.size();i++){
        monomerA.push_back(std::round(1.0*num_of_molecules_in_simulation*concentrationsA[i]/total_monomer_concentration));
    }
    for (int i=0; i<concentrationsB.size();i++){
        monomerB.push_back(std::round(1.0*num_of_molecules_in_simulation*concentrationsB[i]/total_monomer_concentration));
    }

    // for calculation of overall conversion of isocyanate
    const double total_initial_B_functional_groups = 2.0*(std::accumulate(monomerB.begin(), monomerB.end(), decltype(monomerB)::value_type(0))); // calculate the initial total number of B functional groups
    double B_groups_remaining=total_initial_B_functional_groups; // calculate remaining unreacted B functional groups

    // declaring variables
    const double Na = 6.02E23; // avogadro's number: items per mole
    const double R = 0.008314; // ideal gas constant, kJ/(K mol)
    double V = num_of_molecules_in_simulation/(Na*total_monomer_concentration); // Volume in L. For 1 mol/L the volume is L (which is also dm^3)
    double over_x; // overall conversion with respect to A functional groups!
    chain_pool all_chains; // vector of chain objects to store sequences and chain information in
    chain_pool loops; // vector of chain objects which have formed loops
    // Need to keep track of all unreacted functional groups, which are found in either monomers or chain ends
    std::vector<int> chainsA(monomermassA.size(),0); // number of chain ends of each A monomer type, initialize to 0. 
    // chainsA[0] stores # of A1 chain ends, [1] stores # of A2 chain ends, etc. 
    std::vector<int> chainsB(monomermassB.size(),0); // number of chain ends of each B monomer type, initialize to 0. 
    // chainsB[0] stores # of B1 chain ends, [1] stores # of B2 chain ends, etc.

    // For difunctional monomers, the first dyad in a chain will always have frontend = 0 and backend = 1.
    std::ofstream molwt;
    time_t t = time(0);   // get current time
    struct tm * now = localtime( & t );
    char date[80];
    strftime (date,80,"%Y%m%d_%I%M%S%p_%Z",now); 
    std::string str(date);
    std::string mweight = "molecular_weight";
    std::string path = "/Users/mimivirus/Documents/Research/Polymer_Recycling_Project/MiscellaneousResources/KMC_PU/Output/";
    std::string molwtfilename = path+date+mweight+filename;
    molwt.open(molwtfilename); // store time, Mn, Mw, Dispersity
    molwt << "time           Conversion     Mn         Mw         D\n";

    // declare variables needed to track molecular weight
    double sumMiNi=0; // "Mi" refers to the molecular weight of a particular chain, and Ni refers to the number of chains of that particular molecular weight
    double sumNi=0; // Ni refers to the number of chains of molecular weight Mi, so sumNi is equal to the total number of chains present
    double sumMi2Ni=0; // "Mi" refers to the molecular weight of a particular chain, and Ni refers to the number of chains of that particular molecular weight
    double dispersity; // polydispersity index (Mw/Mn)
    double Mn=0; // Number average molecular weight
    double Mw=0; // Weight average molecular weight
    bool isloop = false; // Track whether a given reaction event caused a polymer chain to form a loop
    bool isnewchain=false; // Track whether a given reaction event created a new polymer chain
    bool ismonomerA=false; // Track whether a given reaction event caused an A monomer (diol) to be consumed
    bool ismonomerB=false; // Track whether a given reaction event caused a B monomer (diisocyanate) to be consumed
    double Mi_A; double Mi_B; // When two chains combine, Mi_A is the molecular weight of the chain containing the A-type functional group, while Mi_B is the molecular weight of the chain containing the B functional group
        
    // calculate all rate constants from input kinetic parameters and temperature assuming Arrhenius relationship
    std::vector<double> kf(M);
    for (int i=0;i<M;i++) {
        kf[i]=A_kf[i]*exp(-Ea_kf[i]/(R*T)); // adding a to b
    }
    std::srand(std::time(0)); // seed the random number generator using the current time
    double time = 0; // Elapsed reaction time (seconds)
    // KMC loop
    while ((over_x<target_conversion && useConversion) || (!useConversion && time<simulation_time)) {

        // calculate propensity functions
        double c[M];
        for (int i=0;i<M;i++) {
            c[i] = kf[i]/(V*Na); // 1 / (molecules s)
        }
        // update total rate
        double total_rate=0;
        // note: Entries in R[v] to compute total rate are calculated in the order A1+B1, A1+B2, A1+B3... A1+Bn, A2+B1,A2+B2, A2+B3... etc.
        int v=0;
        double Rv[M];
        for (int j=0;j<monomerA.size();j++) {
            for (int k=0;k<monomerB.size();k++) {
                Rv[v]=c[v]*(2*monomerA[j]+chainsA[j])*(2*monomerB[k]+chainsB[k]); // [=] molecules/s. 2*monomers because each monomer is bifunctional
                total_rate+=Rv[v]; // molecules/s
                v+=1;
            }
        }
        // choose timestep tau
        double r1 = 1.0*std::rand()/RAND_MAX; 
        while (r1==0){
            // this loop is to ensure that r1 is not 0 (on my laptop, less than 1/2,000,000,000 chance per random number generation event)
            r1 = 1.0*std::rand()/RAND_MAX;
        }
        double tau = (1/total_rate)*log(1/r1); // calculate timestep tau
        
        // choose reaction to take place using Gillespie algorithm
        double r2 = 1.0*std::rand()/RAND_MAX;
        int mu=0;
        double sumRv=0; // Lin Wang eqn 1 multiplied by total rate
        while (sumRv<r2*total_rate) { 
            sumRv+=Rv[mu];
            if (sumRv<r2*total_rate) mu++;
        }
        // Translation from reaction channel mu to specific AA monomer type and BB monomer type 
        int A_type=0; int B_type=0; int count=0;
        reactionchannelselector(A_type, B_type, monomerA, monomerB, mu);
        
        // pick which A, B functional group
        double r3 = 1.0*std::rand()/RAND_MAX;
        while (r3==1) {
            // this whole loop is to ensure that r3 excludes 1, want whichA on range [0, # of functional groups) (on my mac, less than 1/2,000,000,000 chance per random number generation event)
            r3 = 1.0*std::rand()/RAND_MAX;
        }
        double whichA = r3*(chainsA[A_type]+2*monomerA[A_type]);
        double r4 = 1.0*std::rand()/RAND_MAX;
        while (r4==1){
            // this whole loop is to ensure that whichB excludes the upper bound, want number on range [0, # of functional groups) (on my mac, less than 1/2,000,000,000 chance per random number generation event)
            r4 = 1.0*std::rand()/RAND_MAX;
        }
        double whichB = r4*(chainsB[B_type]+2*monomerB[B_type]); 
        // update chains
        explicit_sequence_record(whichA,whichB,monomerA,monomerB,A_type,B_type,chainsA,chainsB,all_chains,loops,isloop,isnewchain,ismonomerA,ismonomerB,Mi_A,Mi_B,monomermassA,monomermassB);
        time += tau;
        // END KMC CALCULATIONS. 

        // Update Mn, Mw, and dispersity
        molecular_weight(Mn, Mw, all_chains, loops, monomerA, monomerB, monomermassA, monomermassB,isloop,isnewchain,ismonomerA,ismonomerB,sumNi,sumMiNi,sumMi2Ni,Mi_A,Mi_B);
        dispersity=Mw/Mn; // calculate polydispersity index PDI
        B_groups_remaining = 2.0*(std::accumulate(monomerB.begin(), monomerB.end(), decltype(monomerB)::value_type(0))); // Count unreacted functional groups on monomer
        B_groups_remaining += std::accumulate(chainsB.begin(), chainsB.end(), decltype(chainsB)::value_type(0)); // Count unreacted functional groups on chain ends
        // calculate overall conversion of B functional group
        over_x=1-(B_groups_remaining/total_initial_B_functional_groups); 
        // Output molecular weight at timestep to file
        molwt << std::left << std::setw(10) << time << "     " << std::setw(10) << over_x << "     " << std::setw(6) << Mn << "     " << std::setw(6) << Mw << "     " << std::setw(6) << dispersity << "     " << "\n";
        }

    molwt.close();

    // can easily print all sequences here too if desired
    return 0;
}