# linear PU AB polymerization model

Welcome to documentation for the simple kinetic Monte Carlo (kMC) polyaddition code. This can be expanded to additional cases, as described by the paper: https://doi.org/10.1002/mats.202100058

How to run the model in a bash shell:
1. If you have not already compiled the model, run `make` to compile the model
2. Open the file input.txt, and enter the reaction conditions, kinetic parameters, and simulation size.
3. Run the model: `./model` Output results will go in the folder "./Output".

This code is designed to simulate the one-shot synthesis of linear polyurethanes. These are typically produced from a diisocyanate, a short chain extender diol, and a long soft diol. For generalizability, these diols and diisocyanates monomers are referred to as A or B, respectively. So for the scenario of 1 isocyanate and two diol types (chain extender and soft diol), there would be two A type monomers, and one B type monomer.

Sequences of these monomers are stored inside a polymer chain vector. Each entry in a polymer chain vector contains two numbers which track the identity of the monomer at that position. The first number is a 0 or 1 tracking whether the monomer is an A or B monomer, respectively, while the second number allows differentiation if multiple isocyanate or multiple diols are present. For example, if chain extender diol and soft diol are present and input in that order in input.txt, then they would be numbered 0 and 1 respectively. So the pair [0,0] would refer to diol chain extender, while the pair [0,1] would refer to a soft diol. 

A vector termed all_chains stores all the individual polymer chains in the simulation. Loops are allowed to form with the ends of the polymer chains assumed to react with each other with the same statistical likelihood that they react with any other functional group of that kind--which means that loop formation is generally negligible until full conversion is attained 

Much of the code found in the chain_update.cpp file is straightforward vector manipulation to update these vectors as chains, as illustrated in the figure below:

![image](https://github.com/mwcoile/kMC-linear-PU/blob/master/modelSchematic2.png)

The parts of the chain_update.cpp file that are not vector manipulation are selecting the specific functional groups to react. Essentially, if there are n functional groups of a particular type (e.g. unreacted chain extender alcohols), one of them must be chosen at random to be reacted.

The molecular_weight.cpp file updates the number average and weight average molecular weights after each time step. Monomers are not considered in this calculation, but dimers and all larger oligomers are included in the calculation of molecular weight as written here. 

The code may be readily modified to include reverse reactions, unequal reactivity of 1st and 2nd functional groups to react, or to track additional reactions such as those that might occur in the presence of a small amount of water. 

Sequences can be output at the end of the main function. For example, if 3 monomers are used, one diisocyanate (a B monomer), one chain extender diol (an A monomer), and one soft diol (another type of A monomer):

```
// record all chain sequences explicitly for CLD post processing and analysis
std::ofstream exp_sequences;
std::string seq = "seq";
std::string seqfilename=path+date+seq+s+filename;
exp_sequences.open(seqfilename,std::ios::app);
for (int i=0;i<all_chains.size();i++) {
    for (int j=0;j<all_chains[i].v.size();j++){
        if (j==0){
            // if this is a new sequence, start a new line
            exp_sequences << "\n";
        }
        if (all_chains[i].v[j][0]==0 && all_chains[i].v[j][1]==0) {
            // if it's an A0, label that as a "B"
            exp_sequences << "B";
        }
        else if (all_chains[i].v[j][0]==0 && all_chains[i].v[j][1]==1){
            // if it's an A1, label that as a "C"
            exp_sequences << "C";
        }
        else if (all_chains[i].v[j][0]==1) {
            // if it's a B0, label that as an "A"
            exp_sequences << "A";
        }
    }
}
```
