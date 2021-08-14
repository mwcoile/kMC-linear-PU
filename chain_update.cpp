#include"chain_update.h"
#include"molecular_weight.h"

void explicit_sequence_record(int whichA, int whichB, std::vector<int>& monomerA, std::vector<int>& monomerB, int& A_type, int& B_type, std::vector<int>& chainsA, std::vector<int>& chainsB, chain_pool& all_chains, chain_pool& loops, bool& isloop, bool& isnewchain, bool& ismonomerA,bool& ismonomerB, double& Mi_A, double& Mi_B, std::vector<double>& monomermassA, std::vector<double>& monomermassB) {
    // front and back refer to front and back of vector used to store a chain
    bool front = false; 
    bool back = false;
    Mi_A=0; // reset Mi_A (mass of A monomer or chain that reacts)
    Mi_B=0; // reset Mi_B (mass of B monomer or chain that reacts)
    // case for which 1 A monomer reacts with 1 B monomer
    if (whichA<2*monomerA[A_type] && whichB<2*monomerB[B_type]) {
        chain newchain; // make  a new chain
        newchain.frontend = 0; // i.e. front end is A
        newchain.backend = 1; // backend is defined arbitrarily as species B on new polymer chain
        Mi_A=monomermassA[A_type]; // record mass of reacted "A chain"
        Mi_B=monomermassB[B_type]; // record mass of reacted "B chain"
        newchain.chain_mass=Mi_A+Mi_B; // update total chain mass
        newchain.v.push_back({newchain.frontend,A_type});
        newchain.v.push_back({newchain.backend,B_type});
        all_chains.push_back(newchain);
        monomerA[A_type]--; // update the trackers
        monomerB[B_type]--; // update the trackers
        chainsA[A_type]++; // update the trackers
        chainsB[B_type]++;  // update the trackers
        isnewchain=true; // if two monomers react, then a new chain is formed
    }
    // case for which 1 A monomer reacts with 1 B chain
    else if (whichA<2*monomerA[A_type] && whichB>=2*monomerB[B_type]) {
        // STEP 1: Select B chain
        int B_functional_group = (int) (whichB-2*monomerB[B_type]); // Random int corresponding to a number on range [0, # of B functional groups) 
        int selected_chain = 0; // pick B chain index in all chains
        int countB = 0; // count number of B functional groups matching B_type

        while (selected_chain<all_chains.size()){
            // select which B chain and which B end by looking for the functional group corresponding to B_functional_group
            if (all_chains[selected_chain].backend==1 && all_chains[selected_chain].v.back()[1]==B_type) { // if B functional group on back AND correct type of B monomer
                countB++; // increment countB: another matching B functional group was found
            }
            // if we have now counted B_functional_groups+1, break from the loop, notating the end (back/front) of the chain and the chain itself (selected_chain)
            if (B_functional_group<countB) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==1 && all_chains[selected_chain].v.front()[1]==B_type){ // if B functional group on front AND correct type of B monomer
                countB++; // increment countB: another matching B functional group was found
            }
            // if we have now counted B_functional_groups+1, break from the loop, notating the end (back/front) of the chain and the chain itself (selected_chain)
            if (B_functional_group<countB) {
                front=true;
                break;
            }
            selected_chain++;
        }
        // update Mi_A and Mi_B
        // A monomer weight
        Mi_A=monomermassA[A_type];
        // Weight of the B chain 
        Mi_B=all_chains[selected_chain].chain_mass;
        // update total chain mass
        all_chains[selected_chain].chain_mass+=Mi_A;

        // STEP 2: which end of the chain has B?
        // CASE 1: if the B is at the end of the chain
        if (back) { 
            // then add A to the end
            all_chains[selected_chain].v.push_back({0,A_type}); 
            // and update the trackers accordingly
            all_chains[selected_chain].backend=0;
            chainsB[B_type]--; // there is one fewer chain end terminated in B
            chainsA[A_type]++; // there is one more chain end terminated in A
            monomerA[A_type]--; // there is one less unreacted A monomer
        }
        // CASE 2: if the B is at the front of the chain
        else if (front) { 
            // then add A to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),{0,A_type});
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=0;
            chainsB[B_type]--; // One less B chain end 
            chainsA[A_type]++; // One more A chain end
            monomerA[A_type]--; // One less unreacted A monomer
        }
        ismonomerA=true; // the reaction that happened consumed an A monomer
    }
    // case for which 1 B monomer reacts with 1 A chain
    else if (whichA>=2*monomerA[A_type] && whichB<2*monomerB[B_type]) {
        // STEP 1: Select A chain
        int A_functional_group = (int) (whichA-2*monomerA[A_type]);
        int selected_chain = 0;
        int countA = 0;
        while (selected_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_chain].backend==0 && all_chains[selected_chain].v.back()[1]==A_type) { // if A functional group on back AND correct type of A monomer
                countA++; // another matching A was found, increment the counter
            }
            // if we have now counted A_functional_groups+1, break from the loop, notating the end (back/front) of the chain and the chain itself (selected_chain)
            if (A_functional_group<countA) {
                back=true;
                break;
            }
            if (all_chains[selected_chain].frontend==0 && all_chains[selected_chain].v.front()[1]==A_type) { // if A functional group on front AND correct type of A monomer
                countA++; // another matching A was found, increment the counter
            }
            // if we have now counted A_functional_groups+1, break from the loop, notating the end (back/front) of the chain and the chain itself (selected_chain)
            if (A_functional_group<countA) {
                front=true;
                break;
            }
            selected_chain++; // loop through all chains until the functional group corresponding to A_functional_group is found
        }
        // update Mi_A and Mi_B
        // B monomer weight
        Mi_B=monomermassB[B_type];
        // Weight of the A chain 
        Mi_A=all_chains[selected_chain].chain_mass;
        // update total chain mass
        all_chains[selected_chain].chain_mass+=Mi_B;

        // STEP 2: which end of the chain has A?
        // CASE 1: if the A is at the end of the chain .. but what to do if both ends have A?
        if (back) {
            // then add B to the end
            all_chains[selected_chain].v.push_back({1,B_type}); // can I store a pointer to all_chains[selected_chain] as something shorter?
            // and update the trackers accordingly
            all_chains[selected_chain].backend=1;
            chainsB[B_type]++;
            chainsA[A_type]--;
            monomerB[B_type]--;
        }
        // CASE 2: if the A is at the front of the chain
        else if (front) { 
            // then add B to the front
            all_chains[selected_chain].v.insert(all_chains[selected_chain].v.begin(),{1,B_type});
            // and update the trackers accordingly
            all_chains[selected_chain].frontend=1;
            chainsB[B_type]++;
            chainsA[A_type]--;
            monomerB[B_type]--;
        }
        ismonomerB=true;
    }
    // case for which 1 A chain reacts with 1 B chain
    else {
        if (whichA<2*monomerA[A_type] || whichB<2*monomerB[B_type]) {
            throw std::runtime_error( "Ruh roh! Something bad happened when updating chain sequences" );
        }
        // STEP 1: select A chain
        int A_functional_group = (int) (whichA-2*monomerA[A_type]);
        int selected_A_chain = 0;
        int countA = 0;
        bool frontA=false;
        bool backA=false;

        while (selected_A_chain<all_chains.size()) {
            // select which A chain and which A end
            if (all_chains[selected_A_chain].backend==0 && all_chains[selected_A_chain].v.back()[1]==A_type){ // if A functional group on back AND correct type of A monomer
                countA++;
            }
            if (A_functional_group<countA) {
                backA=true;
                break;
            }
            if (all_chains[selected_A_chain].frontend==0 && all_chains[selected_A_chain].v.front()[1]==A_type){ // if A functional group on back AND correct type of A monomer
                countA++;
            }
            if (A_functional_group<countA) {
                frontA=true;
                break;
            }
            selected_A_chain++;
        }
        // STEP 2: select B chain
        int B_functional_group = (int) (whichB-2*monomerB[B_type]);
        int selected_B_chain = 0;
        int countB = 0;
        bool backB=false;
        bool frontB=false;
        while (selected_B_chain<all_chains.size()){
            // select which B chain and which B end
            if (all_chains[selected_B_chain].backend==1 && all_chains[selected_B_chain].v.back()[1]==B_type){ // if B functional group on back AND correct type of B monomer
                countB++;
            }
            if (B_functional_group<countB) {
                backB=true;
                break;
            }
            if (all_chains[selected_B_chain].frontend==1 && all_chains[selected_B_chain].v.front()[1]==B_type){ // if B functional group on front AND correct type of B monomer
                countB++;
            }
            if (B_functional_group<countB) {
                frontB=true;
                break;
            }
            selected_B_chain++;
        }
        // update Mi_A and Mi_B
        // A chain_weight
        Mi_A=all_chains[selected_A_chain].chain_mass;
        // Weight of the B chain 
        Mi_B=all_chains[selected_B_chain].chain_mass;
        // update total chain mass -- in this case, update both, because one of them will be deleted
        if (selected_A_chain != selected_B_chain) { // if a loop is formed, do nothing
            all_chains[selected_A_chain].chain_mass=Mi_A+Mi_B;
            all_chains[selected_B_chain].chain_mass=Mi_A+Mi_B;
        }

        // which chain is getting added to and which chain is getting deleted?
        bool add_to_chain_A = false;
        if (selected_A_chain<selected_B_chain) {
            add_to_chain_A = true;
        }
        // case 0: loop formation
        if (selected_A_chain==selected_B_chain){
            // add selected chain to chain pool of loops which cannot further react, then delete from the chain pool of all chains which can react
            // Note: This does not currently consider whether sterically it is possible for this to occur (i.e. loops consisting of 2 monomers are permitted)
        
            loops.push_back(all_chains[selected_A_chain]);
            all_chains.erase(all_chains.begin()+selected_A_chain); 
            isloop = true; // Record that the kMC event was loop formation for the molecular weight calculation--molecular weight of sample stays unchanged
        }
        // case 1: front of A chain to front of B chain
        else if (frontA && frontB) {
            if (add_to_chain_A) {
                // reverse A chain, append B chain to end of A chain, update the ends, delete B chain
                std::reverse(all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].frontend=all_chains[selected_A_chain].backend;
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // reverse B chain, append A to end of B chain, update the ends, delete A chain
                std::reverse(all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end()); 
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].frontend=all_chains[selected_B_chain].backend;
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        // case 2: front of A chain to back of B chain
        else if (frontA && backB) {
            if (add_to_chain_A) {
                // add A chain to B chain, then put it in the A chain location, then delete chain B
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                std::iter_swap(all_chains.begin()+selected_A_chain,all_chains.begin()+selected_B_chain); 
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // add A chain to back of B chain, delete A chain
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].backend;
                all_chains.erase(all_chains.begin()+selected_A_chain); 
            }
        }
        // case 3: back of A chain to front of B chain
        else if (backA && frontB) {
            if (add_to_chain_A) {
                // add B chain to back of A chain
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // add B chain to back of A chain, then put it in the B chain location, then
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].backend;
                std::iter_swap(all_chains.begin()+selected_A_chain,all_chains.begin()+selected_B_chain); 
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        // case 4: back of A chain to back of B chain
        else if (backA && backB) {
            if (add_to_chain_A) {
                // reverse chain B, add to end of chain A, delete chain B
                std::reverse(all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].v.insert(all_chains[selected_A_chain].v.end(),all_chains[selected_B_chain].v.begin(),all_chains[selected_B_chain].v.end());
                all_chains[selected_A_chain].backend=all_chains[selected_B_chain].frontend;
                all_chains.erase(all_chains.begin()+selected_B_chain);
            }
            else {
                // reverse chain A, add to end of chain B, delete chain A
                std::reverse(all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].v.insert(all_chains[selected_B_chain].v.end(),all_chains[selected_A_chain].v.begin(),all_chains[selected_A_chain].v.end());
                all_chains[selected_B_chain].backend=all_chains[selected_A_chain].frontend;
                all_chains.erase(all_chains.begin()+selected_A_chain);
            }
        }
        else
        {
            throw std::runtime_error( "Something bad happened when updating chain sequences" );
        }
        chainsA[A_type]--;
        chainsB[B_type]--;
    }
}