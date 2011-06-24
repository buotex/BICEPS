#include "Tuples.h"
#include "bicepsdefinitions.h"
#include <exception>
#include <set>
#define cout cerr
//#define penalty_mutation 0.0

namespace Pepsplice{
    extern Pool<Pepsplice::Tuple, POOLSIZE> __POOL__;    

    
    
    
    Tuples::Tuples(Services *se0, DnaAA *d1)
    {
        se1 = se0;
        dnaAA1 = d1;
        scoring1 = NULL;
        tuples = new vector<Tuple*>();
    
        //  oldtuples = new vector<Tuple*>();
        //outputlevel = se1->outputlevel;	
		
        //automatic buffer size adjustment
        starttrigger = 0; endtrigger = 0; 
        interval = 0;
        starttuples = 0;
        starttime = 0;
        bufferincrease = 0; lastbufferincrease = 0;
        tuples_per_second = 0; last_tuples_per_second = 0;
        performanceincrease = 0; lastperformanceincrease = 0;
        
      //  randomprotein = new Protein("RANDOM");
    }
    
    Tuples::~Tuples()
    {
        
        //Notes:
        //Tuples can have several pointers to them, especially because of the internal modifications.
        //But after going through all lists, no tuple should have a pointercount greater than 0?
        
        
        //delete randomprotein;
        for (vector<Tuple*>::iterator it = tuples->begin(); it != tuples->end(); ++it)
        {
            if ((*it)->pointercount > 0) cout << "still in tuples " << *it << endl;
            //(*it)->decreasePointerCount();
            delete *it;
        }
        delete tuples;

    }
    
    void Tuples::setScoring(Scoring *sc1)
    {
        scoring1 = sc1;
    }
    
	
    void Tuples::addTuple(Tuple* tup)
    {
        enumerateSNPs(tup);
    }    

    void Tuples::enumerateSNPs(Tuple* tup){
        //cout<<"\n starting enumerate";
        setPenalty(tup);
        //cout << "\nO: " << tup->getAASeq() << " " << tup->length;
        if(tup->penalty > se1->penalty_max){
            //changes BYR: added Hightup to keep to high results
            //Tuple *Hightup = new Tuple(tup);
            delete tup;
        }else{
//            if( (se1->doSNPnt || se1->doSNPaa) && tup->penalty <= se1->penalty_max - PEN_SNP){
            if(tup->penalty <= se1->penalty_max - PEN_SNP){
                for(int i = 0; i < tup->length; i++){
                    for(int j = 65; j < 65+26; j++){
                        bool condition = false;
                        //if(se1->doSNPnt == true && tup->ntsequence != NULL && dnaAA1->tripletSNPaa[i*3][i*3+1][i*3+2][j]==true) condition = true;
                        //start changes BYR
                       // if(se1->doSNPaa == true && tup->aasequence != NULL && tup->aasequence[i]!=j) condition = true;
                        if(tup->aasequence[i]!=j) condition = true;
                        if(condition == true){
                            Tuple *newtup = new Tuple(tup);
                            //cout << "\n" << newtup->getAASeq() << " " << newtup->thparentmassMH;
                            newtup->isSNP = true;
                            newtup->isSNP2 = false;
                            newtup->SNPpos = i;
                            newtup->SNPoldAA = tup->aasequence[newtup->SNPpos]; //store old amino acid
                            newtup->SNPnewAA = j;
                            newtup->aasequence[newtup->SNPpos] = newtup->SNPnewAA;
                            newtup->thparentmassMH = newtup->thparentmassMH  + dnaAA1->aa_monomasses256[newtup->SNPnewAA] - dnaAA1->aa_monomasses256[newtup->SNPoldAA]; //adjust parent mass
                            //added BYR for allowing for search for 2nd error
                            //cout<<"\n within loop, see current parametrization";
                            //cout<<"\n mods"<<newtup->modifications*1;
                            //cout<<"\n missCleav"<<newtup->getMissedCleavages()*1;
                            //cout<<"\n oxy"<<newtup->oxidizedmethionines*1;
                            //cout<<"\n current penalty"<<newtup->penalty;
                            //cout<<"\n isSNP"<<newtup->isSNP;
                            //cout<<"\n isSNP2"<<newtup->isSNP2;
                            //cout<<"\n within loop before setPenalty(newtup)";
                            setPenalty(newtup);
                            //cout<<"\n penalty newtup"<< newtup->penalty;
                            //cout<<"\n se1->penalty_max"<< se1->penalty_max;
                            
                            if(newtup->penalty > se1->penalty_max){          
                                //try {
                                delete newtup;                              
                                //} catch (std::exception& e) {
                                //	std::cerr << e.what() << std::endl;
                                //}
                            }else{
                                //if( (se1->doSNPnt || se1->doSNPaa) && newtup->penalty <= se1->penalty_max - PEN_SNP){
                                if(newtup->penalty <= se1->penalty_max - PEN_SNP){

                                    for(int ii = 0; ii < newtup->length;ii++){
                                        for(int jj = 65; jj < 65+26; jj++){
                                            bool condition = false;                    
                                            if(se1->doSNPnt == true && tup->ntsequence != NULL && dnaAA1->tripletSNPaa[ii*3][ii*3+1][ii*3+2][jj]==true) condition = true;
                                            if(se1->doSNPaa == true && newtup->aasequence != NULL && newtup->aasequence[ii]!=jj) condition = true;
                                            if(condition == true){
                                                Tuple *newtup2 = new Tuple(newtup);
                                                newtup2->isSNP = true;
                                                newtup2->isSNP2 = true;
                                                newtup2->SNPpos = ii;
                                                newtup2->SNPoldAA = newtup->aasequence[newtup2->SNPpos]; //store old amino acid
                                                newtup2->SNPnewAA = jj;
                                                newtup2->aasequence[newtup2->SNPpos] = newtup2->SNPnewAA;
                                                newtup2->thparentmassMH = newtup2->thparentmassMH  + dnaAA1->aa_monomasses256[newtup2->SNPnewAA] - dnaAA1->aa_monomasses256[newtup2->SNPoldAA]; //adjust parent mass 
                                                
                                                enumerateTerminalModifications(newtup2);
                                            }
                                        }
                                    }
                                }
                                enumerateTerminalModifications(newtup);	
                            }	
                        } 	
                    }
                }
                
            }
            
            enumerateTerminalModifications(tup);
        }
    }
    
    void Tuples::enumerateTerminalModifications(Tuple* tup){
        
        //if there are general modifications for the protein/peptide N- && C-terminus, create children
        //cout<<"\nenumerateTerminalMods";	
        setPenalty(tup);
        
        if(tup->penalty > se1->penalty_max){
            delete tup;
        }else{
            //if(se1->doModifications && tup->penalty <= se1->penalty_max - PEN_PTM){
            if(tup->penalty <= se1->penalty_max - PEN_PTM){
                
                if(dnaAA1->aamod_nbmodperaa['X'] > 0){
                    //go through all 'X'-modifications, i.e. modifications that are independent of amino acids, e.g. at the terminus
                    int startascii = dnaAA1->aamod_startascii['X'];
                    for(int j = startascii; j < startascii + dnaAA1->aamod_nbmodperaa['X']; j++){
                        signed char t = dnaAA1->aamod_type[j];
                        //general N-terminal modification
                        if( (t == '(') || (t == '[' && tup->protNterm) ){
                            Tuple *newtup = new Tuple(tup);
                            newtup->modifications++;
                            newtup->charNterm = j; //replace amino acid by modified one
                            newtup->thparentmassMH = newtup->thparentmassMH  + dnaAA1->aa_monomasses256[j]; //adjust parent mass
                            enumerateInternalModifications(newtup, 0); //recurse
                        }
                        //general C-terminal modification
                        if( t == ')' || (t == ']' && tup->protCterm) ){
                            Tuple *newtup = new Tuple(tup);
                            newtup->modifications++;
                            newtup->charCterm = j; //replace amino acid by modified one
                            newtup->thparentmassMH = newtup->thparentmassMH  + dnaAA1->aa_monomasses256[j]; //adjust parent mass
                            enumerateInternalModifications(newtup, 0); //recurse					
                        }
                    }
                }
                
            }
            
            //CAREFUL: MUST BE WITHIN ELSE! GIVES A NASTY ERROR HOURS AFTER START IF OUTSIDE
            enumerateInternalModifications(tup, 0);		
        }
        
        
    }
    
    //enumerate modifications in a tree-like fashion
    void Tuples::enumerateInternalModifications(Tuple* tup, int pos)
    {
        //cout<<"\n This is in enumerate Internal Modificatinos";
        setPenalty(tup);
        
        if(tup->penalty > se1->penalty_max){
            delete tup;
        }else{
            
            //enumerate modifications
            
            //if(se1->doModifications)
            //{
                
                //go through AA sequence
                for(int i = pos; i < tup->length; i++){
                    signed char aa_i = tup->aasequence[i];
                    
                    //if there are modifications for a particular amino acid, create children
                    if(dnaAA1->aamod_nbmodperaa[aa_i] > 0){
                        if(aa_i == 'M' || tup->penalty <= se1->penalty_max - PEN_PTM){
                            int startascii = dnaAA1->aamod_startascii[aa_i];
                            
                            //enumerate children for current amino acid
                            for(int j = startascii; j < startascii + dnaAA1->aamod_nbmodperaa[aa_i]; j++){
                                signed char t = dnaAA1->aamod_type[j];
                                if(  (t == '_')  ||  (t == '(' && i == 0)  ||  (t == ')' && i == tup->length - 1)  ||  (t=='[' && tup->protNterm == true)  ||  (t==']' && tup->protCterm == true) ){
                                    Tuple *newtup = new Tuple(tup);
                                    if(aa_i != 'M'){newtup->modifications++;}else{newtup->oxidizedmethionines++;}
                                    newtup->aasequence[i] = j; //replace amino acid by modified one
                                    newtup->thparentmassMH = newtup->thparentmassMH  + dnaAA1->aa_monomasses256[j] - dnaAA1->aa_monomasses256[aa_i]; //adjust parent mass
                                    enumerateInternalModifications(newtup, i + 1); //recurse
                                }
                            }
                        }
                    }
                }
                
            //}
            
            addTuple2(tup); //forward tuple
            
        }
    }
    
    
    
    
    //how exotic, how special is a tuple with respect to
    //tryptic digest, DNA search, splicing, modifications, missed cleavages
    void Tuples::setPenalty(Tuple *tup){
        ///cout<<"\n mods"<<tup->modifications*1;
        //cout<<"\n missCleav"<<tup->getMissedCleavages()*1;
        //cout<<"\n oxy"<<tup->oxidizedmethionines*1;
        //cout<<"\n current penalty"<<tup->penalty;
        //cout<<"\n isSNP"<<tup->isSNP;
        //cout<<"\n isSNP2"<<tup->isSNP2;
        //cout<<"\n sequence"<<tup->penalty;
        float penalty0 = 0;
        
        
//        if(tup->trypticstart == false) penalty0+=se1->penalty_tryptic;
//        if(tup->trypticend == false) penalty0+=se1->penalty_tryptic;
//        if(tup->chromosome != NULL) penalty0+=se1->penalty_genomic;
//        if(tup->isSpliced() == true) penalty0+=se1->penalty_spliced;
        if(tup->trypticstart == false) penalty0+=PEN_TRYPTIC;
        if(tup->trypticend == false) penalty0+=PEN_TRYPTIC;
        //if(tup->chromosome != NULL) penalty0+=PEN_GENOMIC;
        if(tup->isSpliced() == true) penalty0+=PEN_SPLICED;
        
        
        //cout<<"\n isSNPtrue"<<tup->isSNP;
        //cout<<"\n tup-SNP2"<<tup->isSNP2;
        //cout<<"\nright before SNP computation";
        if(tup->isSNP == true) penalty0+=dnaAA1->aa_pampenalty(tup->SNPnewAA,tup->SNPoldAA);
        //cout<<"\n penalty"<<penalty0;
        //cout<<"\n isSNP 2 true,but give SNP1"<<tup->isSNP;
        //cout<<"\n right before SNP2 computation";
        if(tup->isSNP2 == true) penalty0+=dnaAA1->aa_pampenalty(tup->SNPnewAA,tup->SNPoldAA);
        //cout<<"\n right after SNP2";
        penalty0 += tup->modifications*PEN_PTM;
        penalty0 += tup->getMissedCleavages()*PEN_MISSCLEAV;
        penalty0 += tup->oxidizedmethionines*PEN_METHOX;
        tup->penalty = penalty0;
        //cout << "\n now in setPenalty in Tuples" << tup->aasequence;
    }
    
    void Tuples::addTuple2(Tuple* tup)
    {
        if (se1->outputlevel > 2)
        {            
            se1->os << "Tuple added " << tup-> penalty << " " << tup->getAASeq() << " " << tup->aasequenceorig<< '\n';
        }
        tuples->push_back(tup); 
        //cout << "\n" << tup->getAASeq();
        se1->incrementTuples(false);
        //if(se1->randomtuples * se1->realtuples_per_randomtuple < se1->realtuples){
        //CHANGES BYR - removed dummyTupleconstruction
        //Tuple *dummytup = dummyTuple(tup);
        //tuples->push_back(dummytup);
        //se1->incrementTuples(true);
        //}
        if(tuples->size() > TUPLEBUFFER) forwardTuples();
    }
    
    // void Tuples::addTupleDeleted(Tuple* tup) //add Tuples which would otherwise be deleted
    // {
    //     tup->increasePointerCount();
    //     oldtuples->push_back(tup);	
    // }
        
    void Tuples::forwardTuples()
    {
        
        if(tuples->size() <= 0){
            //	cout << "\nTuples::forwardTuples tuples->size() <= 0, forward is thus skipped. This can happen when result files are written." << flush;
            return;
        }   
        
        //forward to scoring
        if(se1->dismisstuples == false){	
            //sortTuples();
            sort(tuples->begin(),tuples->end());
            scoring1->addTuples(tuples);
        }else{
            for(size_t i = 0; i < tuples->size(); i++){
                delete (*tuples)[i];
            }
        }
        tuples->clear();
        if (se1->outputlevel > 0) se1->suggestProgressReport();
        
        //optionally adapt tuple buffer
        //if(se1->adapttuplebuffer == true) adaptTupleBuffer(); //not viable for the current implementation, BX
        
    }
    
    void Tuples::adaptTupleBuffer(){
        
        //start measurement
        if(se1->tuples > starttrigger){
            interval = 100000 + se1->tuplebuffersize * 20 + tuples_per_second * 20;
            endtrigger = se1->tuples + interval;
            starttrigger = endtrigger + interval;
            
            starttime = se1->globaltime->timeSinceStart();
            starttuples = se1->tuples;
        }
        
        //end measurement, analyse, adjust buffersize
        if(se1->tuples > endtrigger){
            starttrigger = se1->tuples + interval;
            endtrigger = starttrigger + interval;
            
            double deltatime = (se1->globaltime->timeSinceStart() - starttime);
            if(deltatime <= 0) deltatime = 1;
            tuples_per_second = (se1->tuples - starttuples) / deltatime;
            se1->tuples_per_second = tuples_per_second;
            
            if(tuples_per_second < 1) tuples_per_second = 0.01;
            performanceincrease = (tuples_per_second - last_tuples_per_second) / last_tuples_per_second;
            
            //decide on increase || decrease
            if(performanceincrease * (lastbufferincrease-1) >= 0){
                bufferincrease = 1.25;
            }else{
                bufferincrease = 0.8;
            }
            
            //adjust buffersize
            //if(se1->force_tb_report == true) se1->doProgressReport();
            se1->tuplebuffersize = (long)((double)se1->tuplebuffersize * (double)bufferincrease);
            if(se1->tuplebuffersize < 1) se1->tuplebuffersize = 1;
            if(se1->tuplebuffersize > 100000) se1->tuplebuffersize = 100000;
            
            //store performance values for next calculation
            last_tuples_per_second = tuples_per_second;
            lastbufferincrease = bufferincrease;
        }
    }
    
    void Tuples::sortTuples(){
        sort(tuples->begin(), tuples->end(), TupleComparator());
        //for(int i = 0; i < 20; i++){cout << "\npointer " << (long)bm[i] << "  score " << bm[i]->score << "  tuple sequence " << bm[i]->tuple->aasequence;}
    }
    
    Tuple *Tuples::dummyTuple(Tuple *tup){
        Tuple *dt = new Tuple(tup); //this creates a new sequence too!
        //dt->chromosome = NULL;
        //dt->protein = randomprotein;
        dt->random = true;
        
        dt->aasequence[tup->length-1] = tup->aasequence[tup->length-1]; //set last amino acid (tryptic) as in original
        
        //either it gets reverted || shifted (always without last amino acid)
        if(se1->randomtype == 0){	
            for(int i = 0; i < tup->length-1; i++){
                dt->aasequence[i] = tup->aasequence[tup->length - 2 - i];
            }
        }else{
            for(int i = 0; i < tup->length-1; i++){
                dt->aasequence[i] = tup->aasequence[(i+se1->randomtype)%(tup->length-1)];
            }
        }
        //cout << "\n" << tup->getAASeq() << " " << dt->getAASeq();
        
        
        //cout << "\n" << dt->getAASeq() << "\t" << tup->getAASeq();
        
        return dt;
    }
    
    
}//namespace

