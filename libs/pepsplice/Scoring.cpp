#include "Scoring.h"
namespace Pepsplice{
    Scoring::Scoring(Services *se0):dnaAA1(DnaAA()),se1(se0) //Changes by BX
    {
        //Begin Changes BX
        hypergeometric1 = se1->hypergeometric1;
        //service objects
        //	se1 = se0;
        //	dnaAA1 = DnaAA();
        //	hypergeometric1 = Hypergeometric();
        
        //End Changes BX
        
        spectra1 = NULL; //spectra is set by main
        vtuples = NULL; //is set by addTuples()
        
        maxparentmass = 100000;
        aamaxlen = maxparentmass / 57;
        
        //cumulatedpeaksmax = se1->processorcacheL2 / 4; //int is 4 bytes
        
        //fragment ions
        //start changes BYR: added double length of b/yions to support handling of charge state 3
        bions3 = new float[2*aamaxlen+1];
        yions3 = new float[2*aamaxlen+1];
        bionsint3 = new int[2*aamaxlen+1];
        yionsint3 = new int[2*aamaxlen+1];
        bions = new float[aamaxlen];
        yions = new float[aamaxlen];
        bionsint = new int[aamaxlen];
        yionsint = new int[aamaxlen];	
        //discretized theoretical spectrum
        thspecdisc2 = new float[maxparentmass];
        for(int i = 0; i < maxparentmass; i++){thspecdisc2[i] = 0;} //initialize
        
        //discretized theoretical spectrum
        thspec = new int[maxparentmass];
        for(int i = 0; i < maxparentmass; i++){thspec[i] = 0;} //initialize
        thspec_prefix = new int[maxparentmass];
        for(int i = 0; i < maxparentmass; i++){thspec_prefix[i] = 61;} //initialize
        thspec_suffix = new int[maxparentmass];
        for(int i = 0; i < maxparentmass; i++){thspec_suffix[i] = 61;} //initialize	
        
        //Hypergeometric score Sadygov/Yates 2003
        //thspec4 = new int[maxparentmass];
        //for(int i = 0; i < maxparentmass; i++){thspec4[i] = 0;} //initialize
        
		
        //parent mass tolerance
        //pmtolneg = se1->masstol_below;  //below means that the theoretical parent mass may be that much below the measured parent mass
        //pmtolpos = se1->masstol_above;
		
    }
    
    Scoring::~Scoring() //Changes by BX
    {
        delete[] bions3;
        delete[] yions3;
        delete[] bionsint3;
        delete[] yionsint3;
        delete[] bions;
        delete[] yions;
        delete[] bionsint;
        delete[] yionsint;
        delete[] thspec;
        delete[] thspec_prefix;
        delete[] thspec_suffix;
        delete[] thspecdisc2;
    }
    
    void Scoring::setSpectra(Spectra *sp1)
    {
        spectra1 = sp1;
    }
    
    void Scoring::addTuples(vector<Tuple*> *vtuples0)
    {
        vtuples = vtuples0;
        //for(int i = 0; i < vtuples->size(); i++){scoreTuple((*vtuples)[i]);}
        //cumulatedpeaksmax is defined in Scoring.h, Services.cpp
        if(se1->tupspec == 0) scoreTuplesAmortizedCacheAware2(true); //bool slicealways
        if(se1->tupspec == 1) scoreTuplesAmortizedCacheAware2(false); //bool slicealways
        //if(se1->tupspec == 2) recursivelyScoreTuples(true); //bool 1by1
        //if(se1->tupspec == 3) recursivelyScoreTuples(false); //bool 1by1
        //scoreTuplesAmortizedCacheAware1();
    }
    
    void Scoring::calculateBYIons(Tuple *tup, int len)
    {
        //ions are shifted by 1 with respect to aasequence
        float pm1 = 0;
        //cout<<"\n scoring " << tup->aasequence;	
        //cout << "\nlen:" << len;
        //cout << tup;
        bions[0] = dnaAA1.monomassH + dnaAA1.aa_monomasses256float[tup->charNterm];
        yions[0] = dnaAA1.monomassHOHH + dnaAA1.aa_monomasses256float[tup->charCterm];
        
        for(int i = 1; i <= len; i++){
            //if(i > aamaxlen - 1) cout << "\nScoring::calculateBYIons: i is too large";
            bions[i] = bions[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[i - 1]];
            yions[i] = yions[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[len - i]];
            //cout << "\n" << bions[i] << " " << yions[i] << " " << aasequence[i - 1] << " " << aasequence[len - i] << " " << (long)aasequence[i - 1];
        }
    }
    void Scoring::calculateBYIons3(Tuple *tup, int len)
    {       //changes BYR function to allow charge state 3
        //ions are shifted by 1 with respect to aasequence
        //cout<<"\n scoring " << tup->aasequence;
        
        float pm1 = 0;
        //cout << "\nlen:" << len;
        bions3[0] = dnaAA1.monomassH + dnaAA1.aa_monomasses256float[tup->charNterm];
        yions3[0] = dnaAA1.monomassHOHH + dnaAA1.aa_monomasses256float[tup->charCterm];
        for(int i = 1; i <= len; i++){
            //if(i > aamaxlen - 1) cout << "\nScoring::calculateBYIons: i is too large";
            bions3[i] = bions3[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[i - 1]];
            yions3[i] = yions3[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[len - i]];
        }	
        
        bions3[len+1] = (dnaAA1.monomassH + dnaAA1.aa_monomasses256float[tup->charNterm])/2;
        yions3[len+1] = (dnaAA1.monomassHOHH + dnaAA1.aa_monomasses256float[tup->charCterm])/2;
        for(int i = len+2; i <= 2*len+2; i++){                                                                                                                                             //if(i > aamaxlen - 1) cout << "\nScoring::calculateBYIons: i is too large";
            bions3[i] = bions3[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[i - 1]];
            yions3[i] = yions3[i - 1] + dnaAA1.aa_monomasses256float[tup->aasequence[2*len +2 - i]]/2;
        }
    }	
    void Scoring::calculateBYIonsInt(Tuple *tup, int len){
        //ions[i] are shifted by 1 with respect to aasequence
        //cout<<"\n scoring " << tup->aasequence;
        //cout << "\nlen:" << len;
        bionsint[0] = dnaAA1.integermassH + dnaAA1.aa_integermasses256[tup->charNterm];
        yionsint[0] = dnaAA1.integermassHOHH + dnaAA1.aa_integermasses256[tup->charCterm];
        
        for(int i = 1; i <= len; i++){
            //if(i > aamaxlen - 1) cout << "\nScoring::calculateBYIons: i is too large";
            bionsint[i] = bionsint[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[i - 1]];
            yionsint[i] = yionsint[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[len - i]];
            //if(bionsint[i] <= 10) cout << "\nScoring::calculateBYIonsInt" << bionsint[i] << " " << yionsint[i] << " " << aasequence[i - 1] << " " << aasequence[len - i] << " " << (long)aasequence[i - 1];
        }
    }
    void Scoring::calculateBYIonsInt3(Tuple *tup, int len){
        // change BYR added function to treat chargestate 3
        //ions[i] are shifted by 1 with respect to aasequence
        //
        //cout << "\nlen:" << len;
        bionsint3[0] = dnaAA1.integermassH + dnaAA1.aa_integermasses256[tup->charNterm];
        yionsint3[0] = dnaAA1.integermassHOHH + dnaAA1.aa_integermasses256[tup->charCterm];
        //
        for(int i = 1; i <= len; i++){
            //if(i > aamaxlen - 1) cout << "\nScoring::calculateBYIons: i is too large";
            bionsint3[i] = bionsint3[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[i - 1]];
            yionsint3[i] = yionsint3[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[len - i]];
        }
        
        bionsint3[len+1] = (dnaAA1.integermassH + dnaAA1.aa_integermasses256[tup->charNterm])/2;
        yionsint3[len+1] = (dnaAA1.integermassHOHH + dnaAA1.aa_integermasses256[tup->charCterm])/2;
        for(int i = len+2; i <= 2*len+1; i++){
            bionsint3[i] = bionsint3[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[i - 2 - len]]/2;
            yionsint3[i] = yionsint3[i - 1] + dnaAA1.aa_integermasses256[tup->aasequence[2*len +2 - i - 1]]/2;
        }								    
    }
    
    
    
    
    //slices either always or only when it encounters a tuple that contains too many spectra
//    void Scoring::scoreTuplesAmortizedCacheAware2(bool slicealways){ //TODO BX change loop here because we only look at 2 spectra
//    	int cumulatedbytesmax = se1->processorcacheL2;
//    	
//    	int lenspectra = spectra1->spectra.size();
//    	int lentuples = vtuples->size();
//    
//    	//prevent tuples from getting deleted while they are still in the scoring process
//    	for(int i = 0; i < lentuples; i++){(*vtuples)[i]->isBeingScored();} 
//    	
//    
//    	int starttuple = 0;
//    	int startspectrum = 0;
//    	int cumulatedbytes = 0;
//    	int stopspectrum = lenspectra - 1;
//    
//    	bool repeatneeded = true;
//    	int k = 0; //repetitions for consecutive tuple-spectrum slices
//    	int i = 0; //tuples
//    	int j = 0; //spectra
//    
//    	//score tuple-spectrum slices which fit into cache, one after the other
//        while(repeatneeded == true && startspectrum < lenspectra){
//    		k++;
//    		repeatneeded = false; //repeat only if requested again
//    		
//    		i = starttuple; //from last slice
//    		if(slicealways == true) cumulatedbytes = 0;
//    		while(i < lentuples){
//    					
//    			pmmin = (*vtuples)[i]->thparentmassMH - pmtolneg;
//    			pmmax = (*vtuples)[i]->thparentmassMH + pmtolpos;
//                if (se1->outputlevel > 2) {
//                    se1->os << (*vtuples)[i]->aasequence << " " << pmmin <<" " << pmmax << '\n';
//                }
//    			//update lower bound; startspectrum is never decreased at any time		
//    			j = startspectrum; //from last tuple
//    			while(j < lenspectra && spectra1->spectra[j]->parentmassMH < pmmin){
//    				j++;
//    			}
//    			startspectrum = j; //save for next tuple
//    			
//    			//if slice is finished, exit i loop
//    			if(startspectrum > stopspectrum){
//    				startspectrum = stopspectrum + 1; //restart right after finished slice; starttuple was defined in if(cumulatedpeaks...)
//    				stopspectrum = lenspectra - 1; //equivalent to no stop
//    				//repeatneeded = true; //continue k loop
//    				break;
//    			}
//    
//    			//if there is/are at least one spectrum within the tolerance, score it or them
//    			if(j < lenspectra && spectra1->spectra[j]->parentmassMH < pmmax){
//    								
//    				//calculate theoretical spectrum
//    				//changes BYR: added treatment of chargestate3
//    				calculateBYIonsInt((*vtuples)[i], (*vtuples)[i]->length);
//    				calculateBYIonsInt3((*vtuples)[i], (*vtuples)[i]->length);
//    
//    				//float thspec_avgint = calculateThSpec3((*vtuples)[i], false, thspec3); //erase == false		
//    				calculateThSpec((*vtuples)[i], false, thspec);
//    				
//    				//spectrum loop
//    				if(slicealways == false) cumulatedbytes = 0;
//    				while(j < lenspectra && j <= stopspectrum && spectra1->spectra[j]->parentmassMH < pmmax){
//    				
//    					//check upper bound for cache limit
//    					if(cumulatedbytes > cumulatedbytesmax && repeatneeded == false){
//    						stopspectrum = j;
//    						starttuple = i;
//    						repeatneeded = true;
//    					}
//    					
//    					//score spectra
//    					calculateScore((*vtuples)[i], spectra1->spectra[j], thspec, 0);
//    					if(se1->scoring < 7){
//    						cumulatedbytes += spectra1->spectra[j]->lengthB*4;
//    					}else{
//    						cumulatedbytes += spectra1->spectra[j]->parentmassMHB;
//    					}
//    				
//                        j++;
//    				}//end while j
//    				
//    				//reset theoretical spectrum
//    				calculateThSpec((*vtuples)[i], true, thspec);
//    				//calculateThSpec4((*vtuples)[i], true, thspec4);
//    				//calculateThSpec3((*vtuples)[i], true, thspec3);  //erase == true
//    				
//    			}//end if stopspectrum and pmmax
//    				
//    		i++;	
//    		}//end while i
//    	}//end while repeatneeded
//    	
//    	for(int i = 0; i < lentuples; i++){(*vtuples)[i]->scoringFinished();} //tuple will delete itself if pointercount is <= 0
//    	
//    }
    

    
    void Scoring::scoreTuplesAmortizedCacheAware2(bool slicealways)
    { 
        //int cumulatedbytesmax = se1->processorcacheL2;
        
        int lenspectra = spectra1->spectra.size();
        int lentuples = vtuples->size();
        
        float thparentmin = spectra1->spectra[0]->parentmassMH / (1. + se1->masstolfactor);
        float thparentmax = spectra1->spectra[lenspectra-1]->parentmassMH / (1. - se1->masstolfactor);
        
        bool charge2 = (spectra1->spectra[0]->chargestate == 2);
        
        
        //        //prevent tuples from getting deleted while they are still in the scoring process
        //
        //for(int i = 0; i < lentuples; i++) (*vtuples)[i]->isBeingScored(); 
        
        for(int i = 0; i < lentuples; (*vtuples)[i]->scoringFinished(),i++)
        {
            
            (*vtuples)[i]->isBeingScored();
            
            if ( (*vtuples)[i]->thparentmassMH < thparentmin) continue;
            if ( (*vtuples)[i]->thparentmassMH > thparentmax) continue;
            //changes BYR: added treatment of chargestate3
            if (charge2) calculateBYIonsInt((*vtuples)[i], (*vtuples)[i]->length);
            else calculateBYIonsInt3((*vtuples)[i], (*vtuples)[i]->length);
            
           // calculateThSpec((*vtuples)[i], false, thspec);
            for (int j = 0; j < lenspectra; ++j)
                if ((*vtuples)[i]->thparentmassMH * (1. - se1->masstolfactor) < spectra1->spectra[j]->parentmassMH &&
                    (*vtuples)[i]->thparentmassMH * (1. + se1->masstolfactor) > spectra1->spectra[j]->parentmassMH )
                    {
                        calculateScore((*vtuples)[i], spectra1->spectra[j], thspec, 0);
                    }
            //calculateThSpec((*vtuples)[i], true, thspec);
            
        }
        //for(int i = 0; i < lentuples; i++) (*vtuples)[i]->scoringFinished();
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
     void Scoring::recursivelyScoreTuples(bool onebyone)
     {
     //cout << "\nScoring::recursively...";
     int lenspectra = spectra1->spectra.size();
     int lentuples = vtuples->size();
     if(lenspectra <= 0 or lentuples <= 0){
     cout << "\nScoring::recursivelyScoreTuples: no spectra or no tuples to score here. Progress report may be a reason." << flush;
     return;
     }
     
     for(int i=0; i < lentuples; i++){(*vtuples)[i]->increasePointerCount();}
     
     //cumulated peaks over all spectra
     spcumpeaks = new int[lenspectra + 2];
     spcumpeaks[0] = 0;
     for(int i = 0; i < lenspectra; i++){
     spcumpeaks[i + 1] = spcumpeaks[i] + spectra1->spectra[i]->lengthB;
     }
     
     //cumulated peaks over all tuples (theoretical spectra)
     tpcumpeaks = new int[lentuples + 2];
     tpcumpeaks[0] = 0;
     for(int i = 0; i < lentuples; i++){
     tpcumpeaks[i + 1] = tpcumpeaks[i] + (int)((*vtuples)[i]->thparentmassMH/dnaAA1.discretization_factor + 0.5);
     }
     
     //coordinates for rectangles in recursion
     int sp1 = 0;
     int sp2 = lenspectra - 1;
     int tp1 = 0;
     int tp2 = lentuples - 1;
     
     //score all tuples vs all spectra
     if(onebyone == true){
     //recurse1by1(sp1, sp2, tp1, tp2);
     }else{
     recurse1byx(sp1, sp2, tp1, tp2);
     }
     
     //delete tuples and dynamic arrays after use
     delete[] spcumpeaks;
     delete[] tpcumpeaks;
     for(int i=0; i < lentuples; i++){(*vtuples)[i]->deleteThSpec();}
     //cout << "\nScoring after deleteThSpec" << flush;
     for(int i=0; i < lentuples; i++){(*vtuples)[i]->decreasePointerCount();}
     //cout << "\nScoring after decreasePointerCount()" << flush;
     }
     
     bool Scoring::recursiveIntersectsArea(int sp1, int sp2, int tp1, int tp2){
     
     bool within = false;
     
     //pm square
     double sp1pmsq = spectra1->spectra[sp1]->parentmassMH;
     double sp2pmsq = spectra1->spectra[sp2]->parentmassMH;
     double tp1pmsq = (*vtuples)[tp1]->thparentmassMH;
     double tp2pmsq = (*vtuples)[tp2]->thparentmassMH;
     
     //pm tolerance boundary
     double sp1pmbd = tp1pmsq - pmtolneg;
     double sp2pmbd = tp2pmsq + pmtolpos;
     
     if(sp1pmsq <= sp2pmbd and sp2pmsq >= sp1pmbd) within = true;
     return within;
     
     }
     
     /*
     void Scoring::recurse1by1(int sp1, int sp2, int tp1, int tp2){
     
     int splen = sp2 - sp1 + 1;
     int tplen = tp2 - tp1 + 1;
     
     //cout << "\nsp1 sp2 tp1 tp2: " << sp1 << " " << sp2 << " " << tp1 << " " << tp2 << " splen: " << splen << " tplen: " << tplen << flush;
     
     if(recursiveIntersectsArea(sp1, sp2, tp1, tp2)){
     if(splen == 1 and tplen == 1){
     //cout << "\nSCORE sp1 tp1: " << sp1 << " " << tp1 << flush;
     scoreTupleVsSpectrum1by1(sp1, tp1);
     }else{
     
     //halve where it is more urgent
     int sppeaks = spcumpeaks[sp2 + 1] - spcumpeaks[sp1]; //spcumpeaks[0] = 0
     int tppeaks = tpcumpeaks[tp2 + 1] - tpcumpeaks[tp1]; //tpcumpeaks[0] = 0
     
     //cout << "\n sppeaks: " << sppeaks << " tppeaks: " << tppeaks << flush;
     
     if(sppeaks > tppeaks or tplen == 1){
     if(splen > 1){
     int sp12 = sp1 + splen/2;
     recurse1by1(sp1, sp12-1, tp1, tp2);
     recurse1by1(sp12, sp2, tp1, tp2);
     }
     }
     if(sppeaks <= tppeaks or splen == 1){
     if(tplen > 1){
     int tp12 = tp1 + tplen/2;
     recurse1by1(sp1, sp2, tp1, tp12 - 1);
     recurse1by1(sp1, sp2, tp12, tp2);
     }
     }
     }
     }
     
     }
     
     
     void Scoring::scoreTupleVsSpectrum1by1(int sp1, int tp1)
     {	//changes BYR added sp1->chargestate2 in call to calculateBYIonsInt
     Spectrum *spec = spectra1->spectra[sp1];
     Tuple *tup = (*vtuples)[tp1];
     //bool chargestate2=se1->chargestate2
     //cout << "\nspec: " << spec->parentmassMH << " tup: " << tup->thparentmassMH;
     //tup->increasePointerCount(); //give right of existence to tuple until all spectra are done
     
     float thspec_avgint = 0;
     
     //calculate theoretical spectrum (each time or only if necessary)
     if(true){
     if(tup->thspec == NULL){
     calculateBYIonsInt(tup->aasequence, tup->length);
     calculateBYIonsInt3(tup->aasequence, tup->length);
     
     int len = (int)(tup->thparentmassMH/dnaAA1.discretization_factor) + se1->safety_margin_disc;
     tup->initializeThSpec(len);
     tup->thspec_avgint = calculateThSpec(tup, false, tup->thspec); //erase == false
     }
     calculateScore(tup, spec, tup->thspec, tup->thspec_avgint);	
     }else{
     //generate theoretical spectrum each time, in cache
     calculateBYIonsInt(tup->aasequence, tup->length);
     calculateBYIonsInt3(tup->aasequence, tup->length);
     
     thspec_avgint = calculateThSpec(tup, false, thspec); //erase == false	
     calculateScore(tup, spec, thspec, thspec_avgint);
     thspec_avgint = calculateThSpec(tup, true, thspec);  //erase == true
     }
     }
     
     
     void Scoring::recurse1byx(int sp1, int sp2, int tp1, int tp2){
     
     int splen = sp2 - sp1 + 1;
     int tplen = tp2 - tp1 + 1;
     
     if(recursiveIntersectsArea(sp1, sp2, tp1, tp2)){
     int sppeaks = spcumpeaks[sp2 + 1] - spcumpeaks[sp1]; //spcumpeaks[0] = 0
     int tppeaks = tpcumpeaks[tp2 + 1] - tpcumpeaks[tp1]; //tpcumpeaks[0] = 0
     
     if(tplen == 1){
     scoreTupleVsSpectrum1byx(tp1, sp1, sp2);
     }else{
     
     //halve where it is more urgent
     if(splen > tplen and sppeaks*4 > se1->processorcacheL2){ //4 bytes per int
     if(splen > 1){
     int sp12 = sp1 + splen/2;
     recurse1byx(sp1, sp12-1, tp1, tp2);
     recurse1byx(sp12, sp2, tp1, tp2);
     }
     }else{
     if(tplen > 1){
     int tp12 = tp1 + tplen/2;
     recurse1byx(sp1, sp2, tp1, tp12 - 1);
     recurse1byx(sp1, sp2, tp12, tp2);
     }
     }
     
     }
     }			
     }
     
     
     void Scoring::scoreTupleVsSpectrum1byx(int tp1, int sp1, int sp2)
     {
     Tuple *tup = (*vtuples)[tp1];
     Spectrum *spec;
     double pmmin = tup->thparentmassMH - pmtolneg;
     double pmmax = tup->thparentmassMH + pmtolpos;
     float thspec_avgint = 0;
     
     //calculate theoretical spectrum
     //changes BYR: added se1->chargestate1 in call of calculateBYIons
     calculateBYIonsInt(tup->aasequence, tup->length);	
     calculateBYIonsInt3(tup->aasequence, tup->length);
     thspec_avgint = calculateThSpec(tup, false, thspec); //erase == false
     for(int i = sp1; i <= sp2; i++){
     spec = spectra1->spectra[i];
     if(spec->parentmassMH > pmmin and spec->parentmassMH < pmmax){
     calculateScore(tup, spec, thspec, thspec_avgint);
     }
     }
     thspec_avgint = calculateThSpec(tup, true, thspec);  //erase == true
     }
     
     */
    
    
    //*************************************************************************************
    
    
    float Scoring::calculateThSpec(Tuple *tup, bool erase, int *thspec){
        float thspec_avgint = 0;
        
		if(se1->scoring == 7){
			return 0;
		}else if(se1->scoring == 3){
			thspec_avgint = calculateThSpec3(tup, erase, thspec);
		}else if(se1->scoring == 4){
			calculateThSpec4(tup, erase, thspec);
		}else if(se1->scoring == 5){
			calculateThSpec5(tup, erase, thspec);
		}else if(se1->scoring == 6){
			calculateThSpec4(tup, erase, thspec);
			calculateThSpec6(tup, erase, thspec);
		}
        return thspec_avgint;
    }
    
    void Scoring::calculateScore(Tuple *tup, Spectrum *spectrum, int *thspec, float thspec_avgint){
        
        float score = 0;
        int sharedpeakcount = 0;
        
        //calculate score once, or multiple times for timing experiments
        for(int i = 0; i < 1 + se1->extraload; i++){		
            if(se1->scoring == 7){
                score = calculateScore7(tup, spectrum, thspec, &sharedpeakcount);
            }else if(se1->scoring == 3){
                score = calculateScore3(spectrum, thspec, thspec_avgint);
            }else if(se1->scoring == 4){
                score = calculateScore4(tup, spectrum, thspec);
            }else if(se1->scoring == 5){
                score = calculateScore5(tup, spectrum, thspec);
            }else if(se1->scoring == 6){
                score = calculateScore6(tup, spectrum, thspec);
                
            }else if(se1->scoring == 9){
                score = 0;
            }
        }
        
        se1->scores++;
        spectrum->scored++;
        
        //add match to best matches of active spectrum if good enough
        if(score >= spectrum->lowerBoundScore && score > 0){
            spectrum->bestmatches1->addMatch(new Match(spectrum, tup, score, sharedpeakcount ));
            se1->incrementMatches();
        }
    }
    
    void Scoring::calculateThSpec2(Tuple *tup){
        int len = tup->length;
        int pm = (int)(tup->thparentmassMH/dnaAA1.discretization_factor);
        int currention = 0;
        int series = 0;
        
        for(int i = 0; i <= pm + 1; i++){
            thspecdisc2[i] = 0;
        }
        
        while(series < 2){
            for(int i = 1; i < len; i++){
                if(series == 0){
                    currention = bionsint[i]; //if series == 0
                }else{
                    currention = yionsint[i]; //if series == 1
                }
                setCurrentIon2(currention-18, 10);
                setCurrentIon2(currention-17, 10);
                setCurrentIon2(currention-1 , 25);
                setCurrentIon2(currention   , 50);
                setCurrentIon2(currention+1 , 25);
            }
            series++;
        }
        
        //integrate over theoretical temppeaks1
        for(int i = 1; i <= pm + 1; i++){
            thspecdisc2[i] += thspecdisc2[i-1];
        }
        
    }
    
    void Scoring::setCurrentIon2(int pos, int weight){
        
        if(pos >= 75){
            thspecdisc2[pos-75] += -weight;
        }else{
            thspecdisc2[     0] += -weight;
        }
        thspecdisc2[pos]        +=  weight*151;
        thspecdisc2[pos+1]      += -weight*151;
        thspecdisc2[pos+76]     += +weight;
        
    }
	
    float Scoring::calculateScore2(Spectrum *spectrum, int seqlen){
        float score = 0;
        /*
         for(int i = 0; i < spectrum->lengthA0; i++){
         score += thspecdisc2[spectrum->mzsA0_disc[i]] * spectrum->intensitiesA0norm[i];
         }
         */
        return score;
    }
    
    //*************************************************************************************
    
    float Scoring::calculateThSpec3(Tuple *tup, bool erase, int *thspec){
        int len = tup->length;
        int pm = (int)(tup->thparentmassMH/dnaAA1.discretization_factor);
        int currention = 0;
        int nbions = 0;
        int series = 0;
		
        while(series < 2){
            for(int i = 1; i < len; i++){
                if(series == 0){
                    currention = bionsint[i]; //if series == 0
                }else{
                    currention = yionsint[i]; //if series == 1
                }
                
                if(currention < 18) cout << "\nScoring::calculateThSpec3: currention: " << currention;
                
                //set or erase spectrum
                if(erase == false){
                    thspec[currention - 28] += 1;
                    thspec[currention - 18] += 1;
                    thspec[currention - 17] += 1;
                    thspec[currention + 0 ] += 4;
                    thspec[currention + 1 ] += 1;
                }else{
                    thspec[currention - 28] = 0;
                    thspec[currention - 18] = 0;
                    thspec[currention - 17] = 0;
                    thspec[currention + 0 ] = 0;
                    thspec[currention + 1 ] = 0;
                }
                
            }
            series++;
        }
        
        //average negative weight
        nbions = (len - 1) * 2;
        float thspec_avgint = (float)nbions * (float)(1+1+1+4+1) / (float)pm;
        return thspec_avgint;
        
    }
	
    float Scoring::calculateScore3(Spectrum *spectrum, int *thspec, float thspec_avgint){
        float score = 0;
        int scorepos = 0;
        float scoreneg = 0;
        
        for(int i = 0; i < spectrum->lengthB; i++){
            scorepos += thspec[spectrum->binsB[i]];
        }
        
        scoreneg = thspec_avgint * spectrum->lengthB;
        score = scorepos - scoreneg;
        return score;
    }
    
    
    //*************************************************************************************
    
    void Scoring::calculateThSpec4(Tuple *tup, bool erase, int *thspec){
        
        int len = tup->length;
        
        if(erase == false){
            for(int i = 1; i < len; i++){
                thspec[bionsint[i]] += 1;
                thspec[yionsint[i]] += 1;
            }
        }else{
            for(int i = 1; i < len; i++){
                thspec[bionsint[i]] = 0;
                thspec[yionsint[i]] = 0;
            }
        }		
    }
	
    float Scoring::calculateScore4(Tuple *tup, Spectrum *spectrum, int *thspec){
        
        float score = 0;
        //int N = spectrum->binsB[spectrum->lengthB - 1] - spectrum->binsB[0]; //area covered by peaks (in experimental spectrum)
        int N = spectrum->parentmassMHB; //area covered by peaks (in experimental spectrum)
        int K = spectrum->lengthB; //number of peaks in experimental spectrum
        int n = (tup->length - 1) * 2; //number of peaks in theoretical spectrum
        int spc = 0; //shared peaks getTotal
        
        //sharedpeakcount only
        if(se1->scorepval == 0){
            for(int i = 0; i < spectrum->lengthB; i++){
                spc += thspec[spectrum->binsB[i]];
            }		
            score = spc;
            
            //ln pdf
        }else if(se1->scorepval == 1){
            for(int i = 0; i < spectrum->lengthB; i++){
                spc += thspec[spectrum->binsB[i]];
            }		
            score = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc);
            
            //ln p-value
        }else if(se1->scorepval == 2){
            for(int i = 0; i < spectrum->lengthB; i++){
                spc += thspec[spectrum->binsB[i]];
            }		
            score = -(hypergeometric1->log10HypergeometricPvalue(N, K, n, spc));
            
            //ln pdf continuity-corrected
        }else if(se1->scorepval == 3){
            for(int i = 0; i < spectrum->lengthB; i++){
                spc += thspec[spectrum->binsB[i]];
            }		
            float y1 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc);
            float y2 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc + 1);
            if(y2 < y1){
                score = 0;
            }else{
                score = y1;
            }
            
            //score highest peaks initially, continue if promising
        }else if(se1->scorepval == 4){	
            int firstpeaks = (int)(N * se1->firstpeaksperhundred / (double)100);
            float neglogpdf = 0;
            if(spectrum->lengthB < firstpeaks) firstpeaks = spectrum->lengthB;
            
            //score highest peaks only
            for(int i = 0; i < firstpeaks; i++){
                //cout << "\n" << spectrum->binsB[i]; 
                spc += thspec[spectrum->binsB[i]];
            }		
            neglogpdf = -hypergeometric1->log10HypergeometricPDFfloat(N, firstpeaks, n, spc);
            
            //score remaining peaks
            if(neglogpdf > se1->firstpdfcutoff){
                se1->postscored++;
                for(int i = firstpeaks; i < spectrum->lengthB; i++){
                    spc += thspec[spectrum->binsB[i]];
                }		
                score = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc);
            }else{
                se1->prescored++;
                score = neglogpdf;
            }
            
            //developer's version, output details	
        }else if(se1->scorepval == 9){
            for(int i = 0; i < spectrum->lengthB; i++){
                spc += thspec[spectrum->binsB[i]];
            }		
            
            cout << "\n" << N << "\t" << K << "\t" << n << "\t" << spc << "\t";
            cout << "\t" << hypergeometric1->log10HypergeometricPvalue(N, K, n, spc);
            cout << "\t" << hypergeometric1->log10HypergeometricPDF(N, K, n, spc);
            cout << "\t" << hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc);
            
            float y1 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc);
            float y2 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, spc + 1);
            if(y2 < y1){
                cout << "\t" << 0;
            }else{
                cout << "\t" << -y1;
            }
        }
        
        return score;
    }
    
    
    // hypergeometric score Sadygov-Yates 2003, including neutral losses
    void Scoring::calculateThSpec5(Tuple *tup, bool erase, int *thspec){
        
        int len = tup->length;
        int lossA = 17;
        int lossB = 18;
        
        if(erase == false){
            for(int i = 1; i < len; i++){
                thspec[bionsint[i]-lossA] += 1;
                thspec[yionsint[i]-lossA] += 1;
                thspec[bionsint[i]-lossB] += 1;
                thspec[yionsint[i]-lossB] += 1;
            }
        }else{
            for(int i = 1; i < len; i++){
                thspec[bionsint[i]-lossA] = 0;
                thspec[yionsint[i]-lossA] = 0;
                thspec[bionsint[i]-lossB] = 0;
                thspec[yionsint[i]-lossB] = 0;
            }
        }
    }
	
    float Scoring::calculateScore5(Tuple *tup, Spectrum *spectrum, int *thspec){
        
        float score = 0;
        //int N = spectrum->binsB[spectrum->lengthB - 1] - spectrum->binsB[0]; //area covered by peaks (in experimental spectrum)
        int N = spectrum->parentmassMHB; //area covered by peaks (in experimental spectrum)
        int K = spectrum->lengthB; //number of peaks in experimental spectrum
        int n = (tup->length - 1) * 4; //number of peaks in theoretical spectrum
        int k = 0; //shared peaks getTotal
        
        for(int i = 0; i < spectrum->lengthB; i++){
            k += thspec[spectrum->binsB[i]];
        }
        
        if(se1->scorepval == 0){
            score = k;
        }else if(se1->scorepval == 1){
            score = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k);
        }else if(se1->scorepval == 2){
            score = -(hypergeometric1->log10HypergeometricPvalue(N, K, n, k));
        }else if(se1->scorepval == 3){
            float y1 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k);
            float y2 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k + 1);
            if(y2 < y1){
                score = 0;
            }else{
                score = y1;
            }
        }else if(se1->scorepval == 9){
            cout << "\n" << N << "\t" << K << "\t" << n << "\t" << k << "\t";
            cout << "\t" << hypergeometric1->log10HypergeometricPvalue(N, K, n, k);
            cout << "\t" << hypergeometric1->log10HypergeometricPDF(N, K, n, k);
            cout << "\t" << hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k);
            
            float y1 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k);
            float y2 = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n, k + 1);
            if(y2 < y1){
                cout << "\t" << 0;
            }else{
                cout << "\t" << -y1;
            }
        }
        
        return score;
    }
    
    // hypergeometric score Sadygov-Yates 2003, including many neutral losses
    void Scoring::calculateThSpec6(Tuple *tup, bool erase, int *thspec){
        
        int len = tup->length;
        int lossA = 17;
        int lossB = 18;
        
        if(erase == false){
            for(int i = 1; i < len; i++){
                for(int j = 0; j < 61; j++){
                    thspec_prefix[bionsint[i] + j - 50] = j;
                    thspec_suffix[yionsint[i] + j - 50] = j;
                }
            }
        }else{
            for(int i = 1; i < len; i++){
                for(int j = 0; j < 61; j++){
                    thspec_prefix[bionsint[i] + j - 50] = 61;
                    thspec_suffix[yionsint[i] + j - 50] = 61;
                }
            }
        }
    }
	
    float Scoring::calculateScore6(Tuple *tup, Spectrum *spectrum, int *thspec){
        
        float score = 0;
        int N = spectrum->parentmassMHB; //area covered by peaks (in experimental spectrum)
        int K = spectrum->lengthB; //number of peaks in experimental spectrum
        int n = (tup->length - 1) * 2; //number of peaks in theoretical spectrum
        int spc = 0; //shared peaks getTotal
        
        int firstpeaks = (int)(N * se1->firstpeaksperhundred / (double)100);
        float neglogpdf = 0;
        if(spectrum->lengthB < firstpeaks) firstpeaks = spectrum->lengthB;
		
        //score highest peaks only
        for(int i = 0; i < firstpeaks; i++){
            spc += thspec[spectrum->binsB[i]];
        }		
        neglogpdf = -hypergeometric1->log10HypergeometricPDFfloat(N, firstpeaks, n, spc);
        
        if(neglogpdf > se1->firstpdfcutoff){
            
            int prefixspcs[62];
            int suffixspcs[62];
            
            //precalculate pdfs
            float neglogpdfs[20];
            int i0 = 20;
            if(i0 > n/2) i0 = n/2;
            for(int i = i0; i >= 0; i--){
                neglogpdfs[i] = -hypergeometric1->log10HypergeometricPDFfloat(N, K, n/2, i);
                //cout << "\ni:" << i << " " << neglogpdfs[i];
                if(i < i0 && (neglogpdfs[i] > neglogpdfs[i+1] || neglogpdfs[i] < 0)){
                    neglogpdfs[i] = 0;
                    //cout << " neglogpdfs2: " << neglogpdfs[i];
                }
            }
            
            float prefixpdfs[62];
            float suffixpdfs[62];
            for(int i = 0; i < 62; i++){
                prefixspcs[i] = 0;
                suffixspcs[i] = 0;
                prefixpdfs[i] = 0;
                suffixpdfs[i] = 0;
            }
            
            for(int i = 0; i < spectrum->lengthB; i++){
                //cout << "\t" << spectrum->binsB[i] << ":" << thspec_prefix[spectrum->binsB[i]];
                prefixspcs[thspec_prefix[spectrum->binsB[i]]]++;
                suffixspcs[thspec_suffix[spectrum->binsB[i]]]++;
            }	
            for(int i = 0; i < 61; i++){
                if(prefixspcs[i] >= i0){
                    prefixpdfs[i] += -hypergeometric1->log10HypergeometricPDFfloat(N, K, n/2, prefixspcs[i]);
                }else{
                    prefixpdfs[i] += neglogpdfs[prefixspcs[i]];
                }
                if(suffixspcs[i] >= i0){
                    suffixpdfs[i] += -hypergeometric1->log10HypergeometricPDFfloat(N, K, n/2, suffixspcs[i]);				
                }else{
                    suffixpdfs[i] += neglogpdfs[suffixspcs[i]];
                }
            }
            vector<float> v;
            for(int i = 0; i < 61; i++){
                if(i == 50 || i == 32 || i == 33){
                    v.push_back(prefixpdfs[i]);
                    v.push_back(suffixpdfs[i]);
                }
            }
            
            sort(v.begin(), v.end());
            cout << "\n" << v.size() << " " << neglogpdf << " " ;
            cout.precision(4);
            for(int i = v.size()-1; i > v.size()-1-se1->bestseries; i--){
                score += v[i];
                cout << "\t" << v[i];
            }
            se1->postscored++;
        }else{
            score = neglogpdf;
            se1->prescored++;
        }
        return score;
    }
    
    
    float Scoring::calculateScore7(Tuple *tup, Spectrum *spectrum, int *thspec, int *sharedpeakcount){
        
        float score = 0;
        int n = (tup->length - 1)*2; //number of peaks in theoretical spectrum in case of charge 2
        if(spectrum->chargestate>2){
        	int n = (tup->length - 1)*4; //number of peaks in theoretical spectrum (since always both charges are allowed in case of charge 3 or up 
        }
        int spc = 0; //shared peaks getTotal
        
        //score b- and y-ions
        //for chargestate 2
        if(spectrum->chargestate==2){
            for(int i = 1; i < tup->length; i++){
                spc += spectrum->binaryB[bionsint[i]]; // += is faster than if (2.60E6 vs 2.45E6 on 5471)
                spc += spectrum->binaryB[yionsint[i]]; // += is faster than if (2.60E6 vs 2.45E6 on 5471)
                //cout<<"\n spc in loop after yions "<<spc;
                //cout<<"and the spectrum->binaryB is "<<spectrum->binaryB[yionsint[i]];
                //cout<<"and yionsint is "<<yionsint[i];
                //cout<<"And i is "<<i;
                //cout<<"size of spec "<<spectrum->parentmassMHB;
                
            }
            if(se1->determineNK){
                for(int i = 1; i < tup->length; i++){
                    spectrum->arrayN[bionsint[i]]=true;
                    spectrum->arrayN[yionsint[i]]=true;
                }
            }
        }
        else{
            //and for higher charges
            for(int i = 1; i <2*tup->length+2; i++){
                //would be better if len+1 case is disregarded
                spc += spectrum->binaryB[bionsint3[i]]; // += is faster than if (2.60E6 vs 2.45E6 on 5471)
		        spc += spectrum->binaryB[yionsint3[i]]; // += is faster than if (2.60E6 vs 2.45E6 on 5471)
            }
     		if(se1->determineNK){
	        	for(int i = 1; i < 2*tup->length+2; i++){
		        	spectrum->arrayN[bionsint3[i]]=true;
                    spectrum->arrayN[yionsint3[i]]=true;
                }
            }
        }	
        //changes BYR output spc
        //cout<<"\n k is currently "<<spc;
        //for( int i=1;i<2*tup->length+1;i++){
        //cout<<"\n bions look like this"<< spectrum->binaryB[bionsint[i]];
        //cout<<"\n yions look like this"<<spectrum->binaryB[yionsint[i]];
        //cout<<"\n bions 3 look like this"<<bionsint3;
        //}
        score = -hypergeometric1->log10HypergeometricPDFfloat(spectrum->N, spectrum->K, n, spc);
        
        //test whether bin is exceptionally bad instead of exceptionally good, but not for low scores
        //ensures continuous scoring function
        if(score > se1->firstpdfcutoff && spc < n){
            if(score > -hypergeometric1->log10HypergeometricPDFfloat(spectrum->N, spectrum->K, n, spc+1)) score = 0;
        }
		
        //return result values
        (*sharedpeakcount) = spc;
        return score;
    }
    
    
    //*************************************************************************************
    
    
    void Scoring::checkParentMass(Tuple *tup){
        double diff = tup->thparentmassMH - dnaAA1.getParentMassMH(tup->aasequence, tup->length);
        if(diff > 10 || diff < -10){
            cout << "\nScoring::scoreTuple: parent mass is not ok! PMfromseq: checkParentMass" << dnaAA1.getParentMassMH(tup->aasequence, tup->length) << " thPM tuple: " << tup->thparentmassMH << "  ps pe ss se " << tup->prefixstart << " " << tup->prefixend << " " << tup->suffixstart << " " << tup->suffixend << " rev:" << tup->chromosome_reverse;	
            tup->getNTSeqFromChromosome();
            cout << " ntseq:" << tup->getNTSeq() << " aaseq:" << tup->getAASeq() << " diff:" << diff;
        }
    }
    
}
