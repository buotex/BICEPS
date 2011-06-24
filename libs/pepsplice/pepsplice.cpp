#include "pepsplice.h"
#include "bicepsdefinitions.h"
namespace Pepsplice{
    
    extern Pool<Pepsplice::Tuple, POOLSIZE> __POOL__;    
    
    Pool<Pepsplice::Tuple, POOLSIZE> __POOL__ = Pool<Pepsplice::Tuple, POOLSIZE>();

    void pepsplice_func(int argc, vector<string> & argv, float penalty_mutation_, std::vector<float> & max_penalties, std::vector<PepspliceResult>& results, const std::vector<std::tuple<unsigned int, std::string, std::string> >& currentfasta){
        try{
            __POOL__.reinit();
            vector<string*> specfiles;	
            vector<string> params;	
            
            //parse arguments (parameters and spectrum files mixed)
            //cout << "\n\nmain.cpp: Reading " << argc - 1 << " arguments\n";
            string paramstring = "";
            for (int i = 1; i < argc; i++)
            {
                if ( (argv[i])[0] == '-'){
                    
                    params.push_back(argv[i]);
                    paramstring += argv[i]; //for output file names later on
                    
                }
                else {
                    
                    specfiles.push_back(new string(argv[i]));
                }
            }
            
            //define THE parameter and service class, available to all classes
            //removed some Clutter, BX
            Services *se1 = new Services();
            //cout<<se1->masstol_belowDa;
            //cout<<"\nposition nach Services";
            se1->parseParameters(params); //vector of parameters
            se1->penalty_mutation = penalty_mutation_;
            if (se1->outputlevel > 2) {
                se1->os.open("__pepspliceDebug__.txt"); 
                if (!se1->os.is_open()) cerr << "Couldn't write debug output to __pepspliceDebug.txt__, please choose a lower debug level" << endl;
            } 
            se1->max_penalties = max_penalties;
            se1->omax_penalty = max_penalties.back();
                        
            //cout<<se1->masstol_belowDa;
            //cout<<"\nposition nach parseParameters";
            //cout<<se1->masstol_below;
            
            //define unique objects
            DnaAA *dnaAA1 = new DnaAA();
            if (se1->outputlevel > 3) dnaAA1->checkInitialization();
            Spectra *spectra1 = new Spectra(se1, dnaAA1);
            SpectrumParser *spectrumparser1 = new SpectrumParser(spectra1, se1);
            Tuples *tuples1 = new Tuples(se1, dnaAA1);
            ProteinParser *proteinparser1 = new ProteinParser(se1, tuples1);
            Chromosomes *chromosomes1 = new Chromosomes(se1);
            SlidingWindow *slidingwindow1 = new SlidingWindow(se1);
            Scoring *scoring1 = new Scoring(se1);
            
            spectra1->currentfasta = &currentfasta;
            
            
            
            
            //connect unique objects, define flow
            chromosomes1->setSlidingWindow(slidingwindow1);
            slidingwindow1->setTuples(tuples1);
            tuples1->setScoring(scoring1);
            scoring1->setSpectra(spectra1);
            
            if(se1->hotspectra == true) spectrumparser1->hotspectra1->parseHotSpectrumFiles("in_hotspectra.param");
            
            //load spectra, learn spectrum intensities, discard spectra
            se1->learning_spectrum_intensities = true;
            spectrumparser1->parseFiles(specfiles);
            se1->invertLearnedPeakDistribution();
            
            //load spectra, preprocess them, keep them
            se1->learning_spectrum_intensities = false;
            se1->spectra = 0;
            spectrumparser1->parseFiles(specfiles);
            
            //add number of spectra to output filename
            paramstring += "_" + se1->intToString(se1->spectra);	
            se1->outfileparams = paramstring;
            
            //last preparations //Changes by BX as these results aren't needed hopefully
            //se1->peakdistribution_rel->writeDistribution(se1->getOutFileName("peakdistribution_rel")); //needs parameter string
            //se1->peakdistribution_abs->writeDistribution(se1->getOutFileName("peakdistribution_abs")); //needs parameter string
            //se1->parentmassdistribution->writeDistribution(se1->getOutFileName("parentmassdistribution")); //needs parameter string
            spectra1->sortAndPrepareSpectra(); //needs parameter string
            
            //slidingwindow1->showParameters(); //show all sliding window parameters //Changed by BX
            
            se1->globaltime->reset();
            
            //PARSE FASTA FILES AND INSTRUCTIONS
            
            //chromprotfile.open("in_fastafiles.param");
            
            //start changes BYR
            //cout << "\n parsedline:" << parsedline;
            
            //	int nbresults = 0;
            
            //bool tag = false;
            //	bool chromosomes = false;
            //	bool hotspots = false;
            //	bool proteins = false;
            //bool write = false;
            
            
            //changes BX 
            //beginning customized pepsplice run
            
            
            se1->penalty_max = max_penalties[0];
            
            // cout << "\nProcessing protein file " << "etsresults.fasta" << flush;
            
            // proteinparser1->parseFASTA("etsresults.fasta");
            proteinparser1->parseFASTA(currentfasta);
            
            //The parseFASTA method reads the Sequences out of the fasta file.
            //Afterwards, "finishProtein" is called which adds the tuples and their mutations while scoring them
            
            tuples1->forwardTuples();
            spectra1->oldSequences = proteinparser1->oldSequences; 

            
            vector<PepspliceResult> samePenaltyResults;
            PepspliceResult dummy;
            samePenaltyResults.push_back(dummy);
            
            spectra1->returnBIC(samePenaltyResults);//spectra has a bestmatch, this is used to calculate the bic to put in res.
            //results.push_back(res);
            stable_sort(samePenaltyResults.begin(), samePenaltyResults.end(), PepspliceResultComparator());
            if (se1->outputlevel > 2)
            {
                for (size_t i = 0; i < samePenaltyResults.size(); i++)
                    se1->os << samePenaltyResults[i];
            }
            
            //Make the tuples ready for writing to a file.
            //catch report here for bic?
            
            //cout << "pass done with penalty_max=" << se1->penalty_max << endl;
            //cout << "Result is \n" << res << endl;
            PepspliceResult & res = samePenaltyResults.front();
            results.push_back(res);
            
            if (se1->outputlevel > 1)
            {
                cout << "pass done with penalty_max=" << se1->penalty_max << endl;
            }
            
            if (se1->outputlevel > 1)
            {
                cout <<'\n' << res << '\n';
            }
            
            
            for (unsigned int i = 1; i < max_penalties.size();++i)
            {
                if (res.bic > 1 - (se1->penalty_max+se1->penalty_mutation+0.3)/2 * log((float)res.n/2)) break;
                
                
                se1->penalty_max = max_penalties[i];
                spectra1->spectra[0]->bestmatches1->reset();
                spectra1->spectra[1]->bestmatches1->reset();
                proteinparser1->doAnotherRun(currentfasta);
                tuples1->forwardTuples();
                spectra1->oldSequences = proteinparser1->oldSequences;
                if (se1->outputlevel > 1)
                {
                    cout << "pass done with penalty_max=" << se1->penalty_max << endl;
                }  
                spectra1->returnBIC(samePenaltyResults);//spectra has a bestmatch, this is used to calculate the bic to put in res.
                
                stable_sort(samePenaltyResults.begin(), samePenaltyResults.end(), PepspliceResultComparator());
                if (se1->outputlevel > 2)
                {
                    for (size_t i = 0; i < samePenaltyResults.size(); i++)
                        se1->os << samePenaltyResults[i];
                }
                PepspliceResult & res = samePenaltyResults.front();
                results.push_back(res);
                
                if (se1->outputlevel > 1)
                {
                    cout << res << '\n';
                }
            }
            
            stable_sort(results.begin(), results.end(), PepspliceResultComparator());
            
            
            //se1->doProgressReport();
            //end Changes
            
            //  stable_sort(results.begin(), results.end(), PepspliceResultComparator());
            
            //results[0].seqIds.assign (currentfasta.find(results[0].OrigSequence)->second.begin(),currentfasta.find(results[0].OrigSequence)->second.end()) ;
            
            for (size_t i = 0; i < specfiles.size(); ++i)
            {
                delete specfiles[i];
            }
            
            
            delete se1;
            delete tuples1;
            delete proteinparser1;
            delete chromosomes1;
            delete slidingwindow1;
            delete spectra1;
            delete spectrumparser1;
            delete scoring1;
            delete dnaAA1;
            //return fastaoutput;
            //end changes BX
        }
        catch(const char* c)
        {
            std::cout << c << endl;
        }
        
        //return 0;

    }
	//cleaning up
	
    
    //while(getline(chromprotfile, parsedline)){
    //		//		//cout << "\nparsedline: " << parsedline;
    //		if(parsedline.size() > 0 && parsedline[0] != '#'){
    //			if(parsedline[0] == '<'){
    //				
    //				//write tag
    //				if(parsedline.substr(0, 15) == "<WRITE RESULTS>"){
    //					tuples1->forwardTuples();
    //					
    //					se1->doProgressReport();
    //					nbresults++;
    //					if(se1->outputlevel > 2) cout << "\n\nmain.cpp Writing result files of round " << nbresults << ".\n\n" << flush;				
    //					se1->outfilenumber = se1->intToString(nbresults);
    //					spectra1->writeResults();
    //					//spectra1->getValues();
    //					cout << "\n";	
    //					
    //					//chromosome tag
    //				}else if(parsedline.substr(0, 19) == "<CHROMOSOMES START>"){
    //					chromosomes = true;
    //					hotspots = false;
    //					proteins = false;
    //					
    //					//hotspots tag
    //				}else if(parsedline.substr(0, 16) == "<HOTSPOTS START>"){
    //					chromosomes = false;
    //					hotspots = true;
    //					proteins = false;
    //					
    //					//proteins tag
    //				}else if(parsedline.substr(0, 16) == "<PROTEINS START>"){
    //					chromosomes = false;
    //					hotspots = false;
    //					proteins = true;
    //				}
    //				
    //			}else{
    //				//IS NO TAG
    //				
    //				//begin changes BXu
    //				if(proteins == true){
    //					
    //					se1->penalty_max = max_penalties[0];
    //					
    //					cout << "\nProcessing protein file " << parsedline << flush;
    //					
    //					proteinparser1->parseFASTA(,false);
    ////					
    ////					
    ////					fastaoutput.push_back(proteinparser1->sequences); 
    ////					fastaoutput.push_back(proteinparser1->sequencesik);
    ////					fastaoutput.push_back(proteinparser1->ids);
    //					
    //					
    //					//The parseFASTA method reads the Sequences out of the fasta file.
    //					//Afterwards, "finishProtein" is called which adds the tuples and their mutations while scoring them
    //					
    //					tuples1->forwardTuples();
    //					
    //					PepspliceResult res;
    //					spectra1->returnBIC(res);
    //					results.push_back(res);
    //					
    //					
    //					//Make the tuples ready for writing to a file.
    //					//catch report here for bic?
    //					
    //					cout << "pass done with penalty_max=" << se1->penalty_max << endl;
    //	
    //					
    //					
    //					
    //					for (unsigned int i = 1; i < max_penalties.size();++i)
    //					{
    //						se1->penalty_max = max_penalties[i];
    //						proteinparser1->doAnotherRun();
    //						tuples1->forwardTuples();
    //						cout << "pass done with penalty_max=" << se1->penalty_max << endl;
    //						spectra1->returnBIC(res);
    //						results.push_back(res);
    //					}
    //					//se1->doProgressReport();
    //					//end Changes
    //					stable_sort(results.begin(), results.end(), PepspliceResultComparator());
    //					
    //					
    //				}
    //				
    //				else if(chromosomes == true){
    //					cout << "\nProcessing chromosome file " << parsedline << flush;
    //					chromosomes1->parseAndFeedChromosome(parsedline); //c_str converts a cpp-string to a c-string
    //					tuples1->forwardTuples();
    //					se1->doProgressReport();
    //				}else if(hotspots == true){
    //					cout << "\nProcessing hot spots file " << parsedline << flush;
    //					chromosomes1->hotspots1->parseHotSpots(parsedline);
    //				}				
    //			}
    //			
    //		} //line without #
    //	} //end while
    //	
    //	delete se1;
    //	delete tuples1;
    //	delete proteinparser1;
    //	delete chromosomes1;
    //	delete slidingwindow1;
    //	delete spectra1;
    //	delete spectrumparser1;
    //	delete scoring1;
    //	
    //	//return fastaoutput;
    //	//end changes BX
    //	return 0;
    //}
    //
    
    
    
    
}




