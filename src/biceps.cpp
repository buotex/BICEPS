#include <tuple>
#include "biceps.h"

const double PI = 3.14159265;



void Biceps::initialize(int argc, char* argv[])
{    
  mgfId = 0;
  int status = parseProgramOptions(
                                   argc, argv, 
                                   this->directOp, this->pepnOp, this->pepsOp, 
                                   this->genOp);
  
  if (status >0){
    return;
  }
  else
  {
    
    
    
    showLicenses();
    //initializing fasta file
    fasta.debuglevel = genOp.debug;
    fasta.loadBigFasta(genOp.fastaname);
    loadAAModifications();
    //going through the mgf
    readMGF();        
  }
}

void Biceps::showLicenses()
{
#ifdef HAVE_PEPNOVO
  if (genOp.tool == PEPNOVO || genOp.tool == PEPNDIREC)
  {
    std::cout  << "\n"<< "PepNovo+, Build " << build_name << "\n"
    << "Copyright 2008, The Regents of the University of California. All Rights Reserved." << "\n"
    << "Created by Ari Frank (arf@cs.ucsd.edu)" << endl;
    
  }
#endif
#ifdef HAVE_DIRECTAG
  if (genOp.tool == DIRECTAG || genOp.tool == PEPNDIREC)
  {
    std::cout << "\n"<< "DirecTag " << DIRECTAG_VERSION_STRING << " (" << DIRECTAG_BUILD_DATE << ")\n" 
    << DIRECTAG_LICENSE << endl; 
  }
  
#endif

    std::cout << "\n" << "Pepsplice " << "Copyright (c) <2007>, <Franz Roos, ETH Zurich> All rights reserved." << endl;


}



bool Biceps::isDosFile(std::string & filename)
{
  std::ifstream file;
  file.open(filename.c_str(),ios::binary);
  std::string linebuffer;
    
    //check  line ending via char comparison  
  if (file.is_open()){
    try{
      while(file.good()){
        getline(file, linebuffer, '\n');
        if (linebuffer.size() > 0) break;
      }
      if(!file.good())
      {
        throw runtime_error("mgf-File not useable.");
      }
      file.close();
    }
    catch(...)
    {
      std::cerr << "Error while reading mgf-File" << std::endl;
      file.close();
    }
    
  }
  else //(file not open)
  {
    throw runtime_error("mgf-File could not be read");
  }
  
  return (linebuffer[linebuffer.size()-1] == '\r');
  
}




void Biceps::readMGF()
{
    std::string filename = genOp.mgfname + string(".mgf");

    std::cout <<"Now reading " << filename << std::endl;
    this->fasta.setNumTags(genOp.numTags);

    std::string linebuffer; //< buffering the textlines in the mgf file
    bool validSpectrum = true; //<So we can just skip if it's invalid.

    //Now checking for the filetype.
    bool dosformat = isDosFile(filename);

    std::ifstream mgfFile;
    mgfFile.open(filename.c_str(), ios::binary); //<The Filename of the spectraFile


    bool inBeginIons = 0; bool inPeakList = 0; //status, where are we in the spectrum?

    //Only for buffer reasons, what's the charge etc. of the current spectrum?
    int charge = 0; 
    float precursormass = 0;
    string title="";

    std::vector<int> indices; //to remember the corresponding spectrum-indices of valid results
    std::vector<string> titles; //remember the spectrum titles
    size_t numGoodscores = 0;

    ofstream bufferfile;
    ofstream results;
    string tempresultname = std::string("Biceps_tempResults") + genOp.mgfname + string(".txt");
    results.open(tempresultname.c_str(), ios::binary); //Biceps_results will handle the final output.

    while(std::getline(mgfFile, linebuffer, '\n') )
    {  
       
       //first case, we're at the beginning of a new spectrum-subsection
        if(linebuffer.substr(0, 8) == "BEGIN IO"){
        
            mgfId++;
            if (inBeginIons)
            {
                //throw runtime_error(("BEGIN IONS tag found without previous BEGIN IONS being closed at" + lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + "\n")); 
                cout << "BEGIN IONS tag found without previous BEGIN IONS being closed at" + lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + "\n" << std::endl;
                cout << "will begin new spectrum here" << std::endl;
                bufferfile.close();
            }
            validSpectrum = true;
            inBeginIons = true;
            if (bufferfile.is_open()) bufferfile.close();
            bufferfile.open(genOp.sSpectrumFN.c_str(), fstream::trunc | ios::binary);
            if (!bufferfile.good()) { cerr << "Problem writing to buffer, check your writing rights please." << std::endl;} 

            if (dosformat)
            {
                bufferfile << linebuffer.substr(0,linebuffer.size()-1).c_str() << "\n";
            }
            else
            {
                bufferfile << linebuffer.c_str() << "\n";
            }
            //continue;
        } //if BEGIN IO - Beginning of a spectrum 
       
       
       //second case
        else if (linebuffer.substr(0, 8) != ("END IONS")) //we are in a spectrum, parse the lines.
        {
            if (!inBeginIons) continue;

            if (dosformat)
            {
                bufferfile << linebuffer.substr(0,linebuffer.size()-1).c_str() << "\n";
            }

            else
            {
                bufferfile << linebuffer.c_str() << "\n";
            }

            try
            {
                if (!inPeakList)
                {

                    if (linebuffer.substr(0, 6) == ("TITLE="))
                    {
                        // if a title is found, use it as the id instead of the index
                        //spectrum.id = lineStr.substr(6);
                        title = linebuffer.substr(6);
                        titles.push_back(title);
                    }
                    else if (linebuffer.substr(0, 8) == ("PEPMASS="))
                    {

                        string pepMassStr = linebuffer.substr(8);
                        //bal::trim(pepMassStr);
                        precursormass= lexical_cast<float>(pepMassStr);
                        //						selectedIon.set(MS_m_z, mz);

                    }
                    else if (linebuffer.substr(0, 7) == ("CHARGE="))
                    {
                        string pepChargeStr = linebuffer.substr(7);
                        //bal::trim_if(pepChargeStr, bal::is_any_of("+- \t\r"));
                        size_t rest = pepChargeStr.find_first_not_of("0123456789",0);
                        if (rest != string::npos)
                            pepChargeStr.erase(rest);
                        charge = lexical_cast<size_t>(pepChargeStr);
                        if (charge > 3) 
                        {
                            std::cout << "chargevalue " << charge << " is unhandled, skipping." << std::endl;                          
                            validSpectrum = false;
                        }
                        //						selectedIon.set(MS_charge_state, charge);

                    }
                    else if(linebuffer.find('=') != string::npos)
                    {
                        continue; // ignored attribute
                    }
                    else
                    {
                        if (inBeginIons)
                            inPeakList = true;
                    }
                } //if (!inPeakList)
            } // if try
            catch(const std::exception& e)
            {
                //throw runtime_error(("[SpectrumList_MGF::parseSpectrum] Error parsing line at offset " +
                //                   lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + ": " + linebuffer + "\n"));
                cout << "[SpectrumList_MGF::parseSpectrum] Error parsing line at offset " +
                    lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + ": " + linebuffer + "\n" << std::endl;
                inBeginIons = false;
                validSpectrum = false;
            }
            //continue;
        }	//else if (not end of spectrum)




        else //END IONS
        {

            //We found the last line of a spectrum, now it's time to finish the mgf buffer file and
            //do pepnovo/directag and pepsplice afterwards.

            if (!inBeginIons)
            { 
                //incorrect mgf file, please fix.

                //throw runtime_error(("END IONS tag found without opening BEGIN IONS at" + lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + "\n")); //
                std::cout << "END IONS tag found without opening BEGIN IONS at" + lexical_cast<string>(size_t(mgfFile.tellg())-linebuffer.length()-1) + "\n" << std::endl;
                validSpectrum = false;
            }
            inBeginIons = false;
            inPeakList = false;
            if (dosformat)
            {
                bufferfile << linebuffer.substr(0,linebuffer.size()-1).c_str() << "\n";
            }
            else
            {
                bufferfile << linebuffer.c_str() << "\n";
            }
            bufferfile.close();


            if (validSpectrum == false) //either charge too high or spectrum invalid
            {
                continue; //just go to the next line
            }

            fasta.initializeMGF(charge, precursormass);

            cout << "Analyzing spectrum " << mgfId << endl;
            
            //Start analyzing by calling run_programs, which will run directag/pepnovo and then pepsplice afterwards
            bool success = run_programs(); //if an actual result was found, return true, else return false;
            if (success){ 
                Pepsplice::PepspliceResult & res = pepResults.back();
                indices.push_back(mgfId);
                titles.push_back(title);
                //if it's good, write the result to a file
                writeResult(results, res, mgfId, title);
                //pepres.pop_back();
                ++numGoodscores;
            }
        } //END IONS, parsed one spectrum, call other programs //end else
    }	//while File.good







    assert(pepResults.size() == numGoodscores);//check if all the results are in pepres

    ofstream resfile;
    string resfileName = string("Biceps_gmm") + genOp.mgfname + string(".txt");
    resfile.open(resfileName.c_str(), ios::trunc | ios::binary);

    //temporary score-output to call bic on.
    for(size_t i = 0; i < pepResults.size(); ++i)
    {
        resfile << pepResults[i].score << "\n";
    }
    resfile << std::flush;
    resfile.close();

    //now beginning BIC part
    std::vector<double> mu,sigma;
    std::vector<int> labels;
    double cutoff;

    //final output
    if (pepResults.size() > 10){
        gmm_bic(2,numGoodscores,resfileName.c_str(), mu, sigma, labels);
        cutoff = findCutoff(mu, sigma, 2); //using 0.05 as cutoff

        Biceps::writeCompleteResult(pepResults, indices, titles, labels, mu, sigma, cutoff);
        Biceps::writeFasta(pepResults);
    }
}


void Biceps::writeCompleteResult
(const std::vector<PepspliceResult> & pepResults, 
 const std::vector<int> & indices, 
 const std::vector<string> & titles, 
 const std::vector<int> & labels, 
 const std::vector<double> & mu, 
 const std::vector<double> & sigma, 
 const double cutoff) const
{
    std::ofstream finalOutput;
    std::string finalOutputName = string("Biceps_Results") + genOp.mgfname+ string(".txt");
    finalOutput.open("Biceps_Results.txt", ios::trunc);
    if (!finalOutput) throw runtime_error(finalOutputName + string(" can't be written, skipping."));
    for (size_t i = 0; i < pepResults.size(); i++)
    {
        if (pepResults[i].score > cutoff)
        {
            writeResult(finalOutput, pepResults[i], indices[i],titles[i]);
            finalOutput << "Label: " << labels[i] << "\n";
            finalOutput << "Confidence: " << returnConfidence(pepResults[i].score, mu[labels[i]-1], sigma[labels[i]-1]) << std::endl;
        }
    }
    finalOutput.close();

}


void Biceps::writeFasta(std::vector<PepspliceResult> & pepResults)
{
    //the parsing/changing isn't commutative, so the order is important - high bic is more important here.
    std::sort(pepResults.begin(),pepResults.end(), PepspliceResultComparator());
    
    std::ofstream fastaoutput;
    string fastaoutputname = string("Biceps_Fasta") + genOp.mgfname + string(".fasta"); 
    fastaoutput.open(fastaoutputname.c_str(), ios::trunc);
    if (!fastaoutput.is_open())
    {
        throw runtime_error(fastaoutputname + " can't be written, skipping fasta output");
    }


    // pick a sequence of the original fasta, check if the pepsplice results can be found, change accordingly
    // write correct output to a file afterwards.
    for (size_t i = 0; i < fasta.bigfasta.size(); ++i)
    {
        string buffer = fasta.bigfasta[i];
        for (size_t j = 0; j < pepResults.size(); ++j)
        {
            size_t pos = 0;
            for(;;)
            {
                pos = buffer.find(pepResults[j].OrigSequence,pos);
                if (pos != std::string::npos)
                    buffer.replace(pos, pepResults[j].OrigSequence.size(), pepResults[j].Sequence);

                else {
                    break;
                }
                ++pos;
            }    
        }
        fastaoutput << fasta.bigfasta_descriptions[i] << "\n";
        size_t limiter = 0;
        while (limiter <buffer.size()) 
        {
            fastaoutput << buffer.substr(limiter, 60) << "\n";
            limiter +=60;
        }

    }

    fastaoutput.close();
}



//This part will call the libraries.
bool Biceps::run_programs(){

    string pepnovoTags; //save the resulting tags there.
    string directagTags;
    pepnovoTags.reserve(genOp.numTags*5 + 1); //hardcoded TagLength = 5 for speedup in fastasearch later
    directagTags.reserve(genOp.numTags*5 + 1); //hardcoded TagLength = 5 for speedup in fastasearch later

    //not only the tags but also their occurence place has to be saved.
    vector<float> plowPeakMzs;
    vector<float> dlowPeakMzs;

    size_t tagcount = 0; //we will save the resulting number of tags here. 


    //the following parts will only be executed if they were activated at the CMake stage. 
#ifdef HAVE_PEPNOVO
    try{
        if (genOp.tool == PEPNOVO || genOp.tool == PEPNDIREC){
            pepnovoFunc(pepnOp.size(),pepnOp, pepnovoTags, plowPeakMzs);
        }
    }
    catch(string & e)
    {
        pepnovoTags.clear();        
        cout << e << endl;

    }

    catch(...)
    {
        pepnovoTags.clear();
    }
#endif    

#ifdef HAVE_DIRECTAG

    try{
        if (genOp.tool == DIRECTAG || genOp.tool == PEPNDIREC){
            std::string cachename = std::string("Biceps_DirectagCache") + genOp.mgfname + std::string(".cache");
            directagFunc(directOp.size(),directOp, directagTags, dlowPeakMzs, cachename); //tags written
        }
    }
    catch(const std::exception & e)
    {
        cerr << "Directag: " << e.what() << "\n";
        directagTags.clear();
    }


    catch(...) //if either of them crashed, try next.
    {
        directagTags.clear();
    }

#endif


    if (directagTags.size() ==0 && pepnovoTags.size() == 0) {
        std::cout << "No tags were found, trying next spectrum" << std::endl; 
        return 0;

    }
    //void createFasta(string & tags,  deque<string> & fastaoutput, deque<size_t> & indices);
    //fasta.createFasta(tags);




    //fasta.loadTags(pepnovoTags, directagTags, plowPeakMzs, dlowPeakMzs);
    fasta.clearTags();
    //save the tags separately
    int pepN = fasta.loadTags(pepnovoTags, plowPeakMzs, fasta.pepntags);
    int dirN = fasta.loadTags(directagTags, dlowPeakMzs, fasta.directags);

    tagcount = fasta.getNumTags();


    if (genOp.debug>0){
        cout << "The following tags were found:" << endl;
        fasta.printTags();
    }

    vector<PepspliceResult> results; //we will pass this by reference and insert good results via pepsplice.
    
    //Pepsplice will run with different penalties according to genOp
    Pepsplice::PepspliceResult bestresult =
        runPepsplice(0, 0,fasta.pepntags.size(), 0,fasta.directags.size() ); //< saving the best result.
    

    bestresult.mutation = false; //add as flag
    results.push_back(bestresult);

    //if mutations should be considered and the result isn't good enough, runPepsplice will create a modified fasta and run again.
    if (genOp.mutation == true)
    {
        if (bestresult.bic < 1-(1.098612)/2*log((float)bestresult.n/2) || bestresult.k == 0)
        {            
            bestresult = runPepsplice(1, 0,pepN, 0, dirN ); //< saving the best result.
        }
        bestresult.mutation = true;
        results.push_back(bestresult);
    }


    stable_sort(results.begin(), results.end(), Pepsplice::PepspliceResultComparator());

    //pick the best result if it's viable.(there will be a dummy in there at the least)
    if (results[0].k > 0)
    {
        pepResults.push_back(results[0]);
        return 1;
    }
    else {
        return 0;
    }

}




//Call pepsplice after actually creating a fasta (this will depend on the mutation switch)
//pepN1 etc. will tell which tags shall be considered.

PepspliceResult Biceps::runPepsplice(bool mutation, int pepN1, int pepN2,  int dirN1, int dirN2 ){
    float penalty_mutation = (mutation)? 2.0:0.0;
    std::vector<PepspliceResult> pepnresults;
    std::vector<PepspliceResult> direcresults;

    PepspliceResult dummy;
    pepnresults.push_back(dummy); //this dummy will ease some if cases.
    direcresults.push_back(dummy);

    std::vector<float> pepsplice_penalties;

    if(mutation == false)
    {
        pepsplice_penalties = this->genOp.max_penalties;
    }
    else 
    {
        pepsplice_penalties = this->genOp.max_penalties_mutated;
    }

//create fasta will use the given tags and indices to create a temporary fasta (memory only), we will pass the results to the pepsplice_func


#ifdef HAVE_PEPNOVO
    if (this->genOp.tool == PEPNOVO || this->genOp.tool == PEPNDIREC)
    {
        std::vector<std::tuple<unsigned int, std::string, std::string> > pepnFasta; //We will use this to keep the fasta created by the Fasta class.
        fasta.createFasta(mutation,pepN1,pepN2, fasta.pepntags, pepnFasta); 
        if ( fasta.getMatches() > 0) 
        {
            Pepsplice::pepsplice_func(pepsOp.size(),pepsOp,penalty_mutation, pepsplice_penalties, pepnresults, pepnFasta); //< this should fill the results and the fastaoutput vector.
            fasta.matchSequences(pepnresults[0], pepnFasta);
        }        
    }
#endif
#ifdef HAVE_DIRECTAG
    if (this->genOp.tool == DIRECTAG || this->genOp.tool == PEPNDIREC)
    {
        std::vector<std::tuple<unsigned int, std::string, std::string> > direcFasta; //We will use this to keep the fasta created by the Fasta class.
        fasta.createFasta(mutation,dirN1, dirN2, fasta.directags, direcFasta); 
        if ( fasta.getMatches() > 0) 
        {
            Pepsplice::pepsplice_func(pepsOp.size(),pepsOp,penalty_mutation, pepsplice_penalties, direcresults, direcFasta); //< this should fill the results and the fastaoutput vector.
            fasta.matchSequences(direcresults[0], direcFasta);
        }
    }
#endif

    //Some final tags for the output, take the best result and return it.
    if (pepnresults[0].bic > direcresults[0].bic)
    {
        pepnresults[0].tool = 0;
        return pepnresults[0];
    }
    else 
    {   
        direcresults[0].tool = 1;
        return direcresults[0];
    }
}


int Biceps::writeResult(std::ostream & os, const Pepsplice::PepspliceResult & res, const int specIndex, const std::string & title) const
{
    os << specIndex << "\n";
    os << "Title=" << title << "\n";
    os << "Sequence: " << convertTupleString(res.Sequence) << "\n";
    os << "OrigSequence: " << res.OrigSequence << "\n";   
    for (std::list<unsigned int>::const_iterator it = res.fastaIds.begin(); it != res.fastaIds.end(); ++it){
      os << "fastaId: " << *it << "\n";
      os << "fastaId: " << fasta.offerDescription(*it) << "\n";
    }


    os << "n: " << res.n << "\n";
    os << "k: " << res.k << "\n";
    os << "bic: " << res.bic << "\n";
    os << "penalty: " << res.penalty << "\n";
    os << "score: " << res.score << "\n";
    os << "Tool: " << res.tool << " pmax: " << res.penalty_max << " mutation: " << res.mutation << std::endl;

    return 1;
}

double Biceps::findCutoff(std::vector<double> & mu, std::vector<double> & sigma, int selector)
{
    double quantile;
    switch (selector){
        case 0: quantile = 2.326348; break;
        case 1: quantile = 2.053749; break;
        default: quantile = 1.644854;
    }
    double cutoff = 999999.0;
    for (size_t i = 0; i < mu.size(); ++i)
    {
        if (mu[i] + quantile * sigma[i] < cutoff)
            cutoff = mu[i] + quantile * sigma[i];
    }



    return cutoff;

}

string Biceps::convertTupleString(const string & sequence)const 
{
    string newSeq = "";
    for (size_t i = 0; i < sequence.size(); i++)
    {
        if ((unsigned char) sequence[i]>91) newSeq += aamod.find((unsigned char)sequence[i])->second;
        else newSeq += sequence[i]; 
    }


    return newSeq;

}

//Source: Slightly Modified Pepsplice Code
void Biceps::loadAAModifications() 
{

    ifstream inFile1;
    string aamodline;

    int aamod_i_ascii = 128;
    inFile1.open("in_AAmodifications.param", ios::binary);	
    while(getline(inFile1, aamodline)){
        if(aamodline[0] != '#'){

            //cout << "\n" << aamodline;

            //parse line field-wise		
            istringstream iss(aamodline);
            string field;
            int i_field = 0;
            unsigned char aa = 0;	
            while( getline(iss, field, ';') ){
                i_field++;
                if(i_field == 1){
                    aa = (unsigned char)field[0];
                    //assign ascii start code for current amino acid
                    aamod[aamod_i_ascii] = aa + 32;
                }else{
                    //assign alternative mass for current amino acid
                    unsigned char c = field[0];
                    string m = field.substr(1);
                    //cout << "\nc: " << c << " m: " << m;
                    //modification type: (pepNterm, )pepCterm, [protNterm, ]protCterm, _internal
                    if(c == '(' || c == ')' || c == '[' || c == ']' || c == '_'){

                        aamod[aamod_i_ascii] = aa + 32;			
                        //aamod_type[aamod_i_ascii] = c;

                        aamod_i_ascii++;

                    }else{
                        cout << "\nDnaAA.cpp line 99: please convert the modification file with dos2unix or else define what type of modification you require: " << aamodline << "\n";
                    }
                }
            }//tag or not tag
        }//not #
    }//while loop
    inFile1.close();
}


double Biceps::returnConfidence(double score, double mu, double sigma) const
{
    return 1/(sigma * std::sqrt(2*PI)) * std::exp( -0.5 * (score-mu)*(score-mu)/(sigma *sigma));
}



void Biceps::checkTupleConversion() const
{
    for (map<unsigned char, char>::const_iterator it = aamod.begin(); it != aamod.end(); ++it)
    {
        cout << "ascii: "<<it->first <<" newascii: " << it->second << endl;
    }


}






