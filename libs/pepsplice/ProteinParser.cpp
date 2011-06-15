#include "ProteinParser.h"
namespace Pepsplice{
  ProteinParser::ProteinParser(Services *se0, Tuples *tp0)
  {
    se1 = se0;
    tuples1 = tp0;
    dnaAA1 = DnaAA();
    //currentprotein = new Protein("initialize");
    oldSequences = 0;
    oldSequenceCount = 0;
    proteincount = 0;
  }
  
  ProteinParser::~ProteinParser()
  {
    if (oldSequences != 0)
      delete oldSequences;
    //delete sequences;
    //delete ids;
    //  delete currentprotein; //Changes by BX
  }
  
  
  
  void ProteinParser::parseFASTA(const map<string,string> & currentfasta)
  {
    
    if (oldSequences != 0) delete oldSequences;
    
    oldSequences = new vector<string>;
    proteincount = 0;
    oldSequenceCount = 0;
    map<string, string >::const_iterator it;
    //BX: second is the original sequence, first the modified sequence
    try
    {
      for (it = currentfasta.begin(); it != currentfasta.end(); ++it)
      {
        currentprotein = new Protein(it->second, it->first);
        if (se1->outputlevel>2) {se1->os << it->second << " " << it->first << '\n';}
        
        finishProtein();
        delete currentprotein;
        proteincount++;
      }
    }
    catch(const runtime_error & e)
    {
      cerr << e.what();
      if (currentprotein != 0)
        delete currentprotein;
    }
    
    
    if (se1->outputlevel>1){ cout << "\n"<< proteincount << " proteins parsed" << endl;}
  }
  
  
  
  inline void ProteinParser::finishProtein()
  {
    generateTuples();
    se1->aminoacids += currentprotein->getLength();
    
  }
  
  void ProteinParser::generateTuples()
  {
    int size = currentprotein->getLength();
    //cout << "\n" << currentprotein->protseq << flush;
    int j = 0;
    int pe = 0;
    int len = 0;
    double monoparentmassMH = 0.0;
    string & seq = currentprotein->protseq;
    string & oldseq = currentprotein->protid;
    bool istrypticstart = false;
    bool istrypticend = false;
    //int counter = 0;
    
    for(int ps = 0; ps < size; ps++){
      
      //cout << seq.substr(ps, 1);
      
      //start (tryptic or not)
      if( ps == 0 || (ps 		> 1 && (seq.at(ps-1) == 'R' || seq.at(ps-1) =='K')) ){istrypticstart = true;}else{istrypticstart = false;}
      //had to use at() instead of [] because the compiler would optimize it into something non-workable
      
      if(istrypticstart >= se1->trypticendsrequired - 1){
        j = ps;
        while(j < size)
        {
          if(dnaAA1.aa_monomasses256[seq[j]] <= 0) break;
          //end (tryptic or not)
          if(seq.at(j) == 'K' || seq.at(j) == 'R'){istrypticend = true;}else{istrypticend = false;}
          if(istrypticstart + istrypticend >= se1->trypticendsrequired){
            pe = j;
            len = pe - ps + 1;
            monoparentmassMH = dnaAA1.getParentMassMH(seq.substr(ps, len));
            
            //start changes BYR: remove ifs
            //if(monoparentmassMH > se1->max_monoparentmassMH) break;
            //changes BYR changed Tuple call to include original
       
            if (len > MAXSEQUENCESIZE) throw runtime_error("Sequence too large to fit in Memorypool, adjust MAXSEQUENCESIZE in definitions.h");
            if (se1->outputlevel > 2)
            {
              se1->os << seq.substr(ps,len) << " " << oldseq.substr(ps,len) << endl;
            }
            //BX: As the original sequences won't change (and several tuples have the same), it's easier/faster to just reference them.
            
            oldSequences->push_back(oldseq.substr(ps, len));
            //changed by BX, using the string sequence now, copying happens in the constructor.
            tuples1->addTuple(new Tuple(monoparentmassMH, len, ps, pe, pe+1, pe, seq, oldSequenceCount,NULL, currentprotein, false, false, istrypticstart, istrypticend));
            oldSequenceCount++;
            //cout << "\n hightup is here"<<hightup;
            //}
            //end changes BYR
            
          }
          
          j++;
        }//end while
      }
    }
    
  }
  
  void ProteinParser::doAnotherRun(const  map<string,string> & currentfasta) //added by BX, using old buffered tuples to improve speed
  { 
    
    se1->tuples = 0;
    se1->randomtuples=0;
    se1->realtuples=0;
    parseFASTA(currentfasta);
    tuples1->forwardTuples();

  }
  
  
  
  
  
  
  
  unsigned char* ProteinParser::getAASeqChar(string seq){
    unsigned char *aaseq = new unsigned char[seq.size()];
    for(int i = 0; i < seq.size(); i++){
      aaseq[i] = seq[i];
    }
    return aaseq;
  }
}
