/*
 *  fasta.cpp
 *  ETS
 *
 *  Created by Buote Xu on 25.09.09.
 *  
 *
 */

#define numAA 20
#include "fasta.h"

ostream & operator<< ( ostream & os, Tag & tag)


{
	os << "tag: " << tag.tag << '\n';
	os << "massTag: " << tag.massTag << '\n';
	os << "massCurrentTag: " << tag.massCurrentTag << '\n';    
    return os;
}

Fasta::Fasta() : tagLength(5),numTags(20), aveAAweight(111.0f){}



inline double Fasta::convAAToMass(const char AA) const
{
    
    switch(AA) {
            
        case 'A' : return 71.03711;
        case 'C' : return (103.00919 + 57.02146);
        case 'D' : return 115.02694;
        case 'E' : return 129.04259;
        case 'F' : return 147.06841;
        case 'G' : return  57.02146;
        case 'H' : return 137.05891;
        case 'I' : return 113.08406;
        case 'K' : return 128.09496;
        case 'L' : return 113.08406;
        case 'M' : return 131.04049;
        case 'N' : return 114.04293;
        case 'P' : return  97.05276;
        case 'Q' : return 128.05858;
        case 'R' : return 156.10111;
        case 'S' : return  87.03203;
        case 'T' : return 101.04768;
        case 'V' : return  99.06841;
        case 'W' : return 186.07931;
        case 'Y' : return 163.06333;
            
        default: return 0.0;
    }
}

inline unsigned int Fasta::convAAToIndex(const char AA) const
{
    
    switch(AA) {
        case 'A' : return 0U;
        case 'C' : return 1U;
        case 'D' : return 2U;
        case 'E' : return 3U;
        case 'F' : return 4U;
        case 'G' : return 5U;
        case 'H' : return 6U;
        case 'I' : return 7U;
        case 'K' : return 8U;
        case 'L' : return 9U; 
        case 'M' : return 10U;
        case 'N' : return 11U;
        case 'P' : return 12U;
        case 'Q' : return 13U;
        case 'R' : return 14U;
        case 'S' : return 15U;
        case 'T' : return 16U;
        case 'V' : return 17U;
        case 'W' : return 18U;
        case 'Y' : return 19U;
            
        default: return 0U;
    }
}


int Fasta::loadTags(const string & tagstring, const vector<float> & lowPeakMzs, vector<Tag> & tags)
{

    if (!tagstring.size()) return 0;    
    unsigned int numoftags =  tagstring.size() / tagLength;
    int firstlimit =0;
    for (unsigned i = 0; i < numoftags; ++i)
    {
        string temptag = tagstring.substr((this->tagLength)*i, this->tagLength);
        //expand, because of Directag/Pepnovo don't create tags with I and J.
        //As we don't consider I and J equal(at least Pepsplice doesn't), we have to do replacements ourselves.
        expandTag(temptag, 0, lowPeakMzs, i, tags);
        if (i < numoftags/3) firstlimit = tags.size();
    }
    
    return ++firstlimit;
}


void Fasta::expandTag(string tag, const size_t pos1, const vector<float> & lowPeakMzs, const size_t pos2, vector<Tag> & tags) 
{
    if (pos1== tag.size()) 
    {
        push_Tag(tag, lowPeakMzs, pos2,tags);
        return ;
    }
        expandTag(tag, pos1+1, lowPeakMzs, pos2, tags);
    
    if (tag[pos1] == 'I')
    {
        
        tag[pos1] = 'L';
        expandTag(tag, pos1+1, lowPeakMzs, pos2, tags);
    }
    else if (tag[pos1] == 'L')
    {
        tag[pos1] = 'I';
        expandTag(tag, pos1+1, lowPeakMzs, pos2, tags);
    }
}


inline void Fasta::push_Tag(const string & tag, const vector<float> & lowPeakMzs, int index, vector<Tag> & tags)
{
    tags.push_back(Tag(tag, this->calcWeight(tag), lowPeakMzs[index]));
}

                                                
void Fasta::printTags()
{
    
    std::cout << "Directag tags:" << '\n';
    for (unsigned int j = 0; j < (this->directags).size(); ++j)
    {
        std::cout << (this->directags)[j] << std::endl;
    }
    
    std::cout << "Pepnovo tags:" << '\n';
    for (unsigned int j = 0; j < (this->pepntags).size(); ++j)
    {
        std::cout << (this->pepntags)[j] << std::endl;
    }
    
}


inline unsigned int Fasta::cTI2(const string & tag, const bool begin) const 
{
    if (begin)
        return convAAToIndex(tag[0])* numAA + convAAToIndex(tag[1]);
    else
        return convAAToIndex(tag[3]) * numAA + convAAToIndex(tag[4]);
}


float Fasta::calcWeight(const string& tag) const
{
    double tagWeight = 0.0f;
    for (unsigned int i = 0; i < tag.size(); ++i)
    {
        tagWeight += this->convAAToMass(tag[i]);
    }
    return tagWeight;
}


void Fasta::createFasta(const bool mutated, const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, std::vector<std::tuple<unsigned int, std::string, std::string> > & currentFasta)
{
    unsigned int numlasttag;
    if (debuglevel > 1) this->outfile.open("etsresults.fasta", fstream::trunc | ios::binary);
    
    this->matches =0;
    //begin parsing of the tags;
    numlasttag = (tagsend < tags.size())? tagsend : tags.size();
    if (!mutated)
        this->findTags2(tagsbegin, tagsend, tags, currentFasta);
    else 
    {
        this->findMutatedTags2(tagsbegin, tagsend, tags, currentFasta);
    }
    
    outfile.close();
    return;
    
}


void Fasta::loadBigFasta(string & file) //initialize sequences 
{
    
	int i = 0;
	string parsedline;
    
	ifstream inFile;
	cout <<"Now indexing the fasta file: "<< file << endl;
	inFile.open(file.c_str());
	if (!inFile.good()) throw runtime_error("FastaFile " + (string) file + " can't be opened, please check");
    
    
    //read file and write sequences into bigfasta
    //convert to uppercase, if needed
    while(inFile.good())
    {
        getline(inFile, parsedline);
        i++;
        if(parsedline.substr(0, 1) == ">")
        {
            bigfasta_descriptions.push_back(parsedline);
            if(this->currentprotein.size() > 0)
            {
                bigfasta.push_back(this->currentprotein);
                this->currentprotein.clear();
            }
        }
        else
        {
            for(unsigned int i = 0; i < parsedline.size(); i++)
            {
                if(90 >= parsedline[i] && parsedline[i] >= 65) 
                {
                    this->currentprotein+=parsedline[i];
                }
                else if(122 >= parsedline[i] && parsedline[i] >= 97) 
                {
                    this->currentprotein+=parsedline[i]-32;
                } //if
            } //for
        } //else
    } //while 
    if(this->currentprotein.size() > 0)
    {
        bigfasta.push_back(this->currentprotein);
        this->currentprotein.clear();
    }
    //^- Adding the last sequence

    this->indexFasta();
    
}             


//the Indexer uses the hash of all possible combinations of 2 AA 
//and searches for them in the fasta, the results are saved in doubleIndices

void Fasta::indexFasta() 
{
    
    doubleIndices = vector<vector<vector<unsigned int> > > (bigfasta.size(), vector<vector<unsigned int> >(numAA*numAA,vector<unsigned int>()));   
    
    char letter1,letter2;
    
    for (unsigned int i = 0; i < bigfasta.size(); ++i)
    {
        string & sequence = bigfasta[i];
        for (unsigned int j = 0; j + 1 < bigfasta[i].size(); ++j)
        {
            letter1 = sequence[j];
            letter2 = sequence[j+1];
            
            doubleIndices[i][this->convAAToIndex(letter1) * numAA + this->convAAToIndex(letter2)].push_back(j); //easy hash function
        }
    }
}





string Fasta::cutSequence(const Tag & tag, const unsigned int sequence_id, const unsigned int posTag, const int mutated, std::vector<std::tuple<unsigned int, std::string, std::string> > & currentFasta)
{
    const string & sequence = bigfasta[sequence_id];
    //as tag.minAvePoss represents the average amount of mass before the tag occurs in a sequence
    //i.e., the starting position of the given tag in a sequence
    //this has to be considered by adding some more AAs before the occurence of the tag
    size_t minPossTag = posTag; //< beginning of the tag
    double buffer = 0.0;
    while(minPossTag > 0)
    {        
        --minPossTag;
        buffer += convAAToMass(sequence[minPossTag]);
        if (buffer > (tag.massCurrentTag * (charge-1) + 2 * this->aveAAweight) )
            break;
    }
    size_t maxPossTag = posTag + tagLength; //< end of the tag, sequence[maxPossTag] is the last AA
    buffer = 0.0;
    while ( maxPossTag < sequence.size()) {
        
        buffer += convAAToMass(sequence[maxPossTag]);
        if (buffer > (this->precursorMassCharged - tag.massCurrentTag * (charge-1) - tag.massTag + 2 * this->aveAAweight) )
            break;
        ++maxPossTag;
    }


    try{         
        if (mutated == -1)
        {  
            //const string AAsequence = this->bigfasta[sequence_id].substr(minPossTag, maxPossTag - minPossTag);
            
            ///This index-list will save the occurence of the tag in the original fasta.
            const string result = this->bigfasta[sequence_id].substr(minPossTag, maxPossTag - minPossTag);
            currentFasta.push_back(std::make_tuple(sequence_id, result, result));
            return result;
        }
        // return sequence.substr(minPossTag,posTag-minPossTag) + tag.tag + sequence.substr(posTag + tagLength, maxPossTag - tagLength - posTag);
        else 
        {
            const string result = this->bigfasta[sequence_id].substr(minPossTag,posTag-minPossTag + mutated) + tag.tag.substr(mutated,1) + this->bigfasta[sequence_id].substr(posTag + mutated+ 1, maxPossTag-posTag-mutated-1);
            currentFasta.push_back(std::make_tuple(sequence_id, result, this->bigfasta[sequence_id].substr(minPossTag, maxPossTag - minPossTag)));
            return result;
        }
    }
    catch(const std::out_of_range & e)
    {
        cerr << e.what() << "\n";
        cout << "posTag = " << posTag << "\n";
        cout << "minPossTag = " << minPossTag << "\n";
        cout << "maxPossTag = " << maxPossTag << "\n";
        cout << "sequencelength = " << sequence.size() << flush;
        
        throw;
    }
    
    //return the sequence with one changed AA because the tag is correct.
    //return sequence.substr(minPossTag,posTag-minPossTag) + tag.tag + sequence.substr(posTag + tagLength, maxPossTag - tagLength - posTag);
}


void Fasta::findTags2(const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, std::vector<std::tuple<unsigned int, std::string, std::string> > & currentFasta)
{
    vector<unsigned int> prefixes;
    vector<vector<unsigned int> >tagnums;
    
    
    for (unsigned int i = tagsbegin; i < ((tagsend<tags.size())?tagsend:tags.size()); ++i)
    {
        bool jump = false;
        unsigned int prefix = cTI2(tags[i].tag,true);
        //save the 2-prefix of the tags, look if some are the same. do some mapping stuff?
        //need: map tagnumber -> prefix. then iterate through the different tag values! (is that possible ?) go for multimap!
        //prefixes[cTI2(tags[i]),1] = i;
        
        for (unsigned int j = 0; j < prefixes.size(); ++j)
            if (prefix == prefixes[j]) { tagnums[j].push_back(i); jump = true; break;}
        
        if (jump) continue;
        prefixes.push_back(prefix);
        tagnums.push_back(vector<unsigned int>(1,i));
        //all prefixes of the tagend-tagsbegin tags should now be in prefixes, while tagnums gives us the according tagnumber in this->tags.
    }
    
    
    
    for (unsigned int i = 0; i < prefixes.size(); ++i) //try the different prefixes
    {    
        for (unsigned int j = 0; j < bigfasta.size(); ++j) //go through all sequences
        {      
            string & sequence = bigfasta[j];
            for (unsigned int k = 0; k < doubleIndices[j][prefixes[i]].size(); ++k) //the indexes help going through the sequence
            {
                unsigned int prefixpos = doubleIndices[j][prefixes[i]][k];
                if (prefixpos + 5 > bigfasta[j].size()) break; //change here for tagsize !=5
                {
                    for (unsigned int l = 0; l < tagnums[i].size(); ++l) //the first two chars of the tag fit, test if it works for 2 of the other 3 chars!
                    {
                        const string & tag = tags[tagnums[i][l]].tag;
                        if (tag.compare(2, 3, sequence, prefixpos+2, 3) ==0)
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, prefixpos, -1, currentFasta);
                            ++(this->matches); //DEBUG
                            if (debuglevel > 1) outfile << "> " << j << "\n" << shortseq << "\n";
                        }
                    } //tagnums[i].size()
                    continue;
                } //doubleIndices[prefixes[i]].size()   //going through the matches of (tag[0]+tag[1]) with the sequence
            }   //bigfasta.size() // going through all saved sequences
        }   //prefixes.size(); going through the different prefixes of the used tags
    }
}







void Fasta::findMutatedTags2(const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, std::vector<std::tuple<unsigned int, std::string, std::string> > & currentFasta)
{
    

    
    
    vector<unsigned int> prefixes;
    vector<vector<unsigned int> >tagnums;
    
    
    for (unsigned int i = tagsbegin; i < ((tagsend<tags.size())?tagsend:tags.size()); ++i)
    {
        bool jump = false;
        unsigned int prefix = cTI2(tags[i].tag,true);
        //save the 2-prefix of the tags, look if some are the same. do some mapping stuff?
        //need: map tagnumber -> prefix. then iterate through the different tag values! (is that possible ?) go for multimap!
        //prefixes[cTI2(tags[i]),1] = i;
        
        for (unsigned int j = 0; j < prefixes.size(); ++j)
            if (prefix == prefixes[j]) { tagnums[j].push_back(i); jump = true; break;}
        
        if (jump) continue;
        prefixes.push_back(prefix);
        tagnums.push_back(vector<unsigned int>(1,i));
        //all prefixes of the tagend-tagsbegin tags should now be in prefixes, while tagnums gives us the according tagnumber in this->tags.
    }
    
    
    
    for (unsigned int i = 0; i < prefixes.size(); ++i) //try the different prefixes
    {    
        for (unsigned int j = 0; j < bigfasta.size(); ++j) //go through all sequences
        {      
            string & sequence = bigfasta[j];
            for (unsigned int k = 0; k < doubleIndices[j][prefixes[i]].size(); ++k) //the indexes help going through the sequence
            {
                unsigned int prefixpos = doubleIndices[j][prefixes[i]][k];
                
                if (prefixpos + 5 > bigfasta[j].size()) break; //change here for tagsize !=5
                for (unsigned int l = 0; l < tagnums[i].size(); ++l) //the first two chars of the tag fit, test if it works for 2 of the other 3 chars!
                {
                    const string & tag = tags[tagnums[i][l]].tag;
                    
                    if ( tag[2] == sequence[prefixpos + 2])
                    {
                        if (tag[3] == sequence[prefixpos + 3]) 
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, prefixpos, 4, currentFasta);
                            ++(this->matches); //DEBUG
                            if (debuglevel > 1) outfile << "> " << j <<"\n" << shortseq << "\n";
                            
                        }
                        
                        
                        else if (tag[4] == sequence [ prefixpos + 4])
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, prefixpos, 3, currentFasta);
                            ++(this->matches); //DEBUG
                            if (debuglevel > 1) outfile << "> " <<  j << "\n" << shortseq << "\n";
                            
                        }
                    }
                    else if (tag[3] == bigfasta[j] [ prefixpos + 3])
                    {	 
                        if (tag[4] == bigfasta[j] [ prefixpos + 4])	 
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, prefixpos, 2, currentFasta);
                            ++(this->matches); //DEBUG
                            if (debuglevel > 1) outfile << "> " << j << "\n" << shortseq << "\n";
                            //outfile <<" " << tag;
                            
                            
                        }
                    }
                } //tagnums[i].size()
                continue;
            } //doubleIndices[prefixes[i]].size()   //going through the matches of (tag[0]+tag[1]) with the sequence
        }   //bigfasta.size() // going through all saved sequences
    }   //prefixes.size(); going through the different prefixes of the used tags
    
    
    
    //now the same algorithm backwards for growing to the front!
    
    vector<unsigned int> suffixes;
    //vector<vector<unsigned int> >tagnums;
    tagnums.clear();
    
    for (unsigned int i = tagsbegin; i < ((tagsend<tags.size())?tagsend:tags.size()); ++i)
    {
        bool jump = false;
        unsigned int suffix = cTI2(tags[i].tag,false); //getting the suffix here / converting it to a unsigned int
        //save the 2-suffix of the tags, look if some are the same. do some mapping stuff?
        //need: map tagnumber -> suffix. then iterate through the different tag values! (is that possible ?) go for multimap!
        
        
        for (unsigned int j = 0; j < suffixes.size(); ++j)
            if (suffix == suffixes[j]) { tagnums[j].push_back(i); jump = true; break; }
        
        if (jump) continue;
        suffixes.push_back(suffix);
        tagnums.push_back(vector<unsigned int>(1,i));
        
    }
    
    
    
    for (unsigned int i = 0; i < suffixes.size(); ++i)//go through suffixes
    {    
        for (unsigned int j = 0; j < bigfasta.size(); ++j)//go through all eqs
        {      
			string & sequence = bigfasta[j];
            for (unsigned int k = 0; k < doubleIndices[j][suffixes[i]].size(); ++k)
            {
                unsigned int suffixpos = doubleIndices[j][suffixes[i]][k]; 
				//check if the suffix is at least 3 after start
                if (suffixpos  < 3) break; //change here for tagsize !=5, this checks for boundaries
                for (unsigned int l = 0; l < tagnums[i].size(); ++l) //the first two chars of the tag fit, test if it works for 2 of the other 3 chars!
                {
					const string & tag = tags[tagnums[i][l]].tag;
                    if (  tag[2] == sequence [ suffixpos -1])
                    {
                        if (tag[1] == sequence [ suffixpos -2]) 
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, suffixpos-3, 0, currentFasta);
                            ++(this->matches); //DEBUG
                            outfile << "> " << j<<tag << "\n" << shortseq << "\n";
                            
                        }
                        
                        
                        else if (tag[0] == sequence [ suffixpos -3])
                        {
                            string shortseq = this->cutSequence(tags[tagnums[i][l]], j, suffixpos-3, 1, currentFasta);
                            ++(this->matches); //DEBUG
                            outfile << "> " << j<<tag  << "\n" << shortseq << "\n";
                            
                        }
                    }
                    
                    
                } //tagnums[i].size()
                continue;
            } //doubleIndices[suffixes[i]].size()   //going through the matches of (tag[0]+tag[1]) with the sequence
        }   //bigfasta.size() // going through all saved sequences
    }   //suffixes.size(); going through the different prefixes of the used tags
    
    
    
    
}


void Fasta::matchSequences(PepspliceResult & pepresult, const std::vector<std::tuple<unsigned int, std::string, std::string> >& currentFasta) 
{
  for (size_t i = 0; i < currentFasta.size(); ++i)
  {
    if (std::get<1>(currentFasta[i]).find(pepresult.OrigSequence) != std::string::npos)
    {
      pepresult.fastaIds.push_back(std::get<0>(currentFasta[i]));    
    }
  }
}




