#ifndef __FASTA_H__
#define __FASTA_H__
#include <string>
#include <vector>
#include <list>
#include <deque>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "pepsplice.h"

using namespace std;


/**
 * @name    Tag
 * @brief   This struct is used to save all necessary information for the directag/pepnovo tags.
 *
 * Example Usage:
 * @code
 * Tag tag("AAIO", 500.0, 200.0);
 * @endcode
 */


struct Tag{
    ///sum(weight(tag.AminoAcids))
    float massTag;
    ///LowPeakMzs, either from pepnovo or directag, represents the probable mass position of the first AA in the tag
    float massCurrentTag;
    ///the actual Tag
    std::string tag;

    ///take the parameters and just save them to the struct without changes.
    Tag(const std::string & tag_, const float massTag_, const float massCurrentTag_):massTag(massTag_), massCurrentTag(massCurrentTag_), tag(tag_){}

};


///overloaded cout << Tag for easier logging output
ostream & operator<< ( ostream & os, Tag & tag);


/** Fasta class
 * 
 * This is a class for dealing with fasta-databases, like ipi.BOVIN.fasta.
 * Searching: Given Tags, which should be of the Tag class (with length of 5), this class is able to search for matching subsequences in the database.
 * While it is possible to only aim for identical sequences, mutations can also be considered.
 * To provide this functionality, the complete Database has to be preloaded by the <tt>loadBigFasta(cstring)</tt> method, which will also parse and hash it.
 * This enables us to search for the beginnings of the tags in the created hash-table, which should speed up the search as most of the sequences can be ignored.
 * Considering mutated tags, valid hits consist of subsequences in the database which differ from the provided tags by less than 2 AAs.
 * As actually creating all different mutations and searching for them in the database is very time-consuming, an approach is chosen where it's possible that 1 compare between the tag and the subsequence may fail,
 * though this is only possible if we use a hash table to actually match either the beginning or the ending of a tag to the subsequence first.
 *
 * Creating: Pepsplice, another library which is used in the ETS program, needs fasta-databases to find AA-sequences.
 * As generic databases are too big to provide results in a reasonable time, a smaller customized one is more viable and will thus be created by the <tt>createFasta</tt>-method.
 * The result, a <tt>map<string, string></tt> containing modified and unmodified subsequences can then be processed by pepsplice.
 *
 *
 *
 * @retval Output will be a fasta map<string, string> in the <mutated tagstring, original tagstring> format.
 *
 * Example Usage:
 * @code
 *
 * vector<Tag> tagbuffer;
 * map<string, string> pepnFasta;
 *
 * Fasta fasta;
 * fasta.loadBigFasta("ipi.CHICK.fasta");
 * fasta.loadTags(string("AAAAA"), lowPeakMzs, tagbuffer);
 * fasta.createFasta(1,0,0, tagbuffer, pepnFasta);
 *
 *
 * @endcode
 */




class Fasta
{
    friend class Biceps;
    public:
    ///The Default Constructor just initializes numOfTags to 20 and
    ///averageAAWeight to 111.0, which is used for heuristic expanding of the tag.

    Fasta();


    ///createFasta will search for (partial, depending on the mutated parameter) matches of the tags compared to the previously loaded fasta.
    ///tagsbegin and tagsend limit the indices of the searched tags, results will be written into the given currentFasta.
    /**
      \param[in] mutated boolean, sets if mutated tags should also be used for the matching algorithm
      \param[in] tagsbegin sets the first tag that should be used to search
      \param[in] tagsend use up to tagsend - 1th tag to search for matches
      \param[in] tags the tag-vector containing the tags to search for
      \param[out] currentFasta this will save all the matching-results
     */
    void createFasta(const bool mutated, const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, map<string, string> & currentFasta);

    ///loadBigFasta needs the filename (as a std::string) of the fastafile, it will then call some private functions.

    ///They will then create an easier to search fastafile with replaced I AAs,
    ///a hashing index containing all 2-combinations of AAs will also be made
    ///\param[in] filename of the fasta-database to index and match/search later, e.g. loadBigFasta(string("ipi.BOVIN.fasta"))

    void loadBigFasta(std::string & filename);


    ///This method takes the unmodified tagstring, changes the L AAs to I for performance reasons and fills the internal tag vectors for further processing
    /**
      \param[in] tagstring string which is made up of the different tags found by pepnovo/directag, no spaces.
      \param[in] lowPeakMzs is also calculated by pepnovo/directag, it has to be multiplicated by the spectrum charge-1 to calculate the mass.
      \retval int number of tags found/created + 1
     */
    int loadTags(const std::string & tagstring, const vector<float> & lowPeakMzs, vector<Tag> & tags);

    /// just print all of the currently saved tags in the fasta class, this will be pepnTags and direcTags
    void printTags();

    ///initializeMGF sets the charge and the precursormass, it's used in the main function for every additional spectrum.
    ///the precursormass will be multiplicated with the charge and stored in the Fasta class
    void initializeMGF(unsigned int charge_, float precursorMass_){this->charge = charge_, this->precursorMassCharged = this->charge * precursorMass_;}

    /// change the limit of saved tags, default is 20 per directag/pepnovo program.
    void setNumTags(unsigned int numTags_){ this->numTags = numTags_;}

    ///counts the number of matches btween tags and fasta
    unsigned int getMatches(){ return this->matches;}

    ///insert id, get the description in the fastafile for that specific id.
    ///\param[in] seqid id of the spectrum you want the information for
    ///\retval string the description of that specific spectrum taken from the fasta file
    string offerDescription(unsigned int seqid) const { return bigfasta_descriptions[seqid]; }

    ///clear all tags
    void clearTags(){ pepntags.clear(); directags.clear();}

    ///public interface
    unsigned int getNumTags(){ return this->pepntags.size() + this->directags.size();}


    private:

    ///calculate the sum of the AA-weights in the tag
    float calcWeight(const string & tag) const;

    ///As we don't want to put whole sequences into the new fasta-database, only short sequences containing the Tag will be added to currentFasta, also return the mutated shortSequence
    /**
      \param[in] tag the tag which is matched against the fasta-entry
      \param[in] sequence_id the number of the sequence in the fasta-database
      \param[in] posTag the beginning of the (mutated) match between entry and tag
      \param[in] mutated the boolean deciding whether mutations should be considered
      \param[out] currentFasta this will save all the matching-results
      \retval string mutated substring which matches the tag, with padding on both sides and with optional mutations
     */
    std::string cutSequence(const Tag & tag, const unsigned int sequence_id, const unsigned int posTag, const int mutated, map<string, string> & currentFasta);

    /**calculating the index of a Tag, if begin ==1 it uses the first two AAs, otherwise it's the last two, helper function to search faster for tags.

      Returns the hash value of the string, using trivial hash with 19 * tag[0] + tag[1] here
      Probably has to be enhanced when looking at tagsize > 5
      Very hardcoded, but fast!
      \param[in] tag the string-portion of the tag.
      \param[in] begin if true, get the hash for the pair of tag[0], tag[1], if false return the hash of tag[3] tag[4]
      \retval hash-value of a pair of chars
     */

    inline unsigned int cTI2(const string & tag,const bool begin) const;

    ///parse the complete fasta-database while indexing all occurences of all AA pairs for a quicker tag-search later
    void indexFasta();


    ///tagsbegin and tagsend limit the index of the used tags - we don't want to use all of them
    ///this find function will use a previously calculated index and the hashed sequences.
    ///the cut sequences will be written to currentFasta via cutSequence
    /**
    \param[in] tagsbegin the index of the first tag to search for
    \param[in] tagsend index+1 of the last tag to search for
    \param[in] tags the tag-vector to search for
    \param[out] currentFasta this will save all the matching-results 
    */
    void findTags2(const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, map<string, string> & currentFasta);


    ///this function searches for partial matches (up to 1 differing AA) between given tags and the previously loaded Fasta.
    ///A speedup was possible via hashing all sequences and searching for parts of a given tag.
    ///very hardcoded, works only for size-5 tags, using the index and 'jumping' over the mutated parts.
    /**
    \param[in] tagsbegin the index of the first tag to search for
    \param[in] tagsend index+1 of the last tag to search for
    \param[in] tags the tag-vector to search for
    \param[out] currentFasta this will save all the matching-results 
    */ 
    void findMutatedTags2(const unsigned int tagsbegin, const unsigned int tagsend, const vector<Tag> & tags, map<string, string> & currentFasta);


    ///switch-case for getting the mass of a specific AA, should be fastest solution with jump-table
    inline double convAAToMass(const char AA) const;

    ///switch-case for getting index of specific AA for indexFasta, should be fastest solution with jump-table
    inline unsigned int convAAToIndex(const char AA) const;

    ///helper function, because pepnovo/directag don't differentiate between I and L, the other tags have to be generated hereby
    /**
    \param[in] tag the string of the expanded tag
    \param[in] pos1 the current letter to be expanded, will only expand if it's I or L
    \param[in] lowPeakMzs a vector with the Peakvalues, gotten from directag or pepnovo
    \param[out] tags the vector to write the results into
    */
    void expandTag(string tag, const size_t pos1, const vector<float> & lowPeakMzs, const size_t pos2, vector<Tag> & tags);
    ///nothing special, just a short way to add to the tags
    /**
    \param[in] tag the string of the to be saved tag
    \param[in] lowPeakMzs a vector with the Peakvalues, gotten from directag or pepnovo
    \param[in] index the index of the corresponding peakvalue in lowPeakMzs
    \param[out] tags the vector to write the results into
    */
    inline void push_Tag(const string & tag, const vector<float> & lowPeakMzs, int index, vector<Tag> & tags);




    public:
    ///if Debugoutput is needed, outfile will be used for saving the resulting fasta-database.
    std::ofstream outfile;
    ///change the debug level, 1 is default. Higher debuglevels will provide file-logging.
    int debuglevel;

    private:
    ///number of matches
    unsigned int matches;
    ///indexed fasta file
    vector<string> bigfasta;
    //vector<string> bigfasta_original;

    ///descriptions for later output
    vector<string> bigfasta_descriptions;
    ///the tags from directag, which are produced by loadTags()
    vector<Tag> directags;
    ///the tags from pepnovo, which are produced by loadTags()
    vector<Tag> pepntags;

    ///not used at the moment, as the hardcoded tagLength is 5 for performance reasons.
    unsigned int tagLength;

    ///number of tags
    unsigned int numTags;
    string currentprotein; //<for memory reasons it's a class member
    //string currentproteinI;
    unsigned int charge; //< charge of the current spectrum
    float precursorMassCharged; //<precursorMass in m, so multiplied by charge
    float aveAAweight; //<Average AA weight, for some calculations in cutting


    ///Indexing all Tags with 2 AAs, for easier/faster search later
    ///It's a pretty ugly triple-vector, access is done via
    ///doubleIndices[SequenceID][hash_of_pair][entry];
    vector<vector<vector<unsigned int> > > doubleIndices;
    
    ///ETS needs some of the private functiosn
    friend class ETS;

};
#endif
