#ifndef __BICEPS_H__
#define __BICEPS_H__
//stl
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <cmath>
#include "boost/lexical_cast.hpp"

//libraries
#include "pepsplice.h"

using Pepsplice::PepspliceResult;
using Pepsplice::PepspliceResultComparator;

#ifdef HAVE_DIRECTAG
	#include "directag.h"
#endif

#ifdef HAVE_PEPNOVO
	#include "Pepnovo.h"
#endif

#include "gmm-bic.h"
#include "fasta.h"

//command-line parser
#include "commandparser.h"


using boost::lexical_cast;
using boost::bad_lexical_cast;

/**
 * @author Buote Xu (Buote.Xu@stud.uni-heidelberg.de)
 * @date   January, 2010
 * @brief  The Biceps-Code calls the different libraries,
 * also parses the commandline for parameters.
 *
 * This class represents the glue to unite the different libraries.
 * It works in several steps:
 * At the beginning, <tt>initialize(argc, argv)</tt> will call the commandparser, which is responsible for parsing the command-line parameters and changing the settings/flags accordingly.
 * During the next step, the mandatory fasta-database will be indexed and be made ready for later parsing.
 * The other mandatory parameter corresponds to the MGF-file, which is parsed line by line via <tt> readMGF </tt>.
 * Every time a complete spectrum is parsed (with BEGIN IONS and END IONS lines), it gBiceps passed in a temporary file to Directag/Pepnovo in <tt> run_programs</tt>.
 * One of these two libraries (or even both, depending on the settings), creates a list of viable tags corresponding to the current spectrum found in the mgf.
 * <tt> runPepsplice </tt> calls the Fasta class to process these tags and matching them to the original fasta, to create a customized fasta to pass to the pepsplice-library along with the mgf.
 * This will happen several times per spectrum, as different parameters to pepsplice have to be processed if the result wasn't good enough yet.
 * When all of it is finished, the best pepsplice-outputs for every spectrum is saved and written to a file at the end, while cutting out some of them according to a quantile calculated by <tt>findCutoff </tt>
 *
 * Example Usage:
 * @code
 * int main(int argc, char *argv[])
 * {
 *   Biceps(argc,argv);
 * }
 * @endcode
 */


class Biceps{
private:
    ///Options for directag, for number of tags, mgf-name etc.
    vector<string> directOp;

    ///Options for pepsplice, has fastaname, mgf-name etc.
    vector<string> pepnOp;

    ///Options for pepnovo, holding mgf-filename etc.
    vector<string> pepsOp;

    ///general_options, for Number of tags etc.
    general_options genOp;

    ///Save all the results of the Pepsplice outputs
    vector<PepspliceResult> pepResults;

    ///See class Fasta documentation
    Fasta fasta;

    ///index of the current sub-mgf
    unsigned int mgfId;

    ///Pepsplice uses a separate in_AAmodifications.param to decide the possible mutations.
    ///This map is needed to re-translate the modifications into ascii for logging reasons.
    map<unsigned char, char> aamod;


private:

    /** call gmm-bic and find the quantile, selector sBiceps the number of gaussians.
        @param[out] mu saves the mu-values of the different gaussians
        @param[out] sigma saves the sigma-values of the different gaussians
    */
    double findCutoff(std::vector<double> & mu, std::vector<double> & sigma, int selector = 2);

    /** create a title-textfile and write all results with penalty, bic etc
        @param[out] os stream to log to, default value: file
        @param[in] res a pepspliceresult from runPrograms
        @param[in] specIndex index of the current spectrum in the mgf
        @param[in] title info for the current spectrum
    */
    int writeResult(std::ostream & os, const Pepsplice::PepspliceResult & res, const int specIndex, const std::string & title) const;


    /**
        similar to writeResult, also add confidence to the result and mgf-index
        @param[in] pepResults a vector with pepsplice-results from different spectra
        @param[in] indices a vector with the indices of the spectra, fitting to pepResults
        @param[in] titles the info-tags of the spectra
        @param[in] labels vector with the labels of the results according to gmm-bic
    */
    void writeCompleteResult(const std::vector<PepspliceResult> & pepResults, const std::vector<int> & indices, const std::vector<string> & titles, const std::vector<int> & labels, const std::vector<double> & mu, const std::vector<double> & sigma, const double cutoff) const;

    ///write the final fasta with modifications.
    void writeFasta(std::vector<PepspliceResult> & pepResults);
    ///Call the commandparser to get all the needed parameters, start the Fasta parsing and go through the .mgf afterwards
    void initialize(int argc, char* argv[]);

    ///show Pepnovo/Directag/Pepsplice Licenses
    void showLicenses();

    ///check for line endings, 1 if '\r\n', 0 otherwise
    bool isDosFile(std::string & filename);

    ///run pepnovo/directag and pepsplice, return 1 if it was successful and a legit Sequence was found, 0 otherwise.
    bool run_programs();

    /**pep and dir N1/N2 are the lower / upper limits for the tagindices which should be considered, tags are saved in the fasta class
    and created by Directag/Pepnovo
    @param[in] mutation bool, use mutated tags or not
    @param[in] pepN1 lower bound for the tagindices for pepnovo
    @param[in] pepN2 upper bound for the tagindices for pepnovo
    @param[in] dirN1 lower bound for the tagindices for directag
    @param[in] dirN2 upper bound for the tagindices for directag
    */
    
    PepspliceResult runPepsplice(bool mutation, int pepN1,int pepN2, int dirN1, int dirN2 );

    ///parse the MGF file, this function will always parse inbetween BEGIN IONS and END IONS and write this to a temporary file which will be parsed by the other programs.
    void readMGF();
    //void addResToFile(std::fstream);


    ///parse the inAAModifications.param file which is also necessary for pepsplice, save internal modifications in AA weights in this->aamod;
    void loadAAModifications();

    ///converts Pepsplice output to readable output, has to be done because of the internal modifications
    string convertTupleString(const string & sequence)const ;

    ///print this->aamod;
    void checkTupleConversion()const;

    ///give a confidence for the result-quality
    double returnConfidence(double score, double mu, double sigma) const;

public:
    ///Standard constructor, call it with argc and argv
    Biceps(int argc, char* argv[]){
        initialize(argc, argv);
    }


};




#endif
