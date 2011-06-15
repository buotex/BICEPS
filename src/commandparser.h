///This File contains the Command-line Parser
///The Implementation is done via boost-program-options
///Author: Buote Xu

/**
 * @Author Buote Xu (Buote.Xu@stud.uni-heidelberg.de)
 * @date   January, 2010
 * @brief  Header for including Commandline-Parsing in Biceps
 *
 */




#ifndef __COMMPARSER__H
#define __COMMPARSER__H

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

using std::string;
using std::vector;
using std::fstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using boost::lexical_cast;

///for easier / clearer options
enum TOOLS {
    PEPNOVO,
    DIRECTAG,
    PEPNDIREC
};

///helper function
template <typename T>
void print_vector(std::vector<T> & vec){
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << std::endl;
    }
}


///Some of the general options
///These will be defined by the Parser and used throughout the program
///Default values are given in the Parser itself

struct general_options{

    /// Chooses the Tool to use.

    /// 0 == Pepnovo, 1 == Directag, 2 == Both.
	/// Default: 2 (Int).

    int tool;


    ///Number of tags used, should be in range 1-50.
	int numTags;

    /// Name of complete Spectrum list, divided by END IONS / BEGIN IONS.

    string mgfname;

    /// Name of Fasta-database.

	string fastaname;

    /// Debuglevel.

    /// Changes amount of Debug output
    /// Default: 1, Values: 1-3
    int debug;


    /// If 0, no mutation of the tags will be done by Biceps.
    /// This will limit the number of matches tremendously,
    /// only deactivate for speed reasons
    /// Default: true bool
    bool mutation;


    /// penalty-limits used in pepsplice, to limit the search-space.

    /// found sequences are mutated in pepsplice according to these limits,
    /// the values used will be elevated slightly because of numerical reasons
    /// Default: 2.3,3.2 (list of floats)
    vector<float> max_penalties;

    ///similar to max_penalties, saves the penalties for mutated tags.
    vector<float> max_penalties_mutated;
    /// Name of the Buffer for a single spectrum found in the mgf.

    /// It will be passed to the libraries.
    /// Default: __buffer__.mgf (string)
    string sSpectrumFN;
};


/**
 * @name    parseProgramOptions
 * @brief   The function takes argc and argv
 * Additionally, You have to pass the three vector<string> & and a general_options object
 * In return, this function will fill the objects accordingly so that you can pass them
 * directly to their respective functions in Biceps.
 *
 * This API provides certain actions as an example.
 *
 *
 * @retval 0  Successfully parsed command line, program can continue
 * @retval 1  Version or Help string are written to command line, stop program
 * @retval 2  Couldn't parse something, stop with error
 *
 * Example Usage:
 * @code
 * vector<string> directag_options, pepnovo_options, pepsplice_options;
 * general_options general_options_object;
 * parseProgramOptions(
 *   argc, argv,
 *    directag_options, pepnovo_options, pepsplice_options,
 *    general_options_object);
 *
 * run_pepsplice(... , pepsplice_options, ...);
 *
 * @endcode
 */



int parseProgramOptions(
    int ac, //< argc from main
    char* av[], //< argv from main
    vector<string> & directOp, //< contains the important options for directag
    vector<string> & pepnOp, //< options for pepnovo
    vector<string>& pepsOp, //< options for pepsplice
    general_options & genOp); //< general options

#endif
