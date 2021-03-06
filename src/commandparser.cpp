#include "commandparser.h"



int parseProgramOptions(
    int ac, 
    char* av[], 
    vector<string> & directOp, 
    vector<string> & pepnOp, 
    vector<string>& pepsOp, 
    general_options & genOp)
{


  //for help, read the boost::program_options documentation.


  namespace po = boost::program_options;
  // Declare a group of options that will be 
  // allowed only on command line


  //Some default command-line strings which are always used.

  genOp.debug = 0;  

  //filling, so that the libraries start parsing from the correct entry.
  directOp.push_back("ets");
  pepnOp.push_back("ets");
  pepsOp.push_back ("ets");

  //The spectrum is split by the main-loop in the Ets Class
  //the libraries need the name of the sole spectrum in progress
  //Edit: The names will be changed by the Parser later, to allow several parallel starts
  directOp.push_back("");
  pepnOp.push_back("-file");
  pepnOp.push_back("");
  pepsOp.push_back("");

  //Choosing the model in Pepnovo by default
  pepnOp.push_back(lexical_cast<string>("-model"));
  pepnOp.push_back(lexical_cast<string>("CID_IT_TRYP"));
  pepnOp.push_back(lexical_cast<string>("-tag_length"));
  pepnOp.push_back(lexical_cast<string>("5"));

  //standard directag options
  directOp.push_back(lexical_cast<string>("-UseChargeStateFromMS"));
  directOp.push_back(lexical_cast<string>("1"));



  //Some Pepsplice parameters which should be used every time.
  pepsOp.push_back(lexical_cast<string>("-sa1"));
  pepsOp.push_back(lexical_cast<string>("-bm20"));
  pepsOp.push_back(lexical_cast<string>("-te1"));
  pepsOp.push_back(lexical_cast<string>("-np0"));

  //The Parser starts working now.
  try
  {

    int numtags; //number of different tags
    int tool; //see enum for tool-association
    //int tagLength;
    std::string mgf; //name of the main mgf
    double tol; //tolerance for the pepsplice mass
    std::string fasta; //fasta
    std::string penaltyvector; //penalties
    std::string penaltyvecmutated; //mutated penalties
    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help", "produce help message")
      ;

    //See Boost::ProgramParser documentation for details 

    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
      ("mgf", po::value<string>(&mgf),"path to the mgf file used as input data")
      ("fasta", po::value<string>(&fasta),
       "path to the fasta file used as reference for the search")
      ("tags", po::value<int>(&numtags)->default_value(20), 
       "number of tags to be generated, more tags increase the search space and run time")
      ("tool", po::value<int>(&tool)->default_value(2),
       "which tools to be used for tag generation (0 - pepnovo, 1 - directag, 2- both")
      //("tagLength", po::value<int>(&tagLength)->default_value(5),
      // "tagLength, should be higher than 2")
      ("tol", po::value<double>(&tol)->default_value(10),
       "massTolerance for running pepsplice in ppm")
      ("debug",po::value<int>(&genOp.debug)->default_value(0), "display debug messages")
      ("mutation",po::value<bool>(&genOp.mutation)->default_value(true), "Should substitutions be considered within the tag when necessary (1=yes) or should all tags be considered to be error-free allowing a faster search (0)")
      ("penaltyvector", po::value< std::string >(&penaltyvector)->default_value("2.3,3.2"), "vector of penalties used in pepsplice indicating the size of the search space, important 0 - not including any changes to sequence,")
      ("penaltyvecmutated", po::value< std::string >(&penaltyvecmutated), "penalties for mutated tags, by default same penalties as non-mutated will be used")
      ;


    //po::options_description hidden("Hidden options");
    //hidden.add_options()
    //("mgf", po::value< vector<string> >(), "mgf")
    //;

    po::options_description cmdline_options;
    //cmdline_options.add(generic).add(config).add(hidden);
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
    //        config_file_options.add(config).add(hidden);
    config_file_options.add(config);


    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    //po::positional_options_description p;
    //p.add("mgf", -1);

    po::variables_map vm;
    store(po::command_line_parser(ac, av).
        //	  options(cmdline_options).positional(p).run(), vm);
      options(cmdline_options).run(),vm);

    ifstream ifs("biceps_config.cfg");
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);

    if (vm.count("help")) 
    {
      std::cout << "Use biceps_config.cfg if you want to pass the parameters in a file" << "\n";
      std::cout << visible << "\n";
      return 1;
    }
    if (vm.count("debug"))
    {
      if (genOp.debug>0){
        std::cout << "DebugLevel " << genOp.debug <<  endl;
        pepsOp.insert(pepsOp.begin()+2,(lexical_cast<string> ("-ol")+ lexical_cast<string>(genOp.debug)) );
      }
    }
    if (vm.count("version")) 
    {
      std::cout << "Biceps, version 1.0\n";
      return 1;
    }

    if (vm.count("mgf"))
    {
      std::cout << "Input files are: " 
        << mgf << "\n";
      genOp.mgfAbsPath = boost::filesystem::absolute(boost::filesystem::path(mgf));
      genOp.mgfShortName = genOp.mgfAbsPath.stem().string(); 
      
      //std::cout << vm["input-file"].as<vector<string> >();
      genOp.sSpectrumFN = string("Biceps.buffer.") + genOp.mgfShortName+ string(".mgf");
      directOp[1] = genOp.sSpectrumFN;
      pepnOp[2] = genOp.sSpectrumFN;
      pepsOp[1] = genOp.sSpectrumFN;
    }
    else
    {
      std::cout << visible << "\n";
      return 1;
    }




    if (vm.count("penaltyvector"))
    {
      // std::cout 
      // << "penalties are: " << penaltyvector << "\n";
      penaltyvector +=',';
      unsigned int oldpos = 0;
      try{
        for (unsigned int i = 0; i < penaltyvector.size(); ++i)
        {
          if (penaltyvector[i] == ',')
          {
            genOp.max_penalties.push_back(lexical_cast<float>(penaltyvector.substr(oldpos,i-oldpos) ) + 0.01f);
            oldpos = i+1;
          }
        }
      }
      catch(boost::bad_lexical_cast &e)
      {
        cerr << "Bad penalties, using defaults: 3.21" << "\n";
        genOp.max_penalties.clear();
        genOp.max_penalties.push_back(3.21f);
      }
      if (vm.count("penaltyvecmutated"))
      {
        penaltyvecmutated +=',';
        unsigned int oldpos = 0;
        try{
          for (unsigned int i = 0; i < penaltyvecmutated.size(); ++i)
          {
            if (penaltyvector[i] == ',')
            {
              genOp.max_penalties_mutated.push_back(lexical_cast<float>(penaltyvecmutated.substr(oldpos,i-oldpos) ) + 0.01f);
              oldpos = i+1;
            }
          }
        }
        catch(boost::bad_lexical_cast &e)
        {
          cerr << "Bad mutated penalties, using defaults: " << "\n";
          genOp.max_penalties_mutated.clear();
          genOp.max_penalties_mutated.push_back(2.31f);
        }
      }
      else //no parameters given for mutated max penalties 
      {
        for (size_t i = 0; i < genOp.max_penalties.size(); i++)
        {
          genOp.max_penalties_mutated.push_back(genOp.max_penalties[i]);
        }

      }



      std::cout
        << "penalties are: " <<"\n";
      print_vector(genOp.max_penalties);
      std::cout
        << "mutated penalties are: " <<"\n";
      print_vector(genOp.max_penalties_mutated);


    }


    if (vm.count("fasta"))
    {
      std::cout
        << fasta << " will be used as a protein database." << "\n";
      genOp.fastaAbsPath = boost::filesystem::absolute(boost::filesystem::path(fasta));
    }
    else
    {
      std::cout << visible << "\n";
      return 1;
    }
    if (vm.count("tags"))
    {
      std::cout
        << numtags <<" Tags will be used" << "\n";
      //pepnOp[2] =(lexical_cast<string> (numtags);
      directOp.push_back((string) "-MaxTagCount");
      directOp.push_back(lexical_cast<string>(numtags)); //Using Tags = 20 as a default value
      pepnOp.push_back((string) "-num_solutions");
      pepnOp.push_back(lexical_cast<string>(numtags));
      genOp.numTags = numtags;
    }
    if (vm.count("mutation"))
    {
      if (genOp.mutation == 1){
        std::cout 
          << "If the results aren't good enough with the normal tags, mutated ones will be used" << std::endl;
      }
      else{
        std::cout
          << "No mutation will be done on the tags" << std::endl;


      }

    }


    if (vm.count("tol"))
    {
      std::cout
        << "MassTolerance is set to: " << tol << "\n";
      pepsOp.push_back(("-mt" + lexical_cast<string>(tol)));
    }

    //depending on the settings in CMakeLists, parts here will be disabled

    if (vm.count("tool"))
    {
      switch(tool){

#ifdef HAVE_PEPNOVO
        case PEPNOVO:
          std::cout << "PepNovo will be used." << "\n";
          break;
#endif

#ifdef HAVE_DIRECTAG
        case DIRECTAG:
          std::cout << "DirecTag will be used." << "\n";
          break;
#endif
#ifdef HAVE_PEPNOVO
#ifdef HAVE_DIRECTAG
        case PEPNDIREC:
          std::cout << "Both PepNovo and DirecTag will be used." << "\n";
          break;
#endif
#endif
        default:
          throw std::runtime_error("Can't use selected tool, check your install.");
      }
      genOp.tool = tool;
    }


  }//try parsing block

  catch(std::exception& e)
  {
    cout << e.what() << "\n";
    cout << "check your parameters again" << std::endl;
    return 2;
  }    
  return 0;

}

