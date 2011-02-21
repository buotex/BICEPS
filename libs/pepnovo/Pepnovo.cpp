
/*
 Copyright 2008, The Regents of the University of California
 All Rights Reserved
 
 Permission to use, copy, modify and distribute any part of this
 program for educational, research and non-profit purposes, without fee,
 and without a written agreement is hereby granted, provided that the
 above copyright notice, this paragraph and the following three paragraphs
 appear in all copies.
 
 Those desiring to incorporate this work into commercial
 products or use for commercial purposes should contact the Technology
 Transfer & Intellectual Property Services, University of California,
 San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
 Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.
 
 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
 IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
 OF SUCH DAMAGE.
 
 THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
 OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
 REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
 EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
 */


#include "AdvancedScoreModel.h"
#include "DeNovoRankScore.h" 
#include "FileManagement.h"
#include "DeNovoDp.h"
#include "DeNovoSolutions.h"
#include "auxfun.h"
#include "Pepnovo.h"

void print_help(const char *message) 
{
	printf("***************************************************************************\n\n%s\n",message);
  
	printf("\nPepNovo+ - de Novo peptide sequencing and\nMS-Filter - spectal quality scoring, precursor mass correction and chage determination.\n");
	printf("Release %s.\nAll rights reserved to the Regents of the University of California.\n\n",build_name.c_str());
  
	
	printf("Required arguments:\n");
	printf("-------------------\n\n");
  
	printf("-model <model name>\n\n");
	printf("-file <path to input file>  - PepNovo can analyze dta,mgf and mzXML files\n");
	printf("   OR\n");
	printf("-list <path to text file listing input files>\n\n");
  
	printf("\nOptional PepNovo arguments: \n");
	printf("----------------------------- \n"); 
	
	printf("-prm		- only print spectrum graph nodes with scores.\n");
	printf("-prm_norm   - prints spectrum graph scores after normalization and removal of negative scores.\n");
	printf("-correct_pm - finds optimal precursor mass and charge values.\n");
	printf("-use_spectrum_charge - does not correct charge.\n");
	printf("-use_spectrum_mz     - does not correct the precursor m/z value that appears in the file.\n");
	printf("-no_quality_filter   - does not remove low quality spectra.\n");
  //	printf("-output_aa_probs	 - calculates the probabilities of individual amino acids.\n");
	printf("-fragment_tolerance < 0-0.75 > - the fragment tolerance (each model has a default setting)\n");
	printf("-pm_tolerance       < 0-5.0 > - the precursor masss tolerance (each model has a default setting)\n");
	printf("-PTMs   <PTM string>    - seprated  by a colons (no spaces) e.g., M+16:S+80:N+1\n");	
	printf("-digest <NON_SPECIFIC,TRYPSIN> - default TRYPSIN\n");
	printf("-num_solutions < 1-2000 > - default 20\n");
	printf("-tag_length < 3-6> - returns peptide sequence of the specified length (only lengths 3-6 are allowed).\n");
	printf("-model_dir  < path > - directory where model files are kept (default ./Models)\n\n");
  
  
	printf("\nOptional MS-Filter arguments:\n");
	printf(	 "-----------------------------\n");
	printf("-min_filter_prob <xx=0-1.0> - filter out spectra from denovo/tag/prm run with a quality probability less than x (e.g., x=0.1)\n");
	printf("-sqs_only      - only print the filtering quality score (faster than pmcsqs_only option)\n");
	printf("-pmcsqs_only   - only output the corrected precursor mass, charge and filtering values\n");
	printf("-filter_spectra <sqs thresh> <out dir>  - outputs MGF files for spectra that have a minimal qulaity score above *thresh* (it is recomended to use a value of 0.05-0.1).");
	printf(" These MGF files will be sent to the directory given in out_dir and have a name with the prefix given in the third argument.\n");
	printf(" NOTE: this option must be used in conjuction with \"-sqs_only\" or \"-pmcsqs_only\" the latter option will also correct the m/z value and assign a charge to the spectrum.\n\n");
	printf("-pmcsqs_and_prm <min prob> - print spectrum graph nodes for spectra that have an SQS probability score of at least <min prob> (typically should have a value 0-0.2)\n\n");
  
	printf("\nTag file for InsPecT:\n");
	printf(  "---------------------\n");
	printf("PepNovo can create a tag file that can be read and used by InsPecT.\n");
	printf("To generate such tags you must use the -file option and supply the following command:\n");
	printf("-inspect_tags len1:num1:len2:num2...\n");
	printf("    For example, -inspect_tags 4:5:5:20:6:50, will generate 5 tags of length 4, 20 tags of length 5 and 50 tags of length 6.\n");
	printf("-tag_suffix <X> - This command will genrate a file with the same name as the input file, but with a suffix \"X.txt\" (default suffix tags.txt).\n\n");
  
	printf("\nRescoring InsPecT or Sequest results:\n");
	printf(  "--------------------------\n");
	printf("PepNovo can rescore an InsPecT raw results file, replacing the MQScore and delta score fields with the scores obtained with the rank modesl.");
	printf("To run in this mode you need to supply a model (\"-model\") a list of PTMs (\"-PTMs\") and the following flag:\n");
	printf("-rescore_inspect <X> <Y> - where <X> is the complete path to the original results file and <Y> is the path to the new score file.\n");
	printf("-rescore_sequest <X> <Y> - where X,Y as above. Note that the Sequest result file should have a special format, see documentation for example.\n");
  
	printf("\nPredicting peptide fragmentation:\n");
	printf(  "---------------------------------\n");
	printf("-predict_fragmentation <X> - where X is the path to a file (see readme for more details.)\n");
	printf("-num_peaks <X>  - number of peaks to predict (default X=20).\n");
  
	printf("\nCitations:\n");
	printf(  "----------\n");
	printf("- Frank, A. and Pevzner, P. \"PepNovo: De Novo Peptide Sequencing via Probabilistic Network Modeling\", Analytical Chemistry 77:964-973, 2005.\n");
	printf("- Frank, A., Tanner, S., Bafna, V. and Pevzner, P. \"Peptide sequence tags for fast database search in mass-spectrometry\", J. Proteome Res. 2005 Jul-Aug;4(4):1287-95.\n");
	printf("- Frank, A.M., Savitski, M.M., Nielsen, L.M., Zubarev, R.A., Pevzner, P.A. \"De Novo Peptide Sequencing and Identification with Precision Mass Spectrometry\", J. Proteome Res. 6:114-123, 2007.\n");
  
	printf("\nPlease send comments and bug reports to Ari Frank (arf@cs.ucsd.edu).\n\n");
	exit(1);
}

int pepnovoFunc(int argc, vector<string> & argv, string & tags, vector<float> & lowPeakMzs)

{ 
	static AdvancedScoreModel model;
  
	int i;
	char ann_file[256];
	char out_file[256];
	char input_file[256];
	char inspect_results_file[256];
	char list_file[256];
	char model_file[256];
	char initial_model[256];
	char model_dir[256];
	char PTM_string[256];
	char mgf_out_dir[256];
	char neg_spec_list[256];
	char tag_string[64];
	char tag_suffix[64];
	
	bool got_input_file=false,got_model_file=false, got_list_file=false;
	bool got_model_dir=false, got_initial_model=false, got_PTM_string = false, got_neg_spec_list=false;
	bool prm_only=false;
	bool prm_norm=false;
	bool pmcsqs_only = false;
	bool sqs_only = false;
	bool got_filter_spectra = false;
	bool pmcsqs_and_prm = false;
	bool train_flag = false;
	bool correct_pm = false;
	bool use_spectrum_charge = false;
	bool use_spectrum_mz     = false;
	bool perform_filter		 = true;
	bool output_aa_probs	 = false;
	bool make_inspect_tags   = false;
	bool make_training_fa	 = false;
	bool test_tags			 = false;
	bool got_make_ann_mgf	 = false;
	bool got_rescore_inspect = false;
	bool got_rescore_sequest = false;
	bool got_recalibrate_inspect = false;
	bool got_make_peak_examples  = false;
	bool got_predict_fragmentation = false;
  
	int start_train_idx=0;
	int end_train_idx = POS_INF;
	int specific_charge=-1;
	int specific_size=-1;
	int specific_region=-1;
  
	int specific_idx = -1;
	
	int file_start_idx =0;
	int tag_length = 0;
	int num_solutions = 20;
	int digest_type = TRYPSIN_DIGEST;
	mass_t train_tolerance;
	float min_pmcsqs_prob = -1.0;
	mass_t fragment_tolerance = -1.0;
	mass_t pm_tolerance = -1.0;
	float sqs_filter_thresh = 0.0;
	float min_filter_prob = 0.0;
	int   num_test_cases=-1;
  
	size_t num_peaks  =20;
  
	rand_seed(112233);
	strcpy(tag_suffix,"tags");
  
	// read command line arguments
	i=1;
  
	while (i<argc)
	{
    
		if (! strcmp(argv[i].c_str(),"-make_ann_mgf"))
		{
			if (++i == argc)
				print_help("Missing file ann file!");
      
			strcpy(ann_file,argv[i].c_str());	
      
			if (++i == argc)
				print_help("Missing file out file!");
      
			strcpy(out_file,argv[i].c_str());	
      
			got_make_ann_mgf=true;
		}
		else if (!strcmp(argv[i].c_str(),"-file"))
		{
			if (++i == argc)
				print_help("Missing file name!");
      
			strcpy(input_file,argv[i].c_str());
			got_input_file=true;
		}
		else
      if (!strcmp(argv[i].c_str(),"-list"))
      {
        if (++i == argc)
          print_help("Missing list name!");
        
        strcpy(list_file,argv[i].c_str());
        got_list_file=true;
      }
      else if  (!strcmp(argv[i].c_str(),"-file_start_idx"))
      {
        if (++i == argc)
          print_help("Missing file start idx!");
        
        file_start_idx = atoi(argv[i].c_str());
      }
      else if (!strcmp(argv[i].c_str(),"-model")) 
      {
        if (++i == argc)
          print_help("Missing model name!");
        
        strcpy(model_file,argv[i].c_str());
        got_model_file=true;
      }
      else if (! strcmp(argv[i].c_str(),"-model_dir"))
      {
        if (++i == argc)
          print_help("Missing model dir name!");
        
        strcpy(model_dir,argv[i].c_str());
        got_model_dir=true;
      }
      else if (! strcmp(argv[i].c_str(),"-fragment_tolerance"))
      {
        if (++i == argc)
          print_help("Missing model dir name!");
        
        fragment_tolerance = atof(argv[i].c_str());
        if (fragment_tolerance<0 || fragment_tolerance>0.75)
          print_help("Error: -fragment_toelerance should be 0-0.75\n");
      }
      else if (! strcmp(argv[i].c_str(),"-pm_tolerance"))
      {
        if (++i == argc)
          print_help("Missing model dir name!");
        
        pm_tolerance = atof(argv[i].c_str());
        if (pm_tolerance<0 || pm_tolerance>5.0)
          print_help("Error: -pm_toelerance should be 0-5.0\n");
      }
      else if  (!strcmp(argv[i].c_str(),"-num_solutions"))
      {
        if (++i == argc)
          print_help("Missing number of solutions!");
        
        num_solutions = atoi(argv[i].c_str());
        if (num_solutions<=0 || num_solutions> 2000)
          print_help("Error: -num_solutions should be 1-2000\n");
      }
      else if (!strcmp(argv[i].c_str(),"-tag_length"))
      {
        if (++i == argc)
          print_help("Missing minimum length parameter!");
        
        tag_length = atoi(argv[i].c_str());
        if (tag_length<3 || tag_length>6)
          print_help("Error: -tag_length value must be 3-6\n");
        
      }
      else if (!strcmp(argv[i].c_str(),"-digest"))
      {
        if (++i == argc)
          print_help("Missing digest type parameter : NON_SPECIFIC, TRYPSIN\n");
        
        if (! strcmp(argv[i].c_str(),"NON_SPECIFIC"))
        {
          digest_type = NON_SPECIFIC_DIGEST;
        }
        else if (! strcmp(argv[i].c_str(),"TRYPSIN"))
        {
          digest_type = TRYPSIN_DIGEST;
        }
        else
        {
          printf("Error: bad digest type: %s\n",argv[i].c_str());
          print_help("Supported digest types: NON_SPECIFIC, TRYPSIN.");
        }
      }
      else if (! strcmp(argv[i].c_str(),"-use_spectrum_charge"))
      {
        use_spectrum_charge = true;
      }
      else if (! strcmp(argv[i].c_str(),"-use_spectrum_mz"))
      {
        use_spectrum_mz = true;
      }
      else if (! strcmp(argv[i].c_str(),"-no_quality_filter"))
      {
        perform_filter = false;
      }
      else if (! strcmp(argv[i].c_str(),"-correct_pm"))
      {
        correct_pm = true;
      }
      else if (! strcmp(argv[i].c_str(),"-prm")) 
      {
        prm_only = true;
      }
      else if (! strcmp(argv[i].c_str(),"-prm_norm")) 
      {
        prm_norm = true;
        prm_only = true;
      }
      else if (! strcmp(argv[i].c_str(),"-output_aa_probs"))
      {
        output_aa_probs=true;
      }
      else if (! strcmp(argv[i].c_str(),"-pmcsqs_only"))
      {
        pmcsqs_only = true;
      }
      else if (! strcmp(argv[i].c_str(),"-sqs_only"))
      {
        sqs_only = true;
      }
      else if (! strcmp(argv[i].c_str(),"-min_filter_prob"))
      {
        if (++i == argc)
          print_help("Missing minimum probability parmater after -min_filter_prob !\n");
        
        min_filter_prob = -1.0;
        min_filter_prob = atof(argv[i].c_str());
        if (min_filter_prob<0.0 || min_filter_prob>=1.0 || argv[i].c_str()[0] != '0')
        {
          print_help("The flag -min_filter_prob should be followed by a minimal probability value [0-1.0]\n");
          exit(1);
        }
      }
      else if ( ! strcmp(argv[i].c_str(),"-filter_spectra"))
      {
        got_filter_spectra = true;
        if (++i == argc)
          print_help("Missing minimum probability parmater after -filter_spectra !\n");
        
        sqs_filter_thresh=atof(argv[i].c_str());
        
        if (sqs_filter_thresh <0 || sqs_filter_thresh>1.0)
          print_help("Error: the sqs threshold should be in the range 0-1 (recommended below 0.1)\n");
        
        if (++i == argc)
          print_help("Missing output directory for MGF files (second argument after -filter_spectra)!\n");
        
        strcpy(mgf_out_dir,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-specific_idx"))
      {
        if (++i == argc)
          print_help("Missing idx!");
        specific_idx=atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-train_model"))
      {
        train_flag = true;
        if (++i == argc)
          print_help("Missing training tolerance!");
        
        train_tolerance = atof(argv[i].c_str());
        if (train_tolerance<0.001 || train_tolerance>1.0)
          print_help("Error: training tolerance should be in the range 0.001 - 1.0\n");
      }
      else if (! strcmp(argv[i].c_str(),"-start_train_idx"))
      {
        if (++i == argc)
          print_help("Missing start_train_idx!");
        
        start_train_idx = atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-end_train_idx"))
      {
        if (++i == argc)
          print_help("end_train_idx!");
        
        end_train_idx = atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-specific_reigon_model"))
      {
        if (++i == argc)
          print_help("specific_reigon_model!");
        
        specific_charge = atoi(argv[i++].c_str());
        specific_size	= atoi(argv[i++].c_str());
        specific_region = atoi(argv[i].c_str());
        
      }
      else if (! strcmp(argv[i].c_str(),"-specific_charge"))
      {
        if (++i == argc)
          print_help("specific_charge!");
        
        specific_charge = atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-specific_size"))
      {
        if (++i == argc)
          print_help("specific_size!");
        
        specific_size = atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-initial_model"))
      {
        got_initial_model = true;
        if (++i == argc)
          print_help("Missing initial model name!");
        strcpy(initial_model,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-neg_spec_list"))
      {
        got_neg_spec_list = true;
        if (++i == argc)
          print_help("Missing neg spec list!");
        strcpy(neg_spec_list,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-PTMs"))
      {
        got_PTM_string = true;
        if (++i == argc)
          print_help("Missing PTM list!");
        strcpy(PTM_string,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-inspect_tags"))
      {
        make_inspect_tags=true;
        if (++i == argc)
          print_help("inspect_tags!");
        
        strcpy(tag_string,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-rescore_inspect"))
      {
        got_rescore_inspect = true;
        if (++i == argc)
          print_help("Missing results file!");
        
        strcpy(inspect_results_file,argv[i].c_str());
        
        if (++i == argc)
          print_help("Missing new results file!");
        
        strcpy(out_file,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-rescore_sequest"))
      {
        got_rescore_sequest = true;
        if (++i == argc)
          print_help("Missing results file!");
        
        strcpy(inspect_results_file,argv[i].c_str());
        
        if (++i == argc)
          print_help("Missing new results file!");
        
        strcpy(out_file,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-recalibrate_inspect"))
      {
        got_recalibrate_inspect = true;
        if (++i == argc)
          print_help("Missing results file!");
        
        strcpy(inspect_results_file,argv[i].c_str());
        
        if (++i == argc)
          print_help("Missing new results file!");
        
        strcpy(out_file,argv[i].c_str()); 		
      }
      else if ( ! strcmp(argv[i].c_str(),"-make_peak_examples"))
      {
        got_make_peak_examples=true;
      }
      else if (! strcmp(argv[i].c_str(),"-make_training_fa"))
      {
        make_training_fa=true;
      }
      else if (! strcmp(argv[i].c_str(),"-test_tags"))
      {
        test_tags=true;
        if (++i == argc)
          print_help("test_tags!");
        
        strcpy(tag_string,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-num_test_cases"))
      {
        if (++i == argc)
          print_help("num_test_cases!");
        
        num_test_cases = atoi(argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-tag_suffix"))
      {
        if (++i == argc)
          print_help("tag suffix!");
        strcpy(tag_suffix,argv[i].c_str());
      }
      else if (! strcmp(argv[i].c_str(),"-predict_fragmentation"))
      {
        if (++i == argc)
          print_help("predict fragmentation, missing file path!");
        strcpy(input_file,argv[i].c_str());
        got_predict_fragmentation = true;
      }
      else if (! strcmp(argv[i].c_str(),"-num_peaks"))
      {
        if (++i == argc)
          print_help("predict fragmentation, missing number of peaks!");
        num_peaks = static_cast<size_t>(atoi(argv[i].c_str()));
        if (num_peaks > 1000)
          print_help("invalid number of peaks for prediction (should be < 100)");
      }
      else
      {
        printf("**********************************************************\n");
        printf("\nError: Unkown command line option: %s\n\n",argv[i].c_str());
        print_help("");
        exit(0); 
      }
		i++;
	}
  
  
	if (! got_model_file) 
		print_help("Error: Missing model name!");
  
  
	if (!got_input_file && ! got_list_file && ! got_predict_fragmentation)
		print_help("Error: missing input file (either -file or -list must be used).");
  
	Config *config = model.get_config();
  
	if (got_model_dir)
	{
		config->set_resource_dir(string(model_dir));
	}
  
	
  
	//////////////////////////////////////////////////////////////////
	// Model Training
	if (train_flag)
	{	
		if (got_initial_model)
		{
			model.read_model(initial_model);
			if (got_PTM_string)
				config->apply_selected_PTMs(PTM_string);
			model.read_rank_models(model_file,true);
		}
		else
		{
			config->init_with_defaults();
			config->set_tolerance(train_tolerance);
			config->set_digest_type(digest_type);
			if (got_PTM_string)
				config->apply_selected_PTMs(PTM_string);
		}
    
		model.set_model_name(string(model_file));
    
		FileManager fm; 
    
		if (! got_list_file)
		{
			if (got_input_file)
			{
				fm.init_from_mgf(config,input_file);
			}
			else
			{
				printf("Must supply a list of annotated spectra for training!\n");
				exit(0);
			}
		}
		else
			fm.init_from_list_file(config,list_file);
    
    
		model.train_model_in_stages(model_file,fm,
                                train_tolerance, 
                                start_train_idx, 
                                end_train_idx,
                                specific_charge, 
                                specific_size, 
                                specific_region,
                                (got_neg_spec_list ? neg_spec_list : NULL));
    
		model.write_model();
		exit(0);
	}
	
	///////////////////////////////////////////////////////////////////
	// Model initializing (running some sort of de novo, need a model)
	// 
	const time_t start_time = time(NULL);
  
	if (model.get_rank_model_ptr(0) == 0)
  {
    
    //cout << "PepNovo+, Build " << build_name << endl;
    //cout << "Copyright 2008, The Regents of the University of California. All Rights Reserved." << endl;
    //cout << "Created by Ari Frank (arf@cs.ucsd.edu)" << endl << endl;
    cout << "Initializing models (this might take a few seconds)... " << flush;
    model.read_model(model_file,true); 
    if (got_PTM_string)
      config->apply_selected_PTMs(PTM_string);
    model.read_rank_models(model_file,true);
    cout << "Done." << endl;
  }
	config = model.get_config();
	config->set_digest_type(digest_type);
  
	if (fragment_tolerance>0)
		config->set_tolerance(fragment_tolerance);
  
	if (pm_tolerance>0)
		config->set_pm_tolerance(pm_tolerance);
  
	if (correct_pm)
		config->set_need_to_estimate_pm(1);
  
	if (use_spectrum_mz)
		config->set_use_spectrum_mz(1);
  
	if (use_spectrum_charge)
		config->set_use_spectrum_charge(1);
  
	if (! perform_filter)
		config->set_filter_flag(0);
  
	if (config->get_pm_tolerance()<0.1)
		config->set_need_to_estimate_pm(0);
  
  // BX    
  //	cout << setprecision(4) << fixed;
  //	cout << "Fragment tolerance : " << config->get_tolerance() << endl;
  //	cout << "PM tolerance       : " << config->get_pm_tolerance() << endl;
  //	cout << "PTMs considered    : " ;
  //	if (got_PTM_string)
  //	{
  //		cout << PTM_string << endl;
  //	}
  //	else
  //	{
  //		cout << "None" << endl;
  //	}
  //	
  
  //less output, BX    
  //	if (fabs(config->get_aa2mass()[Cys]-103.0)<1) 
  //	{
  //		cout << endl <<"*** Warning: searching with unmodified cystine, usually the PTM C+57 should be included ***" << endl << endl;
  //	}
  //	cout << endl;
  
  
	if (got_predict_fragmentation)
	{
		predict_fragmentation(&model, input_file, num_peaks);
		exit(0);
	}
  
	///////////////////////////////////////////////////////////////////
	// Training fa
	if (make_training_fa)
	{
		make_denovo_training_fa(model,input_file);
		exit(0);
	}
  
	///////////////////////////////////////////////////////////////////
	// Inspect tags
  
	if (make_inspect_tags)
	{
		create_tag_file_for_inspect(model,input_file,tag_string,tag_suffix);
		exit(0);
	}
  
	if (test_tags)
	{
		benchmark_tags(model,list_file,tag_string,num_test_cases);
		exit(0);
	}
  
  
	////////////////////////////////////////////////////////////////////
	// Rescore InsPecT or Sequest
	if (got_rescore_inspect)
	{
		DeNovoRankScorer *db_score = (DeNovoRankScorer *)model.get_rank_model_ptr(0);
		db_score->rescore_inspect_results(input_file,inspect_results_file,out_file);
		exit(0);
	}
  
	if (got_rescore_sequest)
	{
		DeNovoRankScorer *db_score = (DeNovoRankScorer *)model.get_rank_model_ptr(0);
		db_score->rescore_sequest_results(input_file,inspect_results_file,out_file);
		exit(0);
	}
  
	if (got_recalibrate_inspect)
	{
		cout << "Recalibrating delta scores in " << input_file << endl;
		DeNovoRankScorer *db_score = (DeNovoRankScorer *)model.get_rank_model_ptr(0);
		db_score->recalibrate_inspect_delta_scores(input_file,inspect_results_file,out_file);
		exit(0);
	}
  
	if (got_make_peak_examples)
	{
		cout << "Making peak examples " << input_file << endl;
		DeNovoRankScorer *db_score = (DeNovoRankScorer *)model.get_rank_model_ptr(0);
		db_score->make_peak_table_examples(input_file);
		exit(0);
	}
  
  
  
	///////////////////////////////////////////////////////////////////
	// Make input file list
	vector<string> list_vector;
	if (got_list_file)
	{
		read_paths_into_list(list_file, list_vector);
	}
	else
		list_vector.push_back(input_file);
  
	int correct_benchmark =0;
	int total_benchmark =0;
	int counter=0;
  
  
	if (got_filter_spectra ||  pmcsqs_only || sqs_only)
	{
		PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model.get_pmcsqs_ptr();
		if (! pmcsqs || ! pmcsqs->get_ind_initialized_pmcr() || ! pmcsqs->get_ind_initialized_sqs())
		{
			cout << "Error: no parent mass correction or quality score for this model!" << endl;
			exit(1);
		}
	}
  
  
  
  
	///////////////////////////////////////////////////////////////////
	// FILTER SPECTRA
	if (got_filter_spectra)
	{
		int num_written =0;
		int num_read = 0;
		PMCSQS_Scorer *pmcsqs = (PMCSQS_Scorer *)model.get_pmcsqs_ptr();
    
		pmcsqs->output_filtered_spectra_to_mgfs(config, list_vector, mgf_out_dir, sqs_filter_thresh, num_written, num_read);
		
		time_t curr_time = time(NULL);
		double elapsed_time = (curr_time - start_time);
		cout << "Processed " << list_vector.size() << " (" << num_read << " spectra)." << endl;
		cout << "Wrote " << num_written << " spectra to mgfs in " << mgf_out_dir << endl;
		cout << "Elapsed time " << fixed << elapsed_time << " seconds." << endl;
		return 0;
	}
  
	//////////////////////////////////////////////////////////////////
	// PRM
	if (prm_only)
	{
		
    
		perform_prm_on_list_of_files(model, list_vector, min_filter_prob, file_start_idx, prm_norm);
    //	prm_benchmark(model, list_vector, min_pmcsqs_prob, file_start_idx);
    
    //	FileManager fm;
    //	fm.init_from_list(config,list_vector);
    //	model.learn_prm_normalizer_values(fm);
    //	model.write_prm_normalizer_values();
		return 0;
	}
  
	//////////////////////////////////////////////////////////////////
	// PMCSQS
	if (pmcsqs_only)
	{
		perform_pmcsqs_on_list_of_files(model, list_vector, file_start_idx);
		return 0;
	}
  
	//////////////////////////////////////////////////////////////////
	// SQS
	if (sqs_only)
	{
		perform_sqs_on_list_of_files(model, list_vector, file_start_idx);
		return 0;
	}  
	
	//////////////////////////////////////////////////////////////////
	// DENOVO AND TAGS
  
	if (tag_length<=0)
	{
		perform_denovo_on_list_of_files(model, list_vector, 
                                    file_start_idx, num_solutions, 7, 16, false, min_filter_prob, output_aa_probs, tags, lowPeakMzs, cout);
	}
	else
	{
		perform_tags_on_list_of_files(model,list_vector,file_start_idx,num_solutions,tag_length,
                                  false, min_filter_prob, output_aa_probs, tags, lowPeakMzs, cout);	
	}
	
  
  
  //    for (int i =0; i < 3; i++) //BX
  //    {
  //        if( model.get_rank_model_ptr(i))
  //            delete (DeNovoRankScorer*) model.get_rank_model_ptr(i);
  //    }
  //    for (int i =0; i < 10; i++) //BX
  //    {
  //        if( model.get_rank_tag_model_ptr(i))
  //            delete (DeNovoRankScorer*) model.get_rank_tag_model_ptr(i);
  //    }
  
  
  
  
  
  
  
	return 0;
}



