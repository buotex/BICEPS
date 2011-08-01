#include "Services.h"
#include "bicepsdefinitions.h"
namespace Pepsplice{
Services::Services()
{
	globaltime = new Timer();
	dnaAA1 = new DnaAA();
	hypergeometric1 = new Hypergeometric();
	peakdistribution_rel = new Distribution(0, 1, 4, 200); //min max columns bins
	peakdistribution_abs = new Distribution(0, 10000, 1, 10000);
	parentmassdistribution = new Distribution(0, 10000, 1, 100);
	progressreport = new Distribution(0, 1, 20, 2000);
    reportnumber = 0;
	
	//timing experiment parameters
	eachXspec = 1;
	tupspec = 1;
	
	//parameters
	//parametersparsed = false;
    changedMass = true;
	trypticendsrequired = 2; //ProteinParser.cpp how many tryptic ends must occur in peptides
	spliced = false; //SlidingWindow.cpp whether spliced peptides are enumerated on genome
	tuplebuffersize = TUPLEBUFFER; //Tuples.cpp initial buffer size
	realtuples_per_randomtuple = 1; //Tuples.cpp minimum 1
	randomtype = 0; //default
	rearrangespectra = true; //Spectra.cpp
	processorcacheL2 = 200000; //Scoring.cpp; in bytes, conservative value please
	discriminationscore = 2; // 0 = -ln(pdf), 1 = deltascore12, 2 = pvalue
	scorepval = 3;
	scoring = 7; // > 4 = hypergeometric   3 = shared peak getTotal minus average
	extraload = 0;
	adapttuplebuffer = false;
	//masstol_belowDa = 0.05;
	//masstol_aboveDa = 0.05;
	outfileprefix = "out";
	outfilenumber = "0";
	outfileparams = "";
	outfilesuffix = ".txt";
	dismisstuples = false;
	force_tb_report = false;
	scorebins = 50;
	firstscore = false;
	firstpeaksperhundred = 3; //-fh
	firstpdfcutoff = 2; //-fc
	writespecpos = false;
	bestseries = 2;
	noprefixend = false;
	nosuffixstart = false;
	//modifications = 1;
	blockedbins = 0;
	determineNK = false;
	writethspec = false;
	bestmatches = 5;
	writebestmatches = false;
	write20bestscores = false;
	resizeeveryX = 1;
	chargestate2 = false;
	minpepconfidence = 0.5;
	hotspottolerance = 10000;
	hotspectra = false;
	doWholeGenome = true;
	writespecsubsetasmgf = false;
	writespecsubsetasdta = false;
	doModifications = true;
	//modificationvalue = 2;
	doSNPaa = false;
	doSNPnt = false;

    penalty_mutation=0.0;
	penalty_max=1.51;
	penalty_ptm=PEN_PTM;
	penalty_snp=PEN_SNP; // shouldn't that be the lower limit instead of the upper? Former value: 2.0, edited by BX
	penalty_methox=PEN_METHOX;
	penalty_tryptic=PEN_TRYPTIC;
	penalty_genomic=PEN_GENOMIC;
	penalty_spliced=PEN_SPLICED;
	penalty_misscleav=PEN_MISSCLEAV;
	
	//modifications
	mod_m16 = true;
	mod_s80 = false;

	//updateParentMassTol();
	
	//parameters mass tolerance
	//masstol_below = masstol_belowDa * dnaAA1->scaling_factor; //below means that the theoretical parent mass may be that much below the measured parent mass (see Scoring::scoreTuple pmmin)
	//masstol_above = masstol_aboveDa * dnaAA1->scaling_factor; //Scoring.cpp
	//masstol = masstol_below + masstol_above; //Scoring.cpp
	//safety_margin_disc = (int)((masstol) / dnaAA1->discretization_factor + 2); //Scoring.cpp
	
	//parameters preprocessing
	learning_spectrum_intensities = true;
	specglobintpeaksper100 = 10; //-aa
	specwinintsize = 0;
	specwinintpeaks = 0;
	specwinisosize = 0;
	specwinisopeaks = 0;
	
	//parameters splicing
	gaplenmin = 10;			//minimum gap length 
	gaplenmax = 3000;		//maximum gap length
	totlenmax = 240;		//maximum (prefix + suffix)
	plenmin = 6;			//minimum prefix length (DNA)
	slenmin = 6;			//minimum suffix length	
		
	//output parameters
	outputlevel = 0;
	nexttime = 1; //must not be 0 since it is multiplied for next value
	
	//parent mass
	max_monoparentmassMH = 0; //will be adjusted after all spectra are loaded
	min_monoparentmassMH = 0;
	
	//initialize counters
	nucleotides = 0;
	triggernucleotide = 1;
	aminoacids = 0;
	spectra = 0;
	prescored = 0;
	postscored = 0;
	
	//Tuples.cpp
	tuples = 0;
	randomtuples = 0;
	realtuples = 0;
	tuples_per_second = 0;
	triggertuple = 1;
	matches = 0;
	triggermatch = 1;
	scores = 0;
	saveOld = 1; //TODO
	
	//Spectrum.cpp
	temppeaks1 = new vector<Peak>; //give the same to each spectrum for preprocessing, prevents heap fragmentation, i.e. bad memory allocation
	temppeaks2 = new vector<Peak>;
}

Services::~Services() //Changes by BX
{    
    
    delete temppeaks1;
    delete temppeaks2;
    delete globaltime;
	delete dnaAA1;
    delete hypergeometric1;
	delete peakdistribution_rel;
	delete peakdistribution_abs;
	delete parentmassdistribution;
	delete progressreport;
 //   delete hypergeometric1;

}


void Services::setMinMaxPM(double min_meas_PM, double max_meas_PM){
	//min_monoparentmassMH = min_meas_PM - masstol_above - masstol_below;	
	//max_monoparentmassMH = max_meas_PM + masstol_above + masstol_below;	
  min_monoparentmassMH = min_meas_PM * (1. - 2. * masstolfactor);
  max_monoparentmassMH = max_meas_PM * (1. + 2. * masstolfactor);
}

//void Services::updateParentMassTol(){
//	masstol_below = masstol_belowDa * dnaAA1->scaling_factor; //below means that the theoretical parent mass may be that much below the measured parent mass (see Scoring::scoreTuple pmmin)
//	masstol_above = masstol_aboveDa * dnaAA1->scaling_factor; //Scoring.cpp
//	masstol = masstol_below + masstol_above; //Scoring.cpp
//	safety_margin_disc = (int)((masstol) / dnaAA1->discretization_factor + 2); //Scoring.cpp
//}


void Services::incrementTuples(bool israndom){
	tuples++;
	if(israndom == true){
		randomtuples++;
	}else{
		realtuples++;
	}
  if (this->outputlevel > 1){
    if(tuples == triggertuple){
      suggestProgressReport();
      triggertuple += 1000;
    }
  }
}


void Services::incrementMatches(){
  matches++;
  if(matches == triggermatch){
    suggestProgressReport();
    triggermatch += 1000;
  }
}

void Services::incrementSpectra(){
  spectra++;
  if(spectra % 20 == 0) if (this->outputlevel > 2) cout << "." << flush;
  if(spectra % 1000 == 0) if (this->outputlevel > 2) cout << "\n" << spectra << "\ttime[s]:" << globaltime->timeSinceStart() << "\t" << flush;
}


void Services::invertLearnedPeakDistribution()
{
  //normalize area of raw getTotal to 1.00
  double sum0 = peakdistribution_rel->sumSeries(0);
  for(int i = 0; i < peakdistribution_rel->length; i++){
    peakdistribution_rel->setBinValue(1, i, peakdistribution_rel->getBinValue(0, i)/sum0);
  }

  //invert normalized distribution
  for(int i = 0; i < peakdistribution_rel->length; i++){
    peakdistribution_rel->setBinValue(2, i, 1/(peakdistribution_rel->getBinValue(1, i) + 0.0005));
  }

  //normalize inverted distribution
  double sum2 = peakdistribution_rel->sumSeries(2);
  for(int i = 0; i < peakdistribution_rel->length; i++){
    peakdistribution_rel->setBinValue(3, i, peakdistribution_rel->getBinValue(2, i)/sum2);
  }

}

void Services::incrementNucleotides(){
  nucleotides++;
  if(nucleotides == triggernucleotide) suggestProgressReport();
  triggernucleotide += 1000;
}

void Services::suggestProgressReport(){ 
  if (this->outputlevel > 1){
  if(globaltime->timeSinceStart() > nexttime && force_tb_report == false){
    doProgressReport();
    //writeLogFiles();
    nexttime = nexttime * 1.05 + 1;
  }
  }
}

void Services::doProgressReport(){

  cout.precision(8);
  cout << "\n";
  cout << "t:" << globaltime->timeSinceStart();
  cout << " N:" << nucleotides << " N/t:" << nucleotides/globaltime->timeSinceStart();
  cout << " AA:" << aminoacids << " AA/t:" << aminoacids/globaltime->timeSinceStart();
  cout << " T:" << tuples << " T/s:" << tuples_per_second << " TBuff:" << tuplebuffersize;
  cout << " M:" << matches << " M/t:" << matches/globaltime->timeSinceStart() << " Sc:" << scores << " Sc/t:" << scores/globaltime->timeSinceStart();
  if(scorepval == 4 || firstscore == 1) cout << " prpo: " << prescored/postscored;
  cout << flush;

  //int n = 0;
  //	reportnumber++;
  //	double gt = globaltime->timeSinceStart();
  //	int diffr = 10; //progress report number
  //	double difft = 0; //time
  //	double ratio = 0; //time
  //	
  //	progressreport->setBinValue(n++, reportnumber, gt);
  //	difft = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr));
  //	progressreport->setBinValue(n++, reportnumber, difft);	
  //	
  //	progressreport->setBinValue(n++, reportnumber, nucleotides);
  //	ratio = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr)) / difft;
  //	progressreport->setBinValue(n++, reportnumber, ratio);	
  //	progressreport->setBinValue(n++, reportnumber, nucleotides/gt);
  //	
  //	progressreport->setBinValue(n++, reportnumber, aminoacids);
  //	ratio = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr)) / difft;
  //	progressreport->setBinValue(n++, reportnumber, ratio);	
  //	progressreport->setBinValue(n++, reportnumber, aminoacids/gt);
  //	
  //	progressreport->setBinValue(n++, reportnumber, tuples);
  //	ratio = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr)) / difft;
  //	progressreport->setBinValue(n++, reportnumber, ratio);	
  //	progressreport->setBinValue(n++, reportnumber, tuples/gt);	
  //	progressreport->setBinValue(n++, reportnumber, tuplebuffersize);
  //	
  //	progressreport->setBinValue(n++, reportnumber, matches);
  //	ratio = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr)) / difft;
  //	progressreport->setBinValue(n++, reportnumber, ratio);	
  //	progressreport->setBinValue(n++, reportnumber, matches/gt);
  //	
  //	progressreport->setBinValue(n++, reportnumber, scores);
  //	ratio = (progressreport->getBinValue(n-1, reportnumber) - progressreport->getBinValue(n-1, reportnumber - diffr)) / difft;
  //	progressreport->setBinValue(n++, reportnumber, ratio);	
  //	progressreport->setBinValue(n++, reportnumber, scores/gt);
  //	
  //	progressreport->writeDistribution(getOutFileName("progressreport"));
  //cout << " rndT: " << randomtuples << " realT:" << realtuples << flush;
}

void Services::writeLogFiles(){
  //distr_specwithin->writeDistribution("distr_specwithin.txt");
  //results_desired = true;
}

string Services::getOutFileName(string x){
  string filename = "";
  filename += outfileprefix;
  filename += outfileparams;
  filename += "_";
  filename += outfilenumber;
  filename += x;
  filename += outfilesuffix;
  return filename;
}

string Services::intToString(int i){
  std::ostringstream o;
  o << i;
  return o.str();
}

double Services::string_to_double(string s)
{
  stringstream stream(s);
  double b = -1.0;
  stream >> b;
  return b;
}

void Services::parseParameters(vector<string> & args)
{
  for(unsigned int i = 0; i < args.size(); i++){
    parseParameter(args[i]);		
  }
}

void Services::parseParameter(string arg){

  string argprefix = "";
  argprefix += arg[1];
  argprefix += arg[2];

  if(argprefix == "ol")
  {
    outputlevel = (int)string_to_double(arg.substr(3));
  }


  //parse parameters
  if (this->outputlevel > 2) cout << "\nProcessing argument " << arg << flush;

  //spectrum preprocessing
  if(argprefix == "aa"){
    specglobintpeaksper100 = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\n" << specglobintpeaksper100 << " peaks per 100 Dalton are parsed on average.";
  }
  if(argprefix == "ab"){
    specwinintsize = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nspecwinintsize: " << specwinintsize << flush;				
  }
  if(argprefix == "ac"){
    specwinintpeaks = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nspecwinintpeaks: " << specwinintpeaks << flush;								
  }
  if(argprefix == "ad"){
    specwinisosize = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nspecwinisosize: " << specwinisosize << flush;								
  }
  if(argprefix == "ae"){
    specwinisopeaks = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nspecwinisopeaks: " << specwinisopeaks << flush;												
  }

  if(argprefix == "bb"){
    blockedbins = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nblockedbins: " << blockedbins << flush;												
  }

  if(argprefix == "bm"){
    bestmatches = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nbestmatches: " << bestmatches << flush;												
  }

  //bins in score distribution
  if(argprefix == "bn"){
    scorebins = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nscorebins: " << scorebins << flush;												
  }

  if(argprefix == "bs"){
    bestseries = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nbestseries: " << bestseries << flush;												
  }

  //writing best matches
  if(argprefix == "bw"){
    writebestmatches = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwritebestmatches: " << writebestmatches << flush;												
  }

  //L2 cache size limitation
  if(argprefix == "ch"){
    processorcacheL2 = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nprocessorcacheL2: " << processorcacheL2 << flush;
  }

  //for changed mass 
  if(argprefix == "cm"){
    changedMass = true;
    if (this->outputlevel > 2) cout << "\n using changed Mass for second pass" << flush;
  }


  //charge state 2
  if(argprefix == "cs"){
    chargestate2 = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nchargestate2: " << chargestate2 << flush;
  }

  //use delta score
  if(argprefix == "dl"){
    discriminationscore = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ndeltascore: " << discriminationscore << flush;
  }

  //dismiss tuples
  if(argprefix == "dt"){
    dismisstuples = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ndismisstuples: " << dismisstuples << flush;
  }

  //do calculations twice or multiple times in scoring (for timing purposes)
  if(argprefix == "el"){
    extraload = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nextraload: " << extraload << flush;
  }

  //probability density cutoff in scoring; filter to save time
  if(argprefix == "fc"){
    firstpdfcutoff = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nfirstpdfcutoff: " << firstpdfcutoff << flush;
  }

  if(argprefix == "fh"){
    firstpeaksperhundred = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nfirstpeaksperhundred: " << firstpeaksperhundred << flush;
  }

  if(argprefix == "fs"){
    firstscore = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nfirstscore: " << firstscore << flush;
  }

  if(argprefix == "gn"){
    gaplenmin = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ngaplenmin: " << gaplenmin << flush;
  }

  if(argprefix == "gx"){
    gaplenmax = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ngaplenmax: " << gaplenmax << flush;
  }

  if(argprefix == "hs"){
    hotspectra = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nhotspectra: " << hotspectra << flush;
  }

  if(argprefix == "ht"){
    hotspottolerance = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nhotspottolerance: " << hotspottolerance << flush;
  }

/* Changes by BX
  if(argprefix == "ma"){
    if (this->outputlevel > 2) cout<< "\nmass_aboveDa as read in:"<<arg.substr(3)<<flush;
    masstol_aboveDa = string_to_double(arg.substr(3));
    //updateParentMassTol();
    if (this->outputlevel > 2) cout << "\nmasstol_aboveDa: " << masstol_aboveDa << flush;

    //changes BYR: added recomputation of masstolerance after reading of the parameters 
    masstol_below = masstol_belowDa * dnaAA1->scaling_factor; //below means that the theoretical parent mass may be that much below the measured parent mass (see Scoring::scoreTuple pmmin)
    masstol_above = masstol_aboveDa * dnaAA1->scaling_factor;//Scoring.cpp
    masstol = masstol_below + masstol_above; //Scoring.
    safety_margin_disc = ((masstol) / dnaAA1->discretization_factor + 2);
    //end changes BYR
  }

  if(argprefix == "mb"){
    masstol_belowDa = string_to_double(arg.substr(3));
    //updateParentMassTol();
    if (this->outputlevel > 2) cout << "\nmasstol_belowDa: " << masstol_belowDa << flush;
    //changes BYR: added recomputation of masstolerance after reading of the parameters
    masstol_below = masstol_belowDa * dnaAA1->scaling_factor; //below means that the theoretical parent mass may be that much below the measured parent mass (see Scoring::scoreTuple pmmin)
    masstol_above = masstol_aboveDa * dnaAA1->scaling_factor;//Scoring.cpp
    masstol = masstol_below + masstol_above; //Scoring.
    safety_margin_disc = ((masstol) / dnaAA1->discretization_factor + 2);
    //end changes BYR
  }
*/

  if(argprefix == "md"){
    doModifications = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ndoModifications: " << doModifications << flush;
  }

  if(argprefix == "mo"){
    mod_m16 = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nmod_m16: " << mod_m16 << flush;
  }

  if(argprefix == "mp"){
    minpepconfidence = 1/1000 * (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nminpepconfidence: " << minpepconfidence << flush;
  }

  if(argprefix == "mt"){
    masstolfactor = 1E-6 * (double)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nmassToleranceFactor: " << masstolfactor << flush;
  }
  if(argprefix == "mu"){
    mod_s80 = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nmod_s80: " << mod_s80 << flush;
  }

  if(argprefix == "nk"){
    determineNK = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ndetermineNK: " << determineNK << flush;
  }

  if(argprefix == "np"){
    noprefixend = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nnoprefixend: " << noprefixend << flush;
  }

  if(argprefix == "ns"){
    nosuffixstart = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nnosuffixstart: " << nosuffixstart << flush;
  }


  //choose probability density function or p-value or continuity-corrected p.d.f.
  if(argprefix == "pv"){
    scorepval = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nscorepval: " << scorepval << flush;
  }

  //rearrange spectra
  if(argprefix == "rs"){
    rearrangespectra = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nrearrangespectra: " << rearrangespectra << flush;
  }

  //rearrange spectra
  if(argprefix == "rt"){
    randomtype = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nrandomtype (0=revert except last; 1=no reversion but shift by 1 except last; 2=shift by 2; etc...): " << rearrangespectra << flush;
  }

  if(argprefix == "rw"){
    writespecpos = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwritespecpos: " << writespecpos << flush;
  }

  if(argprefix == "rz"){
    resizeeveryX = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nresizeeveryX: " << resizeeveryX << flush;
  }

  //choose scoring: 3 is equivalent to simplified SEQUEST, 4 is equivalent to Sadygov/Yates
  //scoring 8 includes a page fault for the first peak, scoring 9 does not use the spectrum
  if(argprefix == "sc"){
    scoring = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nscoring: " << scoring << flush;
  }

  //enumerate SNPs based on amino acids
  if(argprefix == "sa"){
    doSNPaa = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nSNPaa: " << doSNPaa << flush;
  }

  //enumerate SNPs based on nucleotides
  if(argprefix == "sn"){
    doSNPnt = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nSNPnt: " << doSNPnt << flush;
  }

  //splicing on or off
  if(argprefix == "sp"){
    spliced = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nspliced: " << spliced << flush;
  }

  //automatic tuple buffer adaptation
  if(argprefix == "ta"){
    adapttuplebuffer = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nadapttuplebuffer: " << adapttuplebuffer << flush;
  }

  //tryptic ends required (2 = fully tryptic, 1 = semi-tryptic, 0 = non-tryptic)
  if(argprefix == "te"){
    trypticendsrequired = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ntrypticendsrequired: " << trypticendsrequired << flush;
  }

  //force progress report when tuplebuffer is readjusted
  if(argprefix == "tr"){
    force_tb_report = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nforce_tb_report: " << force_tb_report << flush;
  }

  //tuple-spectrum matching (options 0..3)
  if(argprefix == "ts"){
    tupspec = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\ntupspec: " << tupspec << flush;
  }

  if(argprefix == "wd"){
    writespecsubsetasdta = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwrite spectrum subset as dta files: " << writespecsubsetasdta << flush;
  }

  if(argprefix == "wg"){
    doWholeGenome = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwholegenome: " << doWholeGenome << flush;
  }

  if(argprefix == "ws"){
    writespecsubsetasmgf = (bool)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwrite spectrum subset as mgf files: " << writespecsubsetasmgf << flush;
  }

  if(argprefix == "wt"){
    writethspec = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\nwritethspec: " << writethspec << flush;
  }

  //use only subset of parsed spectra
  if(argprefix == "xs"){
    eachXspec = (int)string_to_double(arg.substr(3));
    if (this->outputlevel > 2) cout << "\neachXspec: " << eachXspec << flush;
  }


}

}
