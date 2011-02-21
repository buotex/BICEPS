#include "Results.h"

namespace Pepsplice{
    
ostream& operator << (ostream& os, const PepspliceResult& pep)
{
    //std::cout << "OriginalID: " << pep.seqid << '\n';
   // os << "OriginalIDs:"<< '\n';
//    
//    list<unsigned int>::const_iterator it;
//    for (it = pep.seqIds.begin(); it != pep.seqIds.end();)
//    {
//        os << "FastaOccurence" << (*it) <<':' << '\t';
//        ++it;
//        os << (*it);
//        ++it;
//        os << '-' << (*it) << '\n';
//        ++it;
//    }
    
    os << "Sequence: " << pep.Sequence<< '\n';
	os << "OrigSequence: " << pep.OrigSequence << '\n';
    os << "penalty: " << pep.penalty << '\n';
    os << "bic: " << pep.bic << '\n';
	os << "score: " << pep.score << std::endl;
    
    return os;	
}


Results::Results(Services *se0)
{
	se1 = se0;
	resultspattern1 = new ResultsPattern(se1);
	sf = se1->dnaAA1->scaling_factor;
}

Results::~Results()
{
delete resultspattern1;
}

void Results::writeFiles(vector<BestMatches*> b)
{
	//write peptide and protein files according to requested discrimination score
	
	Distribution *deltascore_performances = new Distribution(0, se1->bestmatches+1, 9, se1->bestmatches+1);	
	
	for(int i = 0; i <= se1->bestmatches; i++){
		int outdegree = 0;
		if(se1->discriminationscore == i || se1->discriminationscore == 99) outdegree = 2;

		//calculating discrimination for -dl" << i << " *** 0=-ln(pdf), 1=-ln(cdf), 2=ln(pdf2)-ln(pdf1), 3=ln(pdf3)-ln(pdf1), 4=ln(pdf4)-ln(pdf1), ...";
		bool writescores = false;
		if(outdegree >= 2) writescores = true;
		calcScoreDistribution(b, i, writescores); //write distribution only if outdegree high enough
		sort(b.begin(), b.end(), PeptideConfidenceComparator());  //greater<Spectrum*>());
		calcLocalRecallPrecision(b, i, deltascore_performances);
		if(outdegree >= 2){

			writeResultsPeptideWise(b, i); //TODO this here is a problem!

			writePepXMLLite(b, i);
			writeResultsProteinWise(b, i, false); //details == false
			writeResultsProteinWise(b, i, true);  //details == true
			resultspattern1->calcMassDeviationAndSeriesPattern(b, i);
		}
		if(outdegree >= 2 && se1->writethspec) resultspattern1->writeBestTheoreticalSpectra(b);
	}
	
	deltascore_performances->writeDistribution(se1->getOutFileName("_deltascore_performances"));
    
    delete deltascore_performances;
    
	//delete d1; causes segmentation fault from time to time
	//delete[] b;

}


void Results::calcScoreDistribution(vector<BestMatches*> b, int discriminate_by, bool writedistribution){
	
	Distribution *d1 = new Distribution(0, 25, 50, se1->scorebins); //min, max, columns, bins
	//string scorefilesuffix = "";
	
	double confidence_imprecision = (double)1/b.size();
			
	//fill in naive distribution of good and random identifications
	for(unsigned int i = 0; i < b.size(); i++){
		
		//calculate p-values
		int N = b[i]->ownerSpectrum->N;
		int K = b[i]->ownerSpectrum->K;
		int n = 0;
		int spc = 0;
		//cout << "\n(long)b[i]->bm[0]" << (long)b[i]->bm[0];
		if(b[i]->bm.size() > 0){
			//cout << "done";	
			n = (*b[i]).bm[0]->tuple->length * 2 - 2;
			spc = (*b[i]).bm[0]->sharedpeakcount;
		}else{
			n = 0;
			spc = 0;
		}
		b[i]->pvalue = (double)-(se1->hypergeometric1->log10HypergeometricPvalue(N, K, n, spc));
			
		//choose discrimination score (i.e. sorting order from best to worst)
		if(discriminate_by == 0){
			b[i]->discriminationscore = b[i]->bestscore;
		}else if(discriminate_by == 1){
			b[i]->discriminationscore = b[i]->pvalue;
		}else if(discriminate_by == 2){
			b[i]->discriminationscore = b[i]->deltascore12;
		}else if(discriminate_by > 2){
			b[i]->discriminationscore = b[i]->bestscore - b[i]->bestscores[discriminate_by];
		}
		
		//choose chromosome_reverse or forward column in table/distribution
		if(b[i]->random == true){
			b[i]->confidence_bin = d1->addElement(0, b[i]->discriminationscore); //add element, get back bin index
		}else{
			b[i]->confidence_bin = d1->addElement(1, b[i]->discriminationscore); //dito
		}
	}
	
	
	//*****************************************
	//calculate ROC, recall/precision
	double bad = 0;
	double good = 0;
	double sumbad = d1->sumSeries(0);
	double sumgood = d1->sumSeries(1);
	double tp = 0;
	double fp = 0;
	double tn = 0;
	double fn = 0;
	double fraction = 0;
	double sensitivity = 0; //sensitivity = recall = TP/(TP + FN) "yield"
	double specificity = 0; //specificity = TN/(TN + FP) != precision
	double recall = 0;
	double precision = 0;   //precision = TP/(TP + FP) "purity"

	double shiftS = 0;
	double badS = 0;
	double goodS = 0;
	double confidenceS = 0;

	int col = 0; //column index
	
	//getTotal downwards!
	for(int i = d1->length-1; i >= 0; i--){
		bad = d1->getBinValue(0, i);
		good = d1->getBinValue(1, i);
		col = 2;
		
		fp += bad;
		tp += good;
		fn = sumgood - tp;
		tn = sumbad - fp;
		
		fraction = (tp + fp)/(tp + fp + tn + fn);
		sensitivity = tp/(tp+fn);
		specificity = tn/(tn+fp);
		recall = sensitivity;
		precision = tp/(tp+fp);
		
		shiftS = bad;
		if(shiftS > good) shiftS = good;
		badS = bad + shiftS;
		goodS = good - shiftS;		
		confidenceS = goodS/(goodS+badS) - confidence_imprecision;
		if(confidenceS < 0) confidenceS = 0;
		if(confidenceS > 1) confidenceS = 1;
			
		//TN, FP, FN, TP
		col++;
		d1->setBinValue(col++, i, tn);
		d1->setBinValue(col++, i, fp);
		d1->setBinValue(col++, i, fn);
		d1->setBinValue(col++, i, tp);

		//recall precision
		col++;
		d1->setBinValue(col++, i, recall);
		d1->setBinValue(col++, i, precision);

		//ROC
		col++;
		d1->setBinValue(col++, i, 1-specificity);
		d1->setBinValue(col++, i, sensitivity);
		
		//fraction, confidence
		col++;
		d1->setBinValue(col++, i, bad);
		d1->setBinValue(col++, i, good);
		d1->setBinValue(col++, i, shiftS);
		d1->setBinValue(col++, i, badS);		
		d1->setBinValue(col++, i, goodS);
		d1->setBinValue(col++, i, fraction);
		d1->setBinValue(col++, i, confidenceS);

	}
	
	string scorefile = se1->getOutFileName("_dl" + se1->intToString(discriminate_by) + "_scores");
	if(writedistribution == true) d1->writeDistribution(scorefile);
		
	//set peptide confidences
	for(int i = 0; i < b.size(); i++){
		b[i]->confidence = d1->getBinValue(20, b[i]->confidence_bin);
	}
	delete d1;	
}


/*********************************************************************
 * WRITE PEPTIDE FILE
 * ******************************************************************/

void Results::writeResultsPeptideWise(vector<BestMatches*> b, int dl){

	ofstream outfile;
	string peptidefile = se1->getOutFileName("_dl" + se1->intToString(dl) + "_peptides");
	outfile.open(peptidefile.c_str());
	outfile.precision(10);

	//HEADERS
		outfile << "\n";
		
		//values of best matches
		outfile << "\treverted";
		outfile << "\ttp";
		outfile << "\tfp";
		outfile << "\tfraction";
		outfile << "\trecall";
		outfile << "\tprecision";
		outfile << "\tconfidence";
		outfile << "\tconfidence_bin";
		outfile << "\tpvalue";
		outfile << "\tbestscore";
		outfile << "\tdiscriminationscore";
		outfile << "\tdeltascore1-2";
		if(se1->write20bestscores == true) for(int j = 2; j < 20; j++){outfile << "\tscore" << j+1;}
		outfile << "\tscored";
		outfile << "\tN";
		outfile << "\tK";
		outfile << "\tEmpiricalN";
		outfile << "\tEmpiricalK";
		outfile << "\tn"; //n
		outfile << "\tk"; //k
		outfile << "\tHITS_total";
		outfile << "\tHITS_protein";
		outfile << "\tHITS_whg";
		outfile << "\tHITS_spliced";
		outfile << "\tthPM";
		
		//parent mass and diff
		outfile << "\texpPM";
		outfile << "\tdiffPM";
		
		//file names
		outfile << "\tspectrumname";
		outfile << "\tmgfname";
		
		outfile << "\treverted";
		outfile << "\ttag";
		outfile << "\tscore";
		outfile << "\taasequence"; //n
		//outfile << "\tpenalty";//n
		outfile << "\tmodifications"; //n
		//changes BYR added original sequence
		outfile << "\taaseqoriginal"; //n
		outfile << "\ttrypticstart";
		outfile << "\ttrypticend";
		outfile << "\tmissedcleavages";
		outfile << "\tmodifications";
		outfile << "\toxidizedM";
		outfile << "\tdna";
		outfile << "\tspliced";

		outfile << "\tpenalty";
		outfile << "\tisSNP";
		outfile << "\tSNPpos_0..";
		outfile << "\tSNPoldAA";
		outfile << "\tSNPnewAA";
		//outfile << "\tpenalty";
		outfile << "\tlength";
		outfile << "\tntsequence";
		outfile << "\tDNA_read_dir";
		outfile << "\tprot";
		outfile << "\tthPM";
		outfile << "\tPS";
		outfile << "\tPlen";
		outfile << "\tPE";
		outfile << "\tSS";		
		outfile << "\tSlen";		
		outfile << "\tSE";		
		
	for(int i = 0; i < b.size(); i++){
		outfile << "\n";
		
		//values of best matches
		outfile << "\t" << b[i]->random;
		outfile << "\t" << b[i]->tp;
		outfile << "\t" << b[i]->fp;
		outfile << "\t" << b[i]->fraction;
		outfile << "\t" << b[i]->recall;
		outfile << "\t" << b[i]->precision;
		outfile << "\t" << b[i]->confidence;
		outfile << "\t" << b[i]->confidence_bin;
		outfile << "\t" << b[i]->pvalue;
		outfile << "\t" << b[i]->bestscore; //score
		outfile << "\t" << b[i]->discriminationscore;
		outfile << "\t" << b[i]->deltascore12;
		if(se1->write20bestscores == true) for(int j = 2; j < 20; j++){outfile << "\t" << b[i]->bestscores[j];}
        
		outfile << "\t" << b[i]->ownerSpectrum->scored;
		outfile << "\t" << b[i]->ownerSpectrum->N;
		outfile << "\t" << b[i]->ownerSpectrum->K;
		outfile << "\t" << b[i]->ownerSpectrum->getEmpiricalN();
		outfile << "\t" << b[i]->ownerSpectrum->getEmpiricalK();
		
		if(b[i]->bm.size() > 0){
			outfile << "\t" << (b[i]->bm[0]->tuple->length - 1) * 2; // n
			outfile << "\t" << b[i]->bm[0]->sharedpeakcount; //k
		}else{
			//outfile << "\t" << 0;
			//outfile << "\t" << 0;
		}		
		outfile << "\t" << b[i]->tot_best_hits;
		outfile << "\t" << b[i]->nb_prot;
		outfile << "\t" << b[i]->nb_chromunspliced;
		outfile << "\t" << b[i]->nb_chromspliced;
		outfile << "\t" << b[i]->thpm / sf;
		
		//parent mass and diff
		outfile << "\t" << b[i]->ownerSpectrum->parentmassMH / sf;
		outfile << "\t" << (b[i]->ownerSpectrum->parentmassMH - b[i]->thpm) / sf;
		

		//file names
		outfile << "\t" << b[i]->ownerSpectrum->getFileName();
		outfile << "\t" << b[i]->ownerSpectrum->getMGFName();
		
		outfile << "\t" << b[i]->random;
		
		//enumerate matches
		string ntseq = "";
		int max = 100;
		if(se1->writebestmatches == false) max = 1; 
		for(int j = 0; j < max && j < b[i]->show_hits && j < b[i]->ownerSpectrum->bestmatches1->bm.size(); j++){
			Match *mt = &(*(b[i])->ownerSpectrum->bestmatches1->bm[j]);
			//b[i]->ownerSpectrum->bestmatches1->bm[j]->tuple->getNTSeqFromChromosome();
			
			outfile << "\t<" << j + 1 << ">"; //e.g. <1> <2> <3>
			outfile << "\t" << mt->score;
			outfile << "\t" << se1->dnaAA1->charseq_to_lowercaseseq(mt->tuple->getAASeq());
			outfile << "\t" << se1->dnaAA1->charseq_to_modnumseq(mt->tuple->getAASeq(), mt->tuple->charNterm, mt->tuple->charCterm);
			outfile << "\t" << (int)mt->tuple->trypticstart;
			outfile << "\t" << (int)mt->tuple->trypticend;
			outfile << "\t" << (int)mt->tuple->getMissedCleavages();
			outfile << "\t" << (int)mt->tuple->modifications;
			outfile << "\t" << (int)mt->tuple->oxidizedmethionines;
			outfile << "\t" << (bool)mt->tuple->isDNA();
			outfile << "\t" << (bool)mt->tuple->isSpliced();
			
			mt->tuple->isSNPok("Results.cpp, line 331");
			
			outfile << "\t" << mt->tuple->penalty;
			outfile << "\t" << (bool)mt->tuple->isSNP;
			outfile << "\t" << mt->tuple->SNPpos;
			outfile << "\t" << mt->tuple->SNPoldAA;
			outfile << "\t" << mt->tuple->SNPnewAA;			
			//outfile << "\t" << mt->tuple->penalty;
			outfile << "\t" << (int)mt->tuple->length;
			outfile << "\t" << mt->tuple->getNTSeq();			
			//outfile << "\t" << b[i]->ownerSpectrum->bestmatches1->bm[j]->tuple->getExcessNTSeq();
			
			if(mt->tuple->chromosome_reverse == 0){
				outfile << "\t" << "fwd";
			}else{
				outfile << "\t" << "rev";
			}
			
			outfile << "\t";
			int stringlen = 15;
			if(j == 0) stringlen = 100;
			if(mt->tuple->chromosome != NULL)	outfile << mt->tuple->chromosome->description_line.substr(0, stringlen);
			if(mt->tuple->protein != NULL) outfile << mt->tuple->protein->protid.substr(0, stringlen);

			outfile << "\t" << mt->tuple->thparentmassMH / sf;
			outfile << "\t" << mt->tuple->getPS();
			outfile << "\t" << mt->tuple->getPrefixLength();
			outfile << "\t" << mt->tuple->getPE();
			outfile << "\t" << mt->tuple->getSS();
			outfile << "\t" << mt->tuple->getSuffixLength();
			outfile << "\t" << mt->tuple->getSE();
		}//next j
        
	} //next i
	outfile.close();
}

/*********************************************************************
 * WRITE PEPXML-LIGHT FILE THAT CAN BE READ BY NOGUI-QUALSCORE
 * ******************************************************************/

void Results::writePepXMLLite(vector<BestMatches*> b, int dl){

	ofstream outfile;
	string peptidefile = se1->getOutFileName("_dl" + se1->intToString(dl) + "_pepxml_light");
	outfile.open(peptidefile.c_str());
	outfile.precision(10);

	outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	outfile << "<?xml-stylesheet type=\"text/xsl\" ?>\n";
	outfile << "<msms_run_summary raw_data=\".mzXML\">\n";
	outfile << "\n";
	
	//OUTPUT SPECTRUM VALUES
	
	for(int i = 0; i < b.size(); i++){
	
		//strip .dta off 
		string filedta = b[i]->ownerSpectrum->getFileName();
		int len = filedta.length();
		len = len - 4;
		if(len < 0) len = 0;
		string file = filedta.substr(0, len);

		//write
		outfile << "\n";				
		outfile << "<spectrum_query spectrum=\"" << file << "\" precursor_neutral_mass=\"" << b[i]->ownerSpectrum->parentmassMH / sf << "\" assumed_charge=\"2\">\n";
		outfile << "<search_result>\n";
		outfile << "<search_hit>\n";
		outfile << "<analysis_result analysis=\"peptideprophet\">\n";
		outfile << "<peptideprophet_result probability=\"" << b[i]->confidence << "\">\n";
		outfile << "</peptideprophet_result>\n";
		outfile << "</analysis_result>\n";
		outfile << "</search_hit>\n";
		outfile << "</search_result>\n";
		outfile << "</spectrum_query>\n";

	} //next i
	
	outfile << "\n";
	outfile << "</msms_run_summary>\n";

	outfile.close();
}


/*********************************************************************
 * WRITE PROTEIN FILE
 * ******************************************************************/

void Results::writeResultsProteinWise(vector<BestMatches*> b, int dl, bool details){
	
	ofstream outfile;
	string proteinfile;
	if(details == true){
		proteinfile = se1->getOutFileName("_dl" + se1->intToString(dl) + "_proteindetails");
	}else{
		proteinfile = se1->getOutFileName("_dl" + se1->intToString(dl) + "_proteinsummary");		
	}
	outfile.open(proteinfile.c_str());
	outfile.precision(10);

	vector<Match*> matchindex;
	vector<Protein*> proteinindex;
	
	Match *mt;
	BestMatches *bmt;
	Spectrum *sp;
	
	//initialize protein confidences
	for(unsigned int i = 0; i < b.size(); i++){
		bmt = b[i];
		
		for(unsigned int j = 0; j < bmt->bm.size(); j++){
			if(bmt->bm[j]->tuple->protein != NULL){
				bmt->bm[j]->tuple->protein->initializeConfidence();
			}
		}
	}

	//calculate protein confidences
	for(unsigned int i = 0; i < b.size(); i++){
		bmt = b[i];
		double pepconfidence = bmt->confidence/bmt->nb_prot;
		for(unsigned int j = 0; j < bmt->tot_best_hits && j < bmt->bm.size(); j++){
			mt = (bmt->bm[j]);
			if(mt->tuple->protein != NULL && pepconfidence > se1->minpepconfidence){
				//add protein to index if it has not yet happened
				if(mt->tuple->protein->matches.size() == 0) proteinindex.push_back((*bmt).bm[j]->tuple->protein);
				//add match to index
				matchindex.push_back(mt);
				//adjust protein confidence (spectrum_confidence divided by competing proteins)
				mt->tuple->protein->adjustProteinConfidence(pepconfidence, matchindex.size()-1);
			}
		}
	}


	sort(proteinindex.begin(), proteinindex.end(), ProteinConfidenceComparator());  //greater<Spectrum*>());

	
	//write proteins
	//identifiedproteins = 0;
	for(unsigned int i = 0; i < proteinindex.size(); i++){
		if(proteinindex[i]->protid != "RANDOM"){
			
			double proteinconfidence = proteinindex[i]->getProteinConfidence();
			//if(proteinconfidence > 0.9) identifiedproteins++;
			
			outfile << "\n" << i;
			outfile << "\t" << proteinconfidence;
			outfile << "\t" << proteinindex[i]->matches.size();
			outfile << "\t " << proteinindex[i]->protid.substr(0, 100);
			
			if(details == true){
				for(unsigned int j = 0; j < proteinindex[i]->matches.size(); j++){
							
					mt = matchindex[proteinindex[i]->matches[j]];
					sp = mt->spectrum;
					bmt = sp->bestmatches1;
						
					//BEST MATCHES SPECIFIC			
					//confidence
					outfile << "\n" << i;
					outfile << "\t\t";
					outfile << "\t" << bmt->confidence;
					outfile << "\t" << bmt->confidence_bin;
							
					//values of best matches
					outfile << "\t" << bmt->bestscore; //score
					outfile << "\t" << bmt->deltascore12;
					outfile << "\t" << bmt->tot_best_hits;
					outfile << "\t" << bmt->nb_prot;
					outfile << "\t" << bmt->nb_chromunspliced;
					outfile << "\t" << bmt->nb_chromspliced;
					outfile << "\t" << bmt->thpm / sf;
					
					//SPECTRUM SPECIFIC
					//parent mass and diff
					outfile << "\t" << sp->parentmassMH / sf;
					outfile << "\t" << (sp->parentmassMH - bmt->thpm) / sf;
				
					//file names
					outfile << "\t" << sp->getFileName();
					outfile << "\t" << sp->getMGFName();	
					
					//MATCH SPECIFIC
					//sequences
					outfile << "\t" << mt->tuple->getAASeq();
					outfile << "\t" << mt->tuple->getNTSeq();
					//outfile << "\t" << mt->tuple->getExcessNTSeq();
					
					if(mt->tuple->chromosome_reverse == 0){outfile << "\t" << "fwd";}else{outfile << "\t" << "rev";}
					
					//ID
					outfile << "\t";
					int stringlen = 15;
					if(j == 0) stringlen = 100;
					if(mt->tuple->chromosome != NULL) outfile << mt->tuple->chromosome->description_line.substr(0, stringlen);
					if(mt->tuple->protein != NULL) outfile << mt->tuple->protein->protid.substr(0, stringlen);
		
					outfile << "\t" << mt->tuple->thparentmassMH / sf;
					//outfile << "\t" << mt->tuple->prefixstart;
					//outfile << "\t" << mt->tuple->prefixend;
					//outfile << "\t" << mt->tuple->suffixstart;
					//outfile << "\t" << mt->tuple->suffixend;
					outfile << "\t" << mt->tuple->getPS();
					outfile << "\t" << mt->tuple->getPE();
					outfile << "\t" << mt->tuple->getSS();
					outfile << "\t" << mt->tuple->getSE();
				}
			}
		}
	}
	//if(se1->outputlevel > 2) cout << "w" << flush;

	outfile.close();
}

void Results::calcLocalRecallPrecision(vector<BestMatches*> b, int delta_i, Distribution *deltascore_performances){
	
	int bsize = b.size();
	double tp = 0;
	double fp = 0;
	double fraction = 0;
	double precision = 0;
	double allpos = 0;
	double recall = 0;
	
	for(int i = 0; i < bsize; i++){
		if(b[i]->random == false) allpos++;
	}
	
	for(int i = 0; i < bsize; i++){
				 
		if(b[i]->random == false){
			tp++;
		}else{
			fp++;
		}
		
		//calculate precision, fraction, recall
		fraction = (tp + fp) / bsize;		
		recall = tp / allpos; //=tp/(tp+fn)
		precision = tp / (tp + fp);
		
		if(precision >= 0.9995){
			deltascore_performances->setBinValue(0, delta_i, fraction);
			deltascore_performances->setBinValue(3, delta_i, recall);
			deltascore_performances->setBinValue(6, delta_i, precision);
							//cout << "\nPrecision:" << lastprecision << "  Recall TP/(TP+FP):" << recall << "  Spectra parsed:" << spectra.size() << "  Spectra above threshold:" << tp + fp << "  " << fraction*100 << "%";
		}
			
		if(precision >= 0.995){
			deltascore_performances->setBinValue(1, delta_i, fraction);
			deltascore_performances->setBinValue(4, delta_i, recall);
			deltascore_performances->setBinValue(7, delta_i, precision);
							//cout << "\nPrecision:" << lastprecision << "  Recall TP/(TP+FP):" << recall << "  Spectra parsed:" << spectra.size() << "  Spectra above threshold:" << tp + fp << "  " << fraction*100 << "%";
		}
		if(precision >= 0.95){
			deltascore_performances->setBinValue(2, delta_i, fraction);
			deltascore_performances->setBinValue(5, delta_i, recall);
			deltascore_performances->setBinValue(8, delta_i, precision);
							//cout << "\nPrecision:" << lastprecision << "  Recall TP/(TP+FP):" << recall << "  Spectra parsed:" << spectra.size() << "  Spectra above threshold:" << tp + fp << "  " << fraction*100 << "%";
		}

		//store precision, recall, fraction for use in proteins
		b[i]->precision = precision;
		b[i]->recall = recall;
		b[i]->fraction = fraction;	
		b[i]->tp = tp;
		b[i]->fp = fp;
	}
}


}