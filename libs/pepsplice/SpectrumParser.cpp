#include "SpectrumParser.h"
namespace Pepsplice{
SpectrumParser::SpectrumParser(Spectra *s1, Services *se0)
{
	spectra1 = s1;
	//dnaAA1 = se1->dnaAA1;
	hotspectra1 = new HotSpectra(se0);
	//sf = dnaAA1->scaling_factor;
	se1 = se0;
	//outputlevel = se1->outputlevel;
}

SpectrumParser::~SpectrumParser()
{
	delete hotspectra1;
}

void SpectrumParser::parseFiles(vector<string*> files)
{
	for(size_t i = 0; i < files.size(); i++){
		parseFile(files[i]);
	}
}

void SpectrumParser::parseFile(string *file)
{
	int suffixstart = file->size() - 3;
	if(suffixstart < 0) suffixstart = 0;
	string filesuffix = file->substr(suffixstart);
	if(filesuffix == "dta"){
		parseDTA(file);		
	}else if(filesuffix == "mgf"){
		parseMGF(file);
	}
}
	
void SpectrumParser::parseDTA(string *file)
{
	ifstream inFile;
	double parentmassMH0 = 0;
	double parentmassMH = 0;
	int chargestate;
	long length = 0;
	double mz0, intensity0;
	vector<double> mz;
	vector<double> intensity;
	
	inFile.open(file->c_str());
	inFile >> parentmassMH0;

	parentmassMH = parentmassMH0 * (se1->dnaAA1->scaling_factor);
	if(se1->outputlevel > 5) cout << "\nfile: " << file << " parentmassMH0: " << parentmassMH0 << " parentmassMH: " << parentmassMH;
	inFile >> chargestate;

	while (inFile.good() == true) //continue until end of file
	{
		inFile >> mz0;
		inFile >> intensity0;
		mz.push_back(mz0 * (se1->dnaAA1->scaling_factor));
		intensity.push_back(intensity0);
		length++;
	}
	inFile.close();

	if(intensity.size() > 0 && parentmassMH0 > 0){
		//write to spectrum
		Spectrum *spectrum1 = new Spectrum(parentmassMH, chargestate, se1->dnaAA1, se1); //dnaAA1 is a helper object with AA masses, translation etc...
		spectrum1->setFileName((*file));
		for(int i = 0; i < length; i++){
			spectrum1->addDataPoint(mz[i], intensity[i]);
		}
		spectrum1->preprocessSpectrum();
		spectra1->addSpectrum(spectrum1);
	}else{
		if(se1->outputlevel > 2) cout << "\nSpectrumParser::parseDTA: Spectrum file empty or corrupt: " << file << flush;
	}
}



void SpectrumParser::parseMGF(string *file) //TODO time to change stuff here
{
	ifstream inFile1, inFile2;
	int nb_line = 0;
	int nb_spec = 0;
	double mz0, intensity0;
    
	//bool filegood = true;
		
	string begin_ions, title, chargestring, pepmass;
	string controlstring;
	string dtafilename;
	int chargestate;
	double parentmassMmeas = 0; //e.g. MH2 / 2
	double parentmassMHX = 0; //eg. MH2
	double parentmassMH = 0; //e.g. MH2 - (2-1)H

	
	int line_begin_ion =0;
	int line_end_ion =0;
	vector<int> spec_lengths;
	//vector<double> mz;
	vector<double> mzintensity; // buffering the intensity
	//indexing file (fast scan)
    if(se1->outputlevel > 2)
    {
		cout << "\nIndexing " << (*file) << flush; //BX, less output!
		cout << "\nSpectra scanned: " << flush;
	}
    inFile1.open(file->c_str(), ios::binary);	
	nb_spec = 0;
	nb_line = 0;
	while(inFile1.good()){ //Comment: parsing spectra, adding number of lines to spec_lengths 
		getline(inFile1, begin_ions, '\n');
		nb_line++;
		if(begin_ions.substr(0, 8) == "BEGIN IO"){
			nb_spec++;
			if(nb_spec%500 == 0) cout << nb_spec << " " << flush;
			line_begin_ion = nb_line;			
		}
		if(begin_ions.substr(0, 8) == "END IONS"){
			line_end_ion = nb_line;
			spec_lengths.push_back(line_end_ion - line_begin_ion - 4);
			
		}
	}
	inFile1.close();
	
	//reading file
	//cout << "\nReading spectra from " << (*file) << flush; //BX
	
	//prepare subset outfile
	inFile2.open(file->c_str(),ios::binary);	
	nb_line = 0;
	for(int i = 0; i < nb_spec; i++){ //now going through the different spectra

        ofstream outFileSubset; // one subset at a time, the sel-> etc. can be changed with arguments.
        if(se1->writespecsubsetasmgf && se1->learning_spectrum_intensities){ // splitting the files
            stringstream osstream;
            osstream << file->substr(0,file->size()-3);
            osstream << i; 
            osstream << ".subset.mgf";
            outFileSubset.open(osstream.str().c_str());
            outFileSubset.precision(10);
        }

        //find spectrum start
        while(begin_ions.substr(0, 8) != "BEGIN IO"){
            getline(inFile2, begin_ions, '\n');
            nb_line++;
        }		
        begin_ions = "X";

        string templine = "";
        chargestring = "none_none";
        pepmass = "none_none";
        title = "none_none";
        for(int j = 0; j < 3; j++){
            getline(inFile2, templine, '\n');
            nb_line++;
            string tempsubstr = templine.substr(0, 5);
            if(tempsubstr == "CHARG") chargestring = templine;
            if(tempsubstr == "PEPMA") pepmass = templine;
            if(tempsubstr == "TITLE") title = templine;
        }


        //parse numbers from strings
        dtafilename = title.substr(6);//normal filename 
       // size_t pos = dtafilename.find("\r");

       //dtafilename.erase(pos); // additional carriage return in parameter list / argument actually messes with the output later
        //future TODO: pass some of the parameters from the main function. 
		parentmassMmeas = string_to_double(pepmass.substr(8)); //Comment: own function, at the bottom of this file (does nothing special, just using stringstream), i can probably tweak it to include two pepmasses 
        if(chargestring != "none_none"){
            chargestate = (int)string_to_double(chargestring.substr(7, 1));
        }else{
            chargestate = 0;
        }


        //if(chargestate == 2) cout << "\n" << chargestring << "  " << chargestate;

        //derive parent masses

        parentmassMHX = (parentmassMmeas * chargestate * se1->dnaAA1->scaling_factor); // - (dnaAA1->monomassH * chargestate) + dnaAA1->monomassH;
        parentmassMH = parentmassMHX - se1->dnaAA1->monomassH*(chargestate-1); 


        double parentmassC = 0.0,parentmassCMHX = 0.0,parentmassCMH = 0.0;
        if (se1->changedMass == true && chargestate != 0){
            parentmassC = parentmassMmeas - 1./chargestate;
            parentmassCMHX = (parentmassC * chargestate * se1->dnaAA1->scaling_factor); // - (dnaAA1->monomassH * chargestate) + dnaAA1->monomassH;
            parentmassCMH = parentmassCMHX - se1->dnaAA1->monomassH*(chargestate-1); 
            //Comment: with changed mass, saves writing the spectrumfile twice
        }




        //show parent mass conversion anomaly for first few spectra (too few digits in float)
        int showXspectra = 3;
        if(se1->outputlevel > 2 && i <= showXspectra){
            cout.precision(10);
            cout << "\n spectrum:" << dtafilename <<" parent mass in MGF:" << parentmassMmeas << "   parentmassMH:" << parentmassMH / se1->dnaAA1->scaling_factor << flush;
            
            //if(se1->changedMass) cout << "using changed mass";
            
            cout << "\n spectrum2:" << dtafilename <<" parent mass in MGF:" << parentmassC << "   parentmassCMH:" << parentmassCMH / se1->dnaAA1->scaling_factor << flush;
            
            if(i == showXspectra) cout << "\n";
        }

        //parse spectrum		
        bool ishotspectrum = hotspectra1->checkIfHot(dtafilename); // only deal with hot spectra
        if(i%se1->eachXspec == 0 && spec_lengths[i] > 0 && parentmassMH > 0 && chargestate > 0 && (se1->chargestate2 == false || chargestate == 2) && (se1->hotspectra == false || ishotspectrum) ){
            Spectrum *spectrum1 = new Spectrum(parentmassMH, chargestate, se1->dnaAA1, se1); //dnaAA1 is a helper object with AA masses, translation etc...
            //entering the spectrum.cpp, using derived parentmass
            spectrum1->setFileName(dtafilename);
            spectrum1->setMGFName(*file);
            Spectrum *spectrum2;
            
            bool wsm = false;	
            if(se1->writespecsubsetasmgf && se1->learning_spectrum_intensities) wsm = true;
            bool wsd = false;	
            if(se1->writespecsubsetasdta && se1->learning_spectrum_intensities) wsd = true;

            if(wsm){
                outFileSubset << "\nBEGIN IONS";
                outFileSubset << "\nTITLE=" << dtafilename;
                outFileSubset << "\nCHARGE=" << chargestate << "+";						
                outFileSubset << "\nPEPMASS=" << parentmassMmeas;					
            }

            ofstream dtafile;
            if(wsd){
                dtafile.open(dtafilename.c_str());
                dtafile.precision(10);
                dtafile << parentmassMH/se1->dnaAA1->scaling_factor << " " << chargestate << "\n";
            }
            if (se1->changedMass){
                for(int j = 0; j < spec_lengths[i]; j++){
                    inFile2 >> mz0;
                    inFile2 >> intensity0;
                    mzintensity.push_back(mz0);
                    mzintensity.push_back(intensity0);
                    if(wsm) outFileSubset << "\n" << mz0 << " " << intensity0;
                    if(wsd) dtafile << "\n" << mz0 << " " << intensity0;
                    spectrum1->addDataPoint(mz0 * se1->dnaAA1->scaling_factor, intensity0);
                     
                    nb_line++;
                }



            }
            else{
                for(int j = 0; j < spec_lengths[i]; j++){
                    inFile2 >> mz0;
                    inFile2 >> intensity0;
                    if(wsm) outFileSubset << "\n" << mz0 << " " << intensity0;
                    if(wsd) dtafile << "\n" << mz0 << " " << intensity0;
                    spectrum1->addDataPoint(mz0 * se1->dnaAA1->scaling_factor, intensity0);
                    nb_line++;
                }
            }
            if(wsm){
                outFileSubset << "\nEND IONS\n\n";
            }

            spectrum1->preprocessSpectrum();
            spectra1->addSpectrum(spectrum1);		
            if (se1->changedMass)
            {
                spectrum2 = new Spectrum(parentmassCMH, chargestate, se1->dnaAA1, se1);
                string newfilename(*file);
                spectrum2->setFileName(dtafilename);
                
                spectrum2->setMGFName(newfilename);
                vector<double>::iterator end = mzintensity.end();
                for (vector<double>::iterator it = mzintensity.begin(); it != end; it++){
                    
                    mz0 = *it;
                    ++it;
                    intensity0 = *it;

                    spectrum2->addDataPoint(mz0 * se1->dnaAA1->scaling_factor,intensity0);
                }
                spectrum2->preprocessSpectrum();
                spectra1->addSpectrum(spectrum2);
            }

        }

        outFileSubset.close();
    }
    inFile2.close();

}
double SpectrumParser::string_to_double(const string & s)
{
    stringstream stream(s);
    double b = -1.0;
    stream >> b;
    if (se1->outputlevel > 2 && b==-1){
        cout << "\nSpectrumParser::string_to_double: " << s << " is no double..." << endl;
    }
    return b;
}
}
