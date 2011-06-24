#include "DnaAA.h"
#include "bicepsdefinitions.h"
static string PARAMPATH = biceps::bicepsconfigpath.append("/in_AAmodifications.param");
namespace Pepsplice{
    
    DnaAA::DnaAA()
    {
        initialize();
        //checkInitialization();
        
    }
    
    void DnaAA::initialize()
    {
        
        scaling_factor = 1000000; //calculate with microDalton
        discretization_factor = scaling_factor * 1.00048; //monoisotopic, Perkins & Cottrell 1999 p. 3555
        
        monomassH = 1.0078250 * scaling_factor;
        monomassO = 15.9949146 * scaling_factor;
        monomassN = 14.0030740 * scaling_factor;
        monomassPhosphorylation = 79.96633 * scaling_factor;
        monomassC13minusC12 = (13.0034 - 12.0000) * scaling_factor;
        monomassHO = monomassH + monomassO;
        monomassHOH = monomassH + monomassO + monomassH;
        monomassHOHH = monomassH + monomassO + monomassH + monomassH;
        monomassNH3 = monomassH * 3 + monomassN;
        
        integermassH = (int)(monomassH / scaling_factor);
        integermassHOHH = (int)(monomassHOHH / scaling_factor);
        
        
        //fill in 64 DNA triplets
        char triplet[3] = {'_', '_', '_'};
        for(char i=0; i<4; i++){
            for(char j=0; j<4; j++){
                for(char k=0; k<4; k++){
                    dna_triplets[i][j][k][0] = num_nt(i);
                    dna_triplets[i][j][k][1] = num_nt(j);
                    dna_triplets[i][j][k][2] = num_nt(k);
                    
                    triplet[0] = num_nt(i);
                    triplet[1] = num_nt(j);
                    triplet[2] = num_nt(k);
                    
                    aa_chars[i][j][k] = nt_aa(triplet);
                    aa_monomasses444[i][j][k] = aa_monoisotopic(nt_aa(triplet));
                }
            }
        }
        
        //initialize SNPs
        initializeSNPs();
        
        //load modifications
        loadAAModifications();
        
        //fill in amino acids
        char aa = '-';
        for(int i = 0; i < 256; i++){
            aa_monomasses256[i] = aa_monoisotopic(i);
            aa_monomasses256float[i] = aa_monoisotopic(i);
            aa_integermasses256[i] = (int)(aa_monoisotopic(i)/scaling_factor + 0.5);
            //cout << "\n" << aa_integermasses256[i] << " " << aa_monoisotopic(i);
        }
        
    }
    
    //This version returns invalid AAs, disabled by BX
    //void DnaAA::initializeSNPs()
    //{
    //	
    //	//enumerate all triplets, initialize array
    //	for(char i=0; i<4; i++){
    //		for(char j=0; j<4; j++){
    //			for(char k=0; k<4; k++){
    //				for(int l=0; l<256; l++){
    //					tripletSNPaa[i][j][k][l] = true; //changes BYR
    //				}				
    //			}
    //		}
    //	}
    //	
    //	//initialize array
    //	for(int i = 0; i < 256; i++){
    //		for(int j = 0; j < 256; j++){
    //			aaSNPaa[i][j] = true;//changed BYR	
    //		}
    //	}    
    //}
    //    
    //    
    
    
    
    
    
    void DnaAA::initializeSNPs()
    {
        
        //enumerate all triplets, initialize array
        for(char i=0; i<4; i++){
            for(char j=0; j<4; j++){
                for(char k=0; k<4; k++){
                    for(int l=0; l<256; l++){
                        tripletSNPaa[i][j][k][l] = false;
                    }				
                }
            }
        }
        
        //initialize array
        for(int i = 0; i < 256; i++){
            for(int j = 0; j < 256; j++){
                aaSNPaa[i][j] = false;	
            }
        }	
        
        //enumerate all triplets
        for(char i=0; i<4; i++){
            for(char j=0; j<4; j++){
                for(char k=0; k<4; k++){
                    
                    //enumerate all triplets
                    for(char l=0; l<4; l++){
                        for(char m=0; m<4; m++){
                            for(char n=0; n<4; n++){
                                
                                //count identities, decide if SNP
                                int identical = 0;
                                if(i==l) identical++;
                                if(j==m) identical++;
                                if(k==n) identical++;
                                
                                char oldaa = aa_chars[i][j][k];
                                char newaa = aa_chars[l][m][n];
                                
                                //SNP true
                                if(identical==2 && oldaa != newaa){
                                    tripletSNPaa[i][j][k][newaa] = true; //nucleotide based
                                    aaSNPaa[oldaa][newaa] = true;        //amino acid based
                                }						
                                
                            }
                        }
                    }
                    
                }
            }
        }
        
        
        
    }    
    
    
    
    
    void DnaAA::loadAAModifications()
    {
        //initialize arrays
        for(int i = 0; i < 256; i++){
            aamod_aa[i] = 0;
            aamod_monomasses[i] = 0;
            aamod_startascii[i] =0;
            aamod_nbmodperaa[i] = 0;
        }
        
        ifstream inFile1;
        string aamodline;
        double mass = 0;
        int aamod_i_ascii = 128;
        bool protNterm = false;
        bool protCterm = false;
        inFile1.open(PARAMPATH.c_str(), ios::binary);	
        while(inFile1.good()){
            getline(inFile1, aamodline);
            if(aamodline[0] != '#'){
              if (aamod_line[ammod_line.size()-1] == '\r') aamod_line = aamod_line.substr(0,aamod_line.size() - 1);
              //cout << "\n" << aamodline;

              //parse line field-wise		
              istringstream iss(aamodline);
              string field;
              int i_field = 0;
              unsigned char aa = 0;	
              while( getline(iss, field, ';') ){
                i_field++;
                if(i_field == 1){
                  aa = (unsigned char)field[0];
                  //assign ascii start code for current amino acid
                  aamod_startascii[aa]=aamod_i_ascii;
                }else{
                  //assign alternative mass for current amino acid
                  unsigned char c = field[0];
                  string m = field.substr(1);
                  //cout << "\nc: " << c << " m: " << m;
                  //modification type: (pepNterm, )pepCterm, [protNterm, ]protCterm, _internal
                  if(c == '(' || c == ')' || c == '[' || c == ']' || c == '_'){
                    if(string_to_double(m) > 0){
                      aamod_aa[aamod_i_ascii] = aa;			
                      aamod_type[aamod_i_ascii] = c;
                      aamod_monomasses[aamod_i_ascii] = string_to_double(m);
                      aamod_i_ascii++;
                    }
                  }else{
                    cout << "\nDnaAA.cpp line 99: please convert the modification file with dos2unix or else define what type of modification you require: " << aamodline << "\n";
                  }
                }

                aamod_nbmodperaa[aa] = i_field - 1;
              }//tag or not tag
            }//not #
        }//while loop
        inFile1.close();
    }

    double DnaAA::aa_monoisotopic(unsigned char aa)
    {
      double rt = 0;

      if(aa=='A' || aa=='a') rt =  71.03711;
      if(aa=='C' || aa=='c') rt = 103.00919 + 57.02146;
      if(aa=='D' || aa=='d') rt = 115.02694;
      if(aa=='E' || aa=='e') rt = 129.04259;
      if(aa=='F' || aa=='f') rt = 147.06841;
      if(aa=='G' || aa=='g') rt =  57.02146;
      if(aa=='H' || aa=='h') rt = 137.05891;
      if(aa=='I' || aa=='i') rt = 113.08406;
      if(aa=='K' || aa=='k') rt = 128.09496;
      if(aa=='L' || aa=='l') rt = 113.08406;
      if(aa=='M' || aa=='m') rt = 131.04049;
      if(aa=='N' || aa=='n') rt = 114.04293;
      if(aa=='P' || aa=='p') rt =  97.05276;
      if(aa=='Q' || aa=='q') rt = 128.05858;
      if(aa=='R' || aa=='r') rt = 156.10111;
      if(aa=='S' || aa=='s') rt =  87.03203;
      if(aa=='T' || aa=='t') rt = 101.04768;
      if(aa=='V' || aa=='v') rt =  99.06841;
      if(aa=='W' || aa=='w') rt = 186.07931;
      if(aa=='Y' || aa=='y') rt = 163.06333;

      if(aa > 127) rt = aamod_monomasses[aa];

      //cout.precision(30);
      //cout << "\nCheck rounding: " << rt * 1000000 << "\tcheck conversion double to long" << (long)(rt * 1000000 + 0.5);
      if(scaling_factor == 0) cout << "scaling_factor is 0!";
      return (double)((long)(rt * scaling_factor + 0.5));

    }

    void DnaAA::checkInitialization()
    {
      for(char i=0; i<4; i++){
        for(char j=0; j<4; j++){
          for(char k=0; k<4; k++){
            cout << "\n" << (long)i << (long)j << (long)k;
            cout << "  triplet:" << dna_triplets[i][j][k][0] << dna_triplets[i][j][k][1] << dna_triplets[i][j][k][2];
            cout << "  amino acid:" << aa_chars[i][j][k];
            cout.precision(12);
            cout << "  mass:" << aa_monomasses444[i][j][k]/scaling_factor;
          }
        }
      }

      for(int i=0; i<256; i++){
        cout << "\nascii: " << i << "  mass: " << aa_monomasses256[i]/scaling_factor;
        if(i >= 65) cout << "  symbol: " << (char)i << "  number of modifications: " << aamod_nbmodperaa[i] << "  startascii:" << aamod_startascii[i];
      }

      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          for(int k=0; k<4; k++){
            cout << "\n" << i << j << k << "\t";
            for(int l=65; l<65+26; l++){
              cout << tripletSNPaa[i][j][k][l];
            }				
          }
        }
      }

      for(char i=65; i<65+26; i++){
        cout << "\n" << i << "\t";
        for(char j=65; j<65+26; j++){
          cout << aaSNPaa[i][j];
        }
      }

      cout.precision(7);
    }	


    char DnaAA::num_nt(char n)
    {
      if(n==0) return 'A';
      if(n==1) return 'C';
      if(n==2) return 'G';
      if(n==3) return 'T';
      return 'N';
    }

    char DnaAA::nt_num(char n)
    {
      if(n=='A') return 0;
      if(n=='C') return 1;
      if(n=='G') return 2;
      if(n=='T') return 3;
      if(n=='a') return 0;
      if(n=='c') return 1;
      if(n=='g') return 2;
      if(n=='t') return 3;
      return -1;
    }

    char DnaAA::nt_revnt(char n)
    {
      if(n=='A') return 'T';
      if(n=='C') return 'G';
      if(n=='G') return 'C';
      if(n=='T') return 'A';
      if(n=='a') return 'T';
      if(n=='c') return 'G';
      if(n=='g') return 'C';
      if(n=='t') return 'A';
      return -1;
    }

    char DnaAA::nt_nt_clean(char n)
    {
      return num_nt(nt_num(n));	
    }


    unsigned char DnaAA::nt_aa(char triplet[3])
    {	
      unsigned char rt = 'X';

      if(tripletis(triplet, "TTT")) rt = 'F';
      if(tripletis(triplet, "TTC")) rt = 'F';

      if(tripletis(triplet, "TTA")) rt = 'L';
      if(tripletis(triplet, "TTG")) rt = 'L';

      if(tripletis(triplet, "TCT")) rt = 'S';
      if(tripletis(triplet, "TCC")) rt = 'S';
      if(tripletis(triplet, "TCA")) rt = 'S';
      if(tripletis(triplet, "TCG")) rt = 'S';

      if(tripletis(triplet, "TAT")) rt = 'Y';
      if(tripletis(triplet, "TAC")) rt = 'Y';

      if(tripletis(triplet, "CCT")) rt = 'P';
      if(tripletis(triplet, "CCC")) rt = 'P';
      if(tripletis(triplet, "CCA")) rt = 'P';
      if(tripletis(triplet, "CCG")) rt = 'P';

      if(tripletis(triplet, "TGT")) rt = 'C';	
      if(tripletis(triplet, "TGC")) rt = 'C';

      if(tripletis(triplet, "TGG")) rt = 'W';

      if(tripletis(triplet, "CTT")) rt = 'L';
      if(tripletis(triplet, "CTC")) rt = 'L';
      if(tripletis(triplet, "CTA")) rt = 'L';
      if(tripletis(triplet, "CTG")) rt = 'L';

      if(tripletis(triplet, "CAT")) rt = 'H';
      if(tripletis(triplet, "CAC")) rt = 'H';	

      if(tripletis(triplet, "CAA")) rt = 'Q';
      if(tripletis(triplet, "CAG")) rt = 'Q';	

      if(tripletis(triplet, "CGT")) rt = 'R';
      if(tripletis(triplet, "CGC")) rt = 'R';
      if(tripletis(triplet, "CGA")) rt = 'R';
      if(tripletis(triplet, "CGG")) rt = 'R';

      if(tripletis(triplet, "ATT")) rt = 'I';
      if(tripletis(triplet, "ATC")) rt = 'I';
      if(tripletis(triplet, "ATA")) rt = 'I';

      if(tripletis(triplet, "TAA")) rt = '-';
      if(tripletis(triplet, "TAG")) rt = '-';	
      if(tripletis(triplet, "TGA")) rt = '-';

      if(tripletis(triplet, "ATT")) rt = 'I';
      if(tripletis(triplet, "ATC")) rt = 'I';	
      if(tripletis(triplet, "ATA")) rt = 'I';	

      if(tripletis(triplet, "ATG")) rt = 'M';

      if(tripletis(triplet, "ACT")) rt = 'T';
      if(tripletis(triplet, "ACC")) rt = 'T';
      if(tripletis(triplet, "ACA")) rt = 'T';
      if(tripletis(triplet, "ACG")) rt = 'T';		

      if(tripletis(triplet, "AAT")) rt = 'N';
      if(tripletis(triplet, "AAC")) rt = 'N';		

      if(tripletis(triplet, "AAA")) rt = 'K';
      if(tripletis(triplet, "AAG")) rt = 'K';

      if(tripletis(triplet, "AGT")) rt = 'S';
      if(tripletis(triplet, "AGC")) rt = 'S';

      if(tripletis(triplet, "AGA")) rt = 'R';
      if(tripletis(triplet, "AGG")) rt = 'R';

      if(tripletis(triplet, "GTT")) rt = 'V';
      if(tripletis(triplet, "GTC")) rt = 'V';
      if(tripletis(triplet, "GTA")) rt = 'V';
      if(tripletis(triplet, "GTG")) rt = 'V';

      if(tripletis(triplet, "GCT")) rt = 'A';
      if(tripletis(triplet, "GCC")) rt = 'A';
      if(tripletis(triplet, "GCA")) rt = 'A';
      if(tripletis(triplet, "GCG")) rt = 'A';

      if(tripletis(triplet, "GAT")) rt = 'D';
      if(tripletis(triplet, "GAC")) rt = 'D';

      if(tripletis(triplet, "GAA")) rt = 'E';
      if(tripletis(triplet, "GAG")) rt = 'E';

      if(tripletis(triplet, "GGT")) rt = 'G';
      if(tripletis(triplet, "GGC")) rt = 'G';
      if(tripletis(triplet, "GGA")) rt = 'G';
      if(tripletis(triplet, "GGG")) rt = 'G';

      return rt;
    }

    bool DnaAA::tripletis(char tripletA[3], const char tripletB[])
    {
      bool rt = false;
      tripletA[0] = toupper(tripletA[0]);
      tripletA[1] = toupper(tripletA[1]);
      tripletA[2] = toupper(tripletA[2]);
      if(tripletA[0] == tripletB[0] && tripletA[1] == tripletB[1] && tripletA[2] == tripletB[2]) rt = true;
      return rt;
    }

    //translate a nucleotide sequence to an amino acid sequence
    string DnaAA::translate(string ntseq)
    {
      //cout << "\nntseq: " << ntseq;
      string aaseq;
      char a, b, c;
      for(int i = 0; i < ntseq.size(); i = i + 3){
        a = nt_num(ntseq[i]);
        b = nt_num(ntseq[i+1]);
        c = nt_num(ntseq[i+2]);
        aaseq += aa_chars[a][b][c];
        //cout << "\nDnaAA::translate  a b c " << a << " " << b << " " << c << " aaseq: " << aaseq;
      }
      return aaseq;
    }


    double DnaAA::getPRMfromAA(string aaseq)
    {
      double prm = 0;
      //cout << "\nDnaAA::getPRMfromAA: aaseq: " << aaseq;
      for(int i = 0; i < aaseq.size(); i++){
        prm += aa_monomasses256[aaseq[i]];
        //cout << "\nDnaAA::getPRMfromAA: aaseq[i]:" << aaseq[i] << " prm:" << prm;
      }
      return prm;
    }

    double DnaAA::getParentMassMH(unsigned char *aaseq, long len)
    {
      double pmmh = 0;
      //cout << "\nDnaAA::getpmmhfromAA: aaseq: " << aaseq;
      for(int i = 0; i < len; i++){
        pmmh += aa_monomasses256[aaseq[i]];
        //cout << "\nDnaAA::getpmmhfromAA: aaseq[i]:" << aaseq[i] << " pmmh:" << pmmh;
      }
      pmmh += monomassHOHH;
      return pmmh;
    }

    double DnaAA::getParentMassMH(string aaseq) //changed by BX, debugging!
    {
      double pmmh = 0.0;
      for(int i = 0; i < aaseq.size(); i++){
        pmmh += aa_monomasses256[aaseq[i]];

      }
      pmmh += monomassHOHH;
      return pmmh;
    }

    double DnaAA::string_to_double(string s)
    {
      stringstream stream(s);
      double b = -1.0;
      stream >> b;
      return b;
    }

    string DnaAA::intToString(int i){
      std::ostringstream o;
      o << i;
      return o.str();
    }

    string DnaAA::charseq_to_modnumseq(string s1, unsigned char charNterm, unsigned char charCterm)
    {
      string s2;
      if(charNterm != 0){
        s2+=intToString((int)aamod_monomasses[charNterm]);
        s2+="-";
      }
      if(charCterm != 0){
        s2+="-";
        s2+=intToString((int)aamod_monomasses[charCterm]);
      }

      //cout << "\n" << s1;
      for(int i = 0; i < s1.size(); i++){
        unsigned char aai = s1[i];
        if(aai >= 128){
          s2 += (unsigned char)aamod_aa[aai];
          s2 += intToString((int)aamod_monomasses[aai]);
        }else{
          s2 += aai;
        }
        //cout << "\n" << s2;
      }
      //cout << " " << s2;
      return s2;
    }

    string DnaAA::charseq_to_lowercaseseq(string s1)
    {
      string s2;
      for(int i = 0; i < s1.size(); i++){
        unsigned char aai = s1[i];
        if(aai >= 128){
          s2 += (unsigned char)(aamod_aa[aai]+32);
        }else{
          s2 += aai;

        }
      }
      return s2;
    }

    double DnaAA::aa_pampenalty(unsigned char aa, unsigned char bb)
    {
      //bb is old amino acid, aa new one


      switch(aa)
      {
        case 'A':
          switch(bb){

            case 'R': return 1.7;
            case 'N': return 1.2;
            case 'D': return 1.2;
            case 'C': return 1.6;
            case 'Q': return 1.3;
            case 'E': return 1;
            case 'G': return 0.9;
            case 'H': return 1.7;
            case 'I': return 1.4;
            case 'L': return 1.5;
            case 'K': return 1.7;
            case 'M': return 1.4;
            case 'F': return 1.7;
            case 'P': return 0.8;
            case 'S': return 0.6;
            case 'T': return 0.7;
            case 'W': return 2;
            case 'Y': return 1.7;
            case 'V': return 1;

            default: return 2.5;
          }

        case 'R':
          switch(bb){
            case 'A': return 1.9;
            case 'N': return 1.9;
            case 'D': return 2;
            case 'C': return 1.9;
            case 'Q': return 1.2;
            case 'E': return 2;
            case 'G': return 2;
            case 'H': return 1.2;
            case 'I': return 1.6;
            case 'L': return 1.9;
            case 'K': return 1;
            case 'M': return 1.5;
            case 'F': return 1.9;
            case 'P': return 1.5;
            case 'S': return 1.4;
            case 'T': return 1.9;
            case 'W': return 1.3;
            case 'Y': return 2;
            case 'V': return 1.9;

            default: return 2.5;


          }
        case 'N':
          switch(bb){
            case 'A': return 1.5;
            case 'R': return 1.9;
            case 'D': return 0.5;
            case 'C': return 2;
            case 'Q': return 1.5;
            case 'E': return 1.4;
            case 'G': return 1.4;
            case 'H': return 0.9;
            case 'I': return 1.6;
            case 'L': return 1.9;
            case 'K': return 1.1;
            case 'M': return 2;
            case 'F': return 1.9;
            case 'P': return 1.7;
            case 'S': return 0.9;
            case 'T': return 1.2;
            case 'W': return 1.9;
            case 'Y': return 1.5;
            case 'V': return 1.9;
            default: return 2.5;
          }

        case 'D':
          switch(bb){
            case 'A': return 1.4;
            case 'R': return 2;
            case 'N': return 0.4;
            case 'C': return 2;
            case 'Q': return 1.4;
            case 'E': return 0.3;
            case 'G': return 1.4;
            case 'H': return 1.5;
            case 'I': return 1.9;
            case 'L': return 2;
            case 'K': return 1.6;
            case 'M': return 2;
            case 'F': return 2;
            case 'P': return 1.9;
            case 'S': return 1.4;
            case 'T': return 1.6;
            case 'W': return 2;
            case 'Y': return 2;
            case 'V': return 1.9;
            default: return 2.5;
          }


        case 'C':
          switch(bb){
            case 'A': return 1.9;
            case 'R': return 1.9;
            case 'N': return 2;
            case 'D': return 2;
            case 'Q': return 2;
            case 'E': return 2;
            case 'G': return 2;
            case 'H': return 1.9;
            case 'I': return 1.9;
            case 'L': return 2;
            case 'K': return 2;
            case 'M': return 2;
            case 'F': return 2;
            case 'P': return 1.9;
            case 'S': return 1.4;
            case 'T': return 1.9;
            case 'W': return 2;
            case 'Y': return 1.6;
            case 'V': return 1.7;
            default: return 2.5;

          }
        case 'Q':
          switch(bb){
            case 'A': return 1.6;
            case 'R': return 1.2;
            case 'N': return 1.5;
            case 'D': return 1.4;
            case 'C': return 2;
            case 'E': return 0.7;
            case 'G': return 1.9;
            case 'H': return 0.8;
            case 'I': return 1.9;
            case 'L': return 1.6;
            case 'K': return 1.4;
            case 'M': return 1.5;
            case 'F': return 2;
            case 'P': return 1.4;
            case 'S': return 1.7;
            case 'T': return 1.7;
            case 'W': return 2;
            case 'Y': return 2;
            case 'V': return 1.9;
            default: return 2.5;
          }
        case 'E':
          switch(bb){
            case 'A': return 1.2;
            case 'R': return 2;
            case 'N': return 1.3;
            case 'D': return 0.2;
            case 'C': return 2;
            case 'Q': return 0.6;
            case 'G': return 1.5;
            case 'H': return 1.7;
            case 'I': return 1.6;
            case 'L': return 1.9;
            case 'K': return 1.5;
            case 'M': return 1.9;
            case 'F': return 2;
            case 'P': return 1.6;
            case 'S': return 1.5;
            case 'T': return 1.7;
            case 'W': return 2;
            case 'Y': return 1.9;
            case 'V': return 1.7;
            default: return 2.5;
          }
        case 'G':
          switch(bb){
            case 'A': return 0.9;
            case 'R': return 1.9;
            case 'N': return 1.1;
            case 'D': return 1.2;
            case 'C': return 1.9;
            case 'Q': return 1.6;
            case 'E': return 1.3;
            case 'H': return 1.9;
            case 'I': return 2;
            case 'L': return 1.9;
            case 'K': return 1.7;
            case 'M': return 1.9;
            case 'F': return 1.9;
            case 'P': return 1.6;
            case 'S': return 0.9;
            case 'T': return 1.6;
            case 'W': return 2;
            case 'Y': return 2;
            case 'V': return 1.4;
            default: return 2.5;
          }
        case 'H':
          switch(bb){
            case 'A': return 1.9;
            case 'R': return 1.3;
            case 'N': return 1;
            case 'D': return 1.6;
            case 'C': return 1.9;
            case 'Q': return 0.9;
            case 'E': return 1.9;
            case 'G': return 2;
            case 'I': return 2;
            case 'L': return 1.9;
            case 'K': return 1.9;
            case 'M': return 2;
            case 'F': return 1.7;
            case 'P': return 1.6;
            case 'S': return 1.9;
            case 'T': return 1.9;
            case 'W': return 1.9;
            case 'Y': return 1.5;
            case 'V': return 1.9;
            default: return 2.5;
          }   
        case 'I':
          switch(bb){
            case 'A': return 1.7;
            case 'R': return 1.7;
            case 'N': return 1.6;
            case 'D': return 1.9;
            case 'C': return 1.7;
            case 'Q': return 1.9;
            case 'E': return 1.7;
            case 'G': return 2;
            case 'H': return 2;
            case 'L': return 1.2;
            case 'K': return 1.7;
            case 'M': return 1.1;
            case 'F': return 1.3;
            case 'P': return 2;
            case 'S': return 1.9;
            case 'T': return 1.3;
            case 'W': return 2;
            case 'Y': return 1.9;
            case 'V': return 0.6;
            default: return 2.5;
          }
        case 'L':
          switch(bb){
            case 'A': return 1.6;
            case 'R': return 1.9;
            case 'N': return 1.6;
            case 'D': return 2;
            case 'C': return 2;
            case 'Q': return 1.4;
            case 'E': return 1.9;
            case 'G': return 1.9;
            case 'H': return 1.5;
            case 'I': return 0.8;
            case 'K': return 1.7;
            case 'M': return 0.4;
            case 'F': return 1.1;
            case 'P': return 1.6;
            case 'S': return 1.9;
            case 'T': return 1.6;
            case 'W': return 1.5;
            case 'Y': return 1.7;
            case 'V': return 1;
            default: return 2.5;
          }
        case 'K':
          switch(bb){
            case 'A': return 1.7;
            case 'R': return 0.5;
            case 'N': return 0.8;
            case 'D': return 1.4;
            case 'C': return 2;
            case 'Q': return 1.1;
            case 'E': return 1.3;
            case 'G': return 1.7;
            case 'H': return 1.7;
            case 'I': return 1.5;
            case 'L': return 1.9;
            case 'M': return 0.9;
            case 'F': return 2;
            case 'P': return 1.6;
            case 'S': return 1.3;
            case 'T': return 1.2;
            case 'W': return 2;
            case 'Y': return 1.9;
            case 'V': return 1.9;
            default: return 2.5;
          }
        case 'M':
          switch(bb){
            case 'A': return 1.9;
            case 'R': return 1.9;
            case 'N': return 2;
            case 'D': return 2;
            case 'C': return 2;
            case 'Q': return 1.7;
            case 'E': return 2;
            case 'G': return 2;
            case 'H': return 2;
            case 'I': return 1.4;
            case 'L': return 1.3;
            case 'K': return 1.5;
            case 'F': return 1.9;
            case 'P': return 2;
            case 'S': return 1.9;
            case 'T': return 1.7;
            case 'W': return 2;
            case 'Y': return 2;
            case 'V': return 1.5;
            default: return 2.5;
          }
        case 'F':
          switch(bb){     
            case 'A': return 1.9;
            case 'R': return 1.9;
            case 'N': return 1.9;
            case 'D': return 2;
            case 'C': return 2;
            case 'Q': return 2;
            case 'E': return 2;
            case 'G': return 1.9;
            case 'H': return 1.7;
            case 'I': return 1.3;
            case 'L': return 1.4;
            case 'K': return 2;
            case 'M': return 1.5;
            case 'P': return 2;
            case 'S': return 1.7;
            case 'T': return 1.9;
            case 'W': return 1.6;
            case 'Y': return 0.7;
            case 'V': return 2;
            default: return 2.5;
          }
        case 'P':      
          switch(bb){     
            case 'A': return 1.1;
            case 'R': return 1.4;
            case 'N': return 1.7;
            case 'D': return 1.9;
            case 'C': return 1.9;
            case 'Q': return 1.3;
            case 'E': return 1.6;
            case 'G': return 1.7;
            case 'H': return 1.4;
            case 'I': return 1.9;
            case 'L': return 1.7;
            case 'K': return 1.7;
            case 'M': return 1.9;
            case 'F': return 1.9;
            case 'S': return 1.1;
            case 'T': return 1.5;
            case 'W': return 2;
            case 'Y': return 2;
            case 'V': return 1.7;
            default: return 2.5;
          }
        case 'S':
          switch(bb){     
            case 'A': return 0.7;
            case 'R': return 1.2;
            case 'N': return 0.6;
            case 'D': return 1.3;
            case 'C': return 1.2;
            case 'Q': return 1.5;
            case 'E': return 1.4;
            case 'G': return 1;
            case 'H': return 1.7;
            case 'I': return 1.7;
            case 'L': return 1.9;
            case 'K': return 1.3;
            case 'M': return 1.5;
            case 'F': return 1.6;
            case 'P': return 1;
            case 'T': return 0.4;
            case 'W': return 1.4;
            case 'Y': return 1.7;
            case 'V': return 1.7;
            default: return 2.5;
          }

        case 'T':
          switch(bb){
            case 'A': return 0.8;
            case 'R': return 1.7;
            case 'N': return 1.1;
            case 'D': return 1.5;
            case 'C': return 1.9;
            case 'Q': return 1.6;
            case 'E': return 1.7;
            case 'G': return 1.7;
            case 'H': return 1.9;
            case 'I': return 1.2;
            case 'L': return 1.7;
            case 'K': return 1.3;
            case 'M': return 1.4;
            case 'F': return 1.9;
            case 'P': return 1.4;
            case 'S': return 0.7;
            case 'W': return 2;
            case 'Y': return 1.7;
            case 'V': return 1.2;
            default: return 2.5;
          }
        case 'W':
          switch(bb){
            case 'A': return 2;
            case 'R': return 1.7;
            case 'N': return 2;
            case 'D': return 2;
            case 'C': return 2;
            case 'Q': return 2;
            case 'E': return 2;
            case 'G': return 2;
            case 'H': return 2;
            case 'I': return 2;
            case 'L': return 2;
            case 'K': return 2;
            case 'M': return 2;
            case 'F': return 1.9;
            case 'P': return 2;
            case 'S': return 1.9;
            case 'T': return 2;
            case 'Y': return 1.9;
            case 'V': return 2;
            default: return 2.5;
          }
        case 'Y':
          switch(bb){
            case 'A': return 1.9;
            case 'R': return 2;
            case 'N': return 1.6;
            case 'D': return 2;
            case 'C': return 1.6;
            case 'Q': return 2;
            case 'E': return 1.9;
            case 'G': return 2;
            case 'H': return 1.5;
            case 'I': return 1.9;
            case 'L': return 1.9;
            case 'K': return 2;
            case 'M': return 2;
            case 'F': return 0.9;
            case 'P': return 2;
            case 'S': return 1.9;
            case 'T': return 1.9;
            case 'W': return 1.7;
            case 'V': return 1.9;
            default: return 2.5;

          }
        case 'V':
          switch(bb){
            case 'A': return 1.1;
            case 'R': return 1.7;
            case 'N': return 1.9;
            case 'D': return 1.9;
            case 'C': return 1.6;
            case 'Q': return 1.7;
            case 'E': return 1.7;
            case 'G': return 1.6;
            case 'H': return 1.6;
            case 'I': return 0.2;
            case 'L': return 1.2;
            case 'K': return 1.9;
            case 'M': return 1;
            case 'F': return 1.9;
            case 'P': return 1.6;
            case 'S': return 1.7;
            case 'T': return 1.2;
            case 'W': return 2;
            case 'Y': return 1.7;
            default: return 2.5;

          }
        default: return 2.5;    

      }
      return 2.5;    

    } //function
} //namespace











//        if (aa=='A' && bb=='R') return 1.7;
//        if (aa=='A' && bb=='N') return 1.2;
//        if (aa=='A' && bb=='D') return 1.2;
//        if (aa=='A' && bb=='C') return 1.6;
//        if (aa=='A' && bb=='Q') return 1.3;
//        if (aa=='A' && bb=='E') return 1;
//        if (aa=='A' && bb=='G') return 0.9;
//        if (aa=='A' && bb=='H') return 1.7;
//        if (aa=='A' && bb=='I') return 1.4;
//        if (aa=='A' && bb=='L') return 1.5;
//        if (aa=='A' && bb=='K') return 1.7;
//        if (aa=='A' && bb=='M') return 1.4;
//        if (aa=='A' && bb=='F') return 1.7;
//        if (aa=='A' && bb=='P') return 0.8;
//        if (aa=='A' && bb=='S') return 0.6;
//        if (aa=='A' && bb=='T') return 0.7;
//        if (aa=='A' && bb=='W') return 2;
//        if (aa=='A' && bb=='Y') return 1.7;
//        if (aa=='A' && bb=='V') return 1;
//   if (aa=='R' && bb=='A') return 1.9;
//    if (aa=='R' && bb=='N') return 1.9;
//    if (aa=='R' && bb=='D') return 2;
//    if (aa=='R' && bb=='C') return 1.9;
//    if (aa=='R' && bb=='Q') return 1.2;
//    if (aa=='R' && bb=='E') return 2;
//    if (aa=='R' && bb=='G') return 2;
//    if (aa=='R' && bb=='H') return 1.2;
//    if (aa=='R' && bb=='I') return 1.6;
//    if (aa=='R' && bb=='L') return 1.9;
//    if (aa=='R' && bb=='K') return 1;
//    if (aa=='R' && bb=='M') return 1.5;
//    if (aa=='R' && bb=='F') return 1.9;
//    if (aa=='R' && bb=='P') return 1.5;
//    if (aa=='R' && bb=='S') return 1.4;
//    if (aa=='R' && bb=='T') return 1.9;
//    if (aa=='R' && bb=='W') return 1.3;
//    if (aa=='R' && bb=='Y') return 2;
//    if (aa=='R' && bb=='V') return 1.9;
//    
//    
//    if (aa=='N' && bb=='A') return 1.5;
//    if (aa=='N' && bb=='R') return 1.9;
//    if (aa=='N' && bb=='D') return 0.5;
//    if (aa=='N' && bb=='C') return 2;
//    if (aa=='N' && bb=='Q') return 1.5;
//    if (aa=='N' && bb=='E') return 1.4;
//    if (aa=='N' && bb=='G') return 1.4;
//    if (aa=='N' && bb=='H') return 0.9;
//    if (aa=='N' && bb=='I') return 1.6;
//    if (aa=='N' && bb=='L') return 1.9;
//    if (aa=='N' && bb=='K') return 1.1;
//    if (aa=='N' && bb=='M') return 2;
//    if (aa=='N' && bb=='F') return 1.9;
//    if (aa=='N' && bb=='P') return 1.7;
//    if (aa=='N' && bb=='S') return 0.9;
//    if (aa=='N' && bb=='T') return 1.2;
//    if (aa=='N' && bb=='W') return 1.9;
//    if (aa=='N' && bb=='Y') return 1.5;
//    if (aa=='N' && bb=='V') return 1.9;
//    if (aa=='D' && bb=='A') return 1.4;
//    if (aa=='D' && bb=='R') return 2;
//    if (aa=='D' && bb=='N') return 0.4;
//    if (aa=='D' && bb=='C') return 2;
//    if (aa=='D' && bb=='Q') return 1.4;
//    if (aa=='D' && bb=='E') return 0.3;
//    if (aa=='D' && bb=='G') return 1.4;
//    if (aa=='D' && bb=='H') return 1.5;
//    if (aa=='D' && bb=='I') return 1.9;
//    if (aa=='D' && bb=='L') return 2;
//    if (aa=='D' && bb=='K') return 1.6;
//    if (aa=='D' && bb=='M') return 2;
//    if (aa=='D' && bb=='F') return 2;
//    if (aa=='D' && bb=='P') return 1.9;
//    if (aa=='D' && bb=='S') return 1.4;
//    if (aa=='D' && bb=='T') return 1.6;
//    if (aa=='D' && bb=='W') return 2;
//    if (aa=='D' && bb=='Y') return 2;
//    if (aa=='D' && bb=='V') return 1.9;
//    if (aa=='C' && bb=='A') return 1.9;
//    if (aa=='C' && bb=='R') return 1.9;
//    if (aa=='C' && bb=='N') return 2;
//    if (aa=='C' && bb=='D') return 2;
//    if (aa=='C' && bb=='Q') return 2;
//    if (aa=='C' && bb=='E') return 2;
//    if (aa=='C' && bb=='G') return 2;
//    if (aa=='C' && bb=='H') return 1.9;
//    if (aa=='C' && bb=='I') return 1.9;
//    if (aa=='C' && bb=='L') return 2;
//    if (aa=='C' && bb=='K') return 2;
//    if (aa=='C' && bb=='M') return 2;
//    if (aa=='C' && bb=='F') return 2;
//    if (aa=='C' && bb=='P') return 1.9;
//    if (aa=='C' && bb=='S') return 1.4;
//    if (aa=='C' && bb=='T') return 1.9;
//    if (aa=='C' && bb=='W') return 2;
//    if (aa=='C' && bb=='Y') return 1.6;
//    if (aa=='C' && bb=='V') return 1.7;
//    if (aa=='Q' && bb=='A') return 1.6;
//    if (aa=='Q' && bb=='R') return 1.2;
//    if (aa=='Q' && bb=='N') return 1.5;
//    if (aa=='Q' && bb=='D') return 1.4;
//    if (aa=='Q' && bb=='C') return 2;
//    if (aa=='Q' && bb=='E') return 0.7;
//    if (aa=='Q' && bb=='G') return 1.9;
//    if (aa=='Q' && bb=='H') return 0.8;
//    if (aa=='Q' && bb=='I') return 1.9;
//    if (aa=='Q' && bb=='L') return 1.6;
//    if (aa=='Q' && bb=='K') return 1.4;
//    if (aa=='Q' && bb=='M') return 1.5;
//    if (aa=='Q' && bb=='F') return 2;
//    if (aa=='Q' && bb=='P') return 1.4;
//    if (aa=='Q' && bb=='S') return 1.7;
//    if (aa=='Q' && bb=='T') return 1.7;
//    if (aa=='Q' && bb=='W') return 2;
//    if (aa=='Q' && bb=='Y') return 2;
//    if (aa=='Q' && bb=='V') return 1.9;
//    if (aa=='E' && bb=='A') return 1.2;
//    if (aa=='E' && bb=='R') return 2;
//    if (aa=='E' && bb=='N') return 1.3;
//    if (aa=='E' && bb=='D') return 0.2;
//    if (aa=='E' && bb=='C') return 2;
//    if (aa=='E' && bb=='Q') return 0.6;
//    if (aa=='E' && bb=='G') return 1.5;
//    if (aa=='E' && bb=='H') return 1.7;
//    if (aa=='E' && bb=='I') return 1.6;
//    if (aa=='E' && bb=='L') return 1.9;
//    if (aa=='E' && bb=='K') return 1.5;
//    if (aa=='E' && bb=='M') return 1.9;
//    if (aa=='E' && bb=='F') return 2;
//    if (aa=='E' && bb=='P') return 1.6;
//    if (aa=='E' && bb=='S') return 1.5;
//    if (aa=='E' && bb=='T') return 1.7;
//    if (aa=='E' && bb=='W') return 2;
//    if (aa=='E' && bb=='Y') return 1.9;
//    if (aa=='E' && bb=='V') return 1.7;
//    if (aa=='G' && bb=='A') return 0.9;
//    if (aa=='G' && bb=='R') return 1.9;
//    if (aa=='G' && bb=='N') return 1.1;
//    if (aa=='G' && bb=='D') return 1.2;
//    if (aa=='G' && bb=='C') return 1.9;
//    if (aa=='G' && bb=='Q') return 1.6;
//    if (aa=='G' && bb=='E') return 1.3;
//    if (aa=='G' && bb=='H') return 1.9;
//    if (aa=='G' && bb=='I') return 2;
//    if (aa=='G' && bb=='L') return 1.9;
//    if (aa=='G' && bb=='K') return 1.7;
//    if (aa=='G' && bb=='M') return 1.9;
//    if (aa=='G' && bb=='F') return 1.9;
//    if (aa=='G' && bb=='P') return 1.6;
//    if (aa=='G' && bb=='S') return 0.9;
//    if (aa=='G' && bb=='T') return 1.6;
//    if (aa=='G' && bb=='W') return 2;
//    if (aa=='G' && bb=='Y') return 2;
//    if (aa=='G' && bb=='V') return 1.4;
//    if (aa=='H' && bb=='A') return 1.9;
//    if (aa=='H' && bb=='R') return 1.3;
//    if (aa=='H' && bb=='N') return 1;
//    if (aa=='H' && bb=='D') return 1.6;
//    if (aa=='H' && bb=='C') return 1.9;
//    if (aa=='H' && bb=='Q') return 0.9;
//    if (aa=='H' && bb=='E') return 1.9;
//    if (aa=='H' && bb=='G') return 2;
//    if (aa=='H' && bb=='I') return 2;
//    if (aa=='H' && bb=='L') return 1.9;
//    if (aa=='H' && bb=='K') return 1.9;
//    if (aa=='H' && bb=='M') return 2;
//    if (aa=='H' && bb=='F') return 1.7;
//    if (aa=='H' && bb=='P') return 1.6;
//    if (aa=='H' && bb=='S') return 1.9;
//    if (aa=='H' && bb=='T') return 1.9;
//    if (aa=='H' && bb=='W') return 1.9;
//    if (aa=='H' && bb=='Y') return 1.5;
//    if (aa=='H' && bb=='V') return 1.9;
//    if (aa=='I' && bb=='A') return 1.7;
//    if (aa=='I' && bb=='R') return 1.7;
//    if (aa=='I' && bb=='N') return 1.6;
//    if (aa=='I' && bb=='D') return 1.9;
//    if (aa=='I' && bb=='C') return 1.7;
//    if (aa=='I' && bb=='Q') return 1.9;
//    if (aa=='I' && bb=='E') return 1.7;
//    if (aa=='I' && bb=='G') return 2;
//    if (aa=='I' && bb=='H') return 2;
//    if (aa=='I' && bb=='L') return 1.2;
//    if (aa=='I' && bb=='K') return 1.7;
//    if (aa=='I' && bb=='M') return 1.1;
//    if (aa=='I' && bb=='F') return 1.3;
//    if (aa=='I' && bb=='P') return 2;
//    if (aa=='I' && bb=='S') return 1.9;
//    if (aa=='I' && bb=='T') return 1.3;
//    if (aa=='I' && bb=='W') return 2;
//    if (aa=='I' && bb=='Y') return 1.9;
//    if (aa=='I' && bb=='V') return 0.6;
//    if (aa=='L' && bb=='A') return 1.6;
//    if (aa=='L' && bb=='R') return 1.9;
//    if (aa=='L' && bb=='N') return 1.6;
//    if (aa=='L' && bb=='D') return 2;
//    if (aa=='L' && bb=='C') return 2;
//    if (aa=='L' && bb=='Q') return 1.4;
//    if (aa=='L' && bb=='E') return 1.9;
//    if (aa=='L' && bb=='G') return 1.9;
//    if (aa=='L' && bb=='H') return 1.5;
//    if (aa=='L' && bb=='I') return 0.8;
//    if (aa=='L' && bb=='K') return 1.7;
//    if (aa=='L' && bb=='M') return 0.4;
//    if (aa=='L' && bb=='F') return 1.1;
//    if (aa=='L' && bb=='P') return 1.6;
//    if (aa=='L' && bb=='S') return 1.9;
//    if (aa=='L' && bb=='T') return 1.6;
//    if (aa=='L' && bb=='W') return 1.5;
//    if (aa=='L' && bb=='Y') return 1.7;
//    if (aa=='L' && bb=='V') return 1;
//    if (aa=='K' && bb=='A') return 1.7;
//    if (aa=='K' && bb=='R') return 0.5;
//    if (aa=='K' && bb=='N') return 0.8;
//    if (aa=='K' && bb=='D') return 1.4;
//    if (aa=='K' && bb=='C') return 2;
//    if (aa=='K' && bb=='Q') return 1.1;
//    if (aa=='K' && bb=='E') return 1.3;
//    if (aa=='K' && bb=='G') return 1.7;
//    if (aa=='K' && bb=='H') return 1.7;
//    if (aa=='K' && bb=='I') return 1.5;
//    if (aa=='K' && bb=='L') return 1.9;
//    if (aa=='K' && bb=='M') return 0.9;
//    if (aa=='K' && bb=='F') return 2;
//    if (aa=='K' && bb=='P') return 1.6;
//    if (aa=='K' && bb=='S') return 1.3;
//    if (aa=='K' && bb=='T') return 1.2;
//    if (aa=='K' && bb=='W') return 2;
//    if (aa=='K' && bb=='Y') return 1.9;
//    if (aa=='K' && bb=='V') return 1.9;
//    if (aa=='M' && bb=='A') return 1.9;
//    if (aa=='M' && bb=='R') return 1.9;
//    if (aa=='M' && bb=='N') return 2;
//    if (aa=='M' && bb=='D') return 2;
//    if (aa=='M' && bb=='C') return 2;
//    if (aa=='M' && bb=='Q') return 1.7;
//    if (aa=='M' && bb=='E') return 2;
//    if (aa=='M' && bb=='G') return 2;
//    if (aa=='M' && bb=='H') return 2;
//    if (aa=='M' && bb=='I') return 1.4;
//    if (aa=='M' && bb=='L') return 1.3;
//    if (aa=='M' && bb=='K') return 1.5;
//    if (aa=='M' && bb=='F') return 1.9;
//    if (aa=='M' && bb=='P') return 2;
//    if (aa=='M' && bb=='S') return 1.9;
//    if (aa=='M' && bb=='T') return 1.7;
//    if (aa=='M' && bb=='W') return 2;
//    if (aa=='M' && bb=='Y') return 2;
//    if (aa=='M' && bb=='V') return 1.5;
//    if (aa=='F' && bb=='A') return 1.9;
//    if (aa=='F' && bb=='R') return 1.9;
//    if (aa=='F' && bb=='N') return 1.9;
//    if (aa=='F' && bb=='D') return 2;
//    if (aa=='F' && bb=='C') return 2;
//    if (aa=='F' && bb=='Q') return 2;
//    if (aa=='F' && bb=='E') return 2;
//    if (aa=='F' && bb=='G') return 1.9;
//    if (aa=='F' && bb=='H') return 1.7;
//    if (aa=='F' && bb=='I') return 1.3;
//    if (aa=='F' && bb=='L') return 1.4;
//    if (aa=='F' && bb=='K') return 2;
//    if (aa=='F' && bb=='M') return 1.5;
//    if (aa=='F' && bb=='P') return 2;
//    if (aa=='F' && bb=='S') return 1.7;
//    if (aa=='F' && bb=='T') return 1.9;
//    if (aa=='F' && bb=='W') return 1.6;
//    if (aa=='F' && bb=='Y') return 0.7;
//    if (aa=='F' && bb=='V') return 2;
//    if (aa=='P' && bb=='A') return 1.1;
//    if (aa=='P' && bb=='R') return 1.4;
//    if (aa=='P' && bb=='N') return 1.7;
//    if (aa=='P' && bb=='D') return 1.9;
//    if (aa=='P' && bb=='C') return 1.9;
//    if (aa=='P' && bb=='Q') return 1.3;
//    if (aa=='P' && bb=='E') return 1.6;
//    if (aa=='P' && bb=='G') return 1.7;
//    if (aa=='P' && bb=='H') return 1.4;
//    if (aa=='P' && bb=='I') return 1.9;
//    if (aa=='P' && bb=='L') return 1.7;
//    if (aa=='P' && bb=='K') return 1.7;
//    if (aa=='P' && bb=='M') return 1.9;
//    if (aa=='P' && bb=='F') return 1.9;
//    if (aa=='P' && bb=='S') return 1.1;
//    if (aa=='P' && bb=='T') return 1.5;
//    if (aa=='P' && bb=='W') return 2;
//    if (aa=='P' && bb=='Y') return 2;
//    if (aa=='P' && bb=='V') return 1.7;
//    if (aa=='S' && bb=='A') return 0.7;
//    if (aa=='S' && bb=='R') return 1.2;
//    if (aa=='S' && bb=='N') return 0.6;
//    if (aa=='S' && bb=='D') return 1.3;
//    if (aa=='S' && bb=='C') return 1.2;
//    if (aa=='S' && bb=='Q') return 1.5;
//    if (aa=='S' && bb=='E') return 1.4;
//    if (aa=='S' && bb=='G') return 1;
//    if (aa=='S' && bb=='H') return 1.7;
//    if (aa=='S' && bb=='I') return 1.7;
//    if (aa=='S' && bb=='L') return 1.9;
//    if (aa=='S' && bb=='K') return 1.3;
//    if (aa=='S' && bb=='M') return 1.5;
//    if (aa=='S' && bb=='F') return 1.6;
//    if (aa=='S' && bb=='P') return 1;
//    if (aa=='S' && bb=='T') return 0.4;
//    if (aa=='S' && bb=='W') return 1.4;
//    if (aa=='S' && bb=='Y') return 1.7;
//    if (aa=='S' && bb=='V') return 1.7;
//    if (aa=='T' && bb=='A') return 0.8;
//    if (aa=='T' && bb=='R') return 1.7;
//    if (aa=='T' && bb=='N') return 1.1;
//    if (aa=='T' && bb=='D') return 1.5;
//    if (aa=='T' && bb=='C') return 1.9;
//    if (aa=='T' && bb=='Q') return 1.6;
//    if (aa=='T' && bb=='E') return 1.7;
//    if (aa=='T' && bb=='G') return 1.7;
//    if (aa=='T' && bb=='H') return 1.9;
//    if (aa=='T' && bb=='I') return 1.2;
//    if (aa=='T' && bb=='L') return 1.7;
//    if (aa=='T' && bb=='K') return 1.3;
//    if (aa=='T' && bb=='M') return 1.4;
//    if (aa=='T' && bb=='F') return 1.9;
//    if (aa=='T' && bb=='P') return 1.4;
//    if (aa=='T' && bb=='S') return 0.7;
//    if (aa=='T' && bb=='W') return 2;
//    if (aa=='T' && bb=='Y') return 1.7;
//    if (aa=='T' && bb=='V') return 1.2;
//    if (aa=='W' && bb=='A') return 2;
//    if (aa=='W' && bb=='R') return 1.7;
//    if (aa=='W' && bb=='N') return 2;
//    if (aa=='W' && bb=='D') return 2;
//    if (aa=='W' && bb=='C') return 2;
//    if (aa=='W' && bb=='Q') return 2;
//    if (aa=='W' && bb=='E') return 2;
//    if (aa=='W' && bb=='G') return 2;
//    if (aa=='W' && bb=='H') return 2;
//    if (aa=='W' && bb=='I') return 2;
//    if (aa=='W' && bb=='L') return 2;
//    if (aa=='W' && bb=='K') return 2;
//    if (aa=='W' && bb=='M') return 2;
//    if (aa=='W' && bb=='F') return 1.9;
//    if (aa=='W' && bb=='P') return 2;
//    if (aa=='W' && bb=='S') return 1.9;
//    if (aa=='W' && bb=='T') return 2;
//    if (aa=='W' && bb=='Y') return 1.9;
//    if (aa=='W' && bb=='V') return 2;
//    if (aa=='Y' && bb=='A') return 1.9;
//    if (aa=='Y' && bb=='R') return 2;
//    if (aa=='Y' && bb=='N') return 1.6;
//    if (aa=='Y' && bb=='D') return 2;
//    if (aa=='Y' && bb=='C') return 1.6;
//    if (aa=='Y' && bb=='Q') return 2;
//    if (aa=='Y' && bb=='E') return 1.9;
//    if (aa=='Y' && bb=='G') return 2;
//    if (aa=='Y' && bb=='H') return 1.5;
//    if (aa=='Y' && bb=='I') return 1.9;
//    if (aa=='Y' && bb=='L') return 1.9;
//    if (aa=='Y' && bb=='K') return 2;
//    if (aa=='Y' && bb=='M') return 2;
//    if (aa=='Y' && bb=='F') return 0.9;
//    if (aa=='Y' && bb=='P') return 2;
//    if (aa=='Y' && bb=='S') return 1.9;
//    if (aa=='Y' && bb=='T') return 1.9;
//    if (aa=='Y' && bb=='W') return 1.7;
//    if (aa=='Y' && bb=='V') return 1.9;
//    if (aa=='V' && bb=='A') return 1.1;
//    if (aa=='V' && bb=='R') return 1.7;
//    if (aa=='V' && bb=='N') return 1.9;
//    if (aa=='V' && bb=='D') return 1.9;
//    if (aa=='V' && bb=='C') return 1.6;
//    if (aa=='V' && bb=='Q') return 1.7;
//    if (aa=='V' && bb=='E') return 1.7;
//    if (aa=='V' && bb=='G') return 1.6;
//    if (aa=='V' && bb=='H') return 1.6;
//    if (aa=='V' && bb=='I') return 0.2;
//    if (aa=='V' && bb=='L') return 1.2;
//    if (aa=='V' && bb=='K') return 1.9;
//    if (aa=='V' && bb=='M') return 1;
//    if (aa=='V' && bb=='F') return 1.9;
//    if (aa=='V' && bb=='P') return 1.6;
//    if (aa=='V' && bb=='S') return 1.7;
//    if (aa=='V' && bb=='T') return 1.2;
//    if (aa=='V' && bb=='W') return 2;
//    if (aa=='V' && bb=='Y') return 1.7;
//    return 2.5;
//}

