#include "stdafx.h"
#include "shared_funcs.h"

using namespace freicore;

namespace std
{
	ostream& operator<< ( ostream& o, const ProteinLocusByIndex& rhs )
	{
		return ( o << "(" << rhs.index << ", " << rhs.offset << ")" );
	}

	ostream& operator<< ( ostream& o, const ProteinLocusByName& rhs )
	{
		return ( o << "(" << rhs.name << ", " << rhs.offset << ")" );
	}

	ostream& operator<< ( ostream& o, const CleavageRule& rhs )
	{
		return ( o << rhs.first << " " << rhs.second );
	}

	istream& operator>> ( istream& i, CleavageRule& rhs )
	{
		static boost::char_separator<char> cleavageCandidatesDelimiter("|");

		rhs.first.clear();
		rhs.second.clear();

		string preCleavageCandidatesString, postCleavageCandidatesString;
		i >> preCleavageCandidatesString >> postCleavageCandidatesString;

		tokenizer preCleavageCandidatesParser(	preCleavageCandidatesString.begin(),
												preCleavageCandidatesString.end(),
												cleavageCandidatesDelimiter );
		for( tokenizer::iterator itr = preCleavageCandidatesParser.begin(); itr != preCleavageCandidatesParser.end(); ++itr )
			rhs.first.insert( *itr );

		tokenizer postCleavageCandidatesParser(	postCleavageCandidatesString.begin(),
												postCleavageCandidatesString.end(),
												cleavageCandidatesDelimiter );
		for( tokenizer::iterator itr = postCleavageCandidatesParser.begin(); itr != postCleavageCandidatesParser.end(); ++itr )
			rhs.second.insert( *itr );

		return i;
	}

	istream& operator>> ( istream& i, CleavageRuleSet& rhs )
	{
		CleavageRule newRule;
		while( i >> newRule )
		{
			rhs.push_back( newRule );
			rhs.longestPreCleavageCandidate = max( newRule.first.longestCleavageCandidate, rhs.longestPreCleavageCandidate );
			rhs.longestPostCleavageCandidate = max( newRule.second.longestCleavageCandidate, rhs.longestPostCleavageCandidate );
		}

		return i;
	}

	ostream& operator<< ( ostream& o, MvIntKey& rhs )
	{
		return o << vector<int>( rhs );
	}

	ostream& operator<< ( ostream& o, MvhTable& rhs )
	{
		for( MvhTable::iterator itr = rhs.begin(); itr != rhs.end(); ++itr )
			o << itr->first << " -> " << itr->second << "\n";
		return o;
	}
}

namespace freicore
{
	void CleavageRuleSet::initialize( const string& cfgStr )
	{
		clear();
		stringstream CleavageRulesStream( cfgStr );
		CleavageRulesStream >> *this;
	}

	bool	TestPtmSite( const string& seq, const DynamicMod& mod, size_t pos )
	{
		if( seq[pos] != mod.unmodChar )
			return false;

		for( size_t i=1; i <= mod.NTerminalFilters.size(); ++i )
		{
			if( pos-i < 0 )
				return false;

			const ResidueFilter& filter = *(mod.NTerminalFilters.rbegin()+(i-1));
			if( !filter.testResidue( seq[pos-i] ) )
				return false;
		}

		for( size_t i=1; i <= mod.CTerminalFilters.size(); ++i )
		{
			if( pos+i >= seq.length() )
				return false;

			const ResidueFilter& filter = mod.CTerminalFilters[i-1];
			if( !filter.testResidue( seq[pos+i] ) )
				return false;
		}
		return true;
	}

	void	MakePtmVariants( const CandidateSequenceInfo& c, CandidateSequenceList& ptmSequences, int maxPtmCount, ResidueMap* residueMap )
	{
		/*cout << "PTM map: ";
		for( int i=0; i < (int) dynamicMods.size(); ++i )
			cout << dynamicMods[i].unmodChar << " -> " << dynamicMods[i].modChar << "; ";
		cout << endl;*/

		int possiblePtmCount = 0;
		for( size_t i=0; i < c.sequence.size(); ++i )
			for( DynamicModSet::iterator itr = residueMap->dynamicMods.begin(); itr != residueMap->dynamicMods.end(); ++itr )
				if( TestPtmSite( c.sequence, *itr, i ) )
					++possiblePtmCount;

		int ptmPermutations = 1 << possiblePtmCount;

		//cout << "Possible PTMs in sequence " << seq << ": " << possiblePtmCount << endl;
		//cout << "Number of permutations: " << ptmPermutations << endl;

		for( int p=0; p < ptmPermutations; ++p )
		{
			int ptmCount = 0;
			/*int ptmInt = p;
			for( int i=0; i < sizeof( ptmPermutations ) * 8; ++i )
			{
				ptmCount += ptmInt & 1;
				ptmInt = ptmInt >> 1;
			}*/
			unsigned int const w = p - ((p >> 1) & 0x55555555);                    // temp
			unsigned int const x = (w & 0x33333333) + ((w >> 2) & 0x33333333);     // temp
			ptmCount = ((x + (x >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count

			if( ptmCount > maxPtmCount )
				continue;

			CandidateSequenceInfo ptmCandidate(c);
			ptmCount = 1;
			for( size_t i=0; i < ptmCandidate.sequence.size(); ++i )
			{
				for( DynamicModSet::iterator itr = residueMap->dynamicMods.begin(); itr != residueMap->dynamicMods.end(); ++itr )
					if( TestPtmSite( c.sequence, *itr, i ) )
					{
						if( p & ptmCount )
						{
							ptmCandidate.sequence[i] = itr->uniqueModChar;
							ptmCandidate.mass += itr->modMass;
						}
						ptmCount = ptmCount << 1;
					}
			}
			ptmSequences.push_back( ptmCandidate );
		}
	}

	string GetUnmodifiedSequence( const string& ptmSeq, ResidueMap* residueMap )
	{
		string seq = ptmSeq;
		for( size_t i=0; i < seq.length(); ++i )
		{
			if( seq[i] == '+' )
			{
				seq.erase( i, 1 );
				--i;
			} else
			{
				for( DynamicModSet::iterator itr = residueMap->dynamicMods.begin(); itr != residueMap->dynamicMods.end(); ++itr )
				{
					if( seq[i] == itr->uniqueModChar )
						seq[i] = itr->unmodChar;
				}
			}
		}

		return seq;
	}

	// strips off terminus symbols
	string GetRawSequence( const string& ptmSeq )
	{
		return ptmSeq.substr( 1, ptmSeq.length()-2 );
	}

	string ConvertFreiPtmToSqtPtm( const string& ptmSeq, ResidueMap* residueMap )
	{
		string seq = ptmSeq;
		if( residueMap )
		{
			for( size_t i=0; i < seq.length(); ++i )
			{
				for( DynamicModSet::iterator itr = residueMap->dynamicMods.begin(); itr != residueMap->dynamicMods.end(); ++itr )
					if( itr->uniqueModChar == seq[i] )
					{
						// put the unmodChar in front of the modChar
						seq[i] = itr->userModChar;
						char buf[2] = { itr->unmodChar, '\0' };
						seq.insert( i, buf );
						++i;

						// consolidate all modifications of the same mass difference to use a single modChar
						/*for( int k=0; k < (int) dynamicMods.size(); ++k )
							if( dynamicMods[k].modMass == dynamicMods[j].modMass )
							{
								seq[i] = dynamicMods[k].modChar;
								break;
							}*/
					}
			}
		}
		return seq;
	}

	string ConvertSqtPtmToFreiPtm( const string& ptmSeq, ResidueMap* residueMap )
	{
		if( !residueMap )
			return "";
		return ConvertSqtPtmToFreiPtm( ptmSeq, residueMap->dynamicMods );
	}

	string ConvertSqtPtmToFreiPtm( const string& ptmSeq, DynamicModSet& dynamicMods )
	{
		string seq = ptmSeq;
		for( size_t i=1; i < seq.length(); ++i )
		{
			if( seq[i] == '+' )
				continue;

			for( DynamicModSet::iterator itr = dynamicMods.begin(); itr != dynamicMods.end(); ++itr )
			{
				if( itr->userModChar == seq[i] && itr->unmodChar == seq[i-1] )
				{
					seq[i] = itr->uniqueModChar;
					seq.erase( i-1, 1 );
					--i;
					break;
				}
			}
		}
		return seq;
	}

	void CommitCommonDatatypes()
	{
	#ifdef USE_MPI
		/*MPI_Aint lbAddress, ubAddress;

		BaseSpectrum s[2];
		MPI_Aint spectrumAddresses[5];
		MPI_Address( s,								&lbAddress );
		MPI_Address( &s[0].id.index,						&spectrumAddresses[0] );	// INTS
		MPI_Address( &s[0].mzOfPrecursor,			&spectrumAddresses[1] );	// FLOATS
		MPI_Address( &s[0].mzDataIsDoublePrecision,	&spectrumAddresses[2] );	// BYTES

		MPI_Address( s+1,							&ubAddress );

		int flatSpectrumTypeBlockLengths[5] = { 1,6,9,2,1 };
		MPI_Aint flatSpectrumTypeDisplacements[5] = { lbAddress, 0,0,0, ubAddress };
		for( int i=0; i < 3; ++i )
			flatSpectrumTypeDisplacements[i+1] = spectrumAddresses[i] - lbAddress;
		MPI_Datatype flatSpectrumTypeList[5] = { MPI_LB, MPI_INT, MPI_FLOAT, MPI_BYTE, MPI_UB };
		MPI_Type_struct( 5, flatSpectrumTypeBlockLengths, flatSpectrumTypeDisplacements, flatSpectrumTypeList, &mpi_flatSpectrum );
		MPI_Type_commit( &mpi_flatSpectrum );*/
	#endif
	}

	double poz( double z )
	{
		double y, x, w;
		double Z_MAX = 6.0; 

		if (z == 0.0) {
			x = 0.0;
		} else {
			y = 0.5 * fabs(z);
			if (y >= (Z_MAX * 0.5)) {
				x = 1.0;
			} else if (y < 1.0) {
				w = y * y;
				x = ((((((((0.000124818987 * w
					- 0.001075204047) * w + 0.005198775019) * w
					- 0.019198292004) * w + 0.059054035642) * w
					- 0.151968751364) * w + 0.319152932694) * w
					- 0.531923007300) * w + 0.797884560593) * y * 2.0;
			} else {
				y -= 2.0;
				x = (((((((((((((-0.000045255659 * y
					+ 0.000152529290) * y - 0.000019538132) * y
					- 0.000676904986) * y + 0.001390604284) * y
					- 0.000794620820) * y - 0.002034254874) * y
					+ 0.006549791214) * y - 0.010557625006) * y
					+ 0.011630447319) * y - 0.009279453341) * y
					+ 0.005353579108) * y - 0.002141268741) * y
					+ 0.000535310849) * y + 0.999936657524;
			}
		}
		return z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5);
	}

	static int BIGX = 20;
	static double LOG_SQRT_PI = 0.5723649429247000870717135; /* log(sqrt(pi)) */
	static double I_SQRT_PI = 0.5641895835477562869480795;   /* 1 / sqrt(pi) */

	double ex( double x )
	{
		return (x < -BIGX) ? 0.0 : exp(x);
	}

	float ChiSquaredToPValue( float x, int df )
	{
		double a, y=0.0, s;
		double e, c, z;
		bool even;


		if (x <= 0.0 || df < 1) {
			return 1.0;
		}

		a = 0.5 * x;
		even = !(df & 1);
		if (df > 1) {
			y = ex(-a);
		}
		s = (even ? y : (2.0 * poz(-sqrt(x))));
		if (df > 2) {
			x = 0.5f * float(df - 1.0);
			z = (even ? 1.0 : 0.5);
			if (a > BIGX) {
				e = (even ? 0.0 : LOG_SQRT_PI);
				c = log(a);
				while (z <= x) {
					e = log(z) + e;
					s += ex(c * z - a - e);
					z += 1.0;
				}
				return (float)s;
			} else {
				e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
				c = 0.0;
				while (z <= x) {
					e = e * (a / z);
					c = c + e;
					z += 1.0;
				}
				return float(c * y + s);
			}
		} else {
			return (float)s;
		}
	}

	void LoadTagInstancesFromIndexFile( const string& tagIndexFilename, const string& tag, tagMetaIndex_t& tagMetaIndex, tagIndex_t& tagIndex )
	{
		if( tagMetaIndex.count( tag ) > 0 )
		{
			ifstream tagIndexFile( tagIndexFilename.c_str(), ios::binary );
			tagIndexFile.seekg( (ios::off_type) tagMetaIndex[tag].offset );

			tagIndexFile.read( (char*) &tagMetaIndex[tag].size, sizeof( tagMetaIndex[tag].size ) );

			ProteinIndex idx;
			ProteinOffset off;
			for( int i=0; i < tagMetaIndex[tag].size; ++i )
			{
				tagIndexFile.read( (char*) &idx, sizeof( idx ) );
				tagIndexFile.read( (char*) &off, sizeof( off ) );
				//cout << string(tag) << " " << (int) instanceCount << " " << idx << " " << off << endl;
				tagIndex[ tag ].push_back( tagInstance_t( idx, off ) );
			}
			tagIndexFile.close();
		} else
		{
			cerr << "Tag \"" << tag << "\" not in index!" << endl;
		}
	}

	void LoadIndexFromFile( const string& tagIndexFilename, tagMetaIndex_t& tagMetaIndex )
	{
		ifstream tagIndexFile;

		tagIndexFile.open( tagIndexFilename.c_str(), ios::in | ios::binary );

		tagIndexFile.seekg( 40 ); // skip the checksum

		int tagCount;
		tagIndexFile.read( (char*) &tagCount, sizeof( tagCount ) );
		tagIndexFile.read( (char*) &tagMetaIndex.totalTagInstances, sizeof( tagMetaIndex.totalTagInstances ) );

		for( int i=0; i < tagCount; ++i )
		{
			/* ASCII-style index
			string tag;
			tagIndexFile >> tag;

			int instanceCount;
			tagIndexFile >> instanceCount;

			for( int i=0; i < instanceCount; ++i )
			{
				int idx, off;
				tagIndexFile >> idx >> off;
				tagIndex[ tag ].push_back( tagInstance_t( idx, off ) );
			}
			*/

			/* Binary-style index */
			char* tag = new char[ 4 ];
			tagIndexFile.read( tag, 3 );
			tag[3] = 0;

			tagIndexFile.read( (char*) &tagMetaIndex[tag].offset, sizeof( tagMetaIndex[tag].offset ) );
			ios::pos_type cur = tagIndexFile.tellg();
			tagIndexFile.seekg( (ios::off_type) tagMetaIndex[tag].offset );
			tagIndexFile.read( (char*) &tagMetaIndex[tag].size, sizeof( tagMetaIndex[tag].size ) );
			tagIndexFile.seekg( cur );
			tagMetaIndex[tag].proportion = (float) tagMetaIndex[tag].size / (float) tagMetaIndex.totalTagInstances;
			delete tag;
		}
		tagIndexFile.close();
	}

	void	CalculateSequenceIons(	const string& seq,
									int maxIonCharge,
									vector< float >* sequenceIonMasses,
									bool useSmartPlusThreeModel,
									vector< string >* sequenceIonLabels,
									float precursorMass, const bool ionTypes[4], ResidueMap* residueMap )
	{
		if( !sequenceIonMasses )
			return;

		int yPeaks = (int) seq.length()-2;
		int bPeaks = (int) seq.length()-2;
		size_t terminusIndex = seq.length()-1;

		sequenceIonMasses->clear();
		sequenceIonMasses->reserve( yPeaks + bPeaks );

		if( sequenceIonLabels )
		{
			sequenceIonLabels->clear();
			sequenceIonLabels->reserve( yPeaks + bPeaks );
		}

		bool nTerminusIsPartial = ( *seq.begin() == '-' );
		bool cTerminusIsPartial = ( *seq.rbegin() == '-' );

		// calculate y ion MZs
		if( maxIonCharge > 2 )
		{
			if( useSmartPlusThreeModel )
			{
				size_t totalStrongBasicCount = 0, totalWeakBasicCount = 0;
				for( size_t i=1; i < terminusIndex; ++i )
					if( seq[i] == 'R' || seq[i] == 'K' || seq[i] == 'H' )
						++totalStrongBasicCount;
					else if( seq[i] == 'Q' || seq[i] == 'N' )
						++totalWeakBasicCount;
				size_t totalBasicity = totalStrongBasicCount * 4 + totalWeakBasicCount * 2 + seq.length()-2;

				map< float, int > basicityThresholds;
				basicityThresholds[ 0.0f ] = 1;
				for( int z = 1; z < maxIonCharge-1; ++z )
					basicityThresholds[ (float) z / (float) (maxIonCharge-1) ] = z+1;

				for( size_t c = 1; c < terminusIndex-1; ++c )
				{
					string bSeq = seq.substr( 1, c+1 );
					string ySeq = seq.substr( c+2, terminusIndex-c-2 );

					float bMass = residueMap->GetMassOfResidues( seq.begin(), seq.begin()+(c+2) ) - HYDROGEN_MONO;
					float yMass = residueMap->GetMassOfResidues( seq.begin()+(c+2), seq.end() ) + HYDROGEN_MONO;

					size_t bStrongBasicCount = 0, bWeakBasicCount = 0;
					for( size_t i=0; i < bSeq.length(); ++i )
						if( bSeq[i] == 'R' || bSeq[i] == 'K' || bSeq[i] == 'H' )
							++bStrongBasicCount;
						else if( bSeq[i] == 'Q' || bSeq[i] == 'N' )
							++bWeakBasicCount;

					size_t bScore = bStrongBasicCount * 4 + bWeakBasicCount * 2 + bSeq.length();

					float basicityRatio = (float) bScore / (float) totalBasicity;
					map< float, int >::iterator itr = basicityThresholds.upper_bound( basicityRatio );
					--itr;
					int bZ = itr->second;
					int yZ = maxIonCharge - bZ;

					//cout << "b" << c+1 << "(+" << bZ << ") <-> y" << seq.length()-(c+3) << "(+" << yZ << ")" << endl;
					//cout << bSeq << " <-> " << ySeq << endl;
					//cout << bMass << " <-> " << yMass << endl << endl;

					if( !bSeq.empty() && ionTypes[0] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "(+" + lexical_cast<string>(bZ) + ")" );
						sequenceIonMasses->push_back( ( bMass + ( bZ * PROTON ) ) / bZ );
					}

					if( !ySeq.empty() && ionTypes[1] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "(+" + lexical_cast<string>(yZ) + ")" );
						sequenceIonMasses->push_back( ( yMass + ( yZ * PROTON ) ) / yZ );
					}

					if( !bSeq.empty() && ionTypes[2] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "-H20 (+" + lexical_cast<string>(bZ) + ")" );
						sequenceIonMasses->push_back( ( bMass - WATER_MONO + ( bZ * PROTON ) ) / bZ );
					}

					if( !ySeq.empty() && ionTypes[2] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "-H20 (+" + lexical_cast<string>(yZ) + ")" );
						sequenceIonMasses->push_back( ( yMass - WATER_MONO + ( yZ * PROTON ) ) / yZ );
					}

					if( !bSeq.empty() && ionTypes[3] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "-NH3 (+" + lexical_cast<string>(bZ) + ")" );
						sequenceIonMasses->push_back( ( bMass - AMMONIA_MONO + ( bZ * PROTON ) ) / bZ );
					}

					if( !ySeq.empty() && ionTypes[3] )
					{
						if( sequenceIonLabels )
							sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "-NH3 (+" + lexical_cast<string>(yZ) + ")" );
						sequenceIonMasses->push_back( ( yMass - AMMONIA_MONO + ( yZ * PROTON ) ) / yZ );
					}
				}
			} else
			{
				for( int z = 1; z < maxIonCharge; ++z )
				{
					int bZ = z;
					int yZ = maxIonCharge - bZ;

					for( size_t c = 1; c < terminusIndex-2; ++c )
					{
						string bSeq = seq.substr( 1, c+1 );
						string ySeq = seq.substr( c+2, terminusIndex-c-2 );

						float bMass = residueMap->GetMassOfResidues( seq.begin(), seq.begin()+(c+2) ) - HYDROGEN_MONO;
						float yMass = residueMap->GetMassOfResidues( seq.begin()+(c+2), seq.end() ) + HYDROGEN_MONO;

						if( !bSeq.empty() && ionTypes[0] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "(+" + lexical_cast<string>(bZ) + ")" );
							sequenceIonMasses->push_back( ( bMass + ( bZ * PROTON ) ) / bZ );
						}

						if( !ySeq.empty() && ionTypes[1] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "(+" + lexical_cast<string>(yZ) + ")" );
							sequenceIonMasses->push_back( ( yMass + ( yZ * PROTON ) ) / yZ );
						}

						if( !bSeq.empty() && ionTypes[2] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "-H20 (+" + lexical_cast<string>(bZ) + ")" );
							sequenceIonMasses->push_back( ( bMass - WATER_MONO + ( bZ * PROTON ) ) / bZ );
						}

						if( !ySeq.empty() && ionTypes[2] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "-H20 (+" + lexical_cast<string>(yZ) + ")" );
							sequenceIonMasses->push_back( ( yMass - WATER_MONO + ( yZ * PROTON ) ) / yZ );
						}

						if( !bSeq.empty() && ionTypes[3] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("b") + lexical_cast<string>(c+1) + "-NH3 (+" + lexical_cast<string>(bZ) + ")" );
							sequenceIonMasses->push_back( ( bMass - AMMONIA_MONO + ( bZ * PROTON ) ) / bZ );
						}

						if( !ySeq.empty() && ionTypes[3] )
						{
							if( sequenceIonLabels )
								sequenceIonLabels->push_back( string("y") + lexical_cast<string>(seq.length()-(c+3)) + "-NH3 (+" + lexical_cast<string>(yZ) + ")" );
							sequenceIonMasses->push_back( ( yMass - AMMONIA_MONO + ( yZ * PROTON ) ) / yZ );
						}
					}
				}
			}
		} else
		{
			// calculate b ion MZs
			float sequenceMass = ( nTerminusIsPartial ? precursorMass : residueMap->getMonoMassByName( *seq.begin() ) ) - HYDROGEN_MONO;
			for( size_t b = 1; b < terminusIndex; ++b )
			{
				sequenceMass += residueMap->getMonoMassByName( seq[b] );

				if( ionTypes[0] )
				{
					sequenceIonMasses->push_back( sequenceMass + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("b") + lexical_cast<string>(b) );
				}

				if( ionTypes[2] )
				{
					sequenceIonMasses->push_back( sequenceMass - WATER_MONO + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("b") + lexical_cast<string>(b) + "-H20" );
				}

				if( ionTypes[3] )
				{
					sequenceIonMasses->push_back( sequenceMass - AMMONIA_MONO + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("b") + lexical_cast<string>(b) + "-NH3" );
				}
			}

			sequenceMass = ( cTerminusIsPartial ? precursorMass : residueMap->getMonoMassByName( *seq.rbegin() ) ) + HYDROGEN_MONO;
			for( size_t y = 1; y < terminusIndex; ++y )
			{
				sequenceMass += residueMap->getMonoMassByName( *(seq.rbegin()+y) );

				if( ionTypes[1] )
				{
					sequenceIonMasses->push_back( sequenceMass + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("y") + lexical_cast<string>(y) );
				}

				if( ionTypes[2] )
				{
					sequenceIonMasses->push_back( sequenceMass - WATER_MONO + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("y") + lexical_cast<string>(y) + "-H20" );
				}

				if( ionTypes[3] )
				{
					sequenceIonMasses->push_back( sequenceMass - AMMONIA_MONO + PROTON );
					if( sequenceIonLabels )
						sequenceIonLabels->push_back( string("y") + lexical_cast<string>(y) + "-NH3" );
				}
			}

		}
	}

	void	CalculateSequenceIons2(	const string& seq,
									int maxIonCharge,
									vector< float >* sequenceIonMasses,
									bool useSmartPlusThreeModel,
									vector< string >* sequenceIonLabels,
									float precursorMass, const bool ionTypes[4], ResidueMap* residueMap )
	{
		int yPeaks = (int) seq.length();
		int bPeaks = (int) seq.length();

		if( !sequenceIonMasses )
			return;

		sequenceIonMasses->clear();
		sequenceIonMasses->reserve( yPeaks + bPeaks );

		if( sequenceIonLabels )
		{
			sequenceIonLabels->clear();
			sequenceIonLabels->reserve( yPeaks + bPeaks );
		}

		bool nTerminusIsPartial = ( *seq.begin() == '-' );
		bool cTerminusIsPartial = ( *seq.rbegin() == '-' );

		// calculate y ion MZs
		if( maxIonCharge > 2 )
		{
			if( useSmartPlusThreeModel )
			{
				size_t totalStrongBasicCount = 0, totalWeakBasicCount = 0;
				for( size_t i=0; i < seq.length(); ++i )
					if( seq[i] == 'R' || seq[i] == 'K' || seq[i] == 'H' )
						++totalStrongBasicCount;
					else if( seq[i] == 'Q' || seq[i] == 'N' )
						++totalWeakBasicCount;
				size_t totalBasicity = totalStrongBasicCount * 4 + totalWeakBasicCount * 2 + seq.length();

				map< float, int > basicityThresholds;
				basicityThresholds[ 0.0f ] = 1;
				for( int z = 1; z < maxIonCharge-1; ++z )
					basicityThresholds[ (float) z / (float) (maxIonCharge-1) ] = z+1;

				for( size_t c = 1; c < seq.length(); ++c )
				{
					string bSeq = seq.substr( 0, c+1 );
					string ySeq = seq.substr( c+1 );
					
					float bMass = residueMap->GetMassOfResidues( bSeq );
					float yMass = residueMap->GetMassOfResidues( ySeq );

					size_t bStrongBasicCount = 0, bWeakBasicCount = 0;
					for( size_t i=0; i < bSeq.length(); ++i )
						if( bSeq[i] == 'R' || bSeq[i] == 'K' || bSeq[i] == 'H' )
							++bStrongBasicCount;
						else if( bSeq[i] == 'Q' || bSeq[i] == 'N' )
							++bWeakBasicCount;

					size_t bScore = bStrongBasicCount * 4 + bWeakBasicCount * 2 + bSeq.length();

					float basicityRatio = (float) bScore / (float) totalBasicity;
					map< float, int >::iterator itr = basicityThresholds.upper_bound( basicityRatio );
					--itr;
					int bZ = itr->second;
					int yZ = maxIonCharge - bZ;

					int bIsotopeCount = 1;
					if( !bSeq.empty() && ionTypes[0] )
						for( int n=0; n < bIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream bLabel;
								bLabel << "b" << c+1 << "(+" << bZ << ")";
								if( n > 0 )
									bLabel << "+" << n;
								sequenceIonLabels->push_back( bLabel.str() );
							}

							sequenceIonMasses->push_back( ( bMass + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
						}

					int yIsotopeCount = 1;
					if( !ySeq.empty() && ionTypes[1] )
						for( int n=0; n < yIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream yLabel;
								yLabel << "y" << seq.length()-(c+1) << "(+" << yZ << ")";
								if( n > 0 )
									yLabel << "+" << n;
								sequenceIonLabels->push_back( yLabel.str() );
							}

							sequenceIonMasses->push_back( ( yMass + WATER_MONO + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
						}

					if( !bSeq.empty() && ionTypes[2] )
						for( int n=0; n < bIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream bLabel;
								bLabel << "b" << c+1 << "-H2O (+" << bZ << ")";
								if( n > 0 )
									bLabel << "+" << n;
								sequenceIonLabels->push_back( bLabel.str() );
							}

							sequenceIonMasses->push_back( ( bMass - WATER_MONO + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
						}

					if( !ySeq.empty() && ionTypes[2] )
						for( int n=0; n < yIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream yLabel;
								yLabel << "y" << seq.length()-(c+1) << "-H2O (+" << yZ << ")";
								if( n > 0 )
									yLabel << "+" << n;
								sequenceIonLabels->push_back( yLabel.str() );
							}

							sequenceIonMasses->push_back( ( yMass + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
						}

					if( !bSeq.empty() && ionTypes[3] )
						for( int n=0; n < bIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream bLabel;
								bLabel << "b" << c+1 << "-NH3 (+" << bZ << ")";
								if( n > 0 )
									bLabel << "+" << n;
								sequenceIonLabels->push_back( bLabel.str() );
							}

							sequenceIonMasses->push_back( ( bMass - AMMONIA_MONO + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
						}

					if( !ySeq.empty() && ionTypes[3] )
						for( int n=0; n < yIsotopeCount; ++n )
						{
							if( sequenceIonLabels )
							{
								stringstream yLabel;
								yLabel << "y" << seq.length()-(c+1) << "-NH3 (+" << yZ << ")";
								if( n > 0 )
									yLabel << "+" << n;
								sequenceIonLabels->push_back( yLabel.str() );
							}

							sequenceIonMasses->push_back( ( yMass + WATER_MONO - AMMONIA_MONO + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
						}
				}
			} else
			{
				for( int z = 1; z < maxIonCharge; ++z )
				{
					int bZ = z;
					int yZ = maxIonCharge - bZ;

					for( size_t c = 1; c < seq.length(); ++c )
					{
						string bSeq = seq.substr( 0, c+1 );
						string ySeq = seq.substr( c+1 );

						float bMass = residueMap->GetMassOfResidues( bSeq );
						float yMass = residueMap->GetMassOfResidues( ySeq );

						int bIsotopeCount = 1;
						if( !bSeq.empty() && ionTypes[0] )
							for( int n=0; n < bIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream bLabel;
									bLabel << "b" << c+1 << "(+" << bZ << ")";
									if( n > 0 )
										bLabel << "+" << n;
									sequenceIonLabels->push_back( bLabel.str() );
								}

								sequenceIonMasses->push_back( ( bMass + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
							}

						int yIsotopeCount = 1;
						if( !ySeq.empty() && ionTypes[1] )
							for( int n=0; n < yIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream yLabel;
									yLabel << "y" << seq.length()-(c+1) << "(+" << yZ << ")";
									if( n > 0 )
										yLabel << "+" << n;
									sequenceIonLabels->push_back( yLabel.str() );
								}

								sequenceIonMasses->push_back( ( yMass + WATER_MONO + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
							}

						if( !bSeq.empty() && ionTypes[2] )
							for( int n=0; n < bIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream bLabel;
									bLabel << "b" << c+1 << "-H2O (+" << bZ << ")";
									if( n > 0 )
										bLabel << "+" << n;
									sequenceIonLabels->push_back( bLabel.str() );
								}

								sequenceIonMasses->push_back( ( bMass - WATER_MONO + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
							}

						if( !ySeq.empty() && ionTypes[2] )
							for( int n=0; n < yIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream yLabel;
									yLabel << "y" << seq.length()-(c+1) << "-H20 (+" << yZ << ")";
									if( n > 0 )
										yLabel << "+" << n;
									sequenceIonLabels->push_back( yLabel.str() );
								}

								sequenceIonMasses->push_back( ( yMass + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
							}

						if( !bSeq.empty() && ionTypes[3] )
							for( int n=0; n < bIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream bLabel;
									bLabel << "b" << c+1 << "-NH3 (+" << bZ << ")";
									if( n > 0 )
										bLabel << "+" << n;
									sequenceIonLabels->push_back( bLabel.str() );
								}

								sequenceIonMasses->push_back( ( bMass - AMMONIA_MONO + ( bZ * PROTON ) + ( n * NEUTRON ) ) / bZ );
							}

						if( !ySeq.empty() && ionTypes[3] )
							for( int n=0; n < yIsotopeCount; ++n )
							{
								if( sequenceIonLabels )
								{
									stringstream yLabel;
									yLabel << "y" << seq.length()-(c+1) << "-NH3 (+" << yZ << ")";
									if( n > 0 )
										yLabel << "+" << n;
									sequenceIonLabels->push_back( yLabel.str() );
								}

								sequenceIonMasses->push_back( ( yMass + WATER_MONO - AMMONIA_MONO + ( yZ * PROTON ) + ( n * NEUTRON ) ) / yZ );
							}
					}
				}
			}
		} else
		{
			// calculate b ion MZs
			for( int b = 1; b <= bPeaks; ++b )
			{
				float sequenceMass = ( nTerminusIsPartial ? precursorMass : 0 );
				if( false )
					for( int j = 0; j < b; ++j )
						sequenceMass += residueMap->getAvgMassByName( seq[j] );
				else
					for( int j = 0; j < b; ++j )
						sequenceMass += residueMap->getMonoMassByName( seq[j] );

				if( ionTypes[0] )
				{
					sequenceIonMasses->push_back( sequenceMass + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "b" << b;
						sequenceIonLabels->push_back( label.str() );
					}
				}

				if( ionTypes[2] )
				{
					sequenceIonMasses->push_back( sequenceMass - WATER_MONO + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "b" << b << "-H20";
						sequenceIonLabels->push_back( label.str() );
					}
				}

				if( ionTypes[3] )
				{
					sequenceIonMasses->push_back( sequenceMass - AMMONIA_MONO + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "b" << b << "-NH3";
						sequenceIonLabels->push_back( label.str() );
					}
				}
			}

			for( int y = 1; y <= yPeaks; ++y )
			{
				float sequenceMass = ( cTerminusIsPartial ? precursorMass : 0 );
				if( false )
					for( int j = (int) yPeaks-1; j >= (int) seq.length() - y; --j )
						sequenceMass += residueMap->getAvgMassByName( seq[j] );
				else
					for( int j = (int) yPeaks-1; j >= (int) seq.length() - y; --j )
						sequenceMass += residueMap->getMonoMassByName( seq[j] );

				if( ionTypes[1] )
				{
					sequenceIonMasses->push_back( sequenceMass + WATER_MONO + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "y" << y;
						sequenceIonLabels->push_back( label.str() );
					}
				}

				if( ionTypes[2] )
				{
					sequenceIonMasses->push_back( sequenceMass + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "y" << y << "-H20";
						sequenceIonLabels->push_back( label.str() );
					}
				}

				if( ionTypes[3] )
				{
					sequenceIonMasses->push_back( sequenceMass + WATER_MONO - AMMONIA_MONO + PROTON );
					if( sequenceIonLabels )
					{
						stringstream label;
						label << "y" << y << "-NH3";
						sequenceIonLabels->push_back( label.str() );
					}
				}
			}

		}
	}

	void	CreateScoringTableMVH_R(	const int minValue,
										const int totalValue,
										const int numClasses,
										const vector< int >& classCounts,
										MvhTable& mvTable,
										MvIntKey& key,
										lnFactorialTable& lnFT,
										const MvIntKey* minKey = NULL )
	{
		// At the highest degree of variability the key is fully set
		// Calculate the MVH score and add it to the mvTable
		if( numClasses == 1 )
		{
			key.front() = totalValue;
			//if( minKey == NULL || !( mvTable.comp( *minKey, key ) ) )
			//	return;

			int totalClasses = (int) key.size();
			double lnP = 0.0f;
			for( int i=0; i < totalClasses; ++i )
				lnP += lnCombin( classCounts[i], key[i], lnFT );
			//float p = 0.0f;
			//for( int i=0; i < totalClasses; ++i )
			//	p += lnCombin( classCounts[i], key[i], lnFT );
			int totalClassCount = accumulate( classCounts.begin(), classCounts.end(), 0 );
			int totalValueCount = accumulate( key.begin(), key.end(), 0 );
			lnP -= lnCombin( totalClassCount, totalValueCount, lnFT );
			START_PROFILER(9);
			mvTable[ key ] = lnP;
			STOP_PROFILER(9);

		// Create another level of variability
		} else
		{
			for( int curValue = minValue; (totalValue - curValue) >= minValue ; ++curValue )
			{
				key[numClasses-1] = curValue;
				CreateScoringTableMVH_R( minValue, totalValue - curValue, numClasses-1, classCounts, mvTable, key, lnFT, minKey );
			}
		}
	}

	void	CreateScoringTableMVH(	const int minValue,
									const int totalValue,
									const int numClasses,
									vector< int > classCounts,
									MvhTable& mvTable,
									lnFactorialTable& lnFT,
									bool normalizeOnMode,
									bool adjustRareOutcomes,
									bool convertToPValues,
									const MvIntKey* minKey )
	{
		// Check to see if all classes have a count of at least 1
		bool allClassesUsed = true;
		for( int i=0; i < numClasses; ++i )
		{
			if( classCounts[i] == 0 )
			{
				allClassesUsed = false;
				break;
			}
		}

		// If any class is not populated, increment each class by one
		if( !allClassesUsed )
			for( int i=0; i < numClasses; ++i )
				++ classCounts[i];

		MvIntKey key;
		key.resize( numClasses, 0 );
		START_PROFILER(10);
		CreateScoringTableMVH_R( minValue, totalValue, numClasses, classCounts, mvTable, key, lnFT, minKey );
		STOP_PROFILER(10);

		if( convertToPValues )
		{
			mvTable.ConvertToPValues();
		} else
		if( normalizeOnMode )
		{
			// Normalize on the mode value if desired
			MvhTable::iterator itr;
			MvIntKey modeKey = mvTable.begin()->first;
			double modeValue = mvTable.begin()->second;

			for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
			{
				if( modeValue < itr->second )
				{
					modeKey = itr->first;
					modeValue = itr->second;
				}
			}

			for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
				itr->second -= modeValue;

			if( adjustRareOutcomes )
			{
				// Prevent rare and undesirable outcomes from having good scores
				for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
				{
					key = itr->first;
					bool IsRareAndUndesirable = true;
					int numVarsToAdjust = (int) key.size() - 1;
					for( int i=0; i < numVarsToAdjust ; ++i )
					{
						if( key[i] > modeKey[i] )
						{
							IsRareAndUndesirable = false;
							break;
						}
					}

					if( IsRareAndUndesirable )
						itr->second = - itr->second;
				}
			}
		}
	}

	void	CreateScoringTableMVB_R(	const int minValue,
										const int totalValue,
										const int numClasses,
										const vector< float >& classProbabilities,
										MvhTable& mvTable,
										MvIntKey& key,
										lnFactorialTable& lnFT )
	{
		// At the highest degree of variability the key is fully set
		// Calculate the MVH score and add it to the mvTable
		if( numClasses == 1 )
		{
			key.front() = totalValue;
			int totalClasses = (int) key.size();
			int N = accumulate( key.begin(), key.end(), 0 );
			double sum1 = 0, sum2 = 0;
			for( int i=0; i < totalClasses; ++i )
			{
				sum1 += log( pow( classProbabilities[i], key[i] ) );
				sum2 += lnFT[ key[i] ];
			}
			mvTable[ key ] = ( lnFT[N] - sum2 ) + sum1;

		// Create another level of variability
		} else
		{
			for( int curValue = minValue; (totalValue - curValue) >= minValue ; ++curValue )
			{
				key[numClasses-1] = curValue;
				CreateScoringTableMVB_R( minValue, totalValue - curValue, numClasses-1, classProbabilities, mvTable, key, lnFT );
			}
		}
	}

	void	CreateScoringTableMVB(	const int minValue,
									const int totalValue,
									const int numClasses,
									const vector< float >& classProbabilities,
									MvhTable& mvTable,
									lnFactorialTable& lnFT,
									bool normalizeOnMode,
									bool adjustRareOutcomes )
	{
		MvIntKey key;
		key.resize( numClasses, 0 );
		CreateScoringTableMVB_R( minValue, totalValue, numClasses, classProbabilities, mvTable, key, lnFT );

		if( normalizeOnMode )
		{
			// Normalize on the mode value if desired
			MvhTable::iterator itr;
			MvIntKey modeKey = mvTable.begin()->first;
			double modeValue = mvTable.begin()->second;

			for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
			{
				if( modeValue < itr->second )
				{
					modeKey = itr->first;
					modeValue = itr->second;
				}
			}

			for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
				itr->second -= modeValue;

			if( adjustRareOutcomes )
			{
				// Prevent rare and undesirable outcomes from having good scores
				for( itr = mvTable.begin(); itr != mvTable.end(); ++itr )
				{
					key = itr->first;
					bool IsRareAndUndesirable = true;
					int numVarsToAdjust = (int) key.size() - 1;
					for( int i=0; i < numVarsToAdjust ; ++i )
					{
						if( key[i] > modeKey[i] )
						{
							IsRareAndUndesirable = false;
							break;
						}
					}

					if( IsRareAndUndesirable )
						itr->second = 0;
				}
			}
		}
	}

	
	int	ClassifyError( float error, vector< float > mzFidelityThresholds )
	{
		for( int i=0; i < (int) mzFidelityThresholds.size(); ++i )
		{
			if( error <= mzFidelityThresholds[i] )
				return i;
		}
		return (int) mzFidelityThresholds.size()-1;

		cout.precision(8);
		cout << "ClassifyError: could not classify error " << error << " in thresholds:\n";
		for( int i=0; i < (int) mzFidelityThresholds.size(); ++i )
			cout << mzFidelityThresholds[i] << " ";
		cout << endl;
		return 0;
	}

	double lnCombin( int n, int k, lnFactorialTable& lnTable )
	{
		if( n < 0 || k < 0 || n < k )
			return -1;

		try
		{
			return lnTable[n] - lnTable[n-k] - lnTable[k];
		} catch( exception& e )
		{
			cerr << "lnCombin(): caught exception with n=" << n << " and k=" << k << endl;
			throw e;
		}
	}

	float lnOdds( float p )
	{
		return log( p / (1 - p) );
	}

	int paramIndex( const string& param, const char** atts, int attsCount )
	{
		for( int i=0; i < attsCount; ++ i )
		{
			if( !strcmp( atts[i], param.c_str() ) )
				return i;
		}
		//cerr << "Attribute \"" << param << "\" required but not specified." << endl;
		return -1;
	}

	void FindFilesByMask( const string& mask, fileList_t& filenames )
	{
	#ifdef WIN32
		string maskPathname = GetPathnameFromFilepath( mask );
		WIN32_FIND_DATA fdata;
		HANDLE srcFile = FindFirstFileEx( mask.c_str(), FindExInfoStandard, &fdata, FindExSearchNameMatch, NULL, 0 );
		if( srcFile == INVALID_HANDLE_VALUE )
			return;

		do
		{
			filenames.insert( maskPathname + fdata.cFileName );
		} while( FindNextFile( srcFile, &fdata ) );

		FindClose( srcFile );

	#else

		glob_t globbuf;
		int rv = glob( mask.c_str(), 0, NULL, &globbuf );
		if( rv > 0 && rv != GLOB_NOMATCH )
			throw runtime_error( "FindFilesByMask(): glob() error" );

		DIR* curDir = opendir( "." );
		struct stat curEntryData;

		for( size_t i=0; i < globbuf.gl_pathc; ++i )
		{
			stat( globbuf.gl_pathv[i], &curEntryData );
			if( S_ISREG( curEntryData.st_mode ) )
				filenames.insert( globbuf.gl_pathv[i] );
		}
		closedir( curDir );

		globfree( &globbuf );

	#endif
	}

	endianType_t GetMachineEndianType()
	{
		int testInt = 127;
		char* testIntP = (char*) &testInt;

		if( testIntP[0] == 127 )
			return SYS_LITTLE_ENDIAN;
		else if( testIntP[ sizeof(int)-1 ] == 127 )
			return SYS_BIG_ENDIAN;
		else
			return SYS_UNKNOWN_ENDIAN;
	}
}
