#include "stdafx.h"
#include "proteinStore.h"
#include "Profiler.h"
#include "simplethreads.h"

using namespace freicore;

namespace freicore
{
	const string& proteinData::getSequence() const
	{
		if( !m_data.empty() )
			return m_data;

		if( !m_owner )
			throw runtime_error( "getting sequence for protein \"" + name + "\" without an owning protein store" );

		if( !m_owner->m_file || !m_owner->m_file->is_open() )
			throw runtime_error( "getting sequence for protein \"" + name + "\" from a protein store with no file" );

		m_owner->m_file->clear(); // clear eof

		m_owner->m_file->seekg( m_offset );

		string& data = const_cast<proteinData*>( this )->m_data;

		string dataLine;
		while( !m_owner->m_file->eof() && m_owner->m_file->peek() != '>' )
		{
			std::getline( *m_owner->m_file, dataLine );
			data += dataLine;
		}

		if( !data.empty() )
		{
			data.erase( remove( data.begin(), data.end(), '\r' ), data.end() );

			// Crop * from the end of the protein sequence, if necessary
			if( *data.rbegin() == '*' )
				data.erase( data.length() - 1 );
		}
		return m_data;
	}

	bool proteinData::testCleavage( const CleavageRule& rule, size_t offset ) const
	{
		size_t n_count = offset+1; // # of residues N terminal to the cleavage site (+1 for the N terminus)
		size_t c_count = m_data.length()-offset+1; // # of residues C terminal to the cleavage site (+1 for the C terminus)

		const CleavageHalfRule& n_rules = rule.first;
		bool n_passes = n_rules.hasWildcard;
		for( CleavageHalfRule::const_iterator itr = n_rules.begin(); !n_passes && itr != n_rules.end(); ++itr )
		{
			const string& n_rule = *itr;

			if( n_rule.length() > n_count ) // rule cannot pass if it needs more residues than there are to test
				continue;

			n_passes = true;
			for( size_t i=0; n_passes && i < n_rule.length(); ++i )
			{
				char protein_residue;
				if( offset-i == 0 ) // at the N terminus, check for PROTEIN_N_TERMINUS_SYMBOL
					protein_residue = PROTEIN_N_TERMINUS_SYMBOL;
				else
					protein_residue = m_data[offset-1-i];

				if( *(n_rule.rbegin()+i) != protein_residue )
					n_passes = false;
			}
		}

		if( !n_passes )
			return false;

		const CleavageHalfRule& c_rules = rule.second;
		bool c_passes = c_rules.hasWildcard;
		for( CleavageHalfRule::const_iterator itr = c_rules.begin(); !c_passes && itr != c_rules.end(); ++itr )
		{
			const string& c_rule = *itr;

			if( c_rule.length() > c_count ) // rule cannot pass if it needs more residues than there are to test
				continue;

			c_passes = true;
			for( size_t i=0; c_passes && i < c_rule.length(); ++i )
			{
				char protein_residue;
				if( offset+i == m_data.length() ) // at the C terminus, check for PROTEIN_C_TERMINUS_SYMBOL
					protein_residue = PROTEIN_C_TERMINUS_SYMBOL;
				else
					protein_residue = m_data[offset+i];

				if( c_rule[i] != protein_residue )
					c_passes = false;
			}
		}

		if( c_passes )
			return true;
		else
			return false;
	}

	CleavageRuleSet::const_iterator proteinData::testCleavage( const CleavageRuleSet& rules, size_t offset ) const
	{
		if( m_data.empty() )
			getSequence();

		if( offset > m_data.length() ) // the offset at m_data.length() will test for cleavage at the protein's C terminus
			throw out_of_range( "offset parameter to ProteinData::testCleavage" );

		CleavageRuleSet::const_iterator itr;
		for( itr = rules.begin(); itr != rules.end(); ++itr )
		{
			if( testCleavage( *itr, offset ) )
				return itr; // this cleavage rule passed the test
		}
		return itr; // no cleavage rule passed the test
	}

	void proteinStore::readFASTA( const string& filename, ProteinIndex startIndex, ProteinIndex endIndex, const string& delimiter, bool storeSequences )
	{
		if( m_file )
			delete m_file;
		m_file = new ifstream( filename.c_str(), ios::binary );
		if( !m_file || !m_file->is_open() )
			throw invalid_argument( "unable to open \"" + filename + "\"" );

		/*size_t filesize = (size_t) GetFileSize( filename );
		string fileStr;
		fileStr.resize( filesize );
		file.read( &fileStr[0], (streamsize) filesize );
		file.close();

		boost::char_separator<char> delim("\r\n");
		tokenizer parser( fileStr.begin(), fileStr.end(), delim );
		tokenizer::iterator itr = parser.begin();*/

		const ProteinOffset MAX_PROTEIN_LENGTH = -1;
		ProteinIndex pIndex = 0;

		char* buf = new char[ MAX_PROTEIN_LENGTH ];
		m_file->getline( buf, MAX_PROTEIN_LENGTH );

		map< char, int > invalidResidueCount;

		string realDelimiter = delimiter;
		if( delimiter.find( '\t' ) == string::npos )
			realDelimiter += '\t';
		if( delimiter.find( '\r' ) == string::npos )
			realDelimiter += '\r';
		if( delimiter.find( '\n' ) == string::npos )
			realDelimiter += '\n';

		while( !m_file->eof() )
		{
			if( buf[0] == '>' ) // signifies a new protein record in a FASTA file
			{
				++pIndex;

				if( pIndex >= startIndex && pIndex <= endIndex )
				{
					push_back( proteinData(this) );
					proteinData& p = back();
					string locusMetaData( buf );
					size_t locusEnd = locusMetaData.find_first_of( realDelimiter );

					p.name = locusMetaData.substr( 1, locusEnd-1 );
					if( p.name.find( decoyPrefix ) == 0 )
					{
						p.isDecoy = true;
						++numDecoys;
					} else
					{
						p.isDecoy = false;
						++numReals;
					}

					if( locusEnd < locusMetaData.length() )
					{
						p.desc = locusMetaData.substr( locusEnd+1 );
						p.desc.erase( remove( p.desc.begin(), p.desc.end(), '\r' ), p.desc.end() );
					}

					p.m_offset = m_file->tellg();

					nameToIndex[ p.name ] = size()-1;
					indexToName[ size()-1 ] = p.name;

					if( storeSequences )
					{
						p.getSequence();

						if( p.m_data.empty() )
						{
							cerr << "Warning: protein \'" << p.name << "\' contains no residues." << endl;
							pop_back();
						} else
						{
							// Verify each residue is in the residue map
							if( g_residueMap )
							{
								for( size_t i=0; i < p.m_data.length(); ++i )
									if( p.m_data[i] != 'X' && !g_residueMap->hasResidue( p.m_data[i] ) )
										//cerr << "Warning: protein \'" << p.name << "\' contains unknown residue \'" << p.data[i] << "\'" << endl;
										++invalidResidueCount[ p.m_data[i] ];
							}
						}
					}
				}
			}

			m_file->getline( buf, MAX_PROTEIN_LENGTH );
		}

		delete [] buf;

		m_file->clear();

		for( map< char, int >::iterator itr = invalidResidueCount.begin(); itr != invalidResidueCount.end(); ++itr )
			cerr << "Warning: unknown residue \'" << itr->first << "\' appears " << itr->second << " times in database." << endl;
	}

	void proteinStore::writeFASTA( const string& filename ) const
	{
		ofstream file( filename.c_str(), ios::binary );
		if( !file.is_open() )
			throw invalid_argument( "unable to open \"" + filename + "\"" );

		for( const_iterator itr = begin(); itr != end(); ++itr )
			file << ">" << itr->name << " " << itr->desc << "\n" << itr->getSequence() << "\n";

		file.close();
	}

	void proteinStore::add( const proteinData& p )
	{
		const_iterator itr = std::find( begin(), end(), p );
		if( itr == end() )
		{
			push_back( p );
			nameToIndex[ p.name ] = size()-1;
			indexToName[ size()-1 ] = p.name;
		}
	}

	void proteinStore::rebuildIndex()
	{
		nameToIndex.clear();
		indexToName.clear();
		for( size_t i=0; i < size(); ++i )
		{
			nameToIndex[at(i).name] = i;
			indexToName[i] = at(i).name;
		}
	}

	void proteinStore::random_shuffle()
	{
		std::random_shuffle( begin(), end() );
		rebuildIndex();
	}
}
