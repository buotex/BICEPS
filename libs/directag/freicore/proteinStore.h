#ifndef _PROTEINSTORE_H
#define _PROTEINSTORE_H

#include "shared_defs.h"
#include "shared_funcs.h"

using namespace freicore;

namespace freicore
{
	struct proteinStore;

	struct proteinData
	{
		proteinData( const string& nm = "", const string& ds = "", const string& dt = "" )
			: name(nm), desc(ds), m_data(dt), m_owner(NULL) {}

		proteinData( const proteinStore* owner, const string& nm = "", const string& ds = "", const string& dt = "" )
			: name(nm), desc(ds), m_data(dt), m_owner(owner) {}

		string			name;					// protein's name
		string			desc;					// description of protein
		const string&	getSequence() const;	// reads amino acid sequence if not already in memory
		bool			isDecoy;				// true if the sequence is a decoy

		// These two functions test for cleavability between offset-1 and offset
		bool testCleavage( const CleavageRule& rule, size_t offset ) const;
		CleavageRuleSet::const_iterator testCleavage( const CleavageRuleSet& rules, size_t offset ) const;

		bool operator== ( const proteinData& rhs ) { return name == rhs.name; }

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & name & desc & m_data & isDecoy;
		}

		friend struct proteinStore;

	private:
		size_t			m_offset;		// offset of protein data in file
		string			m_data;			// protein's full amino acid sequence
		const proteinStore*	m_owner;	// proteinStore instance that owns this proteinData instance
	};

	struct proteinStore : public vector< proteinData >
	{
		proteinStore( const string& decoyPrefix = "rev_" )
			:	decoyPrefix(decoyPrefix), numReals(0), numDecoys(0),
				m_file(NULL)
		{}
		~proteinStore() { if( m_file ) delete m_file; }

		void readFASTA( const string& filename, ProteinIndex startIndex = 0, ProteinIndex endIndex = -1, const string& delimiter = " ", bool storeSequences = true );
		void writeFASTA( const string& filename ) const;
		void add( const proteinData& p );

		const proteinData& operator[]( const string& name ) const	{ return vector< proteinData >::operator[]( nameToIndex.find(name)->second ); }
		proteinData& operator[]( const string& name )				{ return vector< proteinData >::operator[]( nameToIndex.find(name)->second ); }

		const proteinData& operator[]( const size_t index ) const	{ return vector< proteinData >::operator[]( index ); }
		proteinData& operator[]( const size_t index )				{ return vector< proteinData >::operator[]( index ); }

		void rebuildIndex();
		void random_shuffle();

		map< ProteinName, ProteinIndex > nameToIndex;
		map< ProteinIndex, ProteinName > indexToName;

		friend struct proteinData;

		string decoyPrefix;
		int numReals;
		int numDecoys;

	private:
		ifstream* m_file;
	};
}

#endif
