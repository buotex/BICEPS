#ifndef _SHARED_TYPES_H
#define _SHARED_TYPES_H

#include "stdafx.h"

namespace freicore
{
	class ResidueMap;
	class lnFactorialTable;
	//class mzDataReader;
	struct BaseRunTimeConfig;
	struct SqtList;

	typedef enum { SYS_BIG_ENDIAN, SYS_LITTLE_ENDIAN, SYS_UNKNOWN_ENDIAN } endianType_t;

	typedef char											AminoAcidResidue;
	const size_t											AminoAcidResidueSize = sizeof(AminoAcidResidue)*CHAR_BIT;

	typedef unsigned long									ProteinIndex;
	typedef unsigned short									ProteinOffset;
	typedef string											ProteinName;

	struct ProteinLocusByIndex
	{
		ProteinLocusByIndex( ProteinIndex pIndex = 0, ProteinOffset pOffset = 0 ) : index( pIndex ), offset( pOffset ) {}
		ProteinIndex	index;
		ProteinOffset	offset;

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & index & offset;
		}

		bool operator< ( const ProteinLocusByIndex& rhs ) const
		{
			if( index == rhs.index )
				return offset < rhs.offset;
			else
				return index < rhs.index;
		}
	};

	struct ProteinLocusByName
	{
		ProteinLocusByName( ProteinName pName = "", ProteinOffset pOffset = 0 ) : name( pName ), offset( pOffset ) {}
		ProteinName		name;
		ProteinOffset	offset;
		string			desc;

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & name & offset;
		}

		bool operator< ( const ProteinLocusByName& rhs ) const
		{
			if( name == rhs.name )
				return offset < rhs.offset;
			else
				return name < rhs.name;
		}
	};

	typedef set< ProteinLocusByIndex >						ProteinLociByIndex;
	typedef set< ProteinLocusByName >						ProteinLociByName;
	typedef map< ProteinIndex, vector< ProteinOffset > >	ProteinLociMap;

	struct CleavageHalfRule : public set< string >
	{
		CleavageHalfRule() : hasWildcard( false ), longestCleavageCandidate( 0 ) {}
		void clear()
		{
			set< string >::clear();
			hasWildcard = false;
			longestCleavageCandidate = 0;
		}

		SetInsertPair( string ) insert( const string& seq )
		{
			if( seq == "." )
			{
				hasWildcard = true;
				return SetInsertPair( string )( end(), false );
			} else
			{
				longestCleavageCandidate = max( seq.length(), longestCleavageCandidate );
				return set< string >::insert( seq );
			}
		}

		bool hasWildcard;
		size_t longestCleavageCandidate;
	};

	struct CleavageRule : public pair< CleavageHalfRule, CleavageHalfRule >
	{
		CleavageRule() {}
	};

	struct CleavageRuleSet : public vector< CleavageRule >
	{
		CleavageRuleSet( const string& cfgStr = "" ) : longestPreCleavageCandidate( 0 ), longestPostCleavageCandidate( 0 )
		{
			initialize( cfgStr );
		}

		size_t longestPreCleavageCandidate;
		size_t longestPostCleavageCandidate;

		void initialize( const string& cfgStr );
	};

	struct ResidueFilter
	{
		ResidueFilter() {};

		bool testResidue( const AminoAcidResidue& r ) const
		{
			return m_filter[r];
		}

		operator string () const;
		CharIndexedVector<bool> m_filter;
	};

	struct DynamicMod
	{
		DynamicMod( char unmodChar = 0, char userModChar = 0, float modMass = 0 )
			: unmodChar( unmodChar ), userModChar( userModChar ), modMass( modMass ) {}

		vector< ResidueFilter > NTerminalFilters;
		vector< ResidueFilter > CTerminalFilters;

		char unmodChar;
		char userModChar;
		char uniqueModChar;
		float modMass;

		bool operator< ( const DynamicMod& rhs ) const
		{
			return uniqueModChar < rhs.uniqueModChar;
		}
	};

	struct StaticMod
	{
		StaticMod() {}
		StaticMod( char r, float m ) : name(r), mass(m) {}
		char name;
		float mass;

		bool operator< ( const StaticMod& rhs ) const
		{
			if( name == rhs.name )
				return mass < rhs.mass;
			return name < rhs.name;
		}
	};

	struct DynamicModSet : public set< DynamicMod >
	{
		typedef map< char, vector< DynamicMod > >	UserToUniqueMap;
		typedef map< char, DynamicMod >				UniqueToUserMap;

		DynamicModSet( const string& cfgStr = "" )
		{
			initialize( cfgStr );
		}

		UserToUniqueMap	userToUniqueMap;
		UniqueToUserMap	uniqueToUserMap;

		void clear();
		void erase( const DynamicMod& mod );
		SetInsertPair(DynamicMod) insert( const DynamicMod& mod );

		void initialize( const string& cfgStr );
		void parseMotif( const string& motif, char modChar, float modMass );
		operator string () const;
	};

	struct StaticModSet : public set< StaticMod >
	{
		StaticModSet( const string& cfgStr = "" )
		{
			boost::char_separator<char> delim(" ");
			tokenizer parser( cfgStr.begin(), cfgStr.begin() + cfgStr.length(), delim );
			tokenizer::iterator itr = parser.begin();
			while( itr != parser.end() )
			{
				char r = (*itr)[0];
				float m = lexical_cast<float>( *(++itr) );
				insert( StaticMod( r, m ) );
				++itr;
			}
		}

		operator string ()
		{
			stringstream modStr;
			for( iterator itr = begin(); itr != end(); ++itr )
				modStr << ( itr == begin() ? "" : " " ) << itr->name << " " << itr->mass;
			return modStr.str();
		}
	};

	struct SpectrumId
	{
		SpectrumId( const string& s = "" )
			: id(s)
		{
			updateFromString();
		}

		SpectrumId( const string& source, int index, int charge = 0 )
			: source(source), index(index), charge(charge)
		{
			updateFromVars();
		}

		SpectrumId( int index, int charge )
			: source(""), index(index), charge(charge)
		{
			updateFromVars();
		}

		void set( const string& source, int index, int charge = 0 )
		{
			this->source = source; this->index = index; this->charge = charge;
			updateFromVars();
		}

		void setId( const SpectrumId& id )	    { *this = id; }
		void setId( const string& id )		    { this->id = id; updateFromString(); }
		void setSource( const string& source )	{ this->source = source; updateFromVars(); }
		void setIndex( int index )              { this->index = index; updateFromVars(); }
		void setCharge( int charge )            { this->charge = charge; updateFromVars(); }

		bool operator< ( const SpectrumId& rhs ) const
		{
			if( source == rhs.source )
				if( index == rhs.index )
					return charge < rhs.charge;
				else
					return index < rhs.index;
			else
				return source < rhs.source;
		}

		bool operator== ( const SpectrumId& rhs ) const
		{
			return source == rhs.source && charge == rhs.charge && index == rhs.index;
		}

		operator string () { return id; }

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & id & source & index & charge;
		}

		string	id;
		string	source;
		int		index;
		int		charge;

	private:
		void updateFromVars()
		{
			stringstream idStream;
			if( !source.empty() )
				idStream << source;

			if( index > 0 )
			{
				if( !source.empty() )
					idStream << '.';
				idStream << index;
			}

			if( charge > 0 )
				idStream << '.' << charge;

			id = idStream.str();
		}

		void updateFromString()
		{
			size_t firstDot = id.find_first_of( '.' );
			size_t lastDot = id.find_last_of( '.' );

			source = id.substr( 0, firstDot );

			if( firstDot < lastDot )
			{
				index = lexical_cast<int>( id.substr( firstDot+1, lastDot-firstDot-1 ) );
				charge = lexical_cast<int>( id.substr( lastDot+1 ) );
			} else if( firstDot != string::npos )
				index = lexical_cast<int>( id.substr( firstDot+1 ) );
		}
	};

	typedef set< string >			fileList_t;
	typedef vector< string >		argList_t;
}
#endif
