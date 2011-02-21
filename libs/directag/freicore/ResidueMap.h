#ifndef _RESIDUEMAP_H
#define _RESIDUEMAP_H

#include "stdafx.h"
#include "shared_types.h"

namespace freicore
{
	//typedef map< char, float > n2m_t;
	typedef map< float, AminoAcidResidue > m2n_t;
	//typedef vector< float > n2m_t;
	typedef CharIndexedVector<float> n2m_t;

	class ResidueMap
	{
	public:
		ResidueMap()
		{
			m_monoMassesToNamesMap[ 71.03711f ] = m_avgMassesToNamesMap[ 71.07880f ] = 'A';	// Alanine
			m_monoMassesToNamesMap[ 156.1011f ] = m_avgMassesToNamesMap[ 156.1875f ] = 'R';	// Arginine
			m_monoMassesToNamesMap[ 114.0429f ] = m_avgMassesToNamesMap[ 114.1038f ] = 'N';	// Asparagine
			m_monoMassesToNamesMap[ 115.0269f ] = m_avgMassesToNamesMap[ 115.0886f ] = 'D';	// Aspartic acid
			m_monoMassesToNamesMap[ 103.0092f ] = m_avgMassesToNamesMap[ 103.1388f ] = 'C';	// Cysteine
			m_monoMassesToNamesMap[ 129.0426f ] = m_avgMassesToNamesMap[ 129.1155f ] = 'E';	// Glutamic acid
			m_monoMassesToNamesMap[ 128.0586f ] = m_avgMassesToNamesMap[ 128.1307f ] = 'Q';	// Glutamine
			m_monoMassesToNamesMap[ 57.02146f ] = m_avgMassesToNamesMap[ 57.05190f ] = 'G';	// Glycine
			m_monoMassesToNamesMap[ 137.0589f ] = m_avgMassesToNamesMap[ 137.1411f ] = 'H';	// Histidine
			m_monoMassesToNamesMap[ 113.0841f ] = m_avgMassesToNamesMap[ 113.1594f ] = 'I';	// Isoleucine
			m_monoMassesToNamesMap[ 128.0841f ] = m_avgMassesToNamesMap[ 128.1741f ] = 'K';	// Lysine
			m_monoMassesToNamesMap[ 131.0405f ] = m_avgMassesToNamesMap[ 131.1926f ] = 'M';	// Methionine
			m_monoMassesToNamesMap[ 147.0684f ] = m_avgMassesToNamesMap[ 147.1766f ] = 'F';	// Phenylalanine
			m_monoMassesToNamesMap[ 97.05276f ] = m_avgMassesToNamesMap[ 97.11670f ] = 'P';	// Proline
			m_monoMassesToNamesMap[ 87.03203f ] = m_avgMassesToNamesMap[ 87.07820f ] = 'S';	// Serine
			m_monoMassesToNamesMap[ 101.0477f ] = m_avgMassesToNamesMap[ 101.1051f ] = 'T';	// Threonine
			m_monoMassesToNamesMap[ 186.0793f ] = m_avgMassesToNamesMap[ 186.2132f ] = 'W';	// Tryptophan
			m_monoMassesToNamesMap[ 163.0633f ] = m_avgMassesToNamesMap[ 163.1760f ] = 'Y';	// Tyrosine
			m_monoMassesToNamesMap[ 99.06841f ] = m_avgMassesToNamesMap[ 99.13260f ] = 'V';	// Valine

			for( m2n_t::iterator itr = m_monoMassesToNamesMap.begin(); itr != m_monoMassesToNamesMap.end(); ++itr )
			{
				m_defaultNamesToMonoMassesMap[ itr->second ] = itr->first;
				m_defaultResidues.insert( itr->second );
			}

			for( m2n_t::iterator itr = m_avgMassesToNamesMap.begin(); itr != m_avgMassesToNamesMap.end(); ++itr )
				m_defaultNamesToAvgMassesMap[ itr->second ] = itr->first;

			m_monoMassesToNamesMap[ HYDROGEN_MONO + OXYGEN_MONO ] = PEPTIDE_C_TERMINUS_SYMBOL;
			m_avgMassesToNamesMap[ HYDROGEN_AVG + OXYGEN_AVG ] = PEPTIDE_C_TERMINUS_SYMBOL;
			m_monoMassesToNamesMap[ HYDROGEN_MONO ] = PEPTIDE_N_TERMINUS_SYMBOL;
			m_avgMassesToNamesMap[ HYDROGEN_AVG ] = PEPTIDE_N_TERMINUS_SYMBOL;

			finalize();
		}

		bool initialized() { return !cfgStr.empty(); }

		int initializeFromFile( const string& residuesCfgFilename = "residue_masses.cfg" )
		{
			ifstream residuesCfgFile( residuesCfgFilename.c_str() );
			if( residuesCfgFile.is_open() )
			{
				// clear old residues
				m_residues.clear();
				m_monoMassesToNamesMap.clear();
				m_avgMassesToNamesMap.clear();
				m_namesToMonoMassesMap.clear();
				m_namesToAvgMassesMap.clear();

				int cfgSize = (int) GetFileSize( residuesCfgFilename );
				cfgStr.resize( cfgSize );
				residuesCfgFile.read( &cfgStr[0], cfgSize );
				int rv = initializeFromBuffer( cfgStr );
				residuesCfgFile.close();
				return rv;
			} else
				return 1;
		}

		int initializeFromBuffer( const string& cfgStr )
		{
			stringstream cfgStream( cfgStr );

			char r;
			float mono, avg;
			while( cfgStream >> r >> mono >> avg )
			{
				if( !m_defaultResidues.count(r) )
					cerr << "Warning: residue map has been initialized with non-standard residue '" << r << "'!" << endl;
				else
				{
					if( m_defaultNamesToMonoMassesMap[r] != mono )
						cerr << "Warning: residue map has initialized standard residue '" << r << "' with non-standard monoisotopic mass!" << endl;
					if( m_defaultNamesToAvgMassesMap[r] != avg )
						cerr << "Warning: residue map has initialized standard residue '" << r << "' with non-standard average mass!" << endl;
				}

				m_monoMassesToNamesMap[ mono ] = r;
				m_avgMassesToNamesMap[ avg ] = r;
			}
			
			if( m_monoMassesToNamesMap.empty() || m_avgMassesToNamesMap.empty() )
				return 1;

			finalize();

			return 0;
		}

		void finalize()
		{
			for( m2n_t::iterator itr = m_monoMassesToNamesMap.begin(); itr != m_monoMassesToNamesMap.end(); ++itr )
			{
				m_namesToMonoMassesMap[ itr->second ] = itr->first;
				m_residues.insert( itr->second );
			}

			for( m2n_t::iterator itr = m_avgMassesToNamesMap.begin(); itr != m_avgMassesToNamesMap.end(); ++itr )
				m_namesToAvgMassesMap[ itr->second ] = itr->first;

			if( m_namesToMonoMassesMap['I'] > 0 )
			{
				m_namesToMonoMassesMap['L'] = m_namesToMonoMassesMap['I'];
				m_namesToAvgMassesMap['L'] = m_namesToAvgMassesMap['I'];
				m_residues.insert('L');
			} else if( m_namesToMonoMassesMap['L'] > 0 )
			{
				m_namesToMonoMassesMap['I'] = m_namesToMonoMassesMap['L'];
				m_namesToAvgMassesMap['I'] = m_namesToAvgMassesMap['L'];
				m_residues.insert('L');
			}

			m_residues.erase(PEPTIDE_N_TERMINUS_SYMBOL);
			m_residues.erase(PEPTIDE_C_TERMINUS_SYMBOL);
		}

		void addDynamicMod( const DynamicMod& mod )
		{
			if( !dynamicMods.count( mod ) )
			{
				dynamicMods.insert( mod );
				m_monoMassesToNamesMap[ m_namesToMonoMassesMap[ mod.unmodChar ] + mod.modMass ] = mod.uniqueModChar;
				m_avgMassesToNamesMap[ m_namesToAvgMassesMap[ mod.unmodChar ] + mod.modMass ] = mod.uniqueModChar;
				m_namesToMonoMassesMap[ mod.uniqueModChar ] = m_namesToMonoMassesMap[ mod.unmodChar ] + mod.modMass;
				m_namesToAvgMassesMap[ mod.uniqueModChar ] = m_namesToAvgMassesMap[ mod.unmodChar ] + mod.modMass;
			}
		}

		void removeDynamicMod( const DynamicMod& mod )
		{
			if( dynamicMods.count( mod ) )
			{
				dynamicMods.erase( mod );
				forceAddDynamicMod( mod );
			}
		}

		void addStaticMod( const StaticMod& mod )
		{
			if( !staticMods.count( mod ) )
			{
				staticMods.insert( mod );
				forceAddStaticMod( mod );
			}
		}

		void clearDynamicMods()
		{
			for( DynamicModSet::iterator itr = dynamicMods.begin(); itr != dynamicMods.end(); ++itr )
			{
				m_monoMassesToNamesMap.erase( m_namesToMonoMassesMap[ itr->unmodChar ] + itr->modMass );
				m_avgMassesToNamesMap.erase( m_namesToAvgMassesMap[ itr->unmodChar ] + itr->modMass );
				m_namesToMonoMassesMap.erase( itr->uniqueModChar );
				m_namesToAvgMassesMap.erase( itr->uniqueModChar );
			}
			dynamicMods.clear();
		}

		void clearStaticMods()
		{
			for( StaticModSet::iterator itr = staticMods.begin(); itr != staticMods.end(); ++itr )
			{
				float monoMass = m_namesToMonoMassesMap[ itr->name ];
				float avgMass = m_namesToAvgMassesMap[ itr->name ];
				m_monoMassesToNamesMap[ monoMass - itr->mass ] = itr->name;
				m_avgMassesToNamesMap[ avgMass - itr->mass ] = itr->name;
				m_monoMassesToNamesMap.erase( monoMass );
				m_avgMassesToNamesMap.erase( avgMass );
				m_namesToMonoMassesMap[ itr->name ] = monoMass - itr->mass;
				m_namesToAvgMassesMap[ itr->name ] = avgMass - itr->mass;
			}
			staticMods.clear();
		}

		void setDynamicMods( const string& cfgStr )
		{
			clearDynamicMods();
			if( cfgStr.empty() ) return;
			DynamicModSet tmp( cfgStr );
			for( DynamicModSet::iterator itr = tmp.begin(); itr != tmp.end(); ++itr )
				addDynamicMod( *itr );
		}

		void setStaticMods( const string& cfgStr )
		{
			clearStaticMods();
			if( cfgStr.empty() ) return;
			StaticModSet tmp( cfgStr );
			for( StaticModSet::iterator itr = tmp.begin(); itr != tmp.end(); ++itr )
				addStaticMod( *itr );
		}

		const set<char>& getResidues() const	{ return m_residues; }
		bool	hasResidue( char r ) const		{ return m_namesToMonoMassesMap[r] > 0; }
		float	largestDynamicModMass() const	{ return ( dynamicMods.empty() ? 0 : dynamicMods.rbegin()->modMass ); }
		float	smallestDynamicModMass() const	{ return ( dynamicMods.empty() ? 0 : dynamicMods.begin()->modMass ); }

		inline float GetMassOfResidues( const string::const_iterator& seqBegin, const string::const_iterator& seqEnd, bool useAvgMass = false ) const
		{
			float mass = 0.0f;
			if( useAvgMass )
			{
				for( string::const_iterator itr = seqBegin; itr != seqEnd; ++itr )
					mass += getAvgMassByName( *itr );
			} else
			{
				for( string::const_iterator itr = seqBegin; itr != seqEnd; ++itr )
					mass += getMonoMassByName( *itr );
			}
			return mass;
		}

		inline float GetMassOfResidues( const string& residues, bool useAvgMass = false ) const
		{
			return GetMassOfResidues( residues.begin(), residues.end(), useAvgMass );
		}

		void dump() const
		{
			//cout << "Residue map size: " << m_namesToMonoMassesMap.size() << endl;
			//for( n2m_t::iterator itr = m_namesToMonoMassesMap.begin(); itr != m_namesToMonoMassesMap.end(); ++itr )
			//	cout << itr->first << ": " << itr->second << " " << m_namesToAvgMassesMap[ itr->first ] << endl;
			for( size_t i=0; i < 128; ++i )
				if( m_namesToMonoMassesMap[i] > 0 )
					cout << (char) i << ": " << m_namesToMonoMassesMap[i] << " " << m_namesToAvgMassesMap[i] << endl;
			cout << dynamicMods.userToUniqueMap << endl << dynamicMods.uniqueToUserMap << endl;
		}


		size_t size() const								{ return m_namesToMonoMassesMap.size(); }

		n2m_t::const_iterator beginMonoNames() const	{ return m_namesToMonoMassesMap.begin(); }
		n2m_t::const_iterator endMonoNames() const		{ return m_namesToMonoMassesMap.end(); }
		n2m_t::const_iterator beginAvgNames() const		{ return m_namesToAvgMassesMap.begin(); }
		n2m_t::const_iterator endAvgNames() const		{ return m_namesToAvgMassesMap.end(); }

		m2n_t::const_iterator beginMonoMasses() const	{ return m_monoMassesToNamesMap.begin(); }
		m2n_t::const_iterator endMonoMasses() const		{ return m_monoMassesToNamesMap.end(); }
		m2n_t::const_iterator beginAvgMasses() const	{ return m_avgMassesToNamesMap.begin(); }
		m2n_t::const_iterator endAvgMasses() const		{ return m_avgMassesToNamesMap.end(); }

		char getNameByMonoMass( float key ) const		{ return m_monoMassesToNamesMap.find( key )->second; }
		char getNameByAvgMass( float key ) const		{ return m_avgMassesToNamesMap.find( key )->second; }

		char getNameByMonoMass( float key, float epsilon ) const
		{
			if( epsilon > 0.0f )
			{
				m2n_t::const_iterator min, max, cur, best;
				min = m_monoMassesToNamesMap.lower_bound( key - epsilon );
				max = m_monoMassesToNamesMap.lower_bound( key + epsilon );
				if( min == max )
					throw runtime_error( string( "No residue matching mass " ) + lexical_cast<string>( key ) );

				float minDelta = fabs( min->first - key );
				for( best = cur = min; cur != max; ++cur )
				{
					float curDelta = fabs( cur->first - key );
					if( curDelta < minDelta )
					{
						minDelta = curDelta;
						best = cur;
					}
				}
				return best->second;
			} else
				return getNameByMonoMass(key);
		}

		char getNameByAvgMass( float key, float epsilon ) const
		{
			if( epsilon > 0.0f )
			{
				m2n_t::const_iterator min, max, cur, best;
				min = m_avgMassesToNamesMap.lower_bound( key - epsilon );
				max = m_avgMassesToNamesMap.lower_bound( key + epsilon );
				if( min == max )
					throw runtime_error( string( "No residue matching mass " ) + lexical_cast<string>( key ) );

				float minDelta = fabs( min->first - key );
				for( best = cur = min; cur != max; ++cur )
				{
					float curDelta = fabs( cur->first - key );
					if( curDelta < minDelta )
					{
						minDelta = curDelta;
						best = cur;
					}
				}
				return best->second;
			} else
				return getNameByAvgMass(key);
		}

		float getMonoMassByName( char key ) const		{ return m_namesToMonoMassesMap[ key ]; }
		float getAvgMassByName( char key ) const		{ return m_namesToAvgMassesMap[ key ]; }

		string			cfgStr;
		DynamicModSet	dynamicMods;
		StaticModSet	staticMods;

	private:

		void forceAddDynamicMod( const DynamicMod& mod )
		{
			m_monoMassesToNamesMap[ m_namesToMonoMassesMap[ mod.unmodChar ] + mod.modMass ] = mod.uniqueModChar;
			m_avgMassesToNamesMap[ m_namesToAvgMassesMap[ mod.unmodChar ] + mod.modMass ] = mod.uniqueModChar;
			m_namesToMonoMassesMap[ mod.uniqueModChar ] = m_namesToMonoMassesMap[ mod.unmodChar ] + mod.modMass;
			m_namesToAvgMassesMap[ mod.uniqueModChar ] = m_namesToAvgMassesMap[ mod.unmodChar ] + mod.modMass;
		}

		void forceAddStaticMod( const StaticMod& mod )
		{
			float monoMass = m_namesToMonoMassesMap[ mod.name ];
			float avgMass = m_namesToAvgMassesMap[ mod.name ];
			m_monoMassesToNamesMap[ monoMass + mod.mass ] = mod.name;
			m_avgMassesToNamesMap[ avgMass + mod.mass ] = mod.name;
			m_monoMassesToNamesMap.erase( monoMass );
			m_avgMassesToNamesMap.erase( avgMass );
			m_namesToMonoMassesMap[ mod.name ] = monoMass + mod.mass;
			m_namesToAvgMassesMap[ mod.name ] = avgMass + mod.mass;
		}

		n2m_t m_defaultNamesToMonoMassesMap;
		n2m_t m_defaultNamesToAvgMassesMap;
		set<char> m_defaultResidues;

		n2m_t m_namesToMonoMassesMap;
		n2m_t m_namesToAvgMassesMap;
		m2n_t m_monoMassesToNamesMap;
		m2n_t m_avgMassesToNamesMap;
		set<char> m_residues;
	};
}

#endif

