#include "stdafx.h"
#include "BaseRunTimeConfig.h"

namespace freicore
{
	BaseRunTimeConfig::BaseRunTimeConfig()
	{
		// Initialize variables to their default values
		// For each variable name: "VariableName = VariableDefaultValue;"
		//BOOST_PP_SEQ_FOR_EACH( RTCONFIG_INIT_DEFAULT_VAR, ~, BASE_RUNTIME_CONFIG )
	}

	void BaseRunTimeConfig::initializeFromBuffer( string& cfgStr, const string& delim )
	{
		size_t boolIdx;
		while( ( boolIdx = cfgStr.find( "true" ) ) != string::npos )
			cfgStr = cfgStr.replace( boolIdx, 4, "1" );
		while( ( boolIdx = cfgStr.find( "false" ) ) != string::npos )
			cfgStr = cfgStr.replace( boolIdx, 5, "0" );
		string strVal;
		// Find the variable name in the buffer of key-value pairs and read its corresponding value
		//BOOST_PP_SEQ_FOR_EACH( RTCONFIG_PARSE_BUFFER, ~, BASE_RUNTIME_CONFIG )
		finalize();
	}

	RunTimeVariableMap BaseRunTimeConfig::getVariables( bool hideDefaultValues )
	{
		// Update the variable map
		// For each variable name: "m_variables[ "VariableName" ] = VariableName;"
		//BOOST_PP_SEQ_FOR_EACH( RTCONFIG_FILL_MAP, m_variables, BASE_RUNTIME_CONFIG )
		return m_variables;
	}

	void BaseRunTimeConfig::setVariables( RunTimeVariableMap& vars )
	{
		for( RunTimeVariableMap::iterator itr = vars.begin(); itr != vars.end(); ++itr )
		{
			string value = UnquoteString( itr->second );
			if( value == "true" )
				itr->second = "\"1\"";
			else if( value == "false" )
				itr->second = "\"0\"";
		}
		// Update the variable map
		// For each variable name: "m_variables[ "VariableName" ] = VariableName;"
		//BOOST_PP_SEQ_FOR_EACH( RTCONFIG_READ_MAP, vars, BASE_RUNTIME_CONFIG )
		finalize();
	}

	void BaseRunTimeConfig::dump()
	{
		getVariables();
		string::size_type longestName = 0;
		for( RunTimeVariableMap::iterator itr = m_variables.begin(); itr != m_variables.end(); ++itr )
			if( itr->first.length() > longestName )
				longestName = itr->first.length();

		for( RunTimeVariableMap::iterator itr = m_variables.begin(); itr != m_variables.end(); ++itr )
		{
			cout.width( (streamsize) longestName + 2 );
			stringstream s;
			s << right << itr->first << ": ";
			cout << s.str() << boolalpha << "\"" << itr->second << "\"" << endl;
		}
	}

	void BaseRunTimeConfig::finalize()
	{
	}

	int BaseRunTimeConfig::initializeFromFile( const string& rtConfigFilename, const string& delimiters )
	{
		// Abort
		if( rtConfigFilename.empty() )
		{
			finalize();
			return 1;
		}

		// Read settings from file; abort if file does not exist
		else
		{
			ifstream rtConfigFile( rtConfigFilename.c_str(), ios::binary );
			if( rtConfigFile.is_open() )
			{
				//cout << GetHostname() << " is reading its configuration file \"" << rtConfigFilename << "\"" << endl;
				int cfgSize = (int) GetFileSize( rtConfigFilename );
				cfgStr.resize( cfgSize );
				rtConfigFile.read( &cfgStr[0], cfgSize );
				initializeFromBuffer( cfgStr, delimiters );
				rtConfigFile.close();
			} else
			{
				finalize();
				return 1;
			}
		}

		return 0;
	}
}
