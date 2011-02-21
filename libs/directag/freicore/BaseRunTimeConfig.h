#ifndef _BASERUNTIMECONFIG_H
#define _BASERUNTIMECONFIG_H

#include "stdafx.h"

#define RTCONFIG_VARIABLE_EX(varType, varName, varDefaultValue, varInit)	((4, (varType, varName, varDefaultValue, varInit)))
#define RTCONFIG_VARIABLE(varType, varName, varDefaultValue)				RTCONFIG_VARIABLE_EX(varType, varName, varDefaultValue, 1)
#define RTCONFIG_VAR_NAME(var)												BOOST_PP_ARRAY_ELEM(1, var)
#define RTCONFIG_VAR_TYPE(var)												BOOST_PP_ARRAY_ELEM(0, var)
#define RTCONFIG_VAR_DEFAULTVALUE(var)										BOOST_PP_ARRAY_ELEM(2, var)
#define RTCONFIG_VAR_INIT(var)												BOOST_PP_ARRAY_ELEM(3, var)

#define RTCONFIG_VAR_NAME_CAT(var, str)										BOOST_PP_CAT( RTCONFIG_VAR_NAME(var), str )
#define RTCONFIG_VAR_NAME_STR(var)											BOOST_PP_STRINGIZE( RTCONFIG_VAR_NAME(var) )

#define RTCONFIG_DECLARE_VAR(r, n_a, var)						RTCONFIG_VAR_TYPE(var) RTCONFIG_VAR_NAME(var);
#define RTCONFIG_INIT_DEFAULT_VAR_(r, n_a, var)					RTCONFIG_VAR_NAME(var) = RTCONFIG_VAR_DEFAULTVALUE(var);

#define RTCONFIG_INIT_DEFAULT_VAR(r, n_a, var) \
	BOOST_PP_IF( RTCONFIG_VAR_INIT(var), RTCONFIG_INIT_DEFAULT_VAR_(r, n_a, var), 0; )

#define RTCONFIG_FILL_MAP(r, varMap, var) \
		string RTCONFIG_VAR_NAME_CAT( var, Val ); \
		try \
		{ \
			if( !hideDefaultValues || !( RTCONFIG_VAR_NAME(var) == RTCONFIG_VAR_DEFAULTVALUE(var) ) ) \
			{ \
				RTCONFIG_VAR_NAME_CAT( var, Val ) = lexical_cast<string>( RTCONFIG_VAR_NAME(var) ); \
				varMap[ RTCONFIG_VAR_NAME_STR(var) ] = RTCONFIG_VAR_NAME_CAT( var, Val ); \
			} \
		} catch( exception& e ) \
		{ \
			cerr << "FILL_MAP: casting " << RTCONFIG_VAR_NAME_STR(var) << " with value " << RTCONFIG_VAR_NAME(var) << ": " << e.what() << endl; \
		}

#define RTCONFIG_READ_MAP(r, varMap, var) \
		RunTimeVariableMap::const_iterator RTCONFIG_VAR_NAME_CAT( var, Itr ) = varMap.find( RTCONFIG_VAR_NAME_STR(var) ); \
		if( RTCONFIG_VAR_NAME_CAT( var, Itr ) != varMap.end() ) \
		{ \
			string RTCONFIG_VAR_NAME_CAT( var, Str ) = UnquoteString( RTCONFIG_VAR_NAME_CAT( var, Itr )->second ); \
			if( !RTCONFIG_VAR_NAME_CAT( var, Str ).empty() ) \
			{ \
				try \
				{ \
					RTCONFIG_VAR_NAME(var) = lexical_cast< RTCONFIG_VAR_TYPE(var) >( RTCONFIG_VAR_NAME_CAT( var, Str ) ) ; \
				} catch( exception& e ) \
				{ \
					cerr << "READ_MAP: casting " << RTCONFIG_VAR_NAME_STR(var) << " with value " << RTCONFIG_VAR_NAME_CAT( var, Str ) << ": " << e.what() << endl; \
				} \
			} \
		}

#define RTCONFIG_PARSE_BUFFER_(r, n_a, var) \
	if( !( strVal = GetStringByKey( cfgStr, RTCONFIG_VAR_NAME_STR(var), delim ) ).empty() ) \
	{ \
		strVal = UnquoteString( TrimWhitespace( strVal ) ); \
		try \
		{ \
			RTCONFIG_VAR_NAME(var) = lexical_cast< RTCONFIG_VAR_TYPE(var) >( strVal ); \
		} catch( exception& e ) \
		{ \
			cerr << "PARSE_BUFFER_: parsing " << RTCONFIG_VAR_NAME_STR(var) << " with value " << strVal << ": " << e.what() << endl; \
		} \
	}

#define RTCONFIG_PARSE_BUFFER(r, n_a, var) \
	BOOST_PP_IF( RTCONFIG_VAR_INIT(var), RTCONFIG_PARSE_BUFFER_(r, n_a, var), 0; )

#define RTCONFIG_PRINT_VAR(r, n_a, var) \
	stringstream RTCONFIG_VAR_NAME_CAT( var, Stream ); \
	RTCONFIG_VAR_NAME_CAT( var, Stream ) << right << RTCONFIG_VAR_NAME_STR(var) << ": "; \
	cout << RTCONFIG_VAR_NAME_CAT( var, Stream ).str() << RTCONFIG_VAR_NAME(var) << endl;

#define RTCONFIG_DEFINE_MEMBERS(configName, configVariables, configDefaultBufferDelimiters, configDefaultFilename, configDefaultFileDelimiters) \
	BOOST_PP_SEQ_FOR_EACH( RTCONFIG_DECLARE_VAR, ~, configVariables ) \
	configName() : BaseRunTimeConfig() \
	{ \
		BOOST_PP_SEQ_FOR_EACH( RTCONFIG_INIT_DEFAULT_VAR, ~, configVariables ) \
	} \
	void initializeFromBuffer( string& cfgStr, const string& delim = configDefaultBufferDelimiters ) \
	{ \
		BaseRunTimeConfig::initializeFromBuffer( cfgStr, delim ); \
		string strVal; \
		BOOST_PP_SEQ_FOR_EACH( RTCONFIG_PARSE_BUFFER, ~, configVariables ) \
		finalize(); \
	} \
	RunTimeVariableMap getVariables( bool hideDefaultValues = false ) \
	{ \
		BaseRunTimeConfig::getVariables( hideDefaultValues ); \
		BOOST_PP_SEQ_FOR_EACH( RTCONFIG_FILL_MAP, m_variables, configVariables ) \
		return m_variables; \
	} \
	void setVariables( RunTimeVariableMap& vars ) \
	{ \
		BaseRunTimeConfig::setVariables( vars ); \
		BOOST_PP_SEQ_FOR_EACH( RTCONFIG_READ_MAP, vars, configVariables ) \
		finalize(); \
	} \
	int initializeFromFile( const string& rtConfigFilename = configDefaultFilename, const string& delim = configDefaultFileDelimiters ) \
	{ \
		return BaseRunTimeConfig::initializeFromFile( rtConfigFilename, delim ); \
	}

/*#define BASE_RUNTIME_CONFIG	\
		RTCONFIG_VARIABLE( int,				NumChargeStates,			3		) \
		RTCONFIG_VARIABLE( float,			StatusUpdateFrequency,		5		) \
		RTCONFIG_VARIABLE( bool,			UseMultipleProcessors,		true	) \
		RTCONFIG_VARIABLE( int,				ThreadCountMultiplier,		10		) \
		RTCONFIG_VARIABLE( float,			PrecursorMzTolerance,		2.5f	) \
		RTCONFIG_VARIABLE( float,			FragmentMzTolerance,		0.5f	) \
		RTCONFIG_VARIABLE( float,			ComplementMzTolerance,		0.5f	)*/

namespace freicore
{
	struct RunTimeVariableMap : public map< string, string >
	{
		RunTimeVariableMap(	const string& initialVarList = "" )
		{
			static const boost::char_separator<char> delim(" ");
			tokenizer parser( initialVarList.begin(), initialVarList.begin() + initialVarList.length(), delim );

			for( tokenizer::iterator itr = parser.begin(); itr != parser.end(); ++itr )
			{
				operator[]( *itr ) = "";
			}
		}
	};

	struct BaseRunTimeConfig
	{
	protected:
		RunTimeVariableMap m_variables;

	public:
		string cfgStr;
		// Declare variables
		// For each variable name: "VariableType VariableName;"
		//BOOST_PP_SEQ_FOR_EACH( RTCONFIG_DECLARE_VAR, ~, BASE_RUNTIME_CONFIG )

									BaseRunTimeConfig();
		virtual						~BaseRunTimeConfig() {}

		virtual void				initializeFromBuffer( string& cfgStr, const string& delim = "\r\n\t" );
		virtual	RunTimeVariableMap	getVariables( bool hideDefaultValues = false );
		virtual void				setVariables( RunTimeVariableMap& vars );
		virtual void				dump();
		virtual void				finalize();

		bool						initialized() { return !cfgStr.empty(); }
		int							initializeFromFile( const string& rtConfigFilename, const string& delimiters = "\r\n#" );
	};
}

#endif
