#include "stdafx.h"
#include "freicore.h"

namespace freicore
{
	void GetHostname( char* buf, int len )
	{
	#ifdef WIN32
		WSADATA wsaData;
		WSAStartup( MAKEWORD(1, 1), &wsaData );
		gethostname( buf, len );
		WSACleanup();
	#else
		gethostname( buf, len );
	#endif
	}

	string GetHostname()
	{
		char buf[256];
		GetHostname( buf, sizeof(buf) );
		return string(buf);
	}

	long long GetAvailablePhysicalMemory()
	{
	#ifdef WIN32
		MEMORYSTATUSEX memoryStats;
		memoryStats.dwLength = sizeof( memoryStats );
		GlobalMemoryStatusEx( &memoryStats );
		return (long long) memoryStats.ullAvailPhys;
	#else
		long long kbMemFree, kbMemCached;
		string memVar;
		ifstream meminfoFile( "/proc/meminfo" );
		while( meminfoFile >> memVar )
		{
			if( memVar == "Cached:" )
			{
				meminfoFile >> kbMemCached;
			} else if( memVar == "MemFree:" )
			{
				meminfoFile >> kbMemFree;
			}
		}
		meminfoFile.close();
		/*struct sysinfo memoryStats;
		sysinfo( &memoryStats );
		cout << memoryStats.freeram * memoryStats.mem_unit << endl;
		return memoryStats.freeram * memoryStats.mem_unit;*/
		long long memAvailable = (kbMemFree + kbMemCached) * 1024;
		cout << memAvailable << endl;
		return memAvailable;
	#endif
	}

	int GetNumProcessors()
	{
	#ifdef WIN32
		SYSTEM_INFO info;
		GetSystemInfo( &info );
		return info.dwNumberOfProcessors;
	#else
		int numProcessors = 0;
		string cpuVar;
		ifstream cpuinfoFile( "/proc/cpuinfo" );
		while( cpuinfoFile >> cpuVar )
		{
			if( cpuVar == "processor" )
				++numProcessors;
		}
		return numProcessors;
	#endif
	}

	// MM/DD/YYYY
	string GetDateString()
	{
		time_t localTime = time(0);
		struct tm* tmObj = localtime( &localTime );
		stringstream s;

		s.fill('0'); s.width(2);
		s << right << tmObj->tm_mon+1;

		s << '/';

		s.fill('0'); s.width(2);
		s << right << tmObj->tm_mday;

		s << '/' << tmObj->tm_year + 1900;

		return s.str();
	}

	// HH:MM:SS
	string GetTimeString()
	{
		time_t localTime = time(0);
		struct tm* tmObj = localtime( &localTime );
		stringstream s;

		s.fill('0'); s.width(2);
		s << right << tmObj->tm_hour;

		s << ':';

		s.fill('0'); s.width(2);
		s << right << tmObj->tm_min;

		s << ':';

		s.fill('0'); s.width(2);
		s << right << tmObj->tm_sec;

		return s.str();
	}

	// YYYY-MM-DDTHH:MM:SS
	string GetDateTime( bool useLocal )
	{
		time_t now;
		time(&now);
		return GetDateTime( now, useLocal );
	}

	string GetDateTime( const time_t& now, bool useLocal )
	{
		struct tm* tmstruct = NULL;

		if( useLocal )
			tmstruct = localtime(&now);
		else 
			tmstruct = gmtime(&now);

		stringstream d;
		d << tmstruct->tm_year + 1900;
		d << '-';
		d.fill('0');
		d.width(2);
		d << right << tmstruct->tm_mon+1;
		d << '-';
		d.fill('0');
		d.width(2);
		d << right << tmstruct->tm_mday;

		stringstream t;
		t.fill('0');
		t.width(2);
		t << right << tmstruct->tm_hour;
		t << ':';
		t.fill('0');
		t.width(2);
		t << right << tmstruct->tm_min;
		t << ':';
		t.fill('0');
		t.width(2);
		t << right << tmstruct->tm_sec;

		return d.str() + 'T' + t.str();
	}

	string MakeProcessFilename( string str, int processID )
	{
		char buf[80];
		sprintf( buf, "%d", processID );
		str.insert( str.find_last_of( '.' ), buf );
		return str;
	}

	float round( float f, int precision )
	{
		if( f == 0.0f )
			return +0.0f;

		float multiplier = pow( 10.0f, (float) precision ); // moves f over <precision> decimal places
		f *= multiplier;
		f = floor( f + 0.5f );
		return f / multiplier;
	}

	double round( double f, int precision )
	{
		if( f == 0.0f )
			return +0.0f;

		double multiplier = pow( 10.0, (double) precision ); // moves f over <precision> decimal places
		f *= multiplier;
		f = floor( f + 0.5f );
		return f / multiplier;
	}

	string GetFileType( const string& filepath )
	{
		ifstream fileStream( filepath.c_str(), ios::binary );
		if( fileStream.is_open() )
		{
			string line1;
			std::getline( fileStream, line1 );

			// Is this an XML file?
			if( line1.find( "<?xml" ) == 0 )
			{
				// Yes, so the first (root) element gives the type
				size_t rootElIdx;
				rootElIdx = line1.find( '<', 1 );
				while( rootElIdx == string::npos || line1[rootElIdx+1] == '?' || line1[rootElIdx+1] == '!' )
				{
					std::getline( fileStream, line1 );
					rootElIdx = line1.find( '<' );
				}
				string rootEl = line1.substr( rootElIdx+1, line1.find_first_of( " >", rootElIdx+2 )-rootElIdx-1 );
				return to_lower_copy( rootEl );
			}

			if( line1.find( "H\tSQTGenerator" ) == 0 )
				return "sqt";
			else if( line1.find( "H\tTagsGenerator" ) == 0 || line1.find( "GutenTag" ) != string::npos )
				return "tags";
			else if( line1.find( '>' ) == 0 )
				return "fasta";
		}
		return "unknown";
	}

	bool TestFileType( const string& filepath, const string& type, bool printErrorMsg )
	{
		string actualType = GetFileType( filepath );
		if( actualType != type )
		{
			if( printErrorMsg )
				cerr << "Error: expected \"" << type << "\" file; the type of \"" << filepath << "\" is \"" << actualType << "\"" << endl;
			return false;
		}
		return true;
	}

	string GetFilenameFromFilepath( const string& filepath )
	{
		return filepath.substr( filepath.find_last_of( SYS_PATH_DELIMITER ) + 1 );
	}

	string GetPathnameFromFilepath( const string& filepath )
	{
		return filepath.substr( 0, filepath.find_last_of( SYS_PATH_DELIMITER ) + 1 );
	}

	string GetTopLevelOfFilepath( const string& filepath )
	{
		return filepath.substr( filepath.find_last_of( SYS_PATH_DELIMITER ) + 1 );
	}

	string GetFilenameWithoutExtension( const string& filename )
	{
		return basename( MAKE_PATH_FOR_BOOST( filename ) );
	}

	string GetFilenameExtension( const string& filename )
	{
		return extension( MAKE_PATH_FOR_BOOST( filename ) );
	}

	string ChangeFilenameExtension( const string& filename, const string& extension )
	{
		return change_extension( MAKE_PATH_FOR_BOOST( filename ), extension ).string();
	}

	long long GetFileSize( const string& filename )
	{
		return file_size( MAKE_PATH_FOR_BOOST( filename ) );
	}

	string GetFileLastModified( const string& filename )
	{
		return GetDateTime( last_write_time( MAKE_PATH_FOR_BOOST( filename ) ) );
	}

	string TrimWhitespace( const string& str )
	{
		string::size_type startIdx = str.find_first_not_of( "\t " );
		if( startIdx == string::npos )
			startIdx = 0;

		string::size_type endIdx = str.find_last_not_of( "\t " );
		if( endIdx == string::npos )
			endIdx = 0;

		return str.substr( startIdx, endIdx - startIdx + 1 );
	}

	string QuoteString( const string& str )
	{
		return string( "\"" ) + str + string( "\"" );
	}

	string UnquoteString( const string& str )
	{
		size_t firstNonQuote = str.find_first_not_of( '"' );
		size_t lastNonQuote = str.find_last_not_of( '"' ) + 1;
		if( firstNonQuote == string::npos )
			return "";

		return str.substr( firstNonQuote, lastNonQuote - firstNonQuote );
	}

	int GetIntByKey( const string& str, const string& key, const string& lineDelim )
	{
		size_t keyIdx = str.find( key );
		if( keyIdx == string::npos )
			return -1;

		size_t predIdx = str.find_first_of( '=', keyIdx ) + 1;
		size_t valIdx = str.find_first_not_of( "\t ", predIdx );
		size_t valCount = str.find_first_of( lineDelim, valIdx ) - valIdx;
		return atoi( str.substr( valIdx, valCount ).c_str() );
	}

	float GetFloatByKey( const string& str, const string& key, const string& lineDelim )
	{
		size_t keyIdx = str.find( key );
		if( keyIdx == string::npos )
			return -1.0f;

		size_t predIdx = str.find_first_of( '=', keyIdx ) + 1;
		size_t valIdx = str.find_first_not_of( "\t ", predIdx );
		size_t valCount = str.find_first_of( lineDelim, valIdx ) - valIdx;
		return (float) atof( str.substr( valIdx, valCount ).c_str() );
	}

	int GetBoolByKey( const string& str, const string& key, const string& lineDelim )
	{
		size_t keyIdx = str.find( key );
		if( keyIdx == string::npos )
			return -1;

		size_t predIdx = str.find_first_of( '=', keyIdx ) + 1;
		size_t valIdx = str.find_first_not_of( "\t ", predIdx );
		size_t valCount = str.find_first_of( lineDelim, valIdx ) - valIdx;
		return int( !str.substr( valIdx, valCount ).compare( "true" ) || !str.substr( valIdx, valCount ).compare( "1" ) );
	}

	string GetStringByKey( const string& str, const string& key, const string& lineDelim )
	{
		/*regex e( string( "^\\s*(" ) + key + ")\\s*=\\s*?(.*?)?\\s*?(#.*?)?$" );
		smatch KeyValuePair;
		if( !regex_search( str, KeyValuePair, e ) )
			return "";

		return KeyValuePair[2];*/
		size_t keyIdx = str.find( key );
		while( keyIdx != string::npos )
		{
			if( ( keyIdx == 0 || str[keyIdx-1] == ' ' || str[keyIdx-1] == '\t' || lineDelim.find( str[keyIdx-1] ) != string::npos ) &&
				( str[keyIdx + key.length()] == ' ' || str[keyIdx + key.length()] == '\t' || str[keyIdx + key.length()] == '=' ) )
				break;
			keyIdx = str.find( key, keyIdx + key.length() );
		}
		if( keyIdx == string::npos )
			return "";

		size_t predIdx = str.find_first_of( '=', keyIdx ) + 1;
		size_t valIdx = str.find_first_not_of( "\t ", predIdx );
		size_t valCount = str.find_first_of( lineDelim, valIdx ) - valIdx;
		return str.substr( valIdx, valCount );
	}
}