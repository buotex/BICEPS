#ifndef _BASE64_H
#define _BASE64_H

#include "stdafx.h"
#include "Profiler.h"

namespace freicore
{
	extern const signed char base64_to_binary[];
	extern Timer foo;

	template< class outType, class inType >
	void base64_decode_sequence(	const string&		encodedSequenceIn,
									vector< outType >&	decodedSequenceOut,
									bool				swapByteOrder )
	{
		size_t encodedSequenceInLength = encodedSequenceIn.length();
		if( encodedSequenceInLength == 0 )
			return;

		size_t numQuads = encodedSequenceInLength / 4;
		size_t numBytes = numQuads * 3;
		//numBytes = numBytes + ( numBytes%4 ? 4 - numBytes%4 : 0 );
		//printf( "\n(%d digits, %d quads, %d bytes)\n", base64Data->length(), numQuads, numBytes );

		unsigned char* bytes = new unsigned char[ numBytes ];

		string tmp;
		tmp.resize( encodedSequenceInLength );

		for( size_t i = 0; i < encodedSequenceInLength; ++i )
		{
			tmp[i] = base64_to_binary[ (size_t) encodedSequenceIn[i] ];
		}

		for( size_t numHandled=0; numHandled < numQuads; ++numHandled )
		{
			bytes[numHandled*3]		= (char) ( tmp[numHandled*4] << 2 | tmp[numHandled*4+1] >> 4);
			bytes[numHandled*3+1]	= (char) ( tmp[numHandled*4+1] << 4 | tmp[numHandled*4+2] >> 2);
			bytes[numHandled*3+2]	= (char) ((( tmp[numHandled*4+2] << 6) & 0xC0) | tmp[numHandled*4+3]);
		}

		size_t stride = sizeof( inType );

		if( swapByteOrder )
		{
			for( size_t i = 0; i <= numBytes - stride; i += stride )
			{
				for( size_t n = 0; n < stride / 2; ++n )
				{
					/*char tmp = bytes[i+n];
					bytes[i+n] = bytes[i+stride-1-n];
					bytes[i+stride-1-n] = tmp;*/
					std::swap( bytes[i+n], bytes[i+stride-1-n] );
				}
			}
		}

		inType* data = reinterpret_cast< inType* >( bytes );
		size_t numData = numBytes / stride;

		decodedSequenceOut.resize(numData);
		for( size_t i = 0; i < numData; ++i )
		{
				decodedSequenceOut[i] = static_cast< outType >( data[i] );
		}

		delete [] bytes;
		//cout << '\n' << numData << '\t' << foo.End();
	}

	int b64_encode (char *dest,
			const unsigned char *src,
			int len);
	int b64_decode (char *dest,
			const char *src);
}

#endif /* BASE64_H */
