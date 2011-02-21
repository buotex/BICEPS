#ifndef _LNFACTORIALTABLE_H
#define _LNFACTORIALTABLE_H

#include "stdafx.h"

namespace freicore
{
	class lnFactorialTable
	{
	public:
		lnFactorialTable()
		{
			m_table.push_back(0);
			m_table.push_back(0);
		}

		double operator[]( size_t index )
		{
			// Is the table big enough?
			size_t maxIndex = m_table.size() - 1;
			if( index > maxIndex )
			{
				while( index > maxIndex )
				{
					m_table.push_back( m_table[ maxIndex ] + log( (float) m_table.size() ) );
					++maxIndex;
				}
			}

			return m_table[ index ];
		}

		void resize( size_t maxIndex )
		{
			this->operator []( maxIndex );
		}

	private:
		std::vector< double > m_table;
	};
}

#endif
