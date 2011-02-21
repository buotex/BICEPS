#ifndef _TOPSET_H
#define _TOPSET_H

#include "stdafx.h"
#include "Profiler.h"
namespace freicore
{
	template< class T, class ComparePredicate = std::less<T> >
	class topset : public set<T, ComparePredicate>
	{
		typedef set<T, ComparePredicate> MyBase;
	public:
		topset( size_t maxSize = 0 ) : MyBase(), m_maxSize( maxSize ) {}

		//Profiler insertTime;

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & boost::serialization::base_object< MyBase >( *this );
			ar & m_maxSize;// & m_permSet;
		}

		size_t max_size()
		{
			return m_maxSize;
		}

		void max_size( size_t maxSize )
		{
			m_maxSize = maxSize;
			trim();
		}

		void trim()
		{
			if( !m_maxSize )
				return;

			while( MyBase::size() > m_maxSize )
				MyBase::erase( set<T>::begin() );
		}

		void clear()
		{
			MyBase::clear();
			//m_permSet.clear();
		}

		TemplateSetInsertPair(T) insert( const T& value, bool noMatterWhat = false )
		{
			// If m_maxSize is not set (0), do a regular insert
			if( !m_maxSize || noMatterWhat )
				return MyBase::insert( value );
			else
			{
				typename MyBase::iterator itr = MyBase::find( value );

				// If the new value is already in the set, the insert fails
				if( itr != MyBase::end() )
				{
					//insertTime.End();
					return TemplateSetInsertPair(T)( itr, false );
				}

				// If set is not full, add the new value
				if( MyBase::size() < m_maxSize )
				{
					TemplateSetInsertPair(T) result = MyBase::insert( value );
					//insertTime.End();
					return result;
				}

				// If set is full and the new value is better than the worst existing value,
				// add the new value and remove the worst value
				else if( *MyBase::begin() < value )
				{
					MyBase::erase( set<T>::begin() );
					TemplateSetInsertPair(T) result = MyBase::insert( value );
					//insertTime.End();
					return result;
				}
			}

			// Something went wrong...
			//insertTime.End();
			return TemplateSetInsertPair(T)( MyBase::end(), false );
		}

	protected:
		size_t	m_maxSize;
	};
}

#endif
