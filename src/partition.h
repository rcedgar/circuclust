/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#ifndef Part_h
#define Part_h

#include "gobuff.h"

class PartMem
	{
public:
	static const unsigned NVEC = 8;

public:
	unsigned *m_Vecs[NVEC];
	unsigned m_VecPos[NVEC];
	unsigned m_MaxValueCount;

	GoBuff<unsigned> m_Sizes;
	GoBuff<unsigned> m_Offsets;

public:
	PartMem()
		{
		m_MaxValueCount = 0;
		zero(m_Vecs, NVEC);
		}

	void Free()
		{
		for (unsigned i = 0; i < NVEC; ++i)
			{
			myfree(m_Vecs[i]);
			m_Vecs[i] = 0;
			}
		m_MaxValueCount = 0;
		}

	void Alloc(unsigned ValueCount)
		{
		if (ValueCount <= m_MaxValueCount)
			return;

		Free();
		
		m_MaxValueCount = ValueCount;
		for (unsigned i = 0; i < NVEC; ++i)
			m_Vecs[i] = myalloc(unsigned, m_MaxValueCount);
		}
	};

unsigned PartOrderDesc(const unsigned *Values, unsigned ValueCount,
  PartMem &Mem, unsigned *Order);
unsigned PartSubsetDesc(const unsigned *Values, unsigned ValueCount,
  PartMem &Mem, const unsigned *Subset, unsigned *Result);

#endif // Part_h
