#include "myutils.h"

void WriteAln(FILE *f, const byte *A, const byte *B,
  const char *Path, unsigned ColCount)
	{
	unsigned p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'D')
			fprintf(f, "%c", A[p++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, "\n");

	unsigned pa = 0;
	unsigned pb = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M')
			{
			byte a = A[pa];
			byte b = B[pb];
			if (toupper(a) == toupper(b))
				fprintf(f, "|");
			else
				fprintf(f, " ");
			}
		else
			fprintf(f, " ");
		if (c == 'M' || c == 'D')
			++pa;
		if (c == 'M' || c == 'I')
			++pb;
		}
	fprintf(f, "\n");

	p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'I')
			fprintf(f, "%c", B[p++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, "\n");
	}
