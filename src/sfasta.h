/* (C) Copyright 2016 Robert C. Edgar, all rights reserved. */

#ifndef sfasta_h
#define sfasta_h

#include "myutils.h"
#include "seq.h"

typedef void (*fn_OnSeq)(const string &Label,
  const string &Seq, void *ptrUserData);

void ReadFasta(const string &FileName, fn_OnSeq OnSeq,
  void *UserData);

#endif // sfasta_h
