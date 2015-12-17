/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

extern char binaryCheckpointName[1024];
extern char binaryCheckpointBackupName[1024];
extern char binaryCheckpointInputName[1024];

extern double masterTime;
extern double accumulatedTime;

extern partitionLengths pLengths[MAX_MODEL];
extern int Thorough;
extern int optimizeRateCategoryInvocations;

static boolean treeNeedString(const char *fp, char c1, int *position)
{
  char 
    c2 = fp[(*position)++];
  
  if(c2 == c1)  
    return TRUE;
  else  
    {   
      int 
	lower = MAX(0, *position - 20),
	upper = *position + 20;
      
      printf("Tree Parsing ERROR: Expecting '%c', found: '%c'\n", c1, c2); 
      printf("Context: \n");
      
      while(lower < upper && fp[lower])
	printf("%c", fp[lower++]);
      
      printf("\n");

      return FALSE;
  }
} 



static boolean treeLabelEndString (char ch)
{
  switch(ch) 
    {   
    case '\0':  
    case '\t':  
    case '\n':  
    case '\r': 
    case ' ':
    case ':':  
    case ',':   
    case '(':   
    case ')':  
    case ';':
      return TRUE;
    default:
      break;
    }
  
  return FALSE;
} 

static boolean  treeGetLabelString (const char *fp, char *lblPtr, int maxlen, int *position)
{
  char 
    ch;
  
  boolean  
    done, 
    lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *)NULL; 
  else 
    if(lblPtr == NULL) 
      maxlen = 0;

  ch = fp[(*position)++];
  
  done = treeLabelEndString(ch);

  lblfound = !done;  

  while(!done) 
    {      
      if(treeLabelEndString(ch)) 
	break;     

      if(--maxlen >= 0) 
	*lblPtr++ = ch;
      
      ch = fp[(*position)++];      
    }
  
  (*position)--; 

  if (lblPtr != NULL) 
    *lblPtr = '\0';

  return lblfound;
}

static boolean  treeFlushLabelString(const char *fp, int *position)
{ 
  return  treeGetLabelString(fp, (char *) NULL, (int) 0, position);
} 


static boolean treeProcessLengthString (const char *fp, double *dptr, int *position)
{ 
  (*position)++;
  
  if(sscanf(&fp[*position], "%lf", dptr) != 1) 
    {
      printf("ERROR: treeProcessLength: Problem reading branch length\n");     
      assert(0);
    }

  while(fp[*position] != ',' && fp[*position] != ')' && fp[*position] != ';')
    *position = *position + 1;
  
  return  TRUE;
}

static int treeFlushLenString (const char *fp, int *position)
{
  double  
    dummy;  
  
  char     
    ch;

  ch = fp[(*position)++];
 
  if(ch == ':') 
    {     
      if(!treeProcessLengthString(fp, &dummy, position)) 
	return 0;
      return 1;	  
    }
    
  (*position)--;

  return 1;
} 

static int treeFindTipNameString (const char *fp, tree *tr, int *position)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabelString(fp, str, nmlngth+2, position))
    n = treeFindTipByLabelString(str, tr, TRUE);
  else
    n = 0;
   
  return  n;
} 

static boolean addElementLenString(const char *fp, tree *tr, nodeptr p, int *position)
{
  nodeptr  
    q;
  
  int      
    n, 
    fres;

  char 
    ch;
  
  if ((ch = fp[(*position)++]) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return FALSE;
	    }
	  else 
	    {	   
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (!addElementLenString(fp, tr, q->next, position))        
	return FALSE;
      if (!treeNeedString(fp, ',', position))             
	return FALSE;
      if (!addElementLenString(fp, tr, q->next->next, position))  
	return FALSE;
      if (!treeNeedString(fp, ')', position))             
	return FALSE;
      
     
      treeFlushLabelString(fp, position);
    }
  else 
    {   
      (*position)--;
     
      if ((n = treeFindTipNameString(fp, tr, position)) <= 0)          
	return FALSE;
      q = tr->nodep[n];
      
      if (tr->start->number > n)  
	tr->start = q;
      (tr->ntips)++;
    }
  
     
  fres = treeFlushLenString(fp, position);
  if(!fres) 
    return FALSE;
  
  hookupDefault(p, q, tr->numBranches);

  return TRUE;          
}




void treeReadTopologyString(char *treeString, tree *tr)
{ 
  char 
    *fp = treeString;

  nodeptr  
    p;
  
  int
    position = 0, 
    i;
  
  char 
    ch;   
    

  for(i = 1; i <= tr->mxtips; i++)    
    tr->nodep[i]->back = (node *)NULL;      
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;           
    }
      
  tr->start       = tr->nodep[1];
  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;    
  tr->rooted      = FALSE;      
  
  p = tr->nodep[(tr->nextnode)++]; 
   
  assert(fp[position++] == '(');  
    
  if (! addElementLenString(fp, tr, p, &position))                 
    assert(0);
  
  if (! treeNeedString(fp, ',', &position))                
    assert(0);
   
  if (! addElementLenString(fp, tr, p->next, &position))           
    assert(0);

  if(!tr->rooted) 
    {
      if ((ch = fp[position++]) == ',') 
	{ 
	  if (! addElementLenString(fp, tr, p->next->next, &position)) 
	    assert(0);	 
	}
      else 
	assert(0);     
    }
  else
    assert(0);
        
  if (! treeNeedString(fp, ')', &position))                
    assert(0);

  treeFlushLabelString(fp, &position);
  
  if (!treeFlushLenString(fp, &position))                         
    assert(0);
  
  if (!treeNeedString(fp, ';', &position))       
    assert(0);
    
  if(tr->rooted)     
    assert(0);           
  else           
    tr->start = tr->nodep[1];   

  printf("Tree parsed\n");

} 



/*
 * checkpoints BEGIN
 * */

static void writeTree(tree *tr, FILE *f)
{
  int
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;

  myfwrite(&(tr->start->number), sizeof(int), 1, f);
  myfwrite(&base, sizeof(nodeptr), 1, f);
  myfwrite(tr->nodeBaseAddress, sizeof(node), x, f);

}

int ckpCount = 0;

void writeCheckpoint(tree *tr)
{
  int
    model;

  if (ckpCount > 0)
     rename(binaryCheckpointName, binaryCheckpointBackupName);

  ckpCount++;

  FILE
    *f;

  f = myfopen(binaryCheckpointName, "w");

  /* cdta */

  tr->ckp.tr_startLH  = tr->startLH;
  tr->ckp.tr_endLH    = tr->endLH;
  tr->ckp.tr_likelihood = tr->likelihood;
  tr->ckp.tr_bestOfNode = tr->bestOfNode;

  tr->ckp.tr_lhCutoff = tr->lhCutoff;
  tr->ckp.tr_lhAVG    = tr->lhAVG;
  tr->ckp.tr_lhDEC    = tr->lhDEC;
  tr->ckp.tr_itCount  = tr->itCount;
  tr->ckp.tr_doCutoff = tr->doCutoff;

  tr->ckp.accumulatedTime = accumulatedTime + (gettime() - masterTime);

  /* printf("Acc time: %f\n", ckp.accumulatedTime); */

  myfwrite(&tr->ckp, sizeof(checkPointState), 1, f);

  myfwrite(tr->constraintVector, sizeof(int), 2 * tr->mxtips, f);

  myfwrite(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfwrite(tr->tree1, sizeof(char), tr->treeStringLength, f);

  myfwrite(&tr->rateHetModel, sizeof(int), 1, f);

  myfwrite(&tr->cdta->endsite, sizeof(int), 1, f);
  myfwrite(tr->patternPosition, sizeof(int), tr->rdta->sites, f);
  myfwrite(tr->columnPosition, sizeof(int), tr->rdta->sites, f);
  myfwrite(tr->cdta->alias, sizeof(int), tr->rdta->sites + 1, f);
  myfwrite(tr->cdta->aliaswgt, sizeof(int), tr->rdta->sites + 1, f);
  myfwrite(tr->model, sizeof(int), tr->rdta->sites + 1, f);
  myfwrite(tr->dataVector, sizeof(int), tr->rdta->sites + 1, f);

  myfwrite(tr->partitionContributions, sizeof(double), tr->NumberOfModels, f);

  if(tr->rateHetModel == CAT)
    {
      myfwrite(tr->cdta->rateCategory, sizeof(int), tr->rdta->sites + 1, f);
      myfwrite(tr->cdta->patrat, sizeof(double), tr->rdta->sites + 1, f);
      myfwrite(tr->cdta->patratStored, sizeof(double), tr->rdta->sites + 1, f);
    }

  /* pInfo */

  for(model = 0; model < tr->NumberOfModels; model++)
  {
    int
      dataType = tr->partitionData[model].dataType;

    myfwrite(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
    myfwrite(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
    myfwrite(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
    myfwrite(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
    myfwrite(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);

    myfwrite(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
    myfwrite(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
    myfwrite(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
    myfwrite(&(tr->partitionData[model].alpha), sizeof(double), 1, f);
    myfwrite(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);
  }

  writeTree(tr, f);

  fclose(f);

//  printBothOpen("\nCheckpoint written to: %s likelihood: %f\n", extendedName, tr->likelihood);
}

static void readTree(tree *tr, FILE *f, analdef *adef)
{
  int
    nodeNumber,
    x = tr->mxtips + 3 * (tr->mxtips - 1);

  nodeptr
    base = tr->nodeBaseAddress;



  nodeptr
    startAddress;

  myfread(&nodeNumber, sizeof(int), 1, f);

  tr->start = tr->nodep[nodeNumber];

  /*printf("Start: %d %d\n", tr->start->number, nodeNumber);*/

  myfread(&startAddress, sizeof(nodeptr), 1, f);

  /*printf("%u %u\n", (size_t)startAddress, (size_t)tr->nodeBaseAddress);*/



  myfread(tr->nodeBaseAddress, sizeof(node), x, f);

  {
    int i;

    size_t
      offset;

    boolean
      addIt;

    if(startAddress > tr->nodeBaseAddress)
    {
      addIt = FALSE;
      offset = (size_t)startAddress - (size_t)tr->nodeBaseAddress;
    }
    else
    {
      addIt = TRUE;
      offset = (size_t)tr->nodeBaseAddress - (size_t)startAddress;
    }

    for(i = 0; i < x; i++)
    {
      if(addIt)
      {
        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next + offset);
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back + offset);
      }
      else
      {

        tr->nodeBaseAddress[i].next = (nodeptr)((size_t)tr->nodeBaseAddress[i].next - offset);
        tr->nodeBaseAddress[i].back = (nodeptr)((size_t)tr->nodeBaseAddress[i].back - offset);
      }
    }

  }
}

void readCheckpoint(tree *tr, analdef *adef, boolean readModel)
{
  int
    model;

  FILE
    *f = myfopen(binaryCheckpointInputName, "r");

  /* cdta */

  myfread(&tr->ckp, sizeof(checkPointState), 1, f);

  tr->ntips = tr->mxtips;

  tr->startLH    = tr->ckp.tr_startLH;
  tr->endLH      = tr->ckp.tr_endLH;
  tr->likelihood = tr->ckp.tr_likelihood;
  tr->bestOfNode = tr->ckp.tr_bestOfNode;

  tr->lhCutoff   = tr->ckp.tr_lhCutoff;
  tr->lhAVG      = tr->ckp.tr_lhAVG;
  tr->lhDEC      = tr->ckp.tr_lhDEC;
  tr->itCount    = tr->ckp.tr_itCount;
  Thorough       = tr->ckp.Thorough;

  accumulatedTime = tr->ckp.accumulatedTime;

  /* printf("Accumulated time so far: %f\n", accumulatedTime); */

  optimizeRateCategoryInvocations = tr->ckp.optimizeRateCategoryInvocations;

  myfread(tr->constraintVector, sizeof(int), 2 * tr->mxtips, f);

  myfread(tr->tree0, sizeof(char), tr->treeStringLength, f);
  myfread(tr->tree1, sizeof(char), tr->treeStringLength, f);

  int ratehet;
  myfread(&ratehet, sizeof(int), 1, f);

  assert(ratehet == tr->rateHetModel);

  if (adef->newCheckpoint)
  {
      myfread(&tr->cdta->endsite, sizeof(int), 1, f);
      myfread(tr->patternPosition, sizeof(int), tr->rdta->sites, f);
      myfread(tr->columnPosition, sizeof(int), tr->rdta->sites, f);
      myfread(tr->cdta->alias, sizeof(int), tr->rdta->sites + 1, f);
      myfread(tr->cdta->aliaswgt, sizeof(int), tr->rdta->sites + 1, f);
      myfread(tr->model, sizeof(int), tr->rdta->sites + 1, f);
      myfread(tr->dataVector, sizeof(int), tr->rdta->sites + 1, f);

      myfread(tr->partitionContributions, sizeof(double), tr->NumberOfModels, f);
  }

  if(tr->rateHetModel == CAT)
    {
      myfread(tr->cdta->rateCategory, sizeof(int), tr->rdta->sites + 1, f);
      myfread(tr->cdta->patrat, sizeof(double), tr->rdta->sites + 1, f);
      myfread(tr->cdta->patratStored, sizeof(double), tr->rdta->sites + 1, f);
    }

  if(tr->searchConvergenceCriterion && !readModel)
  {
    int bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 0) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 0))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 0\n");
#endif

      treeReadTopologyString(tr->tree0, tr);

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->convHashT, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE);

      assert(bCounter == tr->mxtips - 3);
    }

    bCounter = 0;

    if((tr->ckp.state == FAST_SPRS && tr->ckp.fastIterations > 1) ||
        (tr->ckp.state == SLOW_SPRS && tr->ckp.thoroughIterations > 1))
    {

#ifdef _DEBUG_CHECKPOINTING
      printf("parsing Tree 1\n");
#endif

      treeReadTopologyString(tr->tree1, tr);

      bitVectorInitravSpecial(tr->bitVectors, tr->nodep[1]->back, tr->mxtips, tr->vLength, tr->convHashT, 1, BIPARTITIONS_RF, (branchInfo *)NULL,
          &bCounter, 1, FALSE, FALSE);

      assert(bCounter == tr->mxtips - 3);
    }
  }

  /* pInfo */

  if (readModel)
    {
      for(model = 0; model < tr->NumberOfModels; model++)
      {
        int
          dataType = tr->partitionData[model].dataType;

        myfread(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
        myfread(tr->partitionData[model].perSiteRates, sizeof(double), tr->maxCategories, f);
        myfread(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
        myfread(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
        myfread(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);

        myfread(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
        myfread(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);
        myfread(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);
        myfread(&(tr->partitionData[model].alpha), sizeof(double), 1, f);
        myfread(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);
        makeGammaCats(tr->rateHetModel, tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useGammaMedian, tr->partitionData[model].propInvariant);
      }

    #ifdef _FINE_GRAIN_MPI
      masterBarrierMPI(THREAD_COPY_INIT_MODEL, tr);
    #endif

    #ifdef _USE_PTHREADS
      masterBarrier(THREAD_COPY_INIT_MODEL, tr);
    #endif

      updatePerSiteRates(tr, FALSE);

      readTree(tr, f, adef);
    }

  fclose(f);
}

void checkpointTreeEval(tree *tr)
{
  evaluateGenericInitrav(tr, tr->start);

  printBothOpen("RAxML Restart with likelihood: %1.50f\n", tr->likelihood);
}

void restart(tree *tr, analdef *adef)
{
  switch(tr->ckp.state)
  {
    case REARR_SETTING:
      break;
    case FAST_SPRS:
      break;
    case SLOW_SPRS:
      break;
    case EPA_L1OUT:
      break;
    default:
      assert(0);
  }

  checkpointTreeEval(tr);
}

