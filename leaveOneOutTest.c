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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI) || defined (_SATIVA_MPI))
extern int processes;
#endif

extern int processID;

extern char **globalArgv;
extern int globalArgc;

extern char run_id[128];
extern char workdir[1024];
extern double masterTime;

typedef struct {
    double lh;
    double pendantBranchLen;
    int branchNumber;
} LeaveBranchInfo;

LeaveBranchInfo *brInfo;
nodeptr *branchNodeMap;
int* branchList;

#ifdef _PROFILE_L1OUT
  double
    newviewTime,
    evalTime,
    smoothTime;
#endif

static double getBranch(tree *tr, double *b, double *bb)
{
  double z = 0.0;

  if(!tr->multiBranch)
    {
      assert(b[0] == bb[0]);
      z = b[0];
      if (z < zmin)
      z = zmin;
      if(z > zmax)
      z = zmax;
      z = -log(z);
      return z;
    }
  else
    {
      int i;
      double x = 0;

      for(i = 0; i < tr->numBranches; i++)
      {
	assert(b[i] == bb[i]);
	assert(tr->partitionContributions[i] != -1.0);
	x = b[i];
	if (x < zmin)
	  x = zmin;
	if(x > zmax)
	  x = zmax;
	x = -log(x);

	z += x * tr->partitionContributions[i];
      }

      return z;
    }

}

static void testInsertBranch(tree *tr, nodeptr r, int branchNumber, boolean doSmoothing)
{
  nodeptr
    q = branchNodeMap[brInfo[branchNumber].branchNumber];

  if (!q) {
      printf("Missing node: %d\n", branchNumber);
      return;
  }

  double
    result,
    qz[NUM_BRANCHES],
    z[NUM_BRANCHES];

  nodeptr
    x = q->back,
    s = r->back;

  int
    j;

  for(j = 0; j < tr->numBranches; j++)
    {
      qz[j] = q->z[j];
      z[j] = sqrt(qz[j]);

      if(z[j] < zmin)
    z[j] = zmin;

      if(z[j] > zmax)
    z[j] = zmax;
    }

  hookup(r->next,       q, z, tr->numBranches);
  hookup(r->next->next, x, z, tr->numBranches);
  hookupDefault(r, s, tr->numBranches);

#ifdef _PROFILE_L1OUT
  double
    lastTime = gettime();
#endif

  newviewGeneric(tr, r);

#ifdef _PROFILE_L1OUT
  newviewTime += gettime() - lastTime;
  lastTime = gettime();
#endif

  if (doSmoothing)
      localSmooth(tr, r, smoothings);

#ifdef _PROFILE_L1OUT
  smoothTime += gettime() - lastTime;
  lastTime = gettime();
#endif

  result = evaluateGeneric(tr, r);

#ifdef _PROFILE_L1OUT
  evalTime += gettime() - lastTime;
#endif

  brInfo[branchNumber].lh = result;
  brInfo[branchNumber].pendantBranchLen = getBranch(tr, r->z, r->back->z);

  hookup(q, x, qz, tr->numBranches);

  r->next->next->back = r->next->back = (nodeptr) NULL;
}

int branchInfoCompare(const void *p1, const void *p2)
{
  LeaveBranchInfo *rc1 = (LeaveBranchInfo *)p1;
  LeaveBranchInfo *rc2 = (LeaveBranchInfo *)p2;

  double i = rc1->lh;
  double j = rc2->lh;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}

static void setupNodeBranches(tree *tr, nodeptr p, LeaveBranchInfo *brInfo)
{
    if(!isTip(p->number, tr->mxtips))
      {
        setupNodeBranches(tr, p->next->back, brInfo);
        setupNodeBranches(tr, p->next->next->back, brInfo);
      }

    branchNodeMap[tr->branchCounter] = p;
    brInfo[tr->branchCounter].branchNumber = tr->branchCounter;
    tr->branchCounter++;
}

static void initBranchInfo(tree *tr, LeaveBranchInfo *brInfo)
{
    nodeptr
      originalNode = tr->nodep[tr->mxtips + 1];

    tr->branchCounter = 0;

    setupNodeBranches(tr, originalNode->back, brInfo);
    setupNodeBranches(tr, originalNode->next->back, brInfo);
    setupNodeBranches(tr, originalNode->next->next->back, brInfo);
}

void leaveOneOutTest(tree *tr, analdef *adef)
{
  int
    k,
    i,
    tips,
    numTraversalBranches = (2 * (tr->mxtips - 1)) - 3; /* compute number of branches into which we need to insert once we have removed a taxon */

  const int nodeCount = tr->mxtips * 2 - 2;
  const int branchesCount = tr->mxtips * 2 - 3;

  tr->numberOfBranches = 2 * tr->ntips - 3;

  printBothOpen("RAxML Leave-one-out test: each taxa will be pruned from tree and re-inserted using EPA.\n\n");
  printBothOpen("Reference tree: # alignment patterns: %d, # taxa: %d, # branches %d \n\n", tr->cdta->endsite, tr->ntips, tr->numberOfBranches);

  int
    slowInsertions;

  if (tr->useEpaHeuristics)
  {
    slowInsertions =  MAX(5, (int)(0.5 + (double)(tr->numberOfBranches) * tr->fastEPAthreshold));
    printBothOpen("Using EPA heuristics: %d out of %d branches will be thoroughly evaluated\n\n", slowInsertions, tr->numberOfBranches);
  }
  else
    slowInsertions = tr->numberOfBranches - 2;

  brInfo = (LeaveBranchInfo*) rax_calloc(branchesCount, sizeof(LeaveBranchInfo));
//  LeaveBranchInfo* brInfo2 = (LeaveBranchInfo*) rax_calloc(branchesCount, sizeof(LeaveBranchInfo));

  LeaveBranchInfo* brInfoOrig = (LeaveBranchInfo*) rax_calloc(branchesCount, sizeof(LeaveBranchInfo));

  branchNodeMap = (nodeptr*) rax_calloc(branchesCount, sizeof(nodeptr));

  char
    fileName[1024];

  FILE
    *jsonFile;

  strcpy(fileName,         workdir);
  strcat(fileName,         "RAxML_leaveOneOutResults.");
  strcat(fileName,         run_id);
#ifdef _SATIVA_MPI
  {
    char
      buf[64];

    sprintf(buf, "%d", processID);
    strcat(fileName,       ".PID.");
    strcat(fileName,       buf);
  }
#endif
  strcat(fileName,         ".jplace");

  jsonFile = myfopen(fileName, "wb");

  fprintf(jsonFile, "{\n");
  fprintf(jsonFile, "\t\"placements\": [\n");

  printBothOpen("Likelihood of comprehensive tree %f\n\n", tr->likelihood);

  initBranchInfo(tr, brInfoOrig);

  assert(tr->branchCounter == tr->numberOfBranches);

#ifdef _PROFILE_L1OUT
  double
    lastTime,
    fastInsertTime,
    fastSortTime,
    insertTime,
    sortTime;

  int
    newviewFast,
    newviewSlow;
#endif

int
  startTip,
  endTip;

#ifdef _SATIVA_MPI
  int slice = tr->mxtips / processes + 1;
  startTip = processID * slice + 1;
  endTip = MIN(tr->mxtips, (processID + 1) * slice);

  printBothOpen("Running in MPI mode: each of %d processes will process %d taxa\n\n", processes, slice);
#else
  startTip = 1;
  endTip = tr->mxtips;
#endif

  double
    tipStartTime;

  /* prune and re-insert one tip at a time into all branches of the remaining tree */
  for(tips = startTip; tips <= endTip; tips++)
  {

#ifdef _PROFILE_L1OUT
      lastTime = gettime();
      newviewTime = 0.;
      evalTime = 0.;
      smoothTime = 0.;
      tr->mr_thresh = 0;
#endif

      tipStartTime = gettime();

      nodeptr
        myStart,
        p = tr->nodep[tips]->back, /* this is the node at which we are pruning */
        p1 =  p->next->back,
        p2 =  p->next->next->back;

      double
        pz[NUM_BRANCHES],
        p1z[NUM_BRANCHES],
        p2z[NUM_BRANCHES];

      /* store the three branch lengths adjacent to the position at which we prune */

      for(i = 0; i < tr->numBranches; i++)
      {
          p1z[i] = p1->z[i];
          p2z[i] = p2->z[i];
          pz[i] = p->z[i];
      }

      for(i = 0; i < branchesCount; i++)
      {
          brInfo[i].branchNumber = i;
          brInfo[i].lh = unlikely;
      }

      /* prune the taxon, optimizing the branch between p1 and p2 */
      removeNodeBIG(tr, p,  tr->numBranches);

      if (adef->verbose)
	printBothOpen("Pruning taxon Number %d [%s]\n", tips, tr->nameList[tips]);

      /* first pass - heuristic (no smoothing) */
      for(i = 0; i < branchesCount; i++)
      {
          nodeptr np = branchNodeMap[i];

          if (np == tr->nodep[tips] || np->number == p->number)
          {
//              printf("Ignored: %d %d %d\n", i, np->number, np->back ? np->back->number : -1);
              continue;
          }

          testInsertBranch(tr, p, i, !tr->useEpaHeuristics);
      }

#ifdef _PROFILE_L1OUT
      newviewFast = tr->mr_thresh;
      tr->mr_thresh = 0;
      fastInsertTime = gettime() - lastTime;
      lastTime = gettime();
#endif

//      memcpy(brInfo2, brInfo, sizeof(LeaveBranchInfo) * branchesCount);

      qsort(brInfo, branchesCount, sizeof(LeaveBranchInfo), branchInfoCompare);

#ifdef _PROFILE_L1OUT
      fastSortTime = gettime() - lastTime;
      lastTime = gettime();
#endif

      /* second pass - slow inserts (with smoothing) */
      if (tr->useEpaHeuristics)
      {
          for(i = 0; i < slowInsertions; i++)
          {
              testInsertBranch(tr, p, i, TRUE);
          }

//          for(i = 0; i < branchesCount; i++)
//          {
//              if (brInfo[i].lh > brInfo2[slowInsertions].lh)
//        	testInsertBranch(tr, i, TRUE);
//          }

#ifdef _PROFILE_L1OUT
          newviewSlow = tr->mr_thresh;
          tr->mr_thresh = 0;
          insertTime = gettime() - lastTime;
          lastTime = gettime();
#endif

          qsort(brInfo, slowInsertions, sizeof(LeaveBranchInfo), branchInfoCompare);

#ifdef _PROFILE_L1OUT
          sortTime = gettime() - lastTime;
          lastTime = gettime();
#endif
      }

      // TODO reliably determine original branch #
      int origNode = isTip(p1->number, tr->mxtips) ? p1->number : p2->number;
      int origBranch = -1;

      for(i = 0; i < branchesCount; i++)
      {
          if (branchNodeMap[i]->number == origNode)
          {
              origBranch = i;
              break;
          }
      }
//      if (origBranch != brInfo[0].branchNumber)
//          printf("Nodes: %d %d %d %d\n", p1->number, p2->number, p1->back->number, p2->back->number);

      /* re-connect taxon to its original position  */

      hookup(p->next,       p1,      p1z, tr->numBranches);
      hookup(p->next->next, p2,      p2z, tr->numBranches);
      hookup(p,             p->back, pz, tr->numBranches);

      /* fix likelihood vectors */

      newviewGeneric(tr, p);

      int validEntries = 0;
      int j;
      for(j =  0; j < branchesCount; j++)
        if(brInfo[j].lh <= unlikely)
          break;
        else
          validEntries++;

      assert(validEntries > 0);

      double
	lmax = brInfo[0].lh,
        all = 0.;

      for(j =  0; j < validEntries; j++)
        all += EXP(brInfo[j].lh - lmax);

      fprintf(jsonFile, "\t{\"p\":[");

      double
	acc = 0.,
	maxprob = 0.;

      for(j =  0; j < validEntries; j++)
      {
	  double
	    prob = (EXP(brInfo[j].lh - lmax) / all);

	  if ((tr->useAccumulatedEPACutoff && acc > tr->accumulatedEPACutoff) ||
	      (!tr->useAccumulatedEPACutoff && prob < maxprob * tr->probThresholdEPA))
	    break;

          if (j == 0)
            maxprob = prob;
          else
            fprintf(jsonFile, ",");

	  acc += prob;

          fprintf(jsonFile, "[%d, %f, %f, %f, %f]", brInfo[j].branchNumber, brInfo[j].lh, prob, 0.0, brInfo[j].pendantBranchLen);
      }

      fprintf(jsonFile, "], \"n\":[\"%s\"]}", tr->nameList[tips]);
      if (tips == endTip)
          fprintf(jsonFile, "\n");
      else
          fprintf(jsonFile, ",\n");

      if ((tips - startTip - 1) % 10 == 0)
	fflush(jsonFile);

      if (adef->verbose)
	printBothOpen("[%.3f s] Placement:  %d %d %f\n", gettime() - tipStartTime, brInfo[0].branchNumber, origBranch, brInfo[0].lh);

#ifdef _PROFILE_L1OUT
      double finalizeTime = gettime() - lastTime;
      printBothOpen("newview calls fast/slow: %d %d\n", newviewFast, newviewSlow);
      printBothOpen("Elapsed time: %f %f %f %f %f \n", fastInsertTime, fastSortTime, insertTime, sortTime, finalizeTime);
      printBothOpen("newview/smooth/eval time: %f %f %f\n\n", newviewTime, smoothTime, evalTime);
#endif
    }

  fprintf(jsonFile, "\t ],\n");
  fprintf(jsonFile, "\t\"metadata\": {\"invocation\": ");

  fprintf(jsonFile, "\"");
  {
    int i;

    for(i = 0; i < globalArgc; i++)
      fprintf(jsonFile,"%s ", globalArgv[i]);
  }
  fprintf(jsonFile, "\", \"raxml_version\": \"%s\"", programVersion);
  fprintf(jsonFile,"},\n");

  fprintf(jsonFile, "\t\"version\": 2,\n");
  fprintf(jsonFile, "\t\"fields\": [\n");
  fprintf(jsonFile, "\t\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \n");
  fprintf(jsonFile, "\t\"pendant_length\"\n");
  fprintf(jsonFile, "\t]\n");
  fprintf(jsonFile, "}\n");

  fclose(jsonFile);

  rax_free(branchNodeMap);
  rax_free(brInfoOrig);
  rax_free(brInfo);

  printBothOpen("\nTime for EPA leave-one-out test: %f\n", gettime() - masterTime);
  printBothOpen("Best placements written to file %s\n", fileName);

//  exit(0);
}

