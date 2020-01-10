/* bpp.c
   Markov chain Monte Carlo coalescent program for population genetics and
   phylogeographic data.

   Copyright by Ziheng Yang, since July 2002

   Linux/UNIX gcc/icc:
      icc -o bpp -O3 bpp.c tools.c
      icc -o bpp_sse -O3 -DUSE_SSE -msse3 bpp.c tools.c -lm
      icc -o bpp_avx -O3 -DUSE_AVX -mavx bpp.c tools.c -lm
      icc -o MCcoal -DSIMULATION -O3 bpp.c tools.c -lm

      gcc -o bpp -O3 bpp.c tools.c -lm
      gcc -o bpp_sse -O3 -DUSE_SSE -msse3 bpp.c tools.c -lm
      gcc -o bpp_avx -O3 -DUSE_AVX -mavx bpp.c tools.c -lm
      gcc -o MCcoal -DSIMULATION -O3 bpp.c tools.c -lm

   MAC OSX intel:
      cc -o bpp -O3 bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION -O3 bpp.c tools.c -lm

   Windows MSC++:
      cl -Febpp_NoSSE.exe -O2 -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000
      cl -Febpp_sse.exe -DUSE_SSE -arch:SSE2 -Ox -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000
      cl -Febpp.exe -DUSE_AVX -arch:AVX -Ox -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000
      cl -O2 -FeMCcoal.exe -DSIMULATION bpp.c tools.c /F10000000
*/

/*
#define USE_AVX
*/

#define SIMULATION


#define NSPECIES        32         /* max # of species */
#define NS              300        /* max # of sequences per locus */
#define NBRANCH         NS*2-2     /* max # of branches.       Don't change this */
#define MAXNSONS        40         /* max # of sons per node.  Don't change this */
#define NGENE           12500      /* max # of loci */
#define LSPNAME         4000       /* max # characters in a sequence name */
#define NPopLongNames   20         /* Node numbers are used if more than 20 species */

#include "paml.h"

#if (defined(USE_AVX) && defined(USE_SSE))
#error "Cannot use -DUSE_AVX and -DUSE_SSE together"
#endif

#if (defined(USE_AVX) || defined(USE_SSE))
#if (defined(USE_SSE))
#include <pmmintrin.h>
#else
#include <immintrin.h>
#endif
#endif

/* define alignment size for aligned memory allocation */
#if (defined(USE_AVX))
#define PLL_BYTE_ALIGNMENT 32
#else
#define PLL_BYTE_ALIGNMENT 16
#endif

struct CommonInfo {
   unsigned char *z[2 * NS - 1];
   char *spname[NS], seqf[2048], Imapf[2048], treef[2048], outf[2048], mcmcf[2048];
   char checkpointf[2048], concatf[2048], locusratef[2048], modelparaf[2048], heredityf[2048];
   char oldconP[NS * 2 - 1];
   int  cleandata, checkpoint;
   int  ns, ls, lgene[3], model, clock, simulation, SeqFileFormat;
   int  seqtype, ngene, *pose, npatt, np, readpattern;
   int  ntime, ncatG, ncode, fix_rgene, posG[2];
   double alpha, pi[4], piG[1][4], *rates, rgene[1];
   double *conP;  /* not used */
   int curconP;
   size_t sconP, sizespace;
   double *conPin[2], *fpatt, *space;
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[MAXNSONS], ibranch, ipop;  /* ibranch not used? */
   double branch, age, *conP, label, label2;         /* age is uptodate but not branch */
   char fix_age, *annotation;                        /* not used */
}  *nodes, *gnodes[NGENE], nodes_t[2 * NS - 1], *gnodes_t[NGENE];
/*  nodes is a pointer, gnodes holds the gene trees, nodes_t is temp space
    gnodes is not allocated for the simulation program
*/

/* stree.nseqsp is working space, copied from data.nseqsp[locus], which holds
   the info for each locus.  It seems feasible to remove data.nseqsp[locus].  Its
   major use is calculation of ndesc[] in UpdateGB_SPR, and this can be done by
   checking ipop for tips of gene tree.  The same applies to other uses of nseqsp.
   npop is the number of theta's.
   stree.pops[] holds the node numbers in the species tree of the npop populations.
   It is not always true that younger pops are before ancestral pops since the species
   tree may change in the program.
   stree.pptable[i][j] = 1 if pop j is ancestral to i and = 0 otherwise.
   For modern species i, table[i][i] = 1 if and only if there is a theta for pop i
   (but see below).

   28 Sept 2011: All diagonals of the table are set to 1, with pptable[i][i] = 1 even if
   there are never 2 sequences from pop i.  This made the Split and Join moves easier.
   Watch out for errors in the routines for updating the gene trees.
*/

struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nseqsp[NSPECIES];
   int npop, pops[NSPECIES * 2 - 1];
   signed char pptable[NSPECIES * 2 - 1][NSPECIES * 2 - 1];
   int NoTheta, speciesdelimitation, speciestree, SpeciesModelPrior, analysis;
   int migration;
   double *M;
   double roottau;  /* dist is seq distance */
   int nModels, iModel;  /* for species delimitation */
   char *DelimitationModels;  /* for species delimitation */
   double *pmodel, PriorSA11[NSPECIES];  /* pmodel is used in MCMC and SummarizeA10.  Make it local? */
   struct TREESPN {
      char name[LSPNAME * 2];
      int father, nson, sons[2];
      double age, theta;
   } nodes[2 * NSPECIES - 1];
}  stree;


/* spnames[] has format "125", "23", etc., where the indexes refer to those in control file.
   The idea breaks down if NSPECIES>256.   */
struct SMODEL {
   int  nspecies;
   char spnames[NSPECIES][NSPECIES + 1], splits[(NSPECIES - 2)*(NSPECIES + 1)];
};

struct SDELIMIT {
   int  nspecies;
   char spnames[NSPECIES][NSPECIES + 1];
};



#ifdef SIMULATION

struct DATA {
   double e_seqerr[NSPECIES][4 * 4];
   int nseqerr, iseqerr[NSPECIES], diploid[NSPECIES], ngene;
   int fix_Qrates, fix_pi, fix_alpha;
   double Qrates_Dir[6], Qrates[6], pi_Dir[4], a_siterate[2], a_locusrate, a_branchrate;
   double *siterates;
}  data;

#else

enum { A00, A01, A10, A11 } ANALYSES; /* stree.speciesdelimidtation & stree.speciestree */
enum { p0UniformLH, p1UniformRootedTree, p2UniformSLH, p3UniformSRootedTree } SPECIESMODELPRIORS;

struct DATA { /* locus-specific data and gene tree information, used in lnpG */
   int maxns, ngene, lgene[NGENE];
   int ns[NGENE], nseqsp[NGENE][NSPECIES], ls[NGENE], npatt[NGENE], npattdiploid[NGENE];
   int root[NGENE + 1], conP_offset[NGENE], *pose[NGENE], *Imap;
   int diploid[NSPECIES];
   unsigned char *z[NGENE][NS], *nHdiploid[NGENE], diploidlocus[NGENE];
   char cleandata[NGENE], *Indnames[NS];
   double *fpatt[NGENE], *lnpDi, *lnpGi, *T2h[2];  /* for invG prior for theta's */
   int *ncoal[2];
   double lnpSpeciesModel, *locusrate, *heredity;
   double theta_prior[2], tau_prior[3];
   double a_locusrate, a_heredity, b_heredity;
   double a_seqerr[NSPECIES][4 * 4], e_seqerr[NSPECIES][4 * 4];
   int    est_locusrate, est_heredity, nseqerr, iseqerr[NSPECIES];
}  data;

struct MCMCPARAMETERS {
   int resetFinetune, burnin, nsample, sampfreq, usedata, saveconP;
   int print, printGenetree, printlocusrate, printheredity, moveinnode, RJalgorithm;
   double finetune[7], RJfinetune[2], pSlider, ExpandRatio, ShrinkRatio;
}  mcmc;                /* control parameters */

#endif



double BFbeta = 0;


int GetOptions(char *ctlf);
int GetOptionsSimulation(char *ctlf);
int ReadSpeciesTree(FILE* fctl, char *curline);
int ReadMigrationMatrix(FILE *fctl, char *curline);
int SetupPopPopTable(int PrintTable);
int DownSptreeSetSpnames(int inode, int SetSpNames);
void GetRandomGtree(int locus);
void SimulateData(void);
int MakeDiploid(unsigned char *z[], int *ns, int ls);
int GetMem(int ipop[]);
void FreeMem(void);
int ResetSpeciesGeneTree(int locus);
int Coalescence1Pop(int ispecies);
int CoalescentMigration(void);
int DiploidResolution(FILE *fout, int locus, int ipop[]);
int ReadSeqData(char *seqfile, char *locusratef, char *heredityf, FILE*fout, char cleandata, int ipop[]);
double lnpG_S(int locus);
double lnpG_S_NoTheta(void);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
int UseLocus(int locus, int copytreeconP, int useData, int setSeqName);
int AcceptLocus(int locus, int copyconP);
int GetInitials(void);
int collectx(int mode, FILE* fout, double x[]);

int MCMC(FILE* fout);
void SwitchconPin(void);
void CopyconPin(int locus, int FromCurrent);
double UpdateGB_InternalNode(double* lnL, double* lnpG, double finetune);
double UpdateGB_SPR(double* lnL, double* lnpG, double finetune);
double UpdateTheta(double *lnpG, double finetune, double space[]);
double UpdateTau(double *lnL, double *lnpG, double finetune, double space[]);
double mixing(double* lnL, double* lnpG, double finetune, double space[]);
int    UpdateSpeciesSplit(double *lnL, double* lnpG, double space[], double PrSplit);
int    UpdateSpeciesJoin(double *lnL, double* lnpG, double space[], double PrSplit);
double UpdateLocusrateHeredity(double* lnL, double* lnpG, double finetune);
double UpdateSequenceErrors(double* lnL, double finetune, double space[]);
int    BranchWeights(double weight[]);
int    UpdateSpeciesTreeSPR(int NNIonly, double *lnL, double* lnpG, double space[]);
int    UpdateSpeciesTreeNodeSlider(double *lnL, double* lnpG, double finetune, double space[]);
int    NodeSliderScaleClade(int inode, int RDao[], double tauApath[], double factor);
int    SetSonNodeFlags(int ispecies, double tauU, char flag[]);
int    RubberProportional(int ispecies, double tauU, double tau, double taunew, double *lnproposal, double *space);

int    ReScaleSubTree(int inode, double factor);
int    NodeSlider_NotUsed(double eps);

void GraftNode(int source, int target, double age, int ipop);
int MatchGTree(void);
void copySptree(void);
void printSptreeBPP(FILE* fout);
char *printDelimitationModel(void);
double CountLHsTree2(void);
double lnpriorSpeciesModel(int ndspecies);
int PriorS_IntegerPartitions(int n, int prinT, double PriorS[]);
int InitializeDelimitationModel(int prinT, double priorPmodel[]);
int NumberDelimitationModels(int inode);
int EnumerateDelimitationModels(void);
int GetDmodelIndex(void);
int printGtree(int printBlength);
void checkGtree(void);
int ProcessGtrees(FILE* fout);
int PrintSmodel(FILE *fout, struct SMODEL *model, int printPhylogeny);
int DescriptiveStatisticsSimpleBPP(FILE *fout, char infile[], int SkipColumns);
int SpeciesTreeDelimitationModelRepresentation(struct SMODEL *model, int ndspecies);
int  SummarizeA00_DescriptiveStatisticsSimpleBPP(FILE *fout, char infile[], int SkipColumns);
/* SummarizeA01 uses CladeSupport() */
void SummarizeA10_SpeciesDelimitation(FILE* fout, char mcmcf[]);
void SummarizeA11_SpeciesTreeDelimitation(FILE* fout, char mcmcf[]);
int getab_beta(void);
double(*rndSymmetrical)(void);

extern double PjumpOptimum;

extern int noisy, IncludeNodeLabel;
char timestr[32];
double OLDAGE = 999;
int debug = 0, testlnL = 0, NPMat = 0, LASTROUND = 1;

int times = 0;

#define REALSEQUENCE
#define NODESTRUCTURE
#define BPP
#include "treesub.c"


int main(int argc, char*argv[])
{
#ifndef SIMULATION
   char ctlf[2048] = "bpp.ctl";
#else
   char ctlf[2048] = "MCcoal.ctl";
#endif
   char VerStr[64] = "Version 3.4, March 2018";
   FILE *fout;
   int i, k = 5;

   if (k == 0) { rndSymmetrical = rnduM0V1;             PjumpOptimum = 0.4; }
   if (k == 1) { rndSymmetrical = rndTriangle;          PjumpOptimum = 0.4; }
   if (k == 2) { rndSymmetrical = rndLaplace;           PjumpOptimum = 0.4; }
   if (k == 3) { rndSymmetrical = rndNormal;            PjumpOptimum = 0.4; }
   if (k == 4) { rndSymmetrical = rndBactrian;          PjumpOptimum = 0.3; }
   if (k == 5) { rndSymmetrical = rndBactrianTriangle;  PjumpOptimum = 0.3; }
   if (k == 6) { rndSymmetrical = rndBactrianLaplace;   PjumpOptimum = 0.3; }

   if (argc > 2 && !strcmp(argv[argc - 1], "--stdout-no-buf"))
      setvbuf(stdout, NULL, _IONBF, 0);

   starttimer();
   if (argc > 1) strcpy(ctlf, argv[1]);
   com.cleandata = 0;
   com.clock = 1; com.ncode = 4; com.model = 0;
   for (i = 0; i < 4; i++) com.pi[i] = 0.25;
   noisy = 3;
   com.sizespace = max2(1000000, data.ngene*stree.nspecies * 2) * sizeof(double);  /* check ipop size */
   if ((com.space = (double*)malloc(com.sizespace * sizeof(double))) == NULL)
      error2("oom space");

#ifdef SIMULATION
   printf("MCcoal in bp&p (%s)\n", VerStr);
   com.simulation = 1;
   GetOptionsSimulation(ctlf);
   SimulateData();
#else
   printf("bp&p %s\n", VerStr);
   com.simulation = 0;
   com.ngene = -1;  /* not used */
   data.ngene = 1;
   GetOptions(ctlf);

   fout = gfopen(com.outf, "w");
   fprintf(fout, "bp&p (%s) %s\n", VerStr, com.seqf);
   /* The size of ipop[] is total#sequences*sizeof(char) */

   SetMapAmbiguity(com.seqtype, 0);
   ReadSeqData(com.seqf, com.locusratef, com.heredityf, fout, com.cleandata, (int*)com.space);
   GetMem((int*)com.space);

   if (mcmc.print >= 0) {
      MCMC(fout);
   }

   if (mcmc.print && stree.analysis == A00) {
      SummarizeA00_DescriptiveStatisticsSimpleBPP(fout, com.mcmcf, 1);
   }
   else if (mcmc.print && stree.analysis == A01) {
      com.ns = stree.nspecies;
      for (k = 0; k < com.ns; k++) strcpy(com.spname[k], stree.nodes[k].name);
      CladeSupport(fout, com.mcmcf, 0, NULL, 0);
   }
   else if (mcmc.print && stree.analysis == A10) {
      SummarizeA10_SpeciesDelimitation(fout, com.mcmcf);
   }
   else if (mcmc.print && stree.analysis == A11) {
      if (mcmc.print < 0)
         PriorS_IntegerPartitions(stree.nspecies, (stree.nspecies < 20), stree.PriorSA11);
      SummarizeA11_SpeciesTreeDelimitation(fout, com.mcmcf);
   }
   if (mcmc.printGenetree && stree.nspecies > 1) ProcessGtrees(fout);

   fclose(fout);
   FreeMem();
#endif
   exit(0);
}


#ifdef SIMULATION

int GetOptionsSimulation(char *ctlf)
{
   int nopt = 18, lline = 4096, iopt, i, j, is, ierror;
   char line[4096], *pline, opt[32], *comment = "*#";
   char *optstr[] = { "seed", "noisy", "seqfile", "treefile", "Imapfile", "modelparafile", "concatfile",
                     "species&tree", "diploid", "loci&length", "migration", "sequenceerror", "alpha_locusrate",
                     "model", "Qrates", "basefreqs", "alpha_siterate", "clock" };
   char name[LSPNAME];
   int ifield[NSPECIES];
   double t = 1;
   FILE  *fctl = NULL;

   data.fix_alpha = 1;  com.alpha = 0;  com.ncatG = 5;
   for (i = 0; i < 4; i++) com.pi[i] = 0.25;
   com.ngene = 1;
   com.clock = 1;  data.a_branchrate = data.a_locusrate = 0;

   fctl = (FILE*)gfopen(ctlf, "r");
   if (fctl) {
      if (noisy) printf("\nReading options from %s..\n", ctlf);
      for (; ; ) {
         if (fgets(line, lline, fctl) == NULL)
            break;
         if (line[0] == '/' && line[1] == '/')
            break;
         for (i = 0, t = 0, pline = line; i < lline&&line[i]; i++)
            if (isalnum(line[i])) { t = 1; break; }
            else if (strchr(comment, line[i])) break;
            if (t == 0) continue;
            sscanf(line, "%s%*s%lf", opt, &t);
            if ((pline = strstr(line, "=")) == NULL)
               continue;

            for (iopt = 0; iopt < nopt; iopt++) {
               if (strncmp(opt, optstr[iopt], 8) == 0) {
                  if (noisy >= 9)
                     printf("\n%3d %15s | %-20s %6.2f", iopt + 1, optstr[iopt], opt, t);
                  switch (iopt) {
                  case (0): SetSeed((int)t, 1);                  break;
                  case (1): sscanf(pline + 1, "%d", &noisy);       break;
                  case (2): sscanf(pline + 1, "%s%d", com.seqf, &com.SeqFileFormat);    break;
                  case (3): sscanf(pline + 1, "%s", com.treef);    break;
                  case (4): sscanf(pline + 1, "%s", com.Imapf);    break;
                  case (5): sscanf(pline + 1, "%s", com.modelparaf);    break;
                  case (6): sscanf(pline + 1, "%s", com.concatf);  break;
                  case (7):
                     if ((stree.nspecies = com.ns = (int)t) > NSPECIES)
                        error2("raise NSPECIES in bpp.c & recompile?");
                     ReadSpeciesTree(fctl, pline + 1);
                     break;
                  case (8):  /* data.diploid[] */
                     if (stree.nspecies <= 0)
                        error2("diploid should be after the species&tree block.");
                     splitline(pline + 1, stree.nspecies, ifield);
                     for (is = 0; is < stree.nspecies; is++)
                        sscanf(pline + 1 + ifield[is], "%d", &data.diploid[is]);
                     printf("\ndiploid flags for %d species: ", stree.nspecies);
                     for (is = 0; is < stree.nspecies; is++)
                        printf(" %d", data.diploid[is]);
                     break;
                  case (9): sscanf(pline + 1, "%d%d", &com.ngene, &com.ls);  break;
                  case (10):
                     ReadMigrationMatrix(fctl, pline);
                     break;
                  case (11):               /* sequencing errors */
                     data.nseqerr = 0;
                     if (sscanf(pline + 1, "%d", &data.nseqerr) != 1)
                        error2("error in the sequenceerror line of the control file.");
                     for (ierror = 0; ierror < data.nseqerr; ierror++) {
                        fscanf(fctl, "%s", name);
                        for (is = 0; is < stree.nspecies; is++)
                           if (strcmp(name, stree.nodes[is].name) == 0) break;
                        if (is == stree.nspecies) error2("expecting a species name");
                        data.iseqerr[is] = 1;
                        for (i = 0; i < 16; i++)
                           fscanf(fctl, "%lf", &data.e_seqerr[is][i]);
                        for (i = 0; i < 4; i++)
                           abyx(1 / sum(data.e_seqerr[is] + i * 4, 4), data.e_seqerr[is] + i * 4, 4);
                        printf("\nsequence errors in %s:", stree.nodes[is].name);
                        matout(F0, data.e_seqerr[is], 4, 4);
                        for (i = 0; i < 4; i++) for (j = 1; j < 4; j++)
                           data.e_seqerr[is][i * 4 + j] += data.e_seqerr[is][i * 4 + j - 1];
                        matout(F0, data.e_seqerr[is], 4, 4);
                     }
                     break;
                  case (12): sscanf(pline + 1, "%lf", &data.a_locusrate);  break;  /* gamma rates for loci */
                  case (13): sscanf(pline + 1, "%d", &com.model);  break;  /* mutation model */
                  case (14): sscanf(pline + 1, "%d%lf%lf%lf%lf%lf%lf", &data.fix_Qrates,
                     &data.Qrates[0], &data.Qrates[1], &data.Qrates[2], &data.Qrates[3], &data.Qrates[4], &data.Qrates[5]);  /* fix_Qrates & a_Qrates[] */
                     /* Qrates is reassigned later to Qrates_Dir if model==7 */
                     break;
                  case (15): sscanf(pline + 1, "%d%lf%lf%lf%lf", &data.fix_pi, &data.pi_Dir[0], &data.pi_Dir[1], &data.pi_Dir[2], &data.pi_Dir[3]);  /* base freqs */
                     break;
                  case (16): /* alpha for ASRV */
                     sscanf(pline + 1, "%d%lf%d", &data.fix_alpha, &com.alpha, &com.ncatG);
                     if (data.fix_alpha == 0) sscanf(pline + 1, "%d%lf%lf%d", &data.fix_alpha, &data.a_siterate[0], &data.a_siterate[1], &com.ncatG);
                     break;
                  case (17): sscanf(pline + 1, "%d%lf", &com.clock, &data.a_branchrate);  break;  /* relaxed clock */
                  }
                  break;
               }
            }
            if (iopt == nopt)
            {
               printf("\noption %s in %s\n", opt, ctlf);  exit(-1);
            }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");

   /* Diploid sequences are generated by simulating haploid sequences then makediploid. */
   for (i = 0, com.ns = 0; i < stree.nspecies; i++) {
      if (data.diploid[i])  stree.nseqsp[i] *= 2;
      com.ns += stree.nseqsp[i];
   }
   if (com.ns > NS)  error2("raise NS in bpp.c & recompile");

   if (com.model != 0 && com.model != 7) error2("mutation model not allowed");
   if (com.model == 0) {     /* JC model */
      printf("\nJC model assumed...\n");
      data.fix_Qrates = 1;
      for (i = 0; i < 6; i++) data.Qrates[i] = data.Qrates_Dir[i] = 0;
      data.fix_pi = 1;
   }
   else {                    /* GTR and simpler */
      if (com.modelparaf[0] == '\0') strcpy(com.modelparaf, "modelparas.txt");
      if (data.fix_Qrates) {
         for (i = 0; i < 6; i++)  data.Qrates_Dir[i] = 0;
      }
      else {
         for (i = 0; i < 6; i++)  data.Qrates_Dir[i] = data.Qrates[i];
         for (i = 0; i < 6; i++)  data.Qrates[i] = 0;
      }
      if (data.fix_Qrates) {
         printf("\nGTR Qrates are fixed: ");
         for (i = 0; i < 6; i++)  data.Qrates[i] /= data.Qrates[5];
         for (i = 0; i < 6; i++)  printf("%9.6f", data.Qrates[i]);
         printf("\n");
      }
      else {
         printf("\nGTR Qrates are generated from Dirichlet: ");
         for (i = 0; i < 6; i++)  printf(" %9.6f", data.Qrates_Dir[i]);
         printf("\n");
         if (data.Qrates_Dir[5] <= 0) error2("GTR Qrates..");
      }
      if (data.fix_pi) {
         for (i = 0; i < 4; i++)  com.pi[i] = data.pi_Dir[i];
         for (i = 0; i < 4; i++)  data.pi_Dir[i] = 0;
      }
      else {
         printf("\nbase frequencies are generated from Dirichlet: ");
         for (i = 0; i < 4; i++)  printf(" %9.6f", data.pi_Dir[i]);
         printf("\n");
         if (data.pi_Dir[3] <= 0) error2("pi Dir parameters..");
      }
   }
   if (data.fix_alpha == 0) {
      if (data.a_siterate[0] < 0 || data.a_siterate[1] < 0)  error2("alpha_siterate < 0");
      com.alpha = 0;
   }
   else {
      if (com.alpha) {
         if (com.ncatG)  printf("site rates are from discrete gamma(%9.6f), K = %d.", com.alpha, com.ncatG);
         else            printf("site rates are from continuous gamma(%9.6f), K = inf.", com.alpha);
      }
      data.a_siterate[0] = data.a_siterate[1] = 0;
   }

   if (com.clock < 1 || com.clock > 3)
      error2("clock has to be 1 or 2 or 3.");
   return(0);
}

#else

int GetOptions(char *ctlf)
{
   int nopt = 26, lline = 4096, iopt, i, is, ierror;
   char line[4096], *pline, opt[32], *comment = "*#", *seqerrstr = "0EF";
   char *optstr[] = { "seed", "noisy", "seqfile", "Imapfile", "outfile", "mcmcfile", "checkpoint", "BayesFactorBeta",
      "speciesdelimitation", "speciestree", "speciesmodelprior", "species&tree", "diploid",
      "usedata", "nloci", "cleandata", "thetaprior", "tauprior", "locusrate", "heredity",
      "sequenceerror", "finetune", "print", "burnin", "sampfreq", "nsample" };
   char name[LSPNAME], ch;
   int ifield[NSPECIES];
   double t = 1, *eps = mcmc.finetune;
   FILE  *fctl = gfopen(ctlf, "r");

   strcpy(com.checkpointf, "bpp.bkpt");
   stree.SpeciesModelPrior = p1UniformRootedTree;
   mcmc.pSlider = -1;  mcmc.ExpandRatio = -1;  mcmc.ShrinkRatio = -1;
   data.lnpSpeciesModel = 0;
   mcmc.printGenetree = 0;
   for (is = 0; is < NSPECIES; is++)
      for (i = 0; i < 16; i++)
         data.a_seqerr[is][i] = data.e_seqerr[is][i] = 0;

   if (fctl) {
      if (noisy) printf("\nReading options from %s..\n", ctlf);
      for (; ;) {
         if (fgets(line, lline, fctl) == NULL)
            break;
         if (line[0] == '/' && line[1] == '/')
            break;
         for (i = 0, t = 0, pline = line; i < lline&&line[i]; i++)
            if (isalnum(line[i])) { t = 1; break; }
            else if (strchr(comment, line[i])) break;
            if (t == 0) continue;
            sscanf(line, "%s%*s%lf", opt, &t);
            if ((pline = strstr(line, "=")) == NULL) error2("option file.\nExpecting '=' ");

            for (iopt = 0; iopt < nopt; iopt++) {
               if (strncmp(opt, optstr[iopt], 12) == 0) {
                  if (noisy >= 9)
                     printf("\n%3d %15s | %-20s %6.2f", iopt + 1, optstr[iopt], opt, t);
                  switch (iopt) {
                  case (0): SetSeed((int)t, 1);                 break;
                  case (1): sscanf(pline + 1, "%d", &noisy);      break;
                  case (2): sscanf(pline + 1, "%s", com.seqf);    break;
                  case (3): sscanf(pline + 1, "%s", com.Imapf);   break;
                  case (4): sscanf(pline + 1, "%s", com.outf);    break;
                  case (5): sscanf(pline + 1, "%s", com.mcmcf);   break;
                  case (6): sscanf(pline + 1, "%d", &com.checkpoint); break;
                  case (7): sscanf(pline + 1, "%lf", &BFbeta);        break; /* beta for marginal likelihood */
                  case (8):
                     sscanf(pline + 1, "%d%d%lf%lf", &stree.speciesdelimitation, &mcmc.RJalgorithm, &mcmc.RJfinetune[0], &mcmc.RJfinetune[1]);
                     if (stree.speciesdelimitation == 1) {
                        if (mcmc.RJalgorithm != 0 && mcmc.RJalgorithm != 1)  error2("RJalgorithm should be 0 or 1.");
                        if (mcmc.RJfinetune[0] <= 0 || (mcmc.RJalgorithm == 1 && mcmc.RJfinetune[1] <= 0))
                           error2("RJfinetune <= 0.\nError on the line speciesdelimitation?");
                     }
                     break;
                  case (9):
                     sscanf(pline + 1, "%d%lf%lf%lf", &stree.speciestree, &mcmc.pSlider, &mcmc.ExpandRatio, &mcmc.ShrinkRatio);
                     break;
                  case (10): sscanf(pline + 1, "%d", &stree.SpeciesModelPrior);   break;
                     break;
                  case (11):
                     if ((stree.nspecies = com.ns = (int)t) > NSPECIES)
                        error2("raise NSPECIES in bpp.c & recompile?");
                     ReadSpeciesTree(fctl, pline + 1);
                     break;
                  case (12):  /* data.diploid[] */
                     if (stree.nspecies <= 0)
                        error2("diploid should be after the species&tree block.");
                     splitline(pline + 1, stree.nspecies, ifield);
                     for (is = 0; is < stree.nspecies; is++)
                        sscanf(pline + 1 + ifield[is], "%d", &data.diploid[is]);
                     printf("\ndiploid flags for %d species: ", stree.nspecies);
                     for (is = 0; is < stree.nspecies; is++)
                        printf(" %d", data.diploid[is]);
                     break;
                  case (13): mcmc.usedata = (int)t;    break;
                  case (14): data.ngene = (int)t;      break;
                  case (15): com.cleandata = (char)t;  break;
                  case (16): sscanf(pline + 1, "%lf%lf %c", &data.theta_prior[0], &data.theta_prior[1], &ch);
                     stree.NoTheta = (toupper(ch) == 'E' ? 0 : 1);
                     break;
                  case (17): sscanf(pline + 1, "%lf%lf%lf", &data.tau_prior[0], &data.tau_prior[1], &data.tau_prior[2]);
                     break;
                  case (18):               /* locus rate */
                     sscanf(pline + 1, "%d%lf", &data.est_locusrate, &data.a_locusrate);
                     if (data.est_locusrate == 1) {
                        printf("\nRates vary among loci, with ri/L ~ Dirichlet(%.2f)\n", data.a_locusrate);
                        if (data.a_locusrate < .001) error2("alpha very small?");
                     }
                     else if (data.est_locusrate == 2) {
                        /* check this...? */
                        sscanf(pline + 1, "%d%s", &i, com.locusratef);
                     }
                     break;
                  case (19):
                     sscanf(pline + 1, "%d%lf%lf", &data.est_heredity, &data.a_heredity, &data.b_heredity);
                     if (data.est_heredity == 1) {
                        printf("\n\nheredity multiplier ~ gamma(%.2f, %.2f)\n", data.a_heredity, data.b_heredity);
                        if (data.a_heredity < .001 || data.b_heredity < .001) error2("poor prior?");
                     }
                     else if (data.est_heredity == 2) {
                        /* check this...? */
                        sscanf(pline + 1, "%d%s", &i, com.heredityf);
                     }
                     break;

                  case (20):               /* sequencing errors */
#if (defined(USE_AVX) || defined(USE_SSE))
                     fprintf(stderr, "Accomodating sequencing errors in likelihood "
                        "function is not yet available for SSE/AVX version.\n");
                     exit(1);
#endif
                     /* check this...? */
                     data.nseqerr = 0;
                     if (sscanf(pline + 1, "%d", &data.nseqerr) != 1)
                        error2("error in the sequenceerror line of the control file.");
                     for (ierror = 0; ierror < data.nseqerr; ierror++) {
                        fscanf(fctl, "%s", name);
                        for (is = 0; is < stree.nspecies; is++)
                           if (strcmp(name, stree.nodes[is].name) == 0) break;
                        if (is == stree.nspecies) error2("expecting a species name");
                        data.iseqerr[is] = 1;
                        for (i = 0; i < 16; i++)
                           fscanf(fctl, "%lf", &data.a_seqerr[is][i]);
                        printf("\nalpha matrix for sequence errors in %s:", stree.nodes[is].name);
                        matout(F0, data.a_seqerr[is], 4, 4);
                     }
                     break;

                  case (21):
                     sscanf(pline + 1, "%d:%lf%lf%lf%lf%lf%lf%lf", &mcmc.resetFinetune, eps, eps + 1, eps + 2, eps + 3, eps + 4, eps + 5, eps + 6);
                     break;
                  case (22):
                     sscanf(pline + 1, "%d%d%d%d", &mcmc.print, &mcmc.printlocusrate, &mcmc.printheredity, &mcmc.printGenetree);
                     break;
                  case (23): mcmc.burnin = (int)t;    break;
                  case (24): mcmc.sampfreq = (int)t;  break;
                  case (25): mcmc.nsample = (int)t;   break;
                  }
                  break;
               }
            }
            if (iopt == nopt)
            {
               printf("\noption %s in %s\n", opt, ctlf);  exit(-1);
            }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");


   t = data.theta_prior[1] / (data.theta_prior[0] - 1);
   printf("\nprior: theta ~ IG(%.6f, %.6f), with mean %.6f\n", data.theta_prior[0], data.theta_prior[1], t);
   if (t > 1) printf("\n\aare you sure about this prior? ");
   if (stree.nspecies > 1) {
      t = data.tau_prior[1] / (data.tau_prior[0] - 1);
      printf("\nprior:  tau ~ IG(%.6f, %.6f), with mean %.6f\n", data.tau_prior[0], data.tau_prior[1], t);
      if (t > 1) printf("\aare you sure about this prior?");
   }
   if (data.theta_prior[0] <= 2 || (stree.nspecies > 1 && data.tau_prior[0] <= 2))
      puts("theta or tau prior, use alpha > 2");
   stree.analysis = (stree.speciesdelimitation << 1) + stree.speciestree;
   if (stree.analysis < 0 || stree.analysis>3) error2("speciesdelimitation or speciestree");
   if (stree.SpeciesModelPrior > 3) error2("SpeciesModelPrior unheard of");
   if (stree.SpeciesModelPrior >= 2 && stree.analysis != A11)
      error2("SpeciesModelPrior?");

   if (stree.speciestree && mcmc.pSlider < 0) {
      if (mcmc.usedata) { mcmc.pSlider = 0.4; mcmc.ExpandRatio = 0.15; mcmc.ShrinkRatio = 0.1; }
      else { mcmc.pSlider = 0.4; mcmc.ExpandRatio = 0.40; mcmc.ShrinkRatio = 0.5; }
   }
   if (data.est_locusrate != 1) mcmc.printlocusrate = 0;
   if (data.est_heredity != 1)  mcmc.printheredity = 0;
   if (data.ngene < 1) error2("nloci = 0?");

   /* diploid sequences. */
   for (is = 0, t=0; is < stree.nspecies; is++) t += data.diploid[is];
   if (com.cleandata && t > 0) error2("use cleandata = 0 for the diploid option");
   for (i = 0, data.maxns = 0; i < stree.nspecies; i++) {
      if (data.diploid[i])  stree.nseqsp[i] *= 2;
      data.maxns += stree.nseqsp[i];
   }
   for (i = 0, stree.npop = 0; i < stree.nspecies * 2 - 1; i++)
      if (i >= stree.nspecies || stree.nseqsp[i] > 1)
         stree.pops[stree.npop++] = i;
   printf("\n%d populations for coalescent:\n ", stree.npop);
   for (i = 0; i < stree.npop; i++)
      if (stree.pops[i] < stree.nspecies || stree.nspecies <= NPopLongNames)
         printf(" theta_%d%s", stree.pops[i] + 1, stree.nodes[stree.pops[i]].name);
      else
         printf(" theta_%d", stree.pops[i] + 1);
   if (stree.nspecies > 1) {
      printf("\n%d species divergence times in the order:\n ", stree.nspecies - 1);
      for (i = stree.nspecies; i < stree.nspecies * 2 - 1; i++)
         if (stree.nspecies <= 10) printf(" tau_%d%s", i + 1, stree.nodes[i].name);
         else                    printf(" tau_%d", i + 1);
   }
   printf("\n");

   if (data.maxns > NS)  error2("raise NS in bpp.c & recompile");
   if (data.maxns < 2)  error2("<2 seqs?  Too simple to do anything about.");

   if (mcmc.print != -1 && BFbeta && mcmc.usedata == 0)
      error2("marginal like for prior with usedata =0?");
   else if (BFbeta == 0)
      BFbeta = 1;
   return(0);
}


#endif


int ReadSpeciesTree(FILE* fctl, char *curline)
{
   /* It reads the species tree and initializes stree.
      Note that (stree.nseqsp[j] > 1) adds a population.
   */
   char line[16000] = { '\0' };
   int i, j, lline = 16000, ifields[1000];

   nodes = nodes_t;  /* species tree is read into the temp space nodes_t */
   for (i = 0; i < NS; i++) {  /* for both species tree & gene trees.  This wastes space. */
      if (com.spname[i]) free(com.spname[i]);
      com.spname[i] = (char*)malloc((LSPNAME + 1) * sizeof(char));
   }
   splitline(curline, 1 + stree.nspecies, ifields);
   for (i = 0; i < stree.nspecies; i++) {
      sscanf(curline + ifields[i + 1], "%s", com.spname[i]);
      strcpy(stree.nodes[i].name, com.spname[i]);
   }
   for (i = 0; i < stree.nspecies; i++) {
      fscanf(fctl, "%d", &stree.nseqsp[i]);
      if (stree.nseqsp[i] < 1)  /* 0 would be ok with migration. */
         puts("Could we please have at least 1 seq from each species?");
   }
   printf("\n%d species: ", stree.nspecies);
   for (i = 0; i < stree.nspecies; i++)
      printf(" %s (%d)", com.spname[i], stree.nseqsp[i]);
   putchar('\n');

   fgets(line, lline, fctl);
   if (stree.nspecies == 1) {
      stree.root = stree.nodes[0].nson = 0;  stree.nnode = 1; stree.nbranch = 0;
      stree.nodes[0].age = 0;
      stree.npop = 1;
      stree.pops[0] = 0;
      stree.pptable[0][0] = (char)1;
#ifdef SIMULATION
      fscanf(fctl, "%lf", &stree.nodes[0].theta);
      printf("theta = %9.4f", stree.nodes[0].theta);
#endif
   }
   else {
      ReadTreeN(fctl, &i, 0, 1);
      OutTreeN(F0, 1, com.simulation);
      FPN(F0);
      if (!com.simulation && stree.speciesdelimitation == 0 && stree.speciestree == 0 && stree.nspecies > NPopLongNames) {
         OutTreeN(F0, 1, PrNodeNum); FPN(F0);
      }

      /* copy into stree */
      stree.nnode = tree.nnode;
      stree.nbranch = tree.nbranch;
      stree.root = tree.root;
      for (i = 0; i < stree.nnode; i++) {
         stree.nodes[i].father = nodes[i].father;
         stree.nodes[i].nson = nodes[i].nson;
         for (j = 0; j < stree.nodes[i].nson; j++)
            stree.nodes[i].sons[j] = nodes[i].sons[j];
         stree.nodes[i].age = nodes[i].age = nodes[i].branch;
         stree.nodes[i].theta = nodes[i].label;
      }
      for (i = 0; i < stree.nnode; i++) {
         if (i != stree.root && stree.nodes[stree.nodes[i].father].age < stree.nodes[i].age) {
            printf("\nstree node ages wrong: %.6f < %.6f.", stree.nodes[stree.nodes[i].father].age, stree.nodes[i].age);
            error2("control file error.");
         }
      }
      DownSptreeSetSpnames(stree.root, 1);
      SetupPopPopTable(1);
   }
   stree.M = (double*)malloc(stree.nnode*stree.nnode * sizeof(double));
   if (stree.M == NULL) error2("oom M");
   memset(stree.M, 0, stree.nnode*stree.nnode * sizeof(double));

   return(0);
}


int ReadMigrationMatrix(FILE *fctl, char *pline)
{
   /* Rates for impossible ancestor-descendent migrations are set to 0.
      This does not check node ages to make sure that the populations live
      in the same epoch.
   */
   int s21 = stree.nspecies * 2 - 1, i, j, needtheta, status = 0, nMresets = 0;
   double Agei, Agej, Adadi, Adadj;
   char spname[40];

   sscanf(pline + 1, "%d", &stree.migration);
   if (stree.migration == 0)
      return(0);

   if (stree.migration != s21)
      error2("migration = #species*2 - 1?");
   for (i = 0; i < s21; i++) {
      fscanf(fctl, "%s", spname);
      if (strcmp(stree.nodes[i].name, spname))  error2("Migr: spname mismatch");
   }
   for (i = 0; i < s21; i++) {
      fscanf(fctl, "%s", spname);
      if (strcmp(stree.nodes[i].name, spname))
         error2("error in migration matrix: species name mismatch");
      for (j = 0; j < s21; j++)
         fscanf(fctl, "%lf", &stree.M[i*s21 + j]);
   }

   /* set rates for impossible migrations to 0 or -1 */
   for (i = 0; i < s21; i++) {
      Agei = stree.nodes[i].age;
      Adadi = stree.nodes[stree.nodes[i].father].age;
      for (j = i; j < s21; j++) {
         Agej = stree.nodes[j].age;
         Adadj = stree.nodes[stree.nodes[j].father].age;
         if (i == j) {
            if (stree.M[i*s21 + j] > 0) nMresets++;
            stree.M[i*s21 + j] = stree.M[j*s21 + i] = 0;
         }
         else if (stree.pptable[i][j] || stree.pptable[j][i] || Agei > Adadj || Agej > Adadi) {
            if (stree.M[i*s21 + j] > 0) nMresets++;
            if (stree.M[j*s21 + i] > 0) nMresets++;
            if (stree.M[i*s21 + j] > 0 || stree.M[j*s21 + i] > 0)
               printf("\nmigration between %s and %s impossible...", stree.nodes[i].name, stree.nodes[j].name);
            stree.M[i*s21 + j] = stree.M[j*s21 + i] = -1;
         }
      }
   }
   if (nMresets) printf("\n%d migration rates are reset in the migration matrix", nMresets);
   printf("\nmigration matrix\n%10s", "");
   for (j = 0; j < s21; j++)
      printf(" %9s", stree.nodes[j].name);
   for (i = 0, FPN(F0); i < s21; i++, FPN(F0)) {
      printf("%-10s", stree.nodes[i].name);
      for (j = 0; j < s21; j++) {
         printf(" %9.4f", stree.M[i*s21 + j]);
      }
   }

   for (i = 0; i < s21; i++)
      for (j = 0; j < s21; j++)
         if (stree.M[i*s21 + j] < 0)  stree.M[i*s21 + j] = 0;
   for (i = 0; i < stree.nspecies; i++) {
      for (j = 0, needtheta = 0; j < s21; j++)
         if (stree.M[j*s21 + i]) needtheta = 1;
      if (needtheta && stree.nodes[i].theta <= 0) {
         printf("theta for %s should be > 0\n", stree.nodes[i].name);
         status = -1;
      }
   }
   if (status) error2("fix the migration matrix in the control file.");
   return(0);
}


int DownSptreeSetSpnames(int inode, int SetSpNames)
{
   /* This traverse down the species tree to set stree.nodes[].name.
   */
   int k, ison;

   if (inode < stree.nspecies)
      error2("should not be here?");

   if (SetSpNames && stree.nspecies > NPopLongNames)
      sprintf(stree.nodes[inode].name, "%d", inode + 1);

   for (k = 0; k < stree.nodes[inode].nson; k++) {
      ison = stree.nodes[inode].sons[k];
      if (stree.nodes[ison].nson)
         DownSptreeSetSpnames(ison, SetSpNames);
      if (SetSpNames && stree.nspecies <= NPopLongNames) {
         if (strlen(stree.nodes[inode].name) + strlen(stree.nodes[ison].name) > 2 * LSPNAME - 1)
            error2("we are in trouble.  Increase LSPNAME?");
         strcat(stree.nodes[inode].name, stree.nodes[ison].name);
      }
   }
   return(0);
}

int SetupPopPopTable(int PrintTable)
{
   int i, j, s = stree.nspecies, root = stree.root;

   for (i = 0; i < 2 * s - 1; i++) for (j = 0; j < 2 * s - 1; j++)
      stree.pptable[i][j] = (char)0;
   for (i = 0; i < 2 * s - 1; i++) stree.pptable[i][i] = stree.pptable[i][root] = (char)1;
   for (i = 0; i < 2 * s - 1; i++) {
      for (j = i; j != root; j = stree.nodes[j].father)
         stree.pptable[i][j] = (char)1;
   }

   if (PrintTable) {
      printf("\npop by pop table showing node numbers in species tree\n\n%22s", " ");
      for (i = 0; i < 2 * s - 1; i++)
         printf(" %2d", i + 1);
      FPN(F0);
      for (i = 0; i < 2 * s - 1; i++, FPN(F0)) {
         printf("species %2d %-10s ", i + 1, stree.nodes[i].name);
         for (j = 0; j < 2 * s - 1; j++)
            printf(" %2d", (int)stree.pptable[i][j]);
         if (i < stree.nspecies && stree.nodes[i].age)
            printf("\n\a  <-- age>0??\n");
#ifdef SIMULATION
         printf("  tau =%7.4f theta =%7.4f  ", stree.nodes[i].age, stree.nodes[i].theta);
         if ((i >= s || stree.nseqsp[i] > 1) && stree.nodes[i].theta <= 0)
            printf("  this theta must be > 0!!");
#endif
      }
   }
   return(0);
}


/* used by ResetSpeciesGeneTree() and Coalescence1Pop().  Those globals appear
   necessary as Coalescence1Pop is recursive.
*/
static int cNodeNumber, nin_G[NSPECIES * 2 - 1], ins_G[NSPECIES * 2 - 1][NS];
double mtMRCA = 0, mM[(NSPECIES * 2 - 1)*(NSPECIES * 2 - 1)];


void GetRandomGtree(int locus)
{
   /* This generates a random gene tree in nodes[].  This does not initialize ipop.
      This is used by the simulation program, and also by the inference program to
      generate the starting gene tree with coalescent times.
   */
   int i;

   ResetSpeciesGeneTree(locus);
   cNodeNumber = com.ns;

   if (0 && !stree.migration) {
      Coalescence1Pop(stree.root);
      /* NodeToBranch(); */
      tree.root = cNodeNumber - 1;  /* cNodeNumber = com.ns*2 - 1.  root is last node */
   }
   else {
      CoalescentMigration();
   }

   nodes[tree.root].branch = 0;
   nodes[tree.root].father = -1;
   tree.nnode = com.ns * 2 - 1;
   for (i = 0; i < tree.nnode; i++)
      if (i != tree.root)
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;

   mtMRCA += nodes[tree.root].age;
   if (debug == 9) {
      FPN(F0); OutTreeN(F0, 1, 1); FPN(F0); FPN(F0);
   }
}

int ResetSpeciesGeneTree(int locus)
{
   /* This is called by GetRandomGtree() to reset the species tree before calling Coalescence1Pop().
      It initializes nin_G & ins_G[].
      Also resets cNodeNumber for constructing the gene tree.
      gnodes[][].ipop for tips in the gene tree have been set before this routine already.
   */
   int IS = '_', is, j, tip = 0, resetName = 0;

   for (is = 0; is < 2 * stree.nspecies - 1; is++) nin_G[is] = 0;

#ifdef SIMULATION
   for (is = 0; is < stree.nspecies; is++)
      nin_G[is] = stree.nseqsp[is];
   for (is = 0, tip = 0; is < stree.nspecies; is++) {
      for (j = 0; j < stree.nseqsp[is]; j++, tip++)
         ins_G[is][j] = tip;
   }
#else
   com.ns = data.ns[locus];
   nodes = gnodes[locus];
   for (is = 0; is < stree.nspecies; is++)
      stree.nseqsp[is] = data.nseqsp[locus][is];
   for (j = 0; j < com.ns; j++) {
      is = gnodes[locus][j].ipop;
      ins_G[is][nin_G[is]++] = j;
      if (!resetName) continue;
      if (stree.nseqsp[is] > 1)
         sprintf(com.spname[j], "%s%d", stree.nodes[is].name, nin_G[is]);
      else
         sprintf(com.spname[j], "%s", stree.nodes[is].name);
   }
#endif

   for (j = 0; j < com.ns; j++) ClearNode(j);
   return(0);
}


int Coalescence1Pop(int ispecies)
{
   /* This simulates the coalescent process in population or species ispecies.
      It generates the random genealogy tree (possibly consisting of several
      disconnected subtrees) with waiting times.
      t: waiting time; T: node age
      This is used by GetRandomGtree().
      nin_G[] and ins_G[] are set in ResetSpeciesGeneTree().
   */
   int j, k, k1, k2, father = stree.nodes[ispecies].father;
   double t, T;

   for (k = 0; k < stree.nodes[ispecies].nson; k++)
      Coalescence1Pop(stree.nodes[ispecies].sons[k]);
   if (stree.nodes[ispecies].theta > 0) {
      T = stree.nodes[ispecies].age;
      if (nin_G[ispecies] > 1 && stree.nodes[ispecies].theta <= 0) {
         printf("theta for pop %s is %.6f", stree.nodes[ispecies].name, stree.nodes[ispecies].theta);
         printf(" theta <= 0.\n");
      }
      for (; nin_G[ispecies] > 1; nin_G[ispecies]--, cNodeNumber++) {
         j = nin_G[ispecies];  /* # of lineages */
         t = rndexp(stree.nodes[ispecies].theta / (j*(j - 1.)));
         T += t;
         if (ispecies != stree.root && T > stree.nodes[father].age) break;

         /* k1 & k2 are lineages to coalesce */
         k = (int)(j*rndu());
         k1 = ins_G[ispecies][k]; ins_G[ispecies][k] = ins_G[ispecies][j - 1];
         k = (int)((j - 1)*rndu());
         k2 = ins_G[ispecies][k]; ins_G[ispecies][k] = cNodeNumber;

         nodes[cNodeNumber].nson = 2;
         nodes[cNodeNumber].ipop = ispecies;
         nodes[cNodeNumber].sons[0] = k1;
         nodes[cNodeNumber].sons[1] = k2;
         nodes[k1].father = nodes[k2].father = cNodeNumber;
         nodes[cNodeNumber].age = T;
      }
   }
   /* outgoing lineages added to father's list */
   if (ispecies == stree.root && nin_G[ispecies] > 1)
      error2("nin_Gtree > 1 at root");
   if (ispecies != stree.root) {
      for (k = 0; k < nin_G[ispecies]; k++)
         ins_G[father][nin_G[father]++] = ins_G[ispecies][k];
   }
   return (0);
}


int CoalescentMigration(void)
{
   /* This simulates a gene tree with both coalescent and migration events.
      nspecies epochs: the coalescent and migration rates are used to generate
      exponential waiting times.
      The routine works also when there is no migration.
      The algorithm keeps track of # of populations (npop), updates the list of pops in ipop[],
      the number of individuals in each pop nin_G[], and the individuals in each pop ins_G[][].
      Some pops may be empty, with nin_G[] = 0.  Ci[] holds coalescent rate in pop i and Mi[]
      holds migration rate to pop i, while C & M are total coalescent and migration rates.
      ipop[] keeps the list of pops during the current epoch, and is updated when we move
      to the next epoch.
      ipopE holds the order of the epochs.
      n is the number of lineages ancestral to the sample (current sample size).
      Tmax marks the end of the current epoch.
   */
   int n = com.ns, is, i, j, k, k1, k2, ipopE[NSPECIES - 1]; /* ipop for each epoch node */
   int npop = stree.nspecies, ipop[NSPECIES], *sons, s21 = stree.nspecies * 2 - 1;
   double r, T, Tmax, C, Ci[NSPECIES], M = 0, Mi[NSPECIES];
   double y, ages[NSPECIES - 1], tmp1[NSPECIES];

   /* sort node ages in species tree to work out each epoch. */
   for (i = 0; i < stree.nspecies - 1; i++)
      tmp1[i] = stree.nodes[stree.nspecies + i].age;
   /* The sorting is in increasing order, with ties broken in the original order.
      This way, the algorithm runs when nodes on the species tree are collapsed. */
   indexing(tmp1, stree.nspecies - 1, ipopE, 1, (int*)ages);

   if (debug == 9) {
      printf("\nages in stree: ");
      for (i = stree.nspecies; i < stree.nspecies * 2 - 1; i++)
         printf(" %2d: %7.5f ", i, stree.nodes[i].age);
      printf("\nafter ordering: ");
      for (i = stree.nspecies - 1 - 1; i >= 0; i--) {
         is = stree.nspecies + ipopE[i];
         printf(" %2d: %7.5f ", is, stree.nodes[is].age);
      }
      printf("\n\n");
   }

   /* Initially the tips on the species tree are in the list ipop[] */
   for (i = 0; i < npop; i++)  /* populations */
      ipop[i] = i;

   for (T = 0; ; npop--) {   /*  # of epochs */
      if (npop == 1) {
         Tmax = -1;
      }
      else { /* is: the species in the species tree at the end of this epoch */
         is = stree.nspecies + ipopE[npop - 1 - 1];
         Tmax = stree.nodes[is].age;
      }

      for (; ; ) {
         if (Tmax == 0) break;
         /* calculate poisson rates: Ci[i] is coalescent rate in pop i */
         for (i = 0, C = 0; i < npop; i++) {
            Ci[i] = 0;
            if (nin_G[i] >= 2) {
               if (stree.nodes[ipop[i]].theta <= 0)
                  printf("theta_%s = %.6f <= 0!\n", stree.nodes[ipop[i]].name, stree.nodes[ipop[i]].theta);
               C += Ci[i] = nin_G[i] * (nin_G[i] - 1) / 2 * 2.0 / stree.nodes[ipop[i]].theta;
            }
         }
         if (stree.migration)
            /* Mi[i] is the migration rate to pop i */
            for (i = 0, M = 0; i < npop; i++) {
               for (j = 0, Mi[i] = 0; j < npop; j++)
                  Mi[i] += nin_G[i] * stree.M[ipop[j] * s21 + ipop[i]] / stree.nodes[ipop[i]].theta * 4;
               M += Mi[i];
            }

         if (debug == 9) {
            printf("S%d (%5.3f) %d (", is, Tmax, npop);
            for (i = 0; i < npop; i++) printf(" %d", nin_G[i]);
            printf(" ) rate CM %6.1f %6.1f ", C, M);
         }
         if (C + M < 1e-300) {             /* move to next epoch */
            if (debug == 9) FPN(F0);
            break;
         }
         T += rndexp(1 / (C + M));
         if (debug == 9) printf(" T %6.3f ", T);
         if (T > Tmax && Tmax != -1) { /* move to next epoch */
            if (debug == 9) FPN(F0);
            break;
         }
         r = rndu()*(C + M);
         if (r < C) {  /* coalescent in pop i of lineages k1 & k2 */
            /* r ~ U(0, C) */
            for (i = 0, y = 0; i < npop - 1; i++)
               if (r < (y += Ci[i]))
                  break;

            /* k1 & k2 (with k1 < k2) are the two lineages to coalesce */
            j = (int)(nin_G[i] * (nin_G[i] - 1)*rndu());
            k1 = j / (nin_G[i] - 1);
            k2 = j % (nin_G[i] - 1);
            if (k2 >= k1)
               k2++;
            else {   /* swap */
               j = k1; k1 = k2; k2 = j;
            }

            nodes[cNodeNumber].nson = 2;
            nodes[cNodeNumber].ipop = ipop[i];
            nodes[cNodeNumber].sons[0] = ins_G[i][k1];
            nodes[cNodeNumber].sons[1] = ins_G[i][k2];
            nodes[ins_G[i][k1]].father = nodes[ins_G[i][k2]].father = cNodeNumber;
            nodes[cNodeNumber].age = T;

            /* In pop i, replace k1 by new node, remove k2 */
            ins_G[i][k1] = cNodeNumber++;
            --nin_G[i];
            if (k2 != nin_G[i])
               ins_G[i][k2] = ins_G[i][nin_G[i]];

            if (debug == 9)
               printf("C: %s   ", stree.nodes[ipop[i]].name);

            if (--n == 1)
               break;      /* last coalescent */
         }
         else {            /* migration of lineage k from pop j into pop i */
            r -= C;        /* 0 < r < M */
            for (i = 0, y = 0; i < npop; i++)
               if (r < (y += Mi[i])) break;
            if (i == npop) error2("u01 = 1!");
            y -= Mi[i];

            for (j = 0; j < npop - 1; j++)      /* duplicated calculation! */
               if (r < (y += nin_G[i] * stree.M[ipop[j] * s21 + ipop[i]] / stree.nodes[ipop[i]].theta * 4))
                  break;

            mM[ipop[j] * (stree.nspecies * 2 - 1) + ipop[i]] ++;

            k = (int)(nin_G[i] * rndu());      /* k is migrant from pop j to i */
            /* shift up lineages in pop i */
            ins_G[j][nin_G[j] ++] = ins_G[i][k];
            if (k != --nin_G[i])
               ins_G[i][k] = ins_G[i][nin_G[i]];

            if (debug == 9)
               printf("M: %s < %s ", stree.nodes[ipop[i]].name, stree.nodes[ipop[j]].name);
         }
         if (debug == 9) {
            for (i = 0; i < npop; i++) {
               printf(" %-s:", stree.nodes[ipop[i]].name);
               for (j = 0; j < nin_G[i]; j++)
                  printf("%d ", ins_G[i][j]);
            }
            FPN(F0);
         }
      }  /* forever loop inside for(epoch) */
      T = Tmax;
      /* To move to next epoch, update ipop[] and merge lineages from sons[0] &
         sons[1] into pop is.
       */
      if (npop == 1 || n == 1) break;
      sons = stree.nodes[is].sons;
      /* Replace pop k1 by is, move k2 into k1, delete k2. */
      for (i = 0, k1 = k2 = -1; i < npop; i++) {
         if (ipop[i] == sons[0])       k1 = i;
         else if (ipop[i] == sons[1])  k2 = i;
      }

      if (k1 > k2) {
         i = k1; k1 = k2; k2 = i;
      }
      ipop[k1] = is;
      for (i = 0; i < nin_G[k2]; i++)
         ins_G[k1][nin_G[k1] ++] = ins_G[k2][i];
      if (k2 != npop - 1) {
         ipop[k2] = ipop[npop - 1];
         if ((nin_G[k2] = nin_G[npop - 1]) > 0)
            memmove(ins_G[k2], ins_G[npop - 1], nin_G[k2] * sizeof(int));
      }
   }

   tree.root = cNodeNumber - 1;  /* cNodeNumber = com.ns*2 - 1.  root is last node */
   if (cNodeNumber != com.ns * 2 - 1) error2("cNodeNumber incorrect");

   return 0;
}


int printGtree(int printBlength)
{
   int status = 0, inode, j, ipop, ipopOK;
   double t, tb[2];

   for (inode = 0; inode < tree.nnode; inode++) {
      if (inode != tree.root) {
         nodes[inode].branch = nodes[nodes[inode].father].age - nodes[inode].age;
         if (nodes[inode].branch < 0) {
            printf("blength [%d] = %9.6g\a\n", inode, nodes[inode].branch);
            status = -1;
         }
      }
      nodes[inode].label = (inode < com.ns ? 0 : stree.nodes[nodes[inode].ipop].age);
   }
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s%12s (%s) %7s%7s", "father", "node", "age", "ipop", "nson:", "sons");
   for (inode = 0; inode < tree.nnode; inode++) {
      t = nodes[inode].age;
      ipop = nodes[inode].ipop;
      tb[0] = stree.nodes[ipop].age;
      tb[1] = OLDAGE;
      if (ipop != stree.root)
         tb[1] = stree.nodes[stree.nodes[ipop].father].age;

      ipopOK = (t >= tb[0] && t <= tb[1]);
      if (!ipopOK)
         status = -2;
      printf("\n%7d%7d %11.6f (%2d %s %11.6f):%7d  ",
         nodes[inode].father, inode, t, ipop, (ipopOK ? "OK" : "??"), stree.nodes[ipop].age, nodes[inode].nson);
      for (j = 0; j < nodes[inode].nson; j++) printf(" %2d", nodes[inode].sons[j]);
   }
   FPN(F0); OutTreeN(F0, 0, 0); FPN(F0); OutTreeN(F0, 1, 0);
   if (printBlength) {
      FPN(F0); OutTreeN(F0, 1, 1); FPN(F0);
   }

   /* if(status==-1) exit(status); */
   return(status);
}


int MatchGTree(void)
{
   /* This tests the gene tree topology
   */
   return
      (nodes[0].father == nodes[1].father
         && nodes[3].father == nodes[4].father
         && nodes[nodes[3].father].father == tree.root
         && nodes[nodes[0].father].age < stree.nodes[stree.root].age
         && nodes[nodes[2].father].age < stree.nodes[stree.root].age
         && nodes[nodes[3].father].age < stree.nodes[stree.root].age
         );  /* P(gtree) = 0.237913 */
}



#ifdef SIMULATION

static int n123marks[][3] = { {1,2,3},{1,2,3},{2,3,1},{3,1,2} };

void p0124Fromb0b1(int itree, double p[5], double b[2])
{
   /* This calculates p0,p1,p2,p3,p4 for the 5 site patterns for 3 species,
      given branch lengths b0 and b1.  b0 is the gap, and b1 is the distance
      from the ancestor of 1 and 2 to species 1.
   */
   double e1, e2, e3;

   e1 = exp(-4. / 3 * b[1]);  e2 = exp(-8. / 3 * (b[0] + b[1]));  e3 = e1*e2;  e1 = e1*e1;
   p[0] = (1. + 3 * e1 + 6 * e2 + 6 * e3) / 16;
   p[n123marks[itree][0]] = (3. + 9 * e1 - 6 * e2 - 6 * e3) / 16;
   p[n123marks[itree][1]] = p[n123marks[itree][2]] = (3. - 3 * e1 + 6 * e2 - 6 * e3) / 16;
   p[4] = (6. - 6 * e1 - 12 * e2 + 12 * e3) / 16;
}

void PMismatch3s(void)
{
   /* This calculate the tree mismatch probabilities for figure 3 of Yang
      (2002 Genetics).
   */
   int iHCG = stree.root, iHC = 3 + 4 - iHCG; /* node numbers in stree */
   double theta_HC = stree.nodes[iHC].theta;
   double t_HC = stree.nodes[iHC].age, t_HCG = stree.nodes[iHCG].age;
   /* S: species tree; G: gene tree; E: estimated gene tree
      G[gtree], E[gtree][etree] are counts of resolved (used) loci.
   */
   double SG, SE, GE, G[4], E[4][4], b[2], p0[5], p[5], nused, over, under;
   int ii, nii = 7, ls0[100] = { 200, 400, 500, 1000, 2000, 4000, 10000 };
   int nr = 10000000, ir, i, j, gtree = -1, n[5] = { 0 }, etree, etree2, etree3, every = 100;

   printf("\nPr{S-G mismatch} = %f from equation.\n", 2. / 3 * exp(-2 * (t_HCG - t_HC) / theta_HC));
   puts("Ties in genetree are removed in the following calculation.");
   every = max2(every, nr / 1000);
   for (ii = 0; ii < nii; ii++) {
      com.ls = ls0[ii];
      printf("\n# sites? (Ctrl-C to break) ");
      scanf("%d", &com.ls);
      printf("%d sites, %dK replicates.\n", com.ls, nr / 1000);
      FOR(i, 4) G[i] = 0; FOR(i, 4) FOR(j, 4) E[i][j] = 0;
      for (ir = 0, SG = SE = GE = 0, nused = 0; ir < nr; ir++) {
         GetRandomGtree(-1);
         if (nodes[2].father == tree.root)      gtree = (nodes[0].branch > t_HCG);/* (HC) */
         else if (nodes[0].father == tree.root)  gtree = 2; /* (CG) */
         else if (nodes[1].father == tree.root)  gtree = 3; /* (GH) */
         else    error2("binary tree?");
         b[1] = nodes[nodes[tree.root].sons[0]].age;
         b[0] = nodes[tree.root].age - b[1];
         p0124Fromb0b1(gtree, p0, b);
         p[0] = p0[0] + p0[4];
         p[1] = p[0] + p0[1];
         p[2] = p[1] + p0[2];
         p[3] = p[2] + p0[3];
         p[4] = -1;
         /* matout(F0, p, 1, 4); */
         MultiNomial2(com.ls, 4, p, n, NULL);

         for (i = 2, etree = 1; i <= 3; i++)   if (n[i] > n[etree]) etree = i;
         if (etree == 1) { etree2 = 2; etree3 = 3; }
         else if (etree == 2) { etree2 = 3; etree3 = 1; }
         else { etree2 = 1; etree3 = 2; }

         if (n[etree] != n[etree2] && n[etree] != n[etree3]) { /* exclude ties */
            nused++;
            G[gtree]++;
            E[gtree][etree]++;
            if (gtree >= 2)  SG++;
            if (etree >= 2)  SE++;
            if ((gtree <= 1 && etree != 1) || (gtree >= 2 && etree != gtree)) GE++;
         }
         if ((ir + 1) % every == 0) {
            printf("%4.1f%% (%2d %2d %2d %2d) SG %.4f SE %.4f GE %.4f (+%.4f -%.4f) tie %.4f\r",
               (ir + 1.) / nr * 100, n[0], n[1], n[2], n[3],
               SG / nused, SE / nused, GE / nused,
               (E[0][2] + E[0][3] + E[1][2] + E[1][3]) / nused,  /* T -> F */
               (E[2][1] + E[3][1]) / nused,                    /* T -> F */
               (ir + 1 - nused) / (ir + 1.));
         }
      }  /* for(ir) */
      if (G[0] + G[1] + G[2] + G[3] - nused != 0)  error2("gtree counts incorrect.");
      SG /= nused;
      SE /= nused;
      GE /= nused;
      FOR(i, 4) FOR(j, 4) E[i][j] /= G[i];
      FOR(i, 4) G[i] /= nused;
      printf("\nfrequencies of gene trees 0123 (given ties are removed): ");
      matout2(F0, &G[0], 1, 4, 9, 5);
      printf("transition probability matrix (gene tree 0123 -> MLtree (01)23):");
      matout2(F0, &E[0][0], 4, 4, 9, 5);

      printf("\nThe following three should be equal:\n");
      printf("  (1) SE - SG = %.5f\n", SE - SG);
      printf("  (2) f(T0) * P0 = %.5f\n", G[0] * (E[0][2] + E[0][3]));
      over = G[0] * (E[0][2] + E[0][3]) + G[1] * (E[1][2] + E[1][3]);
      under = G[2] * E[2][1] + G[3] * E[3][1];
      printf("  (3) over - under = %.5f - %.5f = %.5f\n", over, under, over - under);
      printf("\nf(2)*P2/2+f(3)*P3/2 = %.5f, which is GE - (over + under) = under.\n",
         G[2] * E[2][3] + G[3] * E[3][2]);
      printf("\nTime used: %s\n", printtime(timestr));
   }
   exit(0);
}


int MakeDiploid(unsigned char *z[], int *ns, int ls)
{
   /* This takes ns sequences in z[ns][] and produces nsHet sequences, with heterozygotes
      by using data.diploid[].  Memory in z[nsHet, ..., ns-1] is released.
      Sequences are already encoded.
   */
   int i, j, k, h, is = 0, isHet = 0;
   unsigned char b[2], c[2], seqname[64], *pch;

   for (i = 0; i < stree.nspecies; i++) {
      for (j = 0; j < stree.nseqsp[i] / (1 + data.diploid[i]); j++) {
         if (data.diploid[i] == 0) {
            for (h = 0; h < ls; h++)
               z[isHet][h] = z[is][h];
            is++;  isHet++;
         }
         else {
            for (h = 0; h < ls; h++) {
               b[0] = z[is][h];  b[1] = z[is + 1][h];
               if (b[0] == b[1])
                  z[isHet][h] = b[0];
               else {
                  c[0] = min2(b[0], b[1]);  c[1] = max2(b[0], b[1]);
                  if (c[0] == 0 && c[1] == 1)         z[isHet][h] = 5;   /* TC = Y */
                  else if (c[0] == 2 && c[1] == 3)    z[isHet][h] = 6;   /* AG = R */
                  else if (c[0] == 1 && c[1] == 2)    z[isHet][h] = 7;   /* CA = M */
                  else if (c[0] == 0 && c[1] == 3)    z[isHet][h] = 8;   /* TG = K */
                  else if (c[0] == 1 && c[1] == 3)    z[isHet][h] = 9;   /* CG = S */
                  else if (c[0] == 0 && c[1] == 2)    z[isHet][h] = 10;  /* TA = W */
               }
            }
            is += 2;  isHet++;
         }
      }
   }
   *ns = isHet;
   for (i = 0, k = 0; i < stree.nspecies; i++) {
      strcpy(seqname, stree.nodes[i].name);
      pch = seqname;
      while (*pch) { *pch = tolower(*pch);  pch++; }
      if (data.diploid[i]) {
         for (j = 0; j < stree.nseqsp[i] / 2; j++, k++)
            sprintf(com.spname[k], "%s%d_%d^%s", seqname, j * 2 + 1, j * 2 + 2, stree.nodes[i].name);
      }
      else {
         for (j = 0; j < stree.nseqsp[i]; j++, k++)
            sprintf(com.spname[k], "%s%d^%s", seqname, j + 1, stree.nodes[i].name);
      }
   }
   return(0);
}



void MakeSeq(unsigned char *z, int ls, double pi[4])
{
   /* generate a random sequence in z
   */
   int j, h;
   double p[4], r;

   /* Make sequence at the root */
   if (com.model == 0) {
      for (h = 0; h < ls; h++)
         z[h] = (char)(rndu() * 4);
   }
   else {
      for (j = 0; j < 4; j++)  p[j] = com.pi[j];
      for (j = 1; j < 4; j++)  p[j] += p[j - 1];
      for (h = 0; h < ls; h++) {
         for (j = 0, r = rndu(); j < 4 - 1; j++)
            if (r < p[j]) break;
         z[h] = (unsigned char)j;
      }
   }
}


double PMat[4 * 4], U[4 * 4], V[4 * 4], Root[4];
double Qfactor = -1;


int eigenQbase(double Qrates[6], double pi[], double Root[], double U[], double V[], double Q[])
{
   /* Construct the rate matrix Q[] for nucleotide model REV.
   */
   int i, j, k = 0;
   double mr, space[4];

   for (i = 0; i < 16; i++) Q[i] = 0;
   for (i = 0; i < 3; i++) for (j = i + 1; j < 4; j++)
      Q[i * 4 + j] = Q[j * 4 + i] = Qrates[k++];
   for (i = 0; i < 4; i++) for (j = 0; j < 4; j++)  Q[i * 4 + j] *= pi[j];
   for (i = 0, mr = 0; i < 4; i++) {
      Q[i * 4 + i] = 0;
      Q[i * 4 + i] = -sum(Q + i * 4, 4);
      mr -= pi[i] * Q[i * 4 + i];
   }
   abyx(1 / mr, Q, 16);

   eigenQREV(Q, com.pi, 4, Root, U, V, space);
   return (0);
}

void EvolveGTR(int inode)
{
   /* This is modified from Evolve() in evolver.c.  It evolves sequence the com.z[tree.root] along
   the tree to generate com.z[], using nodes[].branch, com.model.  Model should be GTR+G or simpler.
   Needs com.z[0,1,...,nnode-1], while com.z[0] -- com.z[ns-1] constitute the data.
   For codon sequences, data.siterates[] has w's for NSsites and NSbranchsite models.
   */
   int is, h, i, j, ison, from, n = com.ncode, longseq = 100000;
   double t, r;

   for (is = 0; is < nodes[inode].nson; is++) {
      ison = nodes[inode].sons[is];
      memcpy(com.z[ison], com.z[inode], com.ls * sizeof(unsigned char));
      t = nodes[ison].branch;

      for (h = 0; h < com.ls; h++) {
         /* decide whether to recalcualte PMat[]. */
         if (h == 0 || (data.siterates && data.siterates[h] != data.siterates[h - 1])) {
            r = (data.siterates ? data.siterates[h] : 1);
            PMatUVRoot(PMat, t*r, com.ncode, U, V, Root);
            for (i = 0; i < n; i++)
               for (j = 1; j < n; j++)
                  PMat[i*n + j] += PMat[i*n + j - 1];
         }
         for (j = 0, from = com.z[ison][h], r = rndu(); j < n - 1; j++)
            if (r < PMat[from*n + j]) break;
         com.z[ison][h] = j;
      }

      if (com.ls > longseq) printf("\r   nodes %2d -> %2d, evolving . .   ", inode + 1, ison + 1);

      if (nodes[ison].nson) EvolveGTR(ison);
   }  /* for (is) */

   if (inode == tree.root && com.ls > longseq)  printf("\r%s", strc(50, ' '));
}

void EvolveJC(int inode)
{
   /* Special version of Evolve that works with JC69-like (poisson) model only.
      For each branch in the tree, the number of mutations is sampled from the possoin and assigned
      to sites at random.
      When alpha>0, data.siterates[] are accumulative probabilities.
   */
   int is, j, h, nmut, imut, ison;
   double r;

   if (com.alpha && fabs(data.siterates[com.ls - 1] - 1) > 1e-4) {
      printf("rates c.d.f.: 1 = %.6f?\n", data.siterates[com.ls - 1]); exit(-1);
   }
   for (is = 0; is < nodes[inode].nson; is++) {
      ison = nodes[inode].sons[is];
      for (h = 0; h < com.ls; h++) com.z[ison][h] = com.z[inode][h];
      nmut = rndpoisson(nodes[ison].branch*com.ls);
      for (imut = 0; imut < nmut; imut++) {
         if (com.alpha == 0)
            h = (int)(rndu()*com.ls);
         else
            for (h = 0, r = rndu(); h < com.ls; h++)
               if (r < data.siterates[h]) break;
         j = (int)(rndu() * 3);
         if (j >= com.z[ison][h]) j++;
         com.z[ison][h] = (char)j;
      }
      if (nodes[ison].nson) EvolveJC(ison);
   }
}


int Count3sTrees(void)
{
   /* This identifies the gene tree for the case of 3 species & 3 sequences (Dalquen et al. 2017).
      This is for debugging.
   */
   int iGtree = 0, outseq = (nodes[4].sons[0] < 3 ? nodes[4].sons[0] : nodes[4].sons[1]);
   double tau0 = stree.nodes[3].age, tau1 = stree.nodes[4].age;
   double t0 = nodes[4].age, t1 = nodes[3].age;

   if (stree.nspecies != 3 || com.ns != 3) error2("should not be here.");
   if (t0 < tau1) {
      if (outseq == 2)       iGtree = 0;
      else if (outseq == 0)  iGtree = 1;
      else if (outseq == 1)  iGtree = 2;
   }
   else if (t1 < tau1 && tau1 < t0 && t0 < tau0) {
      if (outseq == 2)       iGtree = 3;
      else if (outseq == 0)  iGtree = 4;
      else if (outseq == 1)  iGtree = 5;
   }
   else if (t1 < tau1 && t0 > tau0) {
      if (outseq == 2)       iGtree = 6;
      else if (outseq == 0)  iGtree = 7;
      else if (outseq == 1)  iGtree = 8;
   }
   else if (tau1 < t1 && t1 < tau0 && tau1 < t0 && t0 < tau0) {
      if (outseq == 2)       iGtree = 9;
      else if (outseq == 0)  iGtree = 10;
      else if (outseq == 1)  iGtree = 11;
   }
   else if (tau1 < t1 && t1 < tau0 && t0 > tau0) {
      if (outseq == 2)       iGtree = 12;
      else if (outseq == 0)  iGtree = 13;
      else if (outseq == 1)  iGtree = 14;
   }
   else if (t1 > tau0 && t0 > tau0) {
      if (outseq == 2)       iGtree = 15;
      else if (outseq == 0)  iGtree = 16;
      else if (outseq == 1)  iGtree = 17;
   }
   return(iGtree);
}

void ClockSpeciesNodeRates(int snode, double snoderates[])
{
   int i, ison;

   for (i = 0; i < stree.nodes[snode].nson; i++) {
      ison = stree.nodes[snode].sons[i];
      if (stree.nodes[snode].age == 0)
         snoderates[ison] = snoderates[snode];
      else                /* gamma rate with shape a and mean to be old rate */
         snoderates[ison] = rndgamma(data.a_branchrate) / data.a_branchrate * snoderates[snode];
      ClockSpeciesNodeRates(ison, snoderates);
   }
}

void RelaxedClockBranchLengths(double snoderates[])
{
   /* This generates rates for species-tree branches and update gene-tree branch lengths.
   This right now uses independent gamma rates for species-tree branches.
   */
   int debug = 0, is, ipop, i, dad;
   double tdad, t;

   if (data.a_branchrate < 0) error2("a_branchrat < 1 in RelaxedClockBranchLengths");
   if (com.clock == 2) {  /* independent gamma rates  */
      for (is = 0; is < stree.nnode; is++)
         snoderates[is] = rndgamma(data.a_branchrate) / data.a_branchrate;
   }
   else {                 /* autocorrelated gamma rates  */
      snoderates[stree.root] = rndgamma(data.a_branchrate) / data.a_branchrate;
      ClockSpeciesNodeRates(stree.root, snoderates);
   }
   if (debug) {
      matout(F0, snoderates, 1, stree.nnode);
      printGtree(1);
   }
   for (i = 0; i < tree.nnode; i++) {
      nodes[i].branch = 0;
      if (i == tree.root) continue;
      dad = nodes[i].father;
      t = nodes[i].age;
      tdad = nodes[dad].age;
      ipop = nodes[i].ipop;
      is = stree.nodes[ipop].father;
      if (debug) printf("\nbranch %d (%.6f - %.6f)\n", i, t, tdad);
      for (; is != -1 && tdad > stree.nodes[is].age; is = stree.nodes[is].father) {
         nodes[i].branch += (stree.nodes[is].age - t)*snoderates[ipop];
         if (debug) printf("\tin pop %d, add (%.6f - %.6f)*%.6f\n", ipop, stree.nodes[is].age, t, snoderates[ipop]);
         t = stree.nodes[is].age;
         ipop = is;
      }
      nodes[i].branch += (tdad - t)*snoderates[ipop];
      if (debug) printf("\tin pop %d, add (%.6f - %.6f)*%.6f  = %.6f\n", ipop, tdad, t, snoderates[ipop], nodes[i].branch);
   }
}

void SimulateData(void)
{
   char *zcon[NS], seqname[64], *pch, timestr[64];
   FILE *fseq = NULL, *ftree = NULL, *fImap = NULL, *fconcat = NULL, *fmodelpara = NULL;
   FILE *fseqfull = NULL;
   double rlocus, y, snoderates[NSPECIES * 2 - 1];
   int nloci = com.ngene, locus, i, j, k, h, iGtree, maxls = com.ls;
   int *siteorder = NULL, nsFull = com.ns, nsHet = com.ns;
   char *tmpseq = (char *)com.space, line[96] = "", *p;
   char *G3s[18] = { "G1c", "G1a", "G1b", "G2c", "G2a", "G2b", "G3c", "G3a", "G3b", "G4c", "G4a", "G4b", "G5c", "G5a", "G5b", "G6c", "G6a", "G6b" };
   double PG3s[18] = { 0 };

   /* PMismatch3s(); */
   com.ngene = 1;
   if (com.seqf[0] || com.concatf[0]) {  /* set nsHet for diploid data with heterozygotes */
      /* stree.nseqsp[] has the number of full sequences, set in GetOptionsSimulation. */
      for (i = 0, nsHet = 0; i < stree.nspecies; i++)
         nsHet += stree.nseqsp[i] / (data.diploid[i] ? 2 : 1);
      com.cleandata = 1;
      for (i = 0; i < stree.nspecies; i++)  if (data.diploid[i]) com.cleandata = 0;
   }
   if (com.seqf[0]) {
      if (stree.nspecies > 1) printf("Imap file is %s\n", com.Imapf);
      printf("\nsequence data file is %s\n", com.seqf);
      fseq = gfopen(com.seqf, "w");
      fImap = gfopen(com.Imapf, "w");
      if (nsHet < com.ns) {
         strcpy(line, com.seqf);
         p = strstr(line, ".txt");
         if(p)  sprintf(p, "_full.txt");
         else   strcat(line, "_full");
         printf("Fullly resolved haploid sequences are in %s\n", line);
         fseqfull = gfopen(line, "w");
      }
   }
   if (com.treef[0]) {
      printf("trees go into %s.\n", com.treef);
      ftree = gfopen(com.treef, "w");
   }
   if (com.concatf[0]) {
      printf("concatenated sequence alignment is in %s\n", com.concatf);
      fconcat = gfopen(com.concatf, "w");
      for (i = 0; i < nsHet; i++) {
         if ((zcon[i] = (char*)malloc(com.ls*nloci)) == NULL)
            error2("oom zcon");
      }
   }
   if (com.modelparaf[0]) {
      printf("Model parameters for loci are in %s\n", com.modelparaf);
      fmodelpara = gfopen(com.modelparaf, "w");
      if (com.model == 0 && data.a_locusrate)
         fprintf(fmodelpara, "locus\tlocusrate\n");
      else if (com.model)
         fprintf(fmodelpara, "locus\tQrates_abcdef\tpi_TACG\talpha\tlocusrate\n");
   }
   for (j = 0, k = 0; j < stree.nspecies; j++) {
      if (fImap) fprintf(fImap, "%s\t%s\n", stree.nodes[j].name, stree.nodes[j].name);
      for (i = 0; i < stree.nseqsp[j]; i++, k++)
         nodes[k].ipop = j;
   }
   if (fImap) fclose(fImap);

   if (fseq) {
      if (data.fix_alpha == 0 || com.alpha) {
         if ((data.siterates = (double*)malloc(com.ls * sizeof(double))) == NULL) error2("oom");
         if ((siteorder = (int*)malloc(com.ls * sizeof(int))) == NULL) error2("oom2");
      }
   }

   for (locus = 0; locus < nloci; locus++) {
      if (fmodelpara)
         fprintf(fmodelpara, "%d", locus + 1);
      if (com.model) {
         if (data.fix_Qrates == 0) {
            rnddirichlet(data.Qrates, data.Qrates_Dir, 6);
            for (i = 0; i < 6; i++)  data.Qrates[i] /= data.Qrates[5];
         }
         if (data.fix_pi == 0)
            rnddirichlet(com.pi, data.pi_Dir, 4);

         for (i = 0; i < 6; i++) fprintf(fmodelpara, " %9.6f", data.Qrates[i]);
         for (i = 0; i < 4; i++) fprintf(fmodelpara, " %8.6f", com.pi[i]);
         eigenQbase(data.Qrates, com.pi, Root, U, V, PMat);
      }
      if (data.fix_alpha == 0) {
         com.alpha = rndgamma(data.a_siterate[0]) / data.a_siterate[1];
         if (fmodelpara)  fprintf(fmodelpara, " %9.6f", com.alpha);
      }
      if (data.a_locusrate) {
         rlocus = rndgamma(data.a_locusrate) / data.a_locusrate;   /* G(a, a) */
         if (fmodelpara)  fprintf(fmodelpara, " %9.6f", rlocus);
      }

      if (fseq) { /* this specifies variable sequence length.  edit maxls as well. */
                  /* com.ls = 200 + 800*rndu();  */
      }
      if (fseq || ftree) {
         for (i = 0, k = 0; i < stree.nspecies; i++) {
            strcpy(seqname, stree.nodes[i].name);
            pch = seqname;
            while (*pch) { *pch = tolower(*pch);  pch++; }
            if (data.diploid[i])
               for (j = 0; j < stree.nseqsp[i]; j++, k++)
                  sprintf(com.spname[k], "%s%d%c^%s", seqname, j / 2 + 1, 'a' + j % 2, stree.nodes[i].name);
            else
               for (j = 0; j < stree.nseqsp[i]; j++, k++)
                  sprintf(com.spname[k], "%s%d^%s", seqname, j + 1, stree.nodes[i].name);
         }
      }

      GetRandomGtree(-1);
      if (stree.nspecies == 3 && com.ns == 3)
         PG3s[iGtree = Count3sTrees()] ++;

      if (com.clock == 2 || com.clock == 3) {
         RelaxedClockBranchLengths(snoderates);
         for (j = 0; j < stree.nnode; j++) fprintf(fmodelpara, " %.6f", snoderates[j]);
      }
      if (fmodelpara) fprintf(fmodelpara, "\n");

      if (data.a_locusrate) {
         for (j = 0; j < tree.nnode; j++)
            if (j != tree.root)  nodes[j].branch *= rlocus;
      }

      if (ftree) {   /*  branch lengths in tree incorporates locus rate and branch rates */
         OutTreeN(ftree, 1, 1);
         /* tree height is tau, without using locusrate and branchrate */
         if (stree.nspecies == 3 && com.ns == 3)
            fprintf(ftree, " [TH = %.6f, %s]\n", nodes[tree.root].age, G3s[iGtree]);
         else
            fprintf(ftree, " [TH = %.6f]\n", nodes[tree.root].age);
      }

      if (fseq) {  /* evolve and print sequence alignment at locus */
         if (com.alpha)
            Rates4Sites(data.siterates, com.alpha, com.ncatG, com.ls, (com.model == 0), com.space);
         for (i = 0; i < 2 * com.ns - 1; i++) {
            com.z[i] = (unsigned char*)realloc(com.z[i], maxls * sizeof(unsigned char));
            if (com.z[i] == NULL) error2("memory problem for z[]");
         }
         MakeSeq(com.z[tree.root], com.ls, com.pi);

         if (com.model == 0)
            EvolveJC(tree.root);
         else
            EvolveGTR(tree.root);

         /* randomize sites for discrete-gamma model */
         if (com.alpha && com.ncatG > 1) {
            randorder(siteorder, com.ls, (int*)com.space);
            for (j = 0; j < tree.nnode; j++) {
               memcpy(tmpseq, com.z[j], com.ls * sizeof(char));
               for (h = 0; h < com.ls; h++) com.z[j][h] = tmpseq[siteorder[h]];
            }
         }

         /*
         for (i = 0; i < com.ns; i++) {
            if (data.iseqerr[j = nodes[i].ipop] == 0)
               continue;
            for (h = 0; h < com.ls; h++) {
               for (k = 0, u = rndu(); k < 3; k++)
                  if (u < data.e_seqerr[j][com.z[i][h] * 4 + k])
                     break;
               com.z[i][h] = k;
            }
         }
         */
         if (nsHet < com.ns) {
            if(fseqfull)
               printSeqs(fseqfull, NULL, NULL, com.SeqFileFormat);
            MakeDiploid(com.z, &com.ns, com.ls);
         }
         if (fconcat) {
            for (i = 0; i < com.ns; i++)
               memcpy(zcon[i] + locus*com.ls, com.z[i], com.ls * sizeof(char));
         }

         if (com.SeqFileFormat == 1) {
            PatternWeightSimple();
            /*   if (com.model == 0) PatternWeightJC69like(); */
         }
         printSeqs(fseq, NULL, NULL, com.SeqFileFormat); /* printsma not usable as it codes into 0,1,...,60. */
         com.ns = nsFull;
      }   /* if(fseq) */
      if ((locus + 1) % 1000 == 0 || (nloci > 1000 && locus == nloci - 1))
         printf("\r%10d replicates done... mean tMRCA = %9.6f  %s", locus + 1, mtMRCA / (locus + 1), printtime(timestr));
   }  /* for(locus) */

   if (ftree) fclose(ftree);
   if (fseq) {
      fclose(fseq);
      for (i = 0; i < 2 * com.ns - 1; i++)  free(com.z[i]);
      if (data.fix_alpha == 0 || com.alpha) {
         free(data.siterates);
         free(siteorder);
      }
      if (fconcat) {
         printf("\ngenerating concatenated sequence alignment...\n");
         com.ls *= nloci;
         com.ns = nsHet;
         for (i = 0; i < com.ns; i++)  com.z[i] = zcon[i];
         printSeqs(fconcat, NULL, NULL, 0);
         fclose(fconcat);
         for (i = 0; i < com.ns; i++)  free(zcon[i]);
      }
   }

   for (i = 0; i < (stree.nspecies * 2 - 1)*(stree.nspecies * 2 - 1); i++)
      mM[i] /= nloci;
   if (stree.migration) {
      printf("\nCounts of migration events averaged over replicates: %8.4f\n", mtMRCA / nloci);
      matout(F0, mM, (stree.nspecies * 2 - 1), (stree.nspecies * 2 - 1));
   }

   if (stree.nspecies == 3 && com.ns == 3 && nsFull == nsHet) {
      printf("\nProportion of 18 gene trees as defined in Dalquen, Zhu & Yang (2017)\n");
      for (k = 0; k < 18; k++) {
         y = PG3s[k / 3 * 3] + PG3s[k / 3 * 3 + 1] + PG3s[k / 3 * 3 + 2];
         printf("P[%d%c] = %9.5f %9.0f", k / 3 + 1, 'a' + (k + 2) % 3, PG3s[k] / nloci, PG3s[k]);
         if (y > 0)  printf("%12.6f\n", PG3s[k] / y);
         else       printf("\n");
         if ((k + 1) % 3 == 0) printf("\n");
      }
   }
   printf("\nTime used: %s\n", printtime(timestr));
   exit(0);
}


#else


int GetRootTau(void)
{
   /* This scans the sequence alignments at different loci to calculate the distance
      from of the species tree to the present, for use to specify the upper bound tauU.
   */
   int root = stree.root, *sons = stree.nodes[root].sons, ipopi, ipopj, npair;
   int h, i, j, locus, locusused = 0, diploiddata = 0;
   double md = 0, vd = 0, dlocus, dpair, thetaA;
   int debug = 0;

   for (i = 0; i < stree.nspecies; i++) if (data.diploid[i]) diploiddata = 1;
   if (diploiddata) printf("renovate GetRootTau, especially for diploid data.\n");
   stree.roottau = data.tau_prior[1] / (data.tau_prior[0] - 1) * 4;
   if (mcmc.usedata && !diploiddata) {
      for (locus = 0; locus < data.ngene; locus++) {
         for (i = 0, dlocus = 0, npair = 0; i < data.ns[locus]; i++) {
            ipopi = gnodes[locus][i].ipop;
            for (j = 0; j < data.ns[locus]; j++) {
               ipopj = gnodes[locus][j].ipop;
               if (stree.pptable[ipopi][sons[0]] && stree.pptable[ipopj][sons[1]]) {
                  for (h = 0, dpair = 0; h < data.npatt[locus]; h++)
                     if (data.z[locus][i][h] != data.z[locus][j][h]) dpair += data.fpatt[locus][h];
                  dlocus += dpair /= data.ls[locus];
                  npair++;

                  if (debug)
                     printf("locus %2d S %2d (ipop %d) & %2d (ipop %d): d = %7.5f\n", locus + 1, i + 1, ipopi, j + 1, ipopj, dpair);
               }
            }
         }
         if (npair == 0) continue;
         locusused++;
         dlocus /= (2 * npair);
         vd += square(dlocus - md) * (locusused - 1) / locusused;
         md = (md * (locusused - 1) + dlocus) / locusused;
      }

      vd /= data.ngene;
      if (locusused >= 2) {
         thetaA = 2 * sqrt(vd);
         thetaA = sqrt(vd * 4 + 1) - 1;
         thetaA = (2 * sqrt(vd) + sqrt(vd * 4 + 1) - 1) / 2;
         stree.roottau = (md - thetaA / 2 > 0 ? md - thetaA / 2 : md);
      }
      else
         stree.roottau = md;
   }
   printf("\nroot dist = %7.5f\n", stree.roottau);
   return(0);
}


int ReadSeqData(char *seqfile, char *locusratef, char *heredityf, FILE*fout, char cleandata, int *ipop)
{
   /* Read sequences at each locus. This sets data.nseqsp[ig].  The ipop info for
      tips on gene trees is returned in ipop[], and passed to GetMem(), which allocates
      gnodes[].

      use com.cleandata=1 to delete ambiguities at all loci.
   */
   FILE *fseq = gfopen(seqfile, "r"), *frates = NULL, *fheredity = NULL, *fImap = NULL;
   char line[10000], *pline, iname[100], tag = '^';
   int s = stree.nspecies, locus, i, j, k, is, ind, nsFull, lname = 100;
   int maxnind = data.maxns * 2 * data.ngene, nind = 0;
   int maxns; /* data.maxns is from the control file and may be too large, reset to maxns from seq file.  */
   double mr;
   char cleandata0 = cleandata;

   /* Reading Imap file */
   if (s > 1) {
      printf("\nReading Individual-Species map (Imap) from %s\n", com.Imapf);
      data.Imap = (int*)malloc(maxnind * sizeof(int));
      data.Indnames[0] = (char*)malloc(maxnind * lname * sizeof(char));
      if (data.Imap == NULL || data.Indnames[0] == NULL) error2("oom Imap & Indnames");
      memset(data.Imap, 0, maxnind * sizeof(int));
      memset(data.Indnames[0], 0, maxnind * lname * sizeof(char));
      fImap = gfopen(com.Imapf, "r");
      for (i = 0; i < maxnind; i++) {
         data.Indnames[i] = data.Indnames[0] + i*lname;
         if (fscanf(fImap, "%s%s", data.Indnames[i], line) != 2) break;
         if (strstr(data.Indnames[i], "//")) break;   /* this marks the end of file. */
         for (j = 0; j < s; j++)
            if (strcmp(line, stree.nodes[j].name) == 0) {  /* ind i is from pop j. */
               data.Imap[i] = j;
               break;
            }
         if (j == s) {
            printf("\nspecies %s in map file not found in the control file.\n", line);
            exit(-1);
         }
      }
      printf("Individual -> Species map: ");
      nind = i;
      for (i = 0; i < nind; i++)
         printf(" %d", data.Imap[i] + 1);
      printf("\n# individuals = %2d read from Imap file\n", nind);
   }

   printf("\nReading sequence data..  %d loci\n", data.ngene);
   if (data.ngene > NGENE) error2("raise NGENE in bpp.c & recompile?");
   if (data.est_locusrate && data.ngene == 1)
      error2("ngene = 1 & locus rates");

   /* allocate space for heredity scalars and locus rates.  Space is wasted if only one is in use. */
   if (data.est_heredity || data.est_locusrate) {
      data.heredity = (double*)malloc(2 * data.ngene * sizeof(double));
      if (data.heredity == NULL) error2("oom heredity");
      data.locusrate = data.heredity + data.ngene;
      for (j = 0; j < data.ngene; j++) data.heredity[j] = 1;
      for (j = 0; j < data.ngene; j++) data.locusrate[j] = 1;
   }
   /* Read locus rates from the file locusratef.  File name is read from control file
      but data.ngene was unknown there. */
   if (data.est_locusrate == 2) {
      frates = (FILE*)fopen(locusratef, "r");
      if (frates == NULL)
         error2("locus rate file does not exist");
      for (locus = 0, mr = 0; locus < data.ngene; locus++) {
         if (fscanf(frates, "%lf", &data.locusrate[locus]) != 1) {
            printf("\nEOF when reading rates from %s\n", locusratef);
            error2("");
         }
         mr = (mr*locus + data.locusrate[locus]) / (locus + 1.0);
      }
      fclose(frates);
      for (locus = 0; locus < data.ngene; locus++)
         data.locusrate[locus] /= mr;
      printf("using fixed rates for loci from file %s\n", locusratef);
      if (data.ngene <= 200)
         matout2(F0, data.locusrate, 1, data.ngene, 8, 4);
   }

   /* Read heredity scalars for loci from heredityf.  The file name is read from
      the control file but data.ngene was unknown there. */
   if (data.est_heredity == 2) {
      if ((fheredity = (FILE*)fopen(heredityf, "r")) == NULL)
         error2("heredity scalar file does not exist");
      for (locus = 0; locus < data.ngene; locus++) {
         if (fscanf(fheredity, "%lf", &data.heredity[locus]) != 1) {
            printf("\nEOF when reading rates from %s\n", heredityf);
            error2("");
         }
      }
      fclose(fheredity);
   }

   /*** read and pre-process sequence data at each locus ***/
   for (locus = 0, maxns = 0; locus < data.ngene; ipop += data.ns[locus++]) {
      fprintf(fout, "\n\n*** Locus %d ***\n", locus + 1);
      printf("\n*** Locus %d ***", locus + 1);
      com.cleandata = cleandata0;
      ReadSeq(fout, fseq, cleandata0, locus);

      data.cleandata[locus] = com.cleandata;
      if (data.nseqerr == 0) {
         PatternWeightJC69like();
         fprintf(fout, "\nPrinting out site pattern counts, after collapsing for JC69.\n\n");
         printPatterns(fout);
      }
      if (com.ns > data.maxns) {
         printf("\n%d seqs at locus %d.  More than allowed by the control file.", com.ns, locus + 1);
         exit(-1);
      }
      if (maxns < com.ns)  maxns = com.ns;
      data.ns[locus] = com.ns;
      data.ls[locus] = com.ls;
      data.npatt[locus] = com.npatt;
      data.fpatt[locus] = com.fpatt;
      com.fpatt = NULL;
      for (i = 0; i < com.ns; i++) {
         data.z[locus][i] = com.z[i];  com.z[i] = NULL;
      }

      /* Determine ind ID (j or data.Indnames[j]) and ipop (is) for tips at the locus.
         Count data.nseqsp[locus].  */
      for (i = 0; i < s; i++)  data.nseqsp[locus][i] = 0;
      for (i = 0, nsFull = 0; i < com.ns; i++) {
         if (stree.nspecies == 1)
            is = 0;
         else {
            pline = strchr(com.spname[i], tag);
            if (pline == NULL) error2("sequences must be tagged by population ID");
            sscanf(pline + 1, "%s", iname);
            for (ind = 0; ind < nind; ind++)
               if (strcmp(iname, data.Indnames[ind]) == 0)  break;
            if (ind == nind) {
               printf("Individual label %s not recognised.", iname);
               error2("Please fix the Imap file.");
            }
            else
               is = data.Imap[ind];  /* sequence i is from individual j from species is. */
         }
         ipop[i] = is;  nsFull++;  data.nseqsp[locus][is] ++;
         if (data.diploid[is]) { nsFull++;  data.nseqsp[locus][is] ++; }
         if (noisy >= 3)
            printf("seq %2d %-20s is for indiv %d from species %s  \r", i + 1, com.spname[i], ind + 1, stree.nodes[is].name);
      }  /* for(i, com.ns) */
      if (nsFull <= 1)
         error2("can we have at least two sequences or one diploid sequence at a locus?");

      /* rewrite com.spname[] for the locus if diploid data. */
      if (nsFull > com.ns) {
         for (i = 0; i < stree.nspecies; i++)  data.nseqsp[locus][i] = 0;
         for (i = 0, k = 0; i < com.ns; i++) {
            is = ipop[i];
            for (j = 0; j < (data.diploid[is] ? 2 : 1); j++,k++) {
               data.nseqsp[locus][is] ++;
               sprintf(com.spname[k], "%s.%d", stree.nodes[is].name, data.nseqsp[locus][is]);
            }
         }
      }
      for (is = 0; is < s; is++) {
         if (data.nseqsp[locus][is] > 1 && stree.nseqsp[is] <= 1) {
            printf("\nAt locus %d, species %d has %d seqs.", locus + 1, is + 1, data.nseqsp[locus][is]);
            error2("while control file says no more than one.");
         }
         fprintf(fout, "%s (%2d) ", stree.nodes[is].name, data.nseqsp[locus][is]);
      }
      printf("locus %2d: %2d sequences (", locus + 1, data.ns[locus]);
      for (is = 0; is < s; is++) printf(" %2d", data.nseqsp[locus][is]);
      printf("%-30s", " ), ");
      printf("\n%5d sites, %3d patterns, %s\r", com.ls, com.npatt, (com.cleandata ? "clean" : "messy"));

      /* process diploid data, to change data.ns[locus], data.z[locus], data.npattdiploid[locus] etc. */
      data.diploidlocus[locus] = (nsFull > com.ns);
      if (nsFull > com.ns) {
         /* this changes ns to nsFull and ipop for ns to ipop for nsFull */
         DiploidResolution(fout, locus, ipop); 
         fprintf(fout, "\ncleandata = %d\n", data.cleandata[locus]);
         if (maxns < com.ns)  maxns = com.ns;
      }
   }  /* for(locus) */
   if (noisy >= 3) printf("\n");
   fclose(fseq);

   free(data.Imap);
   free(data.Indnames[0]);

   for (i = 0, com.cleandata = 1; i < data.ngene; i++)
      if (data.cleandata[i] == 0)  com.cleandata = 0;

   data.maxns = maxns;

   fprintf(fout, "\n");
   fflush(fout);
   return(0);
}


/* Ziheng, 14 January 2017.
   Sequences for a species are diploid or haploid depending on data.diploid[ispecies] = 1 or 0.
   The original alignment (A0) at any locus has 1 sequence for each diploid individual and
   uses ambiguity characters YRMKSW to represent heterozygotes (Het).  N-? are allowed and
   treated as ?, but other ambiguities (HBVD) are not allowed.  Ambuiguities are allowed 
   in haploid seqs.
   A0 is compresssed into hetero alignment A1 right after sequences are read.
   DiploidResolution() takes A1 as input.  It copies the number of site patterns in A1 (npatt1)
   to data.npattdiploid[locus] and the pattern frequencies into data.fpatt[locus][].  It then
   resolves some individuals/sites by constructing a H matrix and generates the expanded
   alignment A2, with the number of site patterns to be \sum_i 2^h_i, where i is a site pattern
   in A1 and h_i is its number of Hets.  A2 is then JC-compressed into alignment A3, with the
   number of site patterns in data.npatt[locus].  Columns in A2 are mapped to columns in A3
   through data.pose[locus][]: data.pose[locus][4] = 3 means the fifth pattern in A2 is the
   fourth pattern in A3.
   Calculation of conditional probability vectors uses A3, with data.npatt[locus] site patterns.
   The log probability for the locus sums over site patterns in A1, while the probability for a
   column in A1 averages over 2^h_i resolved site patterns in A2, the probabilities of which are
   retrieved from A3 through data.pose[locus][].
   Alignments A0, A1, A2 are destroyed, although A2 can be reconstructed through A3 and data.pose[locus].
   After this pre-processing, the following are kept in memory:
   (.) data.npattdiploid[locus], data.fpatt[locus][], data.nHdiploid[] refer to A1 (het alignment).
   (.) data.pose[locus][] maps columns in the expanded alignment A2 onto those in A3
   (.) data.z[locus][], com.npatt refer to A3.
*/

int DiploidResolution(FILE *fout, int locus, int ipop[])
{
   /* This function processes alignment A1 if the original data contain diploid sequences, with
      nsFull > data.ns[locus].
      A forever loop resolves one individual per round.  In each one, a loop goes through the
      singleton sites to find the most variable, and given the most variable singleton site,
      find the least variable individual for arbitrary resolution.
      Het alignment A1 in data.z[locus] are encoded before calling this routine.
      nHsite[] = data.nHdiploid[] has Hets for sites in het alignment.
      nHind[] has Hets for individuals at the singleton sites only.
      This assumes nH < 255 at any site, while the number of resolved site patterns at the site is 2^nH.
      isA3: maps sequences in A1 to those in A2.
      z1: holds data (nucleotides) for 1 site in A2.
      ipop is for A1 in inputand is for A2 and A3 in output.
   */
   int ns1 = data.ns[locus], nsFull = 0, npatt1 = data.npatt[locus];  /*ns1 is ns in A1 */
   int nsingle, *singles, indchosen, y;
   int *index, *ispace, i1, i2, is, i, j, k, c, h, ht, site, ip, iround, *het, *isA2;
   double *Hsingles;
   unsigned char *zt[NS], *nHsite, *resolved, *nHind, *z1;
   signed char *H;
   int debug = (noisy <= 3 ? 0 : (data.ns[0] > 10 ? 1 : 9));

   fprintf(fout, "\n\nPreprocessing alignment for diploid data");
   for (is = 0; is < stree.nspecies; is++) nsFull += data.nseqsp[locus][is];
   data.cleandata[locus] = 1;  /* initialized to 1.  may be changed to 0. */
   data.nHdiploid[locus] = (unsigned char*)malloc(npatt1 * sizeof(unsigned char));
   if (data.nHdiploid[locus] == NULL) error2("oom nH");
   H = (signed char*)malloc(ns1*npatt1 * sizeof(signed char));
   Hsingles = (double*)malloc(npatt1 * sizeof(double));
   index = (int*)malloc((3 * npatt1 + 5 * ns1) * sizeof(int));
   if (H == NULL || Hsingles == NULL || index == NULL) error2("oom H ...");
   ispace = index + npatt1;
   singles = ispace + npatt1;
   het = singles + npatt1;
   isA2 = het + ns1;
   resolved = (unsigned char*)(isA2 + ns1);
   nHind = resolved + ns1;
   z1 = nHind + ns1;

   data.npattdiploid[locus] = npatt1;
   nHsite = data.nHdiploid[locus];
   memset(H, 0, ns1*npatt1 * sizeof(signed char));
   memset(nHsite, 0, npatt1 * sizeof(unsigned char));
   for (i = 0; i < ns1; i++)  zt[i] = data.z[locus][i];
   for (i = 0; i < ns1; i++) { resolved[i] = 1;  nHind[i] = 0; }
   for (i = 0, k = 0; i < ns1; i++) { /* seq i in A1 becomes seq k (or k and k+1) in A2. */
      isA2[i] = k++;
      if (data.diploid[ipop[i]] == 1) k++;
   }
   for (i = 0; i < ns1; i++) {
      if (data.diploid[ipop[i]] == 0) {       /* phased sequence */
         if(data.cleandata[locus])
            for (h = 0; h < npatt1; h++) {
               if (zt[i][h] >= 5) { data.cleandata[locus] = 0; break; }
         }
      }
      else
         for (h = 0; h < npatt1; h++) {
            if (zt[i][h] < 5) continue;     /* TCAGU */
            if (zt[i][h] < 11) {            /* YRMKSW */
               H[i*npatt1 + h] = (char)1;
               nHsite[h] ++;
               resolved[i] = 0;
               if (data.fpatt[locus][h] == 1)  nHind[i] ++;   /* only singletons are used  */
            }
            else if (zt[i][h] < 15)         /* HBVD or not N */
               printf("HBVD in sequence %d at locus %d not allowed\n", i + 1, locus + 1);
            else                            /* -N?*/
               data.cleandata[locus] = 0;
         }
   }
   if (debug > 1) {
      printf("\nH matrix, initial");
      for (i1 = 0, FPN(F0); i1 < ns1; i1++, FPN(F0)) {
         for (i2 = 0; i2 < npatt1; i2++)
            printf(" %2d", H[i1*npatt1 + i2]);
         if (resolved[i1]) printf(" *");
      }
   }
   for (h = 0, nsingle = 0; h < npatt1; h++)
      if (nHsite[h] && data.fpatt[locus][h] == 1) {
         singles[nsingle] = h;
         Hsingles[nsingle++] = nHsite[h];
      }

   if (debug > 1) {
      printf("\nH for ind at singleton sites\n");
      for (h = 0; h < ns1; h++) printf(" %4d", nHind[h]);
      printf("\nH for sites\n");
      for (h = 0; h < npatt1; h++) printf(" %4d", nHsite[h]);
   }

   /*  loop to resolve one individual at a time.  No need to resolve nonsingleton sites. */
   for (iround = 0; iround < ns1; iround++) {
      for (i = 0; i < ns1; i++) if (!resolved[i]) break;
      if (i == ns1) break;   /* all indivs are resolved. */
      /* ranking singleton sites by H. */
      indexing(Hsingles, nsingle, index, 1, ispace);
      if (debug > 3) {
         printf("\nRound %2d: singletons & ranks (%2d singletons)", iround + 1, nsingle);
         matIout(F0, singles, 1, nsingle);
         matout2(F0, Hsingles, 1, nsingle, 6, 1);
         matIout(F0, index, 1, nsingle);
      }
      for (i = 0; i < nsingle; i++) { 
         site = singles[index[i]];      /* site is most variable singleton site. */
         for (j = 0, y = npatt1 + 1; j < ns1; j++) {   /* find least variable indiv indchosen at site */
            if (resolved[j] || H[j*npatt1 + site] == 0) continue;
            if (y > nHind[j]) { y = nHind[j];  indchosen = j; }
         }
         if (y < npatt1 + 1) { /* We have found indchosen to resolve a singleton site */
            H[indchosen*npatt1 + site] = -1;
            nHsite[site] --;
            resolved[indchosen] = 1;
            /* update singles list */
            if (--Hsingles[index[i]] == 0) {
               for (j = index[i]; j < nsingle - 1; j++) {
                  singles[j] = singles[j + 1];  Hsingles[j] = Hsingles[j + 1];
               }
               nsingle--;
            }
            break;
         }
      }
      if (debug > 3) {
         printf("\nH matrix");
         for (i1 = 0, FPN(F0); i1 < ns1; i1++, FPN(F0)) {
            for (i2 = 0; i2 < npatt1; i2++)
               printf(" %2d", H[i1*npatt1 + i2]);
            if (resolved[i1]) printf(" *");
         }
      }
      if (i == nsingle) break;  /* found no singletone sites to resolve. */
   }  /* for(iround) */

   /* construt alignment A2 of expanded site patterns */
   for (h = 0; h < npatt1; h++) {
      if (nHsite[h])
         data.npatt[locus] += (1 << nHsite[h]) - 1;
   }
   if (debug) {
      printf("\nH for sites\n");
      for (h = 0; h < npatt1; h++) printf(" %2d", nHsite[h]);
      printf("\nnpatt after expansion = %3d\n", data.npatt[locus]);
   }

   if (fout) {
      fprintf(fout, "\nH matrix\n");
      for (i1 = 0, FPN(F0); i1 < ns1; i1++, FPN(fout)) {
         for (i2 = 0; i2 < npatt1; i2++)
            fprintf(fout, " %2d", H[i1*npatt1 + i2]);
         if (resolved[i1]) fprintf(fout, " *");
      }
      fprintf(fout, "\nHets for sites\n");
      for (h = 0; h < npatt1; h++) fprintf(fout, " %2d", nHsite[h]);
      fprintf(fout, "\nnpatt in A2 after expansion = %3d\n", data.npatt[locus]);
   }

   for (i = 0; i < nsFull; i++) {
      data.z[locus][i] = (unsigned char*)malloc(data.npatt[locus] * sizeof(unsigned char));
      if (data.z[locus][i] == NULL)  error2("oom z ...");
   }
   /* loop through site h in heter alignment (A1) to generate resolved site patterns in A2 */
   for (h = 0, ht = 0; h < npatt1; h++) {
      if (debug > 3) {
         printf("\n\n***Pattern %2d in Het alignment A1: ", h + 1);
         for (i = 0; i < ns1; i++) printf(" %c", BASEs[zt[i][h]]);
      }
      for (i = 0, k = 0; i < ns1; i++) {
         c = zt[i][h];
         if (H[i*npatt1 + h] == 0) {            /* haploid or homo */
            z1[isA2[i]] = c;
            if (data.diploid[ipop[i]] == 1) z1[isA2[i] + 1] = c;
         }
         else if (H[i*npatt1 + h] == -1) {      /* fixed resolution */
            z1[isA2[i]] = (unsigned char)CharaMap[c][0];
            z1[isA2[i] + 1] = (unsigned char)CharaMap[c][1];
         }
         else {
            het[k++] = i;
         }
      }
      if (k != nHsite[h]) error2("jgl");

      /* generate resolved site patterns in A2 */
      for (ip = 0; ip < (1 << nHsite[h]); ip++) {
         if (debug > 3) printf("\nip = %2d: ", ip);
         for (i = 0, j = ip; i < nHsite[h]; i++) {
            k = j % 2;  j /= 2;
            i1 = het[nHsite[h] - 1 - i];
            z1[isA2[i1]] = CharaMap[zt[i1][h]][k];
            z1[isA2[i1] + 1] = CharaMap[zt[i1][h]][1 - k];
            if (debug > 3) printf(" %2d", k);
         }
         for (i = 0; i < nsFull; i++)  data.z[locus][i][ht] = z1[i];
         ht++;

         if (debug > 3) {
            printf("  Site %2d: ", ht);
            for (i = 0; i < nsFull; i++) printf("%2c", BASEs[z1[i]]);
         }
      }  /* for(ip) */
   }     /* for(h) */
   if (ht != data.npatt[locus]) error2("npatt incorrect");
   for (i = 0; i < ns1; i++) free(zt[i]);

   /* set up ipop for A3 */
   for (i = ns1-1, k = 0; i >= 0; i--) {
      ipop[isA2[i]] = ipop[i];
      if (data.diploid[ipop[i]] == 1) ipop[isA2[i]+1] = ipop[i];
   }


   /* setting things up to call PatternWeightJC69like() */
   com.ns = data.ns[locus] = nsFull;
   com.npatt = com.ls = data.npatt[locus];
   for (i = 0; i < com.ns; i++)  com.z[i] = data.z[locus][i];
   com.pose = data.pose[locus] = (int*)malloc(data.npatt[locus] * sizeof(int)); /* this is usedful */
   memset(com.pose, 0, com.npatt * sizeof(int));
   /* com.fpatt is needed by compression routine but counts are not usedful. */
   com.fpatt = (double*)malloc(com.npatt * sizeof(double));
   if (com.fpatt == NULL) error2("oom fpatt");
   for (h = 0; h < com.npatt; h++) com.fpatt[h] = 1;
   if (debug) {
      fprintf(F0, "\n\nAfter expansion, # of site patterns = %d\n", com.npatt);
      if (debug > 3) printPatterns(F0);
   }
   PatternWeightJC69like();
   com.pose = NULL;

   if (fout) {
      fprintf(fout, "\n\nAfter JC compression, # of site patterns in A3 = %d\n", com.npatt);
      printPatterns(fout);
      fprintf(fout, "\nMaps from A2 to A3\n");
      for (h = 0; h < data.npatt[locus]; h++) fprintf(fout, " %2d", data.pose[locus][h]);
   }
   for (i = 0; i < com.ns; i++) {
      data.z[locus][i] = com.z[i];  com.z[i] = NULL;
   }
   data.npatt[locus] = com.npatt;   /* conP calculations will use this npatt. */
   free(com.fpatt);  com.fpatt = NULL;
   free(H);
   free(index);
   free(Hsingles);

   return(0);
}


int GetMem(int ipop[])
{
   /* This routine is called after the species tree and sequence data are read
      to allocate memory for gene trees (gnodes) and conditional probabilities
      (conP).  The size of gnodes[locus] is determined using data.ns[locus],
      no larger than nodes_t (static memory allocated using NS).

      This routine also initializes ipop at tips in gene trees: gnodes[].ipop.
      It is not used by the simulation program.

      The conditional probabilities for internal nodes are com.conPin[2],
      allocated according to data.ns[locus] and data.npatt[locus].
      Two copies of the space are allocated, hence the [2].  The copy used for the
      current gene trees is conPin[com.curconP] while the copy for proposed gene
      trees is the alternative copy conPin[!com.curconP].  Even in the alternative
      copy, the conP space for different loci do not overlap.

      UpdateGB_InternalNode:
      UpdateGB_SPR:
         Those two proposal steps use UseLocus() and AcceptLocus() to copy between the
         two copies of conPin.  They do not change com.curconP.

      UpdateTau:
      mixing:
         Those two steps use SwitchconPin to accept a locus, instead of memcpy to copy
         the whole conPin space.

      Space for heredity scalars and locus rates is allocated in ReadSeqData().
      It would be fine to do it here.
   */
   int s = stree.nspecies, locus, i, j, s1, sizeGtrees;
   double *conP;

   /* get mem for gnodes */
   for (locus = 0, sizeGtrees = 0; locus < data.ngene; locus++)
      sizeGtrees += (data.ns[locus] * 2 - 1) * sizeof(struct TREEN);
   if ((gnodes[0] = (struct TREEN*)malloc(2 * sizeGtrees)) == NULL)
      error2("oom gtree");
   for (locus = 0, nodes = gnodes[0]; locus < data.ngene - 1; locus++)
      gnodes[locus + 1] = (nodes += data.ns[locus] * 2 - 1);
   gnodes_t[0] = gnodes[data.ngene - 1] + data.ns[data.ngene - 1] * 2 - 1;
   for (locus = 0; locus < data.ngene - 1; locus++)
      gnodes_t[locus + 1] = (nodes += data.ns[locus] * 2 - 1);
   for (locus = 0; locus < data.ngene; ipop += data.ns[locus++]) {
      for (j = 0; j < 2 * data.ns[locus] - 1; j++) {
         gnodes[locus][j].nson = 0;
         gnodes[locus][j].age = 0;
         gnodes[locus][j].ibranch = -1;
      }

#ifdef SIMULATION
      for (i = 0, k = 0; i < stree.nspecies; i++) /* ipop for tips in gene tree */
         for (j = 0; j < data.nseqsp[locus][i]; j++)
            gnodes[locus][k++].ipop = i;
#else
      for (i = 0; i < data.ns[locus]; i++)      /* ipop for tips in gene tree */
         gnodes[locus][i].ipop = ipop[i];
#endif
   }
   /* get mem for conP (in nodes) */
   if (mcmc.usedata) {
      data.conP_offset[0] = 0;
      for (locus = 0, com.sconP = 0; locus < data.ngene; locus++) {
         s1 = com.ncode * data.npatt[locus];
         com.sconP += s1*(data.ns[locus] - 1) * sizeof(double);
         if (locus < data.ngene - 1)
            data.conP_offset[locus + 1] = data.conP_offset[locus] + (data.ns[locus] - 1)*s1;
      }
#if (defined(USE_AVX) || defined(USE_SSE))
#if (defined(_WIN32) || defined(_WIN64))
      if ((com.conPin[0] = _aligned_malloc(2 * com.sconP, PLL_BYTE_ALIGNMENT)) == NULL)
         error2("oom conP");
#else
      if (posix_memalign((void **)(&(com.conPin[0])), PLL_BYTE_ALIGNMENT, 2 * com.sconP))
         error2("oom conP");
#endif
#else
      if ((com.conPin[0] = (double*)malloc(2 * com.sconP)) == NULL)
         error2("oom conP");
#endif
      com.conPin[1] = com.conPin[0] + com.sconP / sizeof(double);
      printf("\n%d bytes for conP, %d bytes for gene trees.  \n\n", 2 * (int)com.sconP, sizeGtrees);
      /* set gnodes[locus][].conP for tips and internal nodes */
      com.curconP = 0; conP = com.conPin[0];
      for (locus = 0; locus < data.ngene; locus++) {
         /* Is this still necessary
         if(!data.cleandata[locus])
            UseLocus(locus, 0, 1, 0);
         */
         for (j = data.ns[locus]; j < data.ns[locus] * 2 - 1; j++, conP += com.ncode*data.npatt[locus])
            gnodes[locus][j].conP = conP;
      }
   }

#ifndef SIMULATION
   data.lnpDi = (double*)malloc(data.ngene * sizeof(double));
   data.T2h[0] = (double*)malloc(2 * data.ngene*(s * 2 - 1) * sizeof(double));
   data.ncoal[0] = (int*)malloc(2 * data.ngene*(s * 2 - 1) * sizeof(int));
   if (!stree.NoTheta) {
      data.lnpGi = (double*)malloc(data.ngene * sizeof(double));
   }
   if (data.lnpDi == NULL || data.T2h == NULL || data.ncoal == NULL) error2("oom");

   data.T2h[1] = data.T2h[0] + data.ngene*(s * 2 - 1);
   data.ncoal[1] = data.ncoal[0] + data.ngene*(s * 2 - 1);
#endif

   return(0);
}



void FreeMem(void)
{
   int locus, i, diploiddata = 0;

   free(gnodes[0]);
   if (mcmc.usedata)
   {
#if (defined(USE_AVX) || defined(USE_SSE))
#if (defined(_WIN32) || defined(_WIN64))
      _aligned_free(com.conPin[0]);
#else
      free(com.conPin[0]);
#endif
#else
      free(com.conPin[0]);
#endif
   }
   for (i = 0; i < stree.nspecies; i++) if (data.diploid[i]) diploiddata = 1;
   for (locus = 0; locus < data.ngene; locus++) {
      free(data.fpatt[locus]);
      for (i = 0; i < data.ns[locus]; i++) {
         free(data.z[locus][i]);
      }
      if (diploiddata) {
         free(data.nHdiploid[locus]);
         free(data.pose[locus]);
      }
   }

   free(data.lnpDi);
   free(data.T2h[0]);
   free(data.ncoal[0]);
   if (!stree.NoTheta) free(data.lnpGi);

   if (data.est_heredity || data.est_locusrate)
      free(data.heredity);

   free(stree.M);
}

extern unsigned int z_rndu, w_rndu;


int SaveMCMCstate(char *filename, double lnpG, double lnL)
{
   /* This saves stree, the gene trees, and the MCMC step lengths.
   */
   FILE *f = (FILE*)fopen(filename, "wb");
   int i, sizeGtrees;

   printf("\nSaving MCMC state to file %s for checkpointing\n", filename);
   if (f == NULL) error2("file open error");
   fwrite(&stree, sizeof(stree), 1, f);
   for (i = 0, sizeGtrees = 0; i < data.ngene; i++)
      sizeGtrees += (data.ns[i] * 2 - 1) * sizeof(struct TREEN);
   fwrite(gnodes[0], sizeGtrees, 1, f);
   fwrite(data.root, sizeof(double), data.ngene, f);
   fwrite(mcmc.finetune, sizeof(double), 7, f);
   fwrite(&com.np, sizeof(int), 1, f);
   fwrite(&lnpG, sizeof(double), 1, f);
   fwrite(&lnL, sizeof(double), 1, f);
   fclose(f);
   printf("lnpG = %.6f lnL = %.6f\n", lnpG, lnL);

   return(0);
}

int ReadMCMCstate(char *filename)
{
   /* This reads the MCMC state saved earlier for checkpointing.
   */
   FILE *f = (FILE*)fopen(filename, "rb");
   int i, sizeGtrees;
   double lnpG = 0, lnL = 0;

   printf("\nReading MCMC state from %s\n", filename);
   if (f == NULL) error2("file open error");
   fread(&stree, sizeof(stree), 1, f);
   for (i = 0, sizeGtrees = 0; i < data.ngene; i++)
      sizeGtrees += (data.ns[i] * 2 - 1) * sizeof(struct TREEN);
   fread(gnodes[0], sizeGtrees, 1, f);
   fread(data.root, sizeof(double), data.ngene, f);
   fread(mcmc.finetune, sizeof(double), 7, f);
   fread(&com.np, sizeof(int), 1, f);
   fread(&lnpG, sizeof(double), 1, f);
   fread(&lnL, sizeof(double), 1, f);

   printf("lnpG = %.6f lnL = %.6f\n", lnpG, lnL);

   fclose(f);
   return(0);
}



int copyCoalescent(int direction)
{
   /* direction = 0: 0 -> 1;  1: 1-> 0.
   */
   int s = stree.nspecies;

   memmove(data.T2h[1 - direction], data.T2h[direction], data.ngene*(s * 2 - 1) * sizeof(double));
   memmove(data.ncoal[1 - direction], data.ncoal[direction], data.ngene*(s * 2 - 1) * sizeof(int));
   return(0);
}

int UseLocus(int locus, int copytreeconP, int useData, int setSeqName)
{
   /* This point *nodes to the gene tree at locus gnodes[locus] and set com.z[]
      etc. for likelihood calculation for the locus.
      It also updates stree.nseqsp[], which is used in UpdateGB_SPR() to
      calculate ndesc[].
      If(copytreeconP), the genetree gnodes[locus] is copied to nodes_t.
      If(useData), com. and nodes[].conP are repositioned.
      If (copytreeconP && useData), the conP for internal nodes point to
      the alternative space com.conPin[!com.curconP]+data.conP_offset[locus].
      conP for different loci do not overlap, so that this routine is used by
      all proposal steps, some of which change one locus only while others
      change several or all loci.
   */
   int i, s1 = com.ncode*data.npatt[locus], is, nseqsp[NSPECIES] = { 0 };
   double *conPa = com.conPin[!com.curconP] + data.conP_offset[locus];

   com.ns = data.ns[locus]; tree.root = data.root[locus]; tree.nnode = 2 * com.ns - 1;
   for (i = 0; i < stree.nspecies; i++)
      stree.nseqsp[i] = data.nseqsp[locus][i];

   if (copytreeconP) {
      /* Old genetree is copied into nodes_t.  This is needed if the proposal
      changes genetree topology, but may be unneccesary otherwise.
      Old conP is copied into alternative space.  Think whether this is needed.
      */
      nodes = nodes_t;
      memcpy(nodes_t, gnodes[locus], (com.ns * 2 - 1) * sizeof(struct TREEN));
      if (useData) {
         for (i = com.ns; i < com.ns * 2 - 1; i++)
            nodes_t[i].conP = conPa + (i - com.ns)*s1;
         memcpy(conPa, gnodes[locus][com.ns].conP, s1*(com.ns - 1) * sizeof(double));
      }
   }
   else           /* this works on the old gene tree */
      nodes = gnodes[locus];

   if (useData) {
      com.cleandata = data.cleandata[locus];
      com.npatt = data.npatt[locus];
      com.fpatt = data.fpatt[locus];
      for (i = 0; i < com.ns; i++)
         com.z[i] = data.z[locus][i];
   }
   /*  Check it is o.k. to comment out this block?  14/2/2014.
   for(i=0; i<tree.nnode; i++) {
      if(i!=tree.root) {
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
         if(nodes[i].branch <= 0) {
            printf("father-son age: %.6f %.6f\n", nodes[nodes[i].father].age, nodes[i].age);
            error2("Gtree negative branch length.");
         }
      }
   }
   */
   if (setSeqName) {
      for (i = 0; i < com.ns; i++) {
         is = gnodes[locus][i].ipop;
         nseqsp[is]++;
         if (stree.nseqsp[is] > 1)
            sprintf(com.spname[i], "%s%d", stree.nodes[is].name, nseqsp[is]);
         else
            sprintf(com.spname[i], "%s", stree.nodes[is].name);
      }
   }
   return(0);
}

int AcceptLocus(int locus, int copyconP)
{
   /* If (copytree && usedata), this copies nodes_t[].conP into gnodes[].conP for inner nodes.
      This is used in proposals like UpdateGB_InternalNode() that change one locus only.
      node_t is used for the new gtree, and when the proposal is accepted, node_t is copied over to
      gnodes[locus].  This may not work for proposals that change all loci, such as UpdateTau(),
      UpdateSpeciesSplit(), and UpdateSpeciesJoin(), because the same space of nodes_t may be used
      to store the new gtree for each locus and is overwritten when one cycles through the loci.
   */
   int i, ns = data.ns[locus], s1 = com.ncode*data.npatt[locus];
   double *conPt;

   /* copies node_t into gnodes[locus] */
   data.root[locus] = tree.root;
   memcpy(gnodes[locus], nodes_t, (ns * 2 - 1) * sizeof(struct TREEN));
   if (mcmc.usedata) {  /* reposition gnodes[][].conP */
      conPt = com.conPin[com.curconP] + data.conP_offset[locus];
      for (i = ns; i < ns * 2 - 1; i++)
         gnodes[locus][i].conP = conPt + (i - ns)*s1;
      if (copyconP) {/* used in UpdateGB_InternalNode() and UpdateGBTip() only */
         conPt = com.conPin[!com.curconP] + data.conP_offset[locus];
         memcpy(gnodes[locus][ns].conP, conPt, (ns - 1)*s1 * sizeof(double));
      }
   }
   return(0);
}

void SwitchconPin(void)
{
   /* This resets gnodes[locus].conP to the alternative com.conPin, to avoid
      recalculation of conP, when a proposal is accepted in UpdateTau() or mixing().
   */
   int i, locus;
   double *conP;

   com.curconP = !com.curconP;
   conP = com.conPin[com.curconP];
   for (locus = 0; locus < data.ngene; locus++)
      for (i = data.ns[locus]; i < data.ns[locus] * 2 - 1; i++, conP += com.ncode*data.npatt[locus])
         gnodes[locus][i].conP = conP;
}

void CopyconPin(int locus, int FromCurrent)
{
   /* This copies com.conPin[] from the current place to the alternative place or
      vice versa.
   */
   int size = (data.ns[locus] - 1)*com.ncode*data.npatt[locus] * sizeof(double);
   double *from = com.conPin[com.curconP] + data.conP_offset[locus];
   double *to = com.conPin[!com.curconP] + data.conP_offset[locus], *y;

   if (FromCurrent == 0) {
      y = from; from = to; to = y;
   }
   memcpy(to, from, size);
}


/* use original BPP code if not compiled with SSE/AVX support */
#if (!defined(USE_SSE) && !defined(USE_AVX)) 
static double PMat[4 * 4];

int ConditonalPNode(int inode)
{
   int n = com.ncode, i, j, k, h, ison, ispecies;
   double y;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson > 0 && (!mcmc.saveconP || !com.oldconP[ison]))
         ConditonalPNode(nodes[inode].sons[i]);
   }
   for (h = 0; h < com.npatt*n; h++)
      nodes[inode].conP[h] = 1;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      ispecies = nodes[ison].ipop;
      NPMat++;

      /*
      PMatTN93 (PMat, nodes[ison].branch/3, nodes[ison].branch/3, nodes[ison].branch/3, com.pi);
      */
      pijJC69(PMat, nodes[ison].branch);
      if (nodes[ison].nson > 0) {                          /* internal node */
         for (h = 0; h < com.npatt; h++)
            for (j = 0; j < n; j++) {
               for (k = 0, y = 0; k < n; k++)
                  y += PMat[j != k] * nodes[ison].conP[h*n + k];
               nodes[inode].conP[h*n + j] *= y;
            }
      }
      else if (com.cleandata) {                         /* tip && clean */
         if (data.iseqerr[ispecies] == 0) {
            for (h = 0; h < com.npatt; h++)
               for (j = 0; j < n; j++)
                  nodes[inode].conP[h*n + j] *= PMat[j != com.z[ison][h]];
         }
         else {                                         /* has errors */
            for (h = 0; h < com.npatt; h++)
               for (j = 0; j < n; j++) {
                  for (k = 0, y = 0; k < n; k++)
                     y += PMat[j != k] * data.e_seqerr[ispecies][k * 4 + com.z[ison][h]];
                  nodes[inode].conP[h*n + j] *= y;
               }
         }
      }
      else if (!com.cleandata) {                        /* tip & unclean */
         for (h = 0; h < com.npatt; h++) {
            if (!data.iseqerr[ispecies] || nChara[com.z[ison][h]] > 1) { /* no errors || ambiguity */
               for (j = 0; j < n; j++) {
                  for (k = 0, y = 0; k < nChara[com.z[ison][h]]; k++)
                     y += PMat[j != CharaMap[com.z[ison][h]][k]];
                  nodes[inode].conP[h*n + j] *= y;
               }
            }
            else {                                      /* errors & good base (TCAG) */
               for (j = 0; j < n; j++) {
                  for (k = 0, y = 0; k < n; k++)
                     y += PMat[j != k] * data.e_seqerr[ispecies][k * 4 + CharaMap[com.z[ison][h]][0]];
                  nodes[inode].conP[h*n + j] *= y;
               }
            }
         }
      }
   }   /* for(ison) */
   return (0);
}

/* otherwise, include the corresponding vectorized implementation */
#else
#if (defined(USE_AVX))
#include "pll_avx.c"
#else
#include "pll_sse.c"
#endif
#endif

double lnpData(double lnpDi[])
{
   /* This calculates the log likelihood, the log of the probability of the data
      given gtree[nlocus] and blength[nlocus][2] for each locus.
      This updates gnodes[locus][].conP.
   */
   int j, locus;
   double lnL = 0, y;

   if (mcmc.saveconP)
      for (j = 0; j < data.maxns * 2 - 1; j++)  com.oldconP[j] = 0;
   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 0, 1, 0);
      y = lnpD_locus(locus);
      if (testlnL && fabs(lnpDi[locus] - y) > 1e-10)
         printf("\nlnLi %.6f != %.6f at locus %d\n", lnpDi[locus], y, locus + 1);
      lnpDi[locus] = y;
      lnL += y;
   }
   return(lnL);
}

double lnpD_locus(int locus)
{
   /* This calculates ln{D|G, t} for locus, using nodes[].age and tree.root.
      Note that nodes[].branch may not be kept up to date in the program.
      When rates vary among loci, the relative rate is applied to branches but not ages.
   */
   int  h, h2, h3, npatth2, i, j, haserr = 0;
   double lnL = 0, fh, *fh3, *p;

   if (!mcmc.usedata) return(0);
   for (i = 0; i < tree.nnode; i++) {
      if (i == tree.root) continue;
      nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
      if (data.est_locusrate)
         nodes[i].branch *= data.locusrate[locus];
      if (nodes[i].branch < 0) {
         printf("branch length = %.9g < 0\n", nodes[i].branch);
         if (nodes[i].branch < -1e-5) {
            printGtree(1);
            exit(-1);
         }
      }
   }
   ConditonalPNode(tree.root);
   if (data.diploidlocus[locus]) {
      fh3 = (double*)malloc(data.npatt[locus] * sizeof(double));
      if (fh3 == NULL) error2("oom fh3");
      memset(fh3, 0, data.npatt[locus] * sizeof(double));
      p = nodes[tree.root].conP;
      for (h = 0; h < data.npatt[locus]; h++, p += 4)           /* site columns in A3 */
         fh3[h] = (*p + *(p + 1) + *(p + 2) + *(p + 3)) / 4;    /* average over pi */
      for (h = 0, h2 = 0; h < data.npattdiploid[locus]; h++) {  /* site columns in A1 (Ahet) */
         npatth2 = (1 << data.nHdiploid[locus][h]);
         for (j = 0, fh = 0; j < npatth2; j++, h2++) {
            h3 = data.pose[locus][h2];
            fh += fh3[h3];
            if (debug) {
               for (i = 0; i < data.ns[locus]; i++) printf("%c ", BASEs[data.z[locus][i][h3]]);
               printf(" h123 = %2d %2d %2d %12.9f H = %2d\n", h, h2, h3, fh, data.nHdiploid[locus][h]);
            }
         }
         fh /= npatth2;
         if (fh < 1E-300) {
            printGtree(1);
            if (noisy) {
               printf("\nsmall pr(Data) = %.9f at locus %d site pattern %d\n", fh, locus + 1, h + 1);
               printf("cleandata = %d\n", data.cleandata[locus]);
               print1site(F0, h);
            }
            fh = -500;

            printf("\n");
            for (j = 0; j < data.ngene; j++) {
               printf(" %d", data.cleandata[j]);
               if ((j + 1) % 20 == 0) printf("\n");
            }
            printf("\n");
            exit(-1);
         }
         else
            fh = log(fh);
         lnL += com.fpatt[h] * fh;
         if (debug) printf("%50s h = %2d: ln %12.9f x %3.0f\n", "", h, exp(fh), com.fpatt[h]);
      }
      free(fh3);
      /* Here h2 should equal the number of site patterns in A2 */
   }
   else {
      p = nodes[tree.root].conP;
      for (h = 0; h < data.npatt[locus]; h++, p += 4) {
         fh = (*p + *(p + 1) + *(p + 2) + *(p + 3)) / 4;    /* average over pi */
         if (fh < 1E-300) {
            if (noisy) printf("\nsmall pr(Data) = %.9g at locus %d site pattern %d\n", fh, locus + 1, h + 1);
            fh = -500;
         }
         else
            fh = log(fh);
         lnL += com.fpatt[h] * fh;
      }
   }

   return (lnL*BFbeta);
}


double CountCoalescentLocus(int locus)
{
   /* This counts the number of coalescent events in each population for locus.
      tj[ipop*(ns-1)] has node ages for ns-1 coalescent evens in ipop.
   */
   int s = stree.nspecies, ns = data.ns[locus], i, j, k, inode, ip, ipop, ipopdad, nt;
   double lnp = 0, starttime, endtime, *t, y, H = (data.est_heredity == 0 ? 1 : data.heredity[locus]);
   double T2h, a = data.theta_prior[0], b = data.theta_prior[1];
#if(0)
   int *nin, *ncoal, *nout;
   double *tj;

   k = (s * 2 - 1)*(ns - 1) * sizeof(double) + 3 * (s * 2 - 1) * sizeof(int);
   if ((tj = (double*)malloc(k)) == NULL) error2("oom");
   nin = (int*)(tj + (s * 2 - 1)*(ns - 1));  ncoal = nin + s * 2 - 1;  nout = ncoal + s * 2 - 1;
   memset(nin, 0, 3 * (s * 2 - 1) * sizeof(int));
#else
   int nin[NSPECIES * 2 - 1], ncoal[NSPECIES * 2 - 1], nout[NSPECIES * 2 - 1];
   double tj[(NSPECIES * 2 - 1)*(NS - 1)];

   for (i = 0; i < (s * 2 - 1); i++)  nin[i] = ncoal[i] = nout[i] = 0;
#endif

   if (debug == 1) puts("\n**** In CountCoalescentLocus ****");
   memset(data.T2h[0] + (s * 2 - 1)*locus, 0, (s * 2 - 1) * sizeof(double));
   memset(data.ncoal[0] + (s * 2 - 1)*locus, 0, (s * 2 - 1) * sizeof(int));

   for (i = 0; i < stree.nspecies; i++) nin[i] = data.nseqsp[locus][i];

   for (inode = 0; inode < tree.nnode; inode++) {
      ipop = nodes[inode].ipop;
      if (inode >= ns)    /* internal node */
         tj[ipop*(ns - 1) + ncoal[ipop]++] = nodes[inode].age;
      ipopdad = (inode == tree.root ? -1 : nodes[nodes[inode].father].ipop);
      for (; ipop != ipopdad; ) {
         nout[ipop]++;
         ipop = stree.nodes[ipop].father;
         if (ipop == -1 || inode == tree.root) break;
         nin[ipop] ++;
      }
   }

   for (ip = 0; ip < stree.npop; ip++) { /* sort tj for each pop */
      ipop = stree.pops[ip];
      if (nin[ipop] - nout[ipop] != ncoal[ipop]) {
         printf("\nin out coal %d %d %d in ipop %d\n", nin[ipop], nout[ipop], ncoal[ipop], ipop);
         printGtree(1);
         error2("in out ");
      }
      nt = ncoal[ipop];
      if (nt > 0) t = tj + ipop*(ns - 1);
      if (nt > 10)
         qsort(t, (size_t)nt, sizeof(double), (int(*)(const void *, const void *))comparedouble);
      else
         for (i = 0; i < nt - 1; i++)
            for (j = i + 1; j < nt; j++)
               if (t[j] < t[i]) { y = t[i]; t[i] = t[j]; t[j] = y; }
   }

   if (debug == 1) {
      printf("\ncoalescence info in populations:\nspecies  in coal  out   tau     tj\n");
      for (i = 0; i < stree.nnode; i++, FPN(F0)) {
         printf("%4d: %5d%5d%5d %9.5f", i, nin[i], ncoal[i], nout[i], stree.nodes[i].age);
         for (j = 0; j < ncoal[i]; j++) printf(" %9.5f", tj[i*(ns - 1) + j]);
      }
      FPN(F0);
   }

   for (ip = 0; ip < stree.npop; ip++) {
      ipop = stree.pops[ip];
      if (nin[ipop] <= 1) continue;   /* no coalescent possible in ipop */
      nt = ncoal[ipop];  T2h = 0;  t = tj + ipop*(ns - 1);
      starttime = stree.nodes[ipop].age;
      if (ipop != stree.root) endtime = stree.nodes[stree.nodes[ipop].father].age;
      else                  endtime = -1;

      if (debug == 1) printf("species %d: tau range (%8.5f, %8.5f) ", ipop, starttime, endtime);

      if (ncoal[ipop]) {
         for (j = nin[ipop], k = 0; j > nout[ipop]; j--, k++) {
            y = t[k] - (k == 0 ? starttime : t[k - 1]);
            T2h += j*(j - 1)*y / H;

            if (debug == 1) printf(" j=%d: %7.5f", j, y);

            if (y < 0 && noisy) {
               printGtree(1);
               printf("\n       tj = %.20f < 0 in ipop %d!\n", y, ipop);
               printf(" node age = %.20f\n", t[k]);
               printf("starttime = %.20f\n", starttime);
               error2("negative time in CountCoalescentLocus()");
            }
         }
      }
      if (nout[ipop] > 1) {  /* remaining time in ipop */
         y = (nout[ipop] == nin[ipop] ? endtime - starttime : endtime - t[k - 1]);
         T2h += nout[ipop] * (nout[ipop] - 1)*y / H;

         if (debug == 1) printf(" remaining time %.5f", y);
      }
      if (debug == 1) FPN(F0);

      data.ncoal[0][(s * 2 - 1)*locus + ipop] = nt;
      data.T2h[0][(s * 2 - 1)*locus + ipop] = T2h;
   }   /* for(ip)  */

#if(0)
   free(tj);
#endif
   return(0);
}


double lnpG_S(int locus)
{
   /* this calculates the probability of gene tree and coalescent times (using tree
      and nodes[]), given stree (theta and tau).  This is used if(stree.NoTheta == 0).
   */
   int s = stree.nspecies, is, ip, nt;
   double lnp = 0, T2h, a = data.theta_prior[0], b = data.theta_prior[1];

   if (debug == 1) puts("\n**** In lnpG_S ****");

   for (ip = 0; ip < stree.npop; ip++) {
      is = stree.pops[ip];
      nt = data.ncoal[0][(s * 2 - 1)*locus + is];
      T2h = data.T2h[0][(s * 2 - 1)*locus + is];
      if (nt > 0)   lnp += nt * log(2 / stree.nodes[is].theta);
      if (T2h > 0)  lnp -= T2h / stree.nodes[is].theta;
   }

   return(lnp);
}

double lnpG_S_NoTheta(void)
{
   /* this calculates the probability of gene tree and coalescent times (using tree
      and nodes[]), given stree (theta and tau), with theta integrated out.
   */
   int s = stree.nspecies, is, ip, i, nt;
   double lnp = 0, T2h, a = data.theta_prior[0], b = data.theta_prior[1];
   double alnb, lnGa, y;

   if (debug == 1) puts("\n**** In lnpG_S_NoTheta ****");

   alnb = a*log(b);  lnGa = LnGamma(a);
   for (i = 0; i < data.ngene; i++) lnp += (data.ns[i] - 1)*0.6931471805599453;
   if (data.est_heredity) {
      for (i = 0, y = 1; i < data.ngene; i++) {
         y *= data.heredity[i];
         if ((i + 1) % 100 == 0) { lnp -= log(y);  y = 1; }
      }
      lnp -= log(y);
   }

   for (ip = 0; ip < stree.npop; ip++) {
      is = stree.pops[ip];
      nt = 0;  T2h = 0;
      for (i = 0; i < data.ngene; i++) {
         nt += data.ncoal[0][(s * 2 - 1)*i + is];
         T2h += data.T2h[0][(s * 2 - 1)*i + is];
      }
      if (nt == 0)
         lnp -= a * log(1 + T2h / b);
      else
         lnp += alnb - lnGa - (a + nt)*log(b + T2h) + LnGamma(a + nt);
   }

   return(lnp);
}


double UpdateGB_InternalNode(double* lnL, double* lnpG, double finetune)
{
   /* This slides a node in the gene tree without changing the gene tree topology.  The new
      coalescent time tnew may be in a different population from the current, and nodes[inode].ipop
      may change.  The algorithm loops through all internal nodes in the gene tree at each locus.
      The lower bound tb[0] is determined by age of the two sons and also the age of the species
      that is a common ancestor to both sons.  The upper bound tb[1] is the age of father node.
   */
   int  accepted = 0, locus, i, j, inode, sons[2], pops[2], ninodes;
   int  ipopsource, ipoptarget, copytree;
   double lnacceptance, lnLd, lnpGd, lnpGinew, lnpDinew, t, tnew, tb[2];

   if (debug == 2) puts("\nUpdateGB_InternalNode");
   if (finetune <= 0) error2("steplength = 0 in UpdateGB_InternalNode.  Check finetune");
   for (i = 0, ninodes = 0; i < data.ngene; i++) ninodes += data.ns[i] - 1;

   for (locus = 0; locus < data.ngene; locus++) {
      if (debug == 2) {
         printf("\nlocus %d (ns = %d)\n", locus + 1, data.ns[locus]);
         UseLocus(locus, 1, mcmc.usedata, 1);
         printGtree(1);
      }

      for (inode = data.ns[locus], copytree = 1; inode < data.ns[locus] * 2 - 1; inode++) {
         if (copytree) UseLocus(locus, 1, mcmc.usedata, 0);  /* work on a copy */
         copyCoalescent(0);
         t = nodes[inode].age;
         ipopsource = nodes[inode].ipop;
         for (i = 0; i < 2; i++) {
            sons[i] = nodes[inode].sons[i];
            pops[i] = nodes[sons[i]].ipop;
         }
         tb[0] = max2(nodes[sons[0]].age, nodes[sons[1]].age);
         tb[1] = (inode == tree.root ? OLDAGE : nodes[nodes[inode].father].age);

         /* Find the population ancestral to both sons */
         if (pops[0] != pops[1]) { /* no change if(pops[0] == pops[1]). */
            for (i = pops[0]; i != stree.root; i = stree.nodes[i].father) {
               for (j = pops[1]; j != stree.root; j = stree.nodes[j].father)
                  if (i == j) break;
               if (j != stree.root) break;
            }
            tb[0] = max2(tb[0], stree.nodes[i].age);
            if (tb[0] > tb[1])
               printf("t bound wrong: %9.5f > %9.5f", tb[0], tb[1]);
         }
         tnew = t + finetune*rndSymmetrical();
         tnew = reflect(tnew, tb[0], tb[1]);

         /* determine ipoptarget by going up stree */
         for (i = ipoptarget = pops[0]; ; i = stree.nodes[i].father) {
            if (stree.nodes[i].age > tnew) break;
            ipoptarget = i;
            if (i == stree.root) break;
         }

         /* two things are updated in the loop: node age and ipop */
         nodes[inode].age = tnew;
         nodes[inode].ipop = ipoptarget;

         if (mcmc.saveconP) {
            for (i = com.ns; i < com.ns * 2 - 1; i++) com.oldconP[i] = 1;
            for (i = inode; ; i = nodes[i].father)
            {
               com.oldconP[i] = 0; if (i == tree.root) break;
            }
         }
         CountCoalescentLocus(locus);
         if (stree.NoTheta)
            lnpGd = lnpG_S_NoTheta() - *lnpG;
         else {
            lnpGinew = lnpG_S(locus);
            lnpGd = lnpGinew - data.lnpGi[locus];
         }

         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance = lnpGd + lnLd;

         if (debug == 2) {
            printf("node %2d tb (%8.5f %8.5f) t %7.5f -> %7.5f lnLd %9.5f ", inode, tb[0], tb[1], t, tnew, lnLd);
            for (i = com.ns; i < com.ns * 2 - 1; i++) printf("%d", com.oldconP[i]);
         }

         if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
            accepted++;
            *lnL += lnLd;
            data.lnpDi[locus] = lnpDinew;
            *lnpG += lnpGd;
            if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew;
            AcceptLocus(locus, 1);
            copytree = 0;
            if (debug == 2) printf(" Y\n");
         }
         else {
            nodes[inode].age = t;
            nodes[inode].ipop = ipopsource;
            copyCoalescent(1);
            copytree = 1;
            if (debug == 2) printf(" N\n");
         }
      }
   }
   return((double)accepted / ninodes);
}


void GraftNode(int source, int target, double age, int ipop)
{
   /* This prunes off source branch on the gene tree (branch ancestral to node source) and
      regrafts it to branch target.  It assumes binary rooted tree.
      age & ipop are for the new node (father) to be created on the target branch.
      After the operation, the following branches are generated:
         grandpa->sib (at the source)
         father->source, father->target, targetfather->father (at the target)

      IMPORTANT: This keeps the node number for tree.root unchanged.
   */
   int i, k, father, grandpa = -1, sib, targetfather, oldroot = tree.root;
   double y;

   father = nodes[source].father;
   nodes[father].age = age;
   nodes[father].ipop = ipop;
   if (father != tree.root) grandpa = nodes[father].father;
   sib = nodes[father].sons[0] + nodes[father].sons[1] - source;
   targetfather = nodes[target].father;

   if (target != sib && target != father) { /* change to tree topology */
      /* cut at source */
      if (father == tree.root) {         /* old root is lost */
         tree.root = sib;  nodes[sib].father = -1;
      }
      else {
         if (nodes[grandpa].sons[i = 0] != father)
            i = 1;
         nodes[grandpa].sons[i] = sib;
         nodes[sib].father = grandpa;
      }

      /* regraft */
      nodes[father].father = nodes[target].father;
      nodes[target].father = father;
      nodes[father].sons[0] = target;
      nodes[father].sons[1] = source;
      if (target == tree.root) {  /* new root is born */
         tree.root = father;
      }
      else
         for (i = 0; i < 2; i++) {
            if (nodes[targetfather].sons[i] == target) {
               nodes[targetfather].sons[i] = father;
               break;
            }
         }
   }
   if (tree.root != oldroot) { /* swap everything except .conP */
      /* if(debug==3)  printGtree(1); */
      swap2(nodes[tree.root].ipop, nodes[oldroot].ipop, k);
      swap2(nodes[tree.root].age, nodes[oldroot].age, y);
      for (i = 0; i < 2; i++) {
         nodes[nodes[tree.root].sons[i]].father = oldroot;
         nodes[nodes[oldroot].sons[i]].father = tree.root;
      }
      for (i = 0; i < 2; i++)
         swap2(nodes[tree.root].sons[i], nodes[oldroot].sons[i], k);
      for (i = 0; i < 2; i++)
         if (nodes[oldroot].sons[i] == oldroot)
            nodes[oldroot].sons[i] = tree.root;
      k = nodes[oldroot].father;
      for (i = 0; i < 2; i++)
         if (nodes[k].sons[i] == oldroot)
            nodes[k].sons[i] = tree.root;
      if ((nodes[tree.root].father = k) == tree.root)
         nodes[tree.root].father = oldroot;
      nodes[oldroot].father = -1;

      father = nodes[source].father;
      if (grandpa == tree.root)    grandpa = oldroot;
      else if (grandpa == oldroot) grandpa = tree.root;
      tree.root = oldroot;
   }
   if (mcmc.usedata && mcmc.saveconP) {
      for (i = com.ns; i < com.ns * 2 - 1; i++) com.oldconP[i] = 1;
      for (i = father; ; i = nodes[i].father) {
         com.oldconP[i] = 0; if (i == tree.root) break;
      }
      for (i = grandpa; i != -1; i = nodes[i].father) {
         com.oldconP[i] = 0; if (i == tree.root) break;
      }
   }
}

double UpdateGB_SPR(double* lnL, double* lnpG, double finetune)
{
   /* This removes a node (inode and its father) in the gene tree and regrafts it to a feasible branch
      (fig. A1A in RY2003).  There is no change to the node times inside the subtree, but the age of
      the father node (t) is changed.  The gene tree topology may not change (for example, if the
      pruned subtree is a daughter of root).  The algorithm cycles through all loci and for each locus
      all nodes except the root.  sizeSnode[is] keeps the number of descendent tips of population is
      (internal node on the species tree).  For example, sizeSnode[stree.root]==com.ns.  This is used to
      deal with loci at which some species are missing to identify feasible populations (species) to
      determine the min bound tb[0] for tnew.

      Note that inode is the root node for the subtree being moved, species is pop id for inode, and
      t and tnew are the time at which the subtree joins the rest of the gtree, that is, the age of
      inode's father.
      The lower bound tb[0] for tnew is determined by the following requirements:
      (a) tnew must be older than the age of the subtree: nodes[inode].age.
      (b) The population of the target node (ipop) must be ancestral to or the same as the population
          of inode: (stree.pptable[species][ipop]==1).
      (c) There must be other lineages outside the subtree passing that target population to which the
          subtree can join.
   */
   int accepted = 0, nproposal = 0, nsource, ntarget, targets[NS * 2 - 1], target, copytree;
   int locus, is, i, j, k, inode, sizeinode, species, ipop, ipopsource, ipoptarget, father, sib;
   int s = stree.nspecies, sizeSnode[2 * NSPECIES - 1];  /* #descendents for ancestral species */
   double lnacceptance = 0, lnLd, lnpGd, lnpGinew, lnpDinew;
   double tb[2], t, tnew;

   if (debug == 3) puts("\nUpdateGB_SPR ");
   if (finetune <= 0) error2("steplength = 0 in UpdateGB_InternalNode.  Check finetune");

   for (locus = 0; locus < data.ngene; locus++) {
      //if(data.ns[locus]==2) continue;  /* no need to update the gene-tree topology */
      nproposal += (mcmc.moveinnode ? data.ns[locus] * 2 - 2 : data.ns[locus] - 1);
      UseLocus(locus, 0, 0, 1);

      if (debug == 3) {
         printf("\nlocus %d (ns = %d)\n", locus + 1, data.ns[locus]);
         for (is = s; is < 2 * s - 1; is++)
            printf("sp.time[%d] = %.6f\n", is, stree.nodes[is].age);
         printGtree(1);
      }
      for (i = 0; i < 2 * s - 1; i++) sizeSnode[i] = 0;
      for (species = 0; species < s; species++) { /* calculate sizeSnode[] for the locus */
         k = data.nseqsp[locus][species];
         for (i = species; ; i = stree.nodes[i].father) {
            sizeSnode[i] += k;
            if (i == stree.root) break;
         }
      }
      if (sizeSnode[stree.root] != com.ns) error2("sizeSnode");

      for (inode = 0, copytree = 1; inode < (mcmc.moveinnode ? com.ns * 2 - 1 : com.ns); inode++) {
         /* NOTE: tree.root should not be changed in GraftNode() */
         if (inode == tree.root) continue;
         copyCoalescent(0);
         if (copytree) UseLocus(locus, 1, mcmc.usedata, 0);
         species = nodes[inode].ipop;
         father = nodes[inode].father;
         sib = nodes[father].sons[0] + nodes[father].sons[1] - inode;
         ipopsource = nodes[father].ipop;
         t = nodes[father].age;

         /* find lower bound tb[0] for tnew. */
         tb[0] = nodes[inode].age;
         tb[1] = OLDAGE;
         for (i = 0, sizeinode = 0; i < com.ns; i++) {
            for (j = i; j != tree.root; j = nodes[j].father)
               if (j == inode) { sizeinode++;  break; }
         }
         for (is = species; ; is = stree.nodes[is].father)
            if (sizeSnode[is] > sizeinode) break;
         tb[0] = max2(tb[0], stree.nodes[is].age);

         tnew = t + finetune*rndSymmetrical();
         tnew = reflect(tnew, tb[0], tb[1]);

         /* identify the target pop in which tnew is, by going up stree */
         for (ipop = ipoptarget = species; ; ipop = stree.nodes[ipop].father) {
            if (stree.nodes[ipop].age > tnew) break;
            ipoptarget = ipop;
            if (ipop == stree.root) break;
         }

         if (debug == 3) printf("inode %2d father %2d tb (%6.4f %7.4f) time %7.5f (%2d) ->%7.5f (%2d) ",
            inode, father, tb[0], tb[1], t, ipopsource, tnew, ipoptarget);

         /* count and identify feasible target branches. Watch out for root */
         ntarget = 0;  nsource = 1;
         if (tnew >= nodes[tree.root].age)
            targets[ntarget++] = tree.root;
         else
            for (i = 0; i < tree.nnode; i++) { /* Is node i a possible target? */
               if (i != inode && i != tree.root && nodes[i].age <= tnew && nodes[nodes[i].father].age > tnew
                  && stree.pptable[nodes[i].ipop][ipoptarget])
                  targets[ntarget++] = (i == father ? sib : i);
            }
         if (father != tree.root)
            for (i = 0; i < tree.nnode; i++) { /* Is node i a possible source? */
               if (i != inode && i != tree.root && i != sib && i != father
                  && nodes[i].age <= t && nodes[nodes[i].father].age > t
                  && stree.pptable[nodes[i].ipop][ipopsource])
                  nsource++;
            }

         if (nsource < 1 || ntarget < 1) {
            printGtree(1);
            printf("\nlocus %d node %d", locus, inode);
            printf("\nerror: nsource = %d ntarget = %d", nsource, ntarget);
            exit(-1);
         }
         target = targets[(int)(ntarget*rndu())];
         GraftNode(inode, target, tnew, ipoptarget);

         if (debug == 3) {
            printf("line %d > %2d (", nsource, ntarget);
            for (i = 0; i < ntarget; i++) printf(" %2d", targets[i]);
            printf(") > %2d ", target);
            printGtree(1);
         }
         lnacceptance = log((double)ntarget / nsource);
         CountCoalescentLocus(locus);
         if (stree.NoTheta)
            lnpGd = lnpG_S_NoTheta() - *lnpG;
         else {
            lnpGinew = lnpG_S(locus);
            lnpGd = lnpGinew - data.lnpGi[locus];
         }

         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance += lnLd + lnpGd;

         if (debug == 3) {
            printf(" lnLd %9.5f ", lnLd);
            for (i = com.ns; i < com.ns * 2 - 1; i++) printf("%d", com.oldconP[i]);
         }

         if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
            accepted++;
            *lnL += lnLd;
            data.lnpDi[locus] = lnpDinew;
            *lnpG += lnpGd;
            if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew;
            AcceptLocus(locus, 1);
            copytree = 0;
            if (debug == 3) printf(" Y (%4d)\n", NPMat);
         }
         else {
            copyCoalescent(1);
            copytree = 1;
            if (debug == 3) printf(" N (%4d)\n", NPMat);
         }
      }
   }
   return(nproposal ? (double)accepted / nproposal : 0);
}



int ReScaleSubTree(int inode, double factor)
{
   /* this multiplies all node ages inside the subtree at inode by factor.
   */
   int nnodesScaled = 0, i;

   if (nodes[inode].age == 0) return nnodesScaled;
   nnodesScaled++;
   nodes[inode].age *= factor;
   for (i = 0; i < nodes[inode].nson; i++)
      if (nodes[inode].sons[i] >= com.ns)
         nnodesScaled += ReScaleSubTree(nodes[inode].sons[i], factor);
   return (nnodesScaled);
}

int NodeSlider_NotUsed(double eps)
{
   /* This is not yet used or tested.  The code for generating the candidate tree (reflecting)
      is tested and seemed correct, but one has to add the code for the proposal ratio etc.
      to use it as a MCMC proposal.
      This algorithm perturbs the tree by sliding a node.  Right now it is written
      for rooted trees only, and assumes that the tree is binary,  Some changes (not
      done) are necessary to deal with unrooted trees.
      Right now this does not consider the constraints placed on the gene tree by the
      species tree.
      It does not seem a good idea to cycle through the nodes on the tree, because the
      node numbers are changed during each cycle.  Instead it may be better to pick
      up a node at random for sliding.
   */
   int j, inode, dad0, dad, target, iround;
   double offset = 0, tnew, factor;

   debug = 11;
   if (nodes[tree.root].nson != 2)
      error2("NodeSlider: this is not used or tested.  This is for rooted tree only.");
   for (inode = 0; inode < 2 * com.ns - 1; inode++) {
      if (inode == tree.root) continue;
      offset = eps * rndSymmetrical();

      if (debug) printf("\n\nnode %d offset %9.5f\n", inode, offset);

      dad0 = dad = nodes[inode].father;
      target = nodes[dad0].sons[0] + nodes[dad0].sons[1] - inode;  /* sib */
      offset += nodes[dad0].age - nodes[target].age;
      for (iround = 0; ; iround++) {  /* forever reflecting */
         if (offset > 0) {
            dad = nodes[target].father;
            if (dad == dad0)
               dad = nodes[dad].father;
            if (target == tree.root)
               break;
            tnew = nodes[target].age + offset;
            if (dad == -1 || tnew < nodes[dad].age)
               break;
            else if (rndu() < 0.5) {  /* up to dad */
               offset -= nodes[dad].age - nodes[target].age;
               target = dad;
            }
            else {                 /* down to sib */
               j = nodes[dad].sons[0];
               if (j == target || j == nodes[target].father)
                  j = nodes[dad].sons[1];
               target = j;
               offset = 2 * nodes[dad].age - nodes[target].age - tnew;
            }
         }
         else {
            if (target < com.ns) /* tip */
               offset *= -1;
            else {
               target = nodes[target].sons[rndu() < 0.5];
               offset += nodes[nodes[target].father].age - nodes[target].age;
            }
         }

         if (debug) printf("round %d target %2d offset %9.5f\n", iround + 1, target, offset);

      }  /* forever loop */

      factor = (nodes[target].age + offset) / nodes[dad0].age;
      nodes[dad0].age = nodes[target].age + offset;
      /* rescale the branch lengths */
      if (inode >= com.ns) {
         ReScaleSubTree(inode, factor);
         if (debug) printf("factor = %9.5f\n", factor);
      }
      if (nodes[target].father != dad0) {
         /* ipop is not updated and the algorithm does not work for the model anyway. */
         GraftNode(inode, target, nodes[target].age + offset, 0);
      }

      if (debug) printGtree(1);
   }  /*  for(inode) */
   return(0);
}


double UpdateTheta(double *lnpG, double finetune, double space[])
{
   /* This updates theta's one by one, using proportional expansion or shrinkage.
      Perhaps change to sliding window to avoid getting stuck at 0.
      This step does not require calculation of the likelihood.
   */
   int i, ipop, locus, accepted = 0, ntheta = 0;
   double thetaold, thetanew, lnpGd, lnacceptance;
   double *lnpGinew = space, a = data.theta_prior[0], b = data.theta_prior[1];

   if (finetune <= 0) error2("steplength = 0 in UpdateTheta.  Check finetune");
   for (i = 0; i < stree.npop; i++) {
      /* prior and proposal ratios */
      ipop = stree.pops[i];
      if (stree.nodes[ipop].theta < 0) continue;
      copyCoalescent(0);
      ntheta++;
      thetaold = stree.nodes[ipop].theta;

      thetanew = thetaold + finetune*rndSymmetrical();
      if (thetanew < 0) thetanew = -thetanew;
      stree.nodes[ipop].theta = thetanew;

      lnacceptance = (-a - 1)*log(thetanew / thetaold) - b*(1 / thetanew - 1 / thetaold);

      for (locus = 0, lnpGd = 0; locus < data.ngene; locus++) {
         UseLocus(locus, 1, 0, 0);
         CountCoalescentLocus(locus);
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
      lnacceptance += lnpGd;
      if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
         for (locus = 0; locus < data.ngene; locus++)
            data.lnpGi[locus] = lnpGinew[locus];
         *lnpG += lnpGd;
         accepted++;
      }
      else {
         stree.nodes[ipop].theta = thetaold;
         copyCoalescent(1);
      }
   }
   return((double)accepted / ntheta);
}


double UpdateTau (double *lnL, double *lnpG, double finetune, double space[])
{
   /* This updates speciation times tau using the rubber-band algorithm.
      ntj[] are counts of nodes below and above tau (m and n in the paper).
      changetheta=1 and 2 are two algorithms changing theta at the same time.
   */
   int changetheta = 2, k, is, inode, locus, accepted = 0, ntau = 0;
   int ntj[2], ntj_locus[2];
   double lnacceptance = 0, lnLd, lnpGd, *lnpGinew = space, *lnpDinew = space + data.ngene;
   double a = data.tau_prior[0], b = data.tau_prior[1], tauold, taunew, taub[2], t, taufactor[2];
   double atheta = data.theta_prior[0], btheta = data.theta_prior[1], thetaold, thetanew;

   if (debug == 5) puts("\nUpdateTau...");
   if (finetune <= 0) error2("steplength = 0 in UpdateTimes.  Check finetune");
   if (stree.NoTheta) changetheta = 0;
   for (is = stree.nspecies; is < stree.nnode; is++)
      if (stree.nodes[is].age > 0) ntau++;

   for (is = stree.nspecies; is < stree.nnode; is++) {
      if (stree.nodes[is].age == 0) continue;
      copyCoalescent(0);
      tauold = stree.nodes[is].age;
      lnLd = 0;  lnpGd = 0;
      taub[0] = 0;
      taub[1] = OLDAGE;
      if (stree.nodes[is].nson)
         taub[0] = max2(stree.nodes[stree.nodes[is].sons[0]].age,
            stree.nodes[stree.nodes[is].sons[1]].age);
      if (is != stree.root)
         taub[1] = stree.nodes[stree.nodes[is].father].age;

      taunew = tauold + finetune*rndSymmetrical();
      taunew = stree.nodes[is].age = reflect(taunew, taub[0], taub[1]);
      for (k = 0; k < 2; k++)
         taufactor[k] = (taunew - taub[k]) / (tauold - taub[k]);
      lnacceptance = 0;
      if (is == stree.root) /* tau_root ~ invgamma(a, b) */
         lnacceptance = (-a - 1 - ntau + 1)*log(taunew / tauold) - b*(1 / taunew - 1 / tauold);

      if (changetheta) {
         thetaold = stree.nodes[is].theta;
         if (changetheta == 1)       t = taunew / tauold;
         else if (changetheta == 2)  t = (taunew - taub[0]) / (tauold - taub[0]);
         stree.nodes[is].theta = thetanew = thetaold / t;
         lnacceptance += -log(t) + (-atheta - 1)*log(thetanew / thetaold) - btheta*(1 / thetanew - 1 / thetaold);
      }
      if (debug == 5) printf("species %d taub: %8.5f %8.5f tau: %8.5f%8.5f", is, taub[0], taub[1], tauold, taunew);

      ntj[0] = ntj[1] = 0;
      for (locus = 0; locus < data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 1);  /* copy gtree & use alternative space for conP */
         ntj_locus[0] = ntj_locus[1] = 0;
         if (mcmc.saveconP)
            for (k = com.ns; k < 2 * com.ns - 1; k++)
               com.oldconP[k] = 1;

         for (inode = com.ns; inode < tree.nnode; inode++) {
            t = nodes[inode].age;
            if (t >= taub[0] && t < taub[1]
               && (nodes[inode].ipop == is || stree.nodes[nodes[inode].ipop].father == is)) {
               k = (t >= tauold && is != stree.root); /* k=0: below; 1: above */
               nodes[inode].age = t = taub[k] + taufactor[k] * (t - taub[k]);
               ntj_locus[k]++;
               if (t < taub[0] || (t<stree.nodes[stree.root].age && t>taub[1]))
                  error2("tj out of tb");

               if (mcmc.saveconP)
                  for (k = inode; ; k = nodes[k].father) {
                     com.oldconP[k] = 0;
                     if (k == tree.root) break;
                  }
            }
         }

         ntj[0] += ntj_locus[0];
         ntj[1] += ntj_locus[1];
         CountCoalescentLocus(locus);
         if (!stree.NoTheta) {
            lnpGinew[locus] = lnpG_S(locus);
            lnpGd += lnpGinew[locus] - data.lnpGi[locus];
         }
         if (ntj_locus[0] + ntj_locus[1]) {
            lnpDinew[locus] = lnpD_locus(locus);
            lnLd += lnpDinew[locus] - data.lnpDi[locus];
         }
         else
            lnpDinew[locus] = data.lnpDi[locus];

         if (debug == 5) printf(" (%d %d) ", ntj_locus[0], ntj_locus[1]);

      }  /* for(locus) */
      if (stree.NoTheta)
         lnpGd = lnpG_S_NoTheta() - *lnpG;
      lnacceptance += lnLd + lnpGd + ntj[0] * log(taufactor[0]) + ntj[1] * log(taufactor[1]);


      if (debug == 5) {
         printf(" ntj: %d %d", ntj[0], ntj[1]);
         printf(" lnLd %9.5f ", lnLd);
         for (k = com.ns; k < com.ns * 2 - 1; k++)  printf("%d", com.oldconP[k]);
      }

      if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
         accepted++;
         for (locus = 0; locus < data.ngene; locus++) {
            data.lnpDi[locus] = lnpDinew[locus];
            if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
            UseLocus(locus, 0, 0, 0);
            for (inode = com.ns; inode < tree.nnode; inode++) { /* redo changes */
               t = nodes[inode].age;
               if (t >= taub[0] && t < taub[1]
                  && (nodes[inode].ipop == is || stree.nodes[nodes[inode].ipop].father == is)) {
                  k = (t >= tauold && is != stree.root);
                  nodes[inode].age = taub[k] + taufactor[k] * (t - taub[k]);
               }
            }
         }
         *lnL += lnLd;
         *lnpG += lnpGd;
         if (mcmc.usedata) SwitchconPin();
         if (debug == 5) printf(" Y (%4d)\n", NPMat);
      }
      else {
         stree.nodes[is].age = tauold;
         if (changetheta) stree.nodes[is].theta = thetaold;
         copyCoalescent(1);
         if (debug == 5) printf(" N (%4d)\n", NPMat);
      }
   }

   return((double)accepted / ntau);
}

int PriorS_IntegerPartitions(int n, int prinT, double PriorS[])
{
   /* Jerome Kelleher's algorithm for generating ascending partitions.
      see table 1 column "Number of Delimitations" in Yang & Rannala (2014 MBE 12: 3125-3135).
      Number of delimitations should be equal to n!/\prod_k #permutations
   */
   int a[256] = { 0 }, npartitions = 0, x, y, k, i, nsame, start;
   double nfactorial = factorial(n), nD, nTree, nGuide, wLH;

   if (n > 255) puts("n too large here");
   if (prinT && n < 64) printf("\nCounting delimitations\n");

   zero(PriorS, n);
   a[k = 1] = n;
   while (k) {
      y = a[k] - 1;
      k--;
      x = a[k] + 1;
      while (x <= y) {
         a[k] = x;
         y -= x;
         k++;
      }
      a[k] = x + y;
      a[k + 1] = 0;
      npartitions++;

      if (prinT && n < 64) {
         printf("[%3d ]  ", k + 1);
         for (i = 0; i < k + 1; i++) printf(" %d", a[i]);
      }

      /* a has the partition, with k+1 clusters (species). */
      for (i = 0, nD = nfactorial; i < k + 1; i++)
         if (a[i] > 1) nD /= factorial(a[i]);

      start = 0; nsame = 1;
      for (i = 1; i < k + 1; i++) {
         if (a[i] != a[start]) {
            if (nsame > 1) nD /= factorial(nsame);
            start = i;  nsame = 1;
         }
         else nsame++;
      }
      if (nsame > 1) nD /= factorial(nsame);
      nTree = CountTrees(k + 1, 1);

      for (i = 0, nGuide = 1; i < k + 1; i++)
         if (a[i] > 2) nGuide *= CountTrees(a[i], 1);

      wLH = 1;
      if (k + 1 > 3 && (stree.SpeciesModelPrior & 1) == 0)   /* uniform S and uniform labeled histories  */
         wLH = CountLHs(k + 1) / CountTrees(k + 1, 1);

      if (prinT)
         printf("%*s nD = %3.0f nTree = %4.0f nGuide = %4.0f  pro=%4.0f\n", n * 5 - k * 2, "", nD, nTree, nGuide, nD*nTree*nGuide*wLH);

      PriorS[k] += nD * nTree * nGuide * wLH;
   }

   if (prinT) {
      printf("\n# delimitations for %d populations is %d", n, npartitions);
      matout2(F0, stree.PriorSA11, 1, n, 9, 0);
      printf(" sum = %.0f", sum(stree.PriorSA11, n));
   }
   abyx(1 / sum(stree.PriorSA11, n), stree.PriorSA11, n);
   if (prinT)
      matout2(F0, stree.PriorSA11, 1, n, 9, 4);

   return(npartitions);
}

double lnpriorSpeciesModel(int ndspecies)
{
   double p = 1;

   if ((stree.SpeciesModelPrior & 1) == 1)                    /* prior 1 or 3: uniform rooted trees */
      p = 1 / CountLHsTree2();

   if (stree.analysis == A11 && stree.SpeciesModelPrior >= 2)  /* prior 2 or 3: uniform #species for A11 */
      p /= stree.PriorSA11[ndspecies - 1];

   return p < 1e-300 ? -500 : log(p);
}


int SetSonNodeFlags(int ispecies, double tauU, char flag[])
{
   /* This sets the two-bits flags for each node in the gene tree.
      The first bit is set if the node has the age < tauU and has a descendent tip
      in population son0 or its descendent.
      The second bit is set if the node has the age < tauU and has a descendent tip
      in population son1 or its descendent.
   */
   int i, j, *sons = stree.nodes[ispecies].sons;
   char debug = 0;

   for (i = 0; i < com.ns; i++) {
      if (stree.pptable[nodes[i].ipop][sons[0]]) {
         flag[i] = 1;
         for (j = nodes[i].father; flag[j] == 0; j = nodes[j].father) {
            if (nodes[j].age > tauU)  break;
            flag[j] = 1;        /* 1st bit for son0 */
            if (j == tree.root) break;
         }
      }
   }
   for (i = 0; i < com.ns; i++) {
      if (stree.pptable[nodes[i].ipop][sons[1]]) {
         flag[i] = 2;
         for (j = nodes[i].father; flag[j] < 2; j = nodes[j].father) {
            if (nodes[j].age > tauU)  break;
            flag[j] += 2;       /* 2nd bit for son1 */
            if (j == tree.root) break;
         }
      }
   }
   /*
   for(i=0; i<tree.nnode; i++)
      printf("node %3d  flag  %d%d\n", i, flag[i]&1, flag[i]>>1);
   */
   return(0);
}


int RubberProportional(int ispecies, double tauU, double tau, double taunew, double *lnproposal, double *space)
{
   /* This modifies the gene tree (nodes) for a locus for the split and join moves.
      The two-bits flags for each node are set to indicate whether the node is ancestral to sons[0] and sons[1].
      When ispecies is the species tree root, tauU is somewhat arbitrary.
      This adds to lnproposal.
      Split: tau=0; taunew>0;
      Join:  tau>0; taunew=0;
      returns 1 if naffected>0, in which case lnL for the locus needs to be recalculated.
   */
   int i, j, *sons = stree.nodes[ispecies].sons, dad = stree.nodes[ispecies].father;
   int  naffected = 0, nwithin;
   double Jrubber = (tauU - taunew) / (tauU - tau), tnew, y, *c = space;  /* c[i-com.ns] is c for proporitional expansion */
   char *flag = (char*)(space + com.ns - 1);
   char debug = 0;

   if (debug) {
      printf("\ngene tree before RubberProportional\n");
      printGtree(1);
   }
   memset(space, 0, (com.ns - 1) * sizeof(double) + (com.ns * 2 - 1) * sizeof(char));
   SetSonNodeFlags(ispecies, tauU, flag);

   for (i = com.ns; i < com.ns * 2 - 1; i++) {
      if (flag[i] != 3) continue;  /* affected nodes */
      naffected++;
      tnew = tauU - Jrubber*(tauU - nodes[i].age);
      c[i - com.ns] = tnew / nodes[i].age;
      nodes[i].age = tnew;
   }
   if (debug) {
      for (i = com.ns; i < com.ns * 2 - 1; i++)
         printf("node %2d (ipop %2d): flag %d %d t* =%9.6f c =%9.6f\n", i, nodes[i].ipop, flag[i] & 1, flag[i] >> 1, nodes[i].age, c[i - com.ns]);
   }

   if (naffected) {   /* proportional changes to ages of within-move nodes */
      for (i = com.ns, nwithin = 0, y = 1; i < com.ns * 2 - 1; i++) {
         if (flag[i] != 1 && flag[i] != 2) continue;
         if (stree.pptable[nodes[i].ipop][ispecies] == 0) continue;
         for (j = nodes[i].father; ; j = nodes[j].father) {
            if (flag[j] == 3) {
               nwithin++;
               nodes[i].age *= c[j - com.ns];
               y *= c[j - com.ns];
               if (debug)
                  printf("node changed prop. %2d: flag %d %d t* = %8.6f (c = %8.6f from node %d) \n",
                     i, flag[i] & 1, flag[i] >> 1, nodes[i].age, c[j - com.ns], j);
               break;
            }
            if (j == tree.root) break;
         }
      }
      if (nwithin) *lnproposal += log(y);
      *lnproposal += naffected*log(Jrubber);
   }

   if (taunew > 0)  /* for split, change some ipop from i to j or k. */
      for (i = com.ns; i < com.ns * 2 - 1; i++) {
         if (nodes[i].ipop == ispecies && nodes[i].age < taunew)
            nodes[i].ipop = sons[flag[i] - 1];  /* flag[i] must be 1 or 2. */
      }
   else          /* for join, change ipop j & k to i. */
      for (i = com.ns; i < com.ns * 2 - 1; i++) {
         if (nodes[i].ipop == sons[0] || nodes[i].ipop == sons[1])
            nodes[i].ipop = ispecies;
      }

   if (debug) {
      printf("\ngene tree after RubberProportional (nwithin & between = %d %d)\n", nwithin, naffected);
      printGtree(1);
   }

   /* if(naffected==0), lnL for the locus does not change.  */
   return(naffected > 0);
}


/* LnGamma(a) is constant throughout the chain */
#define lnPDFGamma(x, a, b)  ( (a)*log(b) - LnGamma(a) + ((a)-1)*log(x) - (b)*(x) )
#define lnPDFinvGamma(x, a, b)  ( (a)*log(b) - LnGamma(a) - (a+1)*log(x) - (b)/(x) )
#define lnPDFBeta(x, p, q, b)  (-LnBeta(p,q) + (p-1)*log(x/b) + (q-1)*log(1-x/b) - log(b) )

int UpdateSpeciesSplit(double *lnL, double* lnpG, double space[], double PrSplit)
{
   /* This splits a node in the guide species tree.  A node is feasible for splitting
      if it is joined and if its mother node is not joined.  The move changes node ages
      in gene trees, but keep the gene tree topology and ipop unchanged.
   */
   int nfeasible[2] = { 0 }, feasibles[NSPECIES] = { 0 }, locus, is, i, j, *sons, accepted = 0;
   int i1, i2, s = stree.nspecies, ntau, ithetajk[2] = { 0 }, iModelNew = -1;
   double atheta = data.theta_prior[0], btheta = data.theta_prior[1];  /* G(a,b) */
   double atau = data.tau_prior[0], btau = data.tau_prior[1];
   double lnacceptance = log((1 - PrSplit) / PrSplit), taunew, tauU;
   double lnLd, lnpGd, *lnpGinew = space, *lnpDinew = space + data.ngene, lnpSpeciesModelnew = 0;
   double thetafactor = 1, y, thetai, aRJ = mcmc.RJfinetune[0], mRJ = mcmc.RJfinetune[1];
   double pbetatau = 2, qbetatau = 8;

   if (debug == 6) {
      printf("\nUpdateSpeciesSplit\nSpecies tree: %s\n", printDelimitationModel());
      for (i = 0; i < stree.nnode; i++)
         printf("node %2d theta = %9.6f tau = %9.6f\n", i, stree.nodes[i].theta, stree.nodes[i].age);
   }

   /* count nfeasible for splitting at source species tree */
   for (i = s, ntau = 0; i < stree.nnode; i++) {
      if (stree.nodes[i].age > 0)
         ntau++;
      else {
         if (i == stree.root || stree.nodes[stree.nodes[i].father].age > 0) /* i is joined node */
            feasibles[nfeasible[0] ++] = i;
      }
   }
   if (nfeasible[0] == 0)
      return(2);

   copyCoalescent(0);

   if (mcmc.saveconP) /* change this to try to save */
      for (i = data.maxns; i < data.maxns * 2 - 1; i++)
         com.oldconP[i] = 0;


   /* split species is, generate theta_j & theta_k, if they can exist, and
      calculate their prior.
    */
   is = feasibles[(int)(nfeasible[0] * rndu())];
   thetai = stree.nodes[is].theta;  /* not used if(stree.NoTheta) */
   sons = stree.nodes[is].sons;
   if (is == stree.root)
      tauU = stree.roottau*0.6;
   else
      tauU = stree.nodes[stree.nodes[is].father].age;

   /*
  stree.nodes[is].age = taunew = tauU*rndu();
  */
   stree.nodes[is].age = taunew = tauU*rndbeta(pbetatau, qbetatau);
   lnacceptance -= lnPDFBeta(taunew, pbetatau, qbetatau, tauU);

   if (debug == 6)
      printf("Spnode %d (sons: %d %d)\ntU = %.6f taunew = %.6f\n", is, sons[0], sons[1], tauU, taunew);

   /* new theta_j & theta_k */
   if (!stree.NoTheta) {
      for (j = 0; j < 2; j++) {  /* theta_j & theta_k */
         for (i = 0; i < stree.npop; i++)
            if (sons[j] == stree.pops[i]) break;
         if (i == stree.npop) continue;
         ithetajk[j] = 1;
         if (mcmc.RJalgorithm == 0) {
            stree.nodes[sons[j]].theta = thetai * exp(aRJ*(rndu() - 0.5));
            thetafactor *= aRJ * stree.nodes[sons[j]].theta;
         }
         else {
            stree.nodes[sons[j]].theta = rndgamma(aRJ) / (aRJ / (mRJ*thetai));
            thetafactor /= PDFGamma(stree.nodes[sons[j]].theta, aRJ, aRJ / (mRJ*thetai));
         }
         /* prior on theta */
         lnacceptance += lnPDFinvGamma(stree.nodes[sons[j]].theta, atheta, btheta);
      }
   }

   /* prior on species demilitation model and on tau */
   lnpSpeciesModelnew = lnpriorSpeciesModel(ntau + 2);  /* new model has ntau+2 species */
   lnacceptance += lnpSpeciesModelnew - data.lnpSpeciesModel;
   if (is == stree.root)   /* ntau changed from 0 to 1.  */
      lnacceptance += lnPDFinvGamma(taunew, atau, btau);
   else
      lnacceptance += log(ntau / stree.nodes[stree.root].age);  /* Eq 2 in YR2010 */

   /* count nfeasible for joining at target.  This cannot be calculated before we know
      which node to split (is).  A node is feasible for joining if each son is either
      a tip or is already joined.
   */
   for (i = s; i < stree.nnode; i++) {
      if (stree.nodes[i].age == 0) continue;
      i1 = stree.nodes[i].sons[0];
      i2 = stree.nodes[i].sons[1];
      if ((i1 < s || stree.nodes[i1].age == 0)
         && (i2 < s || stree.nodes[i2].age == 0))
         nfeasible[1] ++;
   }

   lnacceptance += log((double)nfeasible[0] / nfeasible[1] * thetafactor);

   /* rubber band to modify node ages in gene trees */
   for (locus = 0, lnLd = lnpGd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & use alternative space for conP */
      if (debug == 6) {
         printf("\ngene tree at locus %d before rubber band\n", locus + 1);
         printGtree(1);
      }
      j = RubberProportional(is, tauU, 0, taunew, &lnacceptance, space + data.ngene * 2);
      if (debug == 6) {
         printf("\ngene tree at locus %d after rubber band\n", locus + 1);
         printGtree(1);
      }

      lnpDinew[locus] = (j == 0 ? data.lnpDi[locus] : lnpD_locus(locus));
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnpGd + lnLd;

   if (stree.speciesdelimitation == 1 && stree.speciestree == 0)  iModelNew = GetDmodelIndex();

   if (debug == 6)
      printf("lnacceptance = %.6f\n", lnacceptance);

   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      accepted = 1;
      if (stree.speciesdelimitation == 1 && stree.speciestree == 0)  stree.iModel = iModelNew;
      for (locus = 0; locus < data.ngene; locus++) {
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
         data.lnpDi[locus] = lnpDinew[locus];
         /* re-do the changes to gtrees (age) */
         UseLocus(locus, 0, 0, 0);
         RubberProportional(is, tauU, 0, taunew, &y, space + data.ngene * 2);
         if (debug == 6) {
            printf("\nsplit success locus %d\n", locus);
            printGtree(1);
         }
      }
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      if (stree.NoTheta)  com.np++;
      else                com.np += 1 + (stree.nodes[sons[0]].theta > 0) + (stree.nodes[sons[1]].theta > 0);
      data.lnpSpeciesModel = lnpSpeciesModelnew;

      if (debug == 6) printf("split move accepted");
   }
   else {                                           /* reject */
      stree.nodes[is].age = 0;
      stree.nodes[sons[0]].theta = -1;
      stree.nodes[sons[1]].theta = -1;
      copyCoalescent(1);
      if (debug == 6) printf("split move rejected");
   }
   return(accepted);
}



int UpdateSpeciesJoin(double *lnL, double* lnpG, double space[], double PrSplit)
{
   /* This joins a node in the guide species tree.  Note that the join move may be disallowed
      if the current tau and theta's are unreachable from the split move.  In this case, the
      function returns 2.  Otherwise it returns 0 for reject and 1 for acceptance.
      See notes in UpdateSpeciesSplit().
   */
   int nfeasible[2] = { 0 }, feasibles[NSPECIES] = { 0 }, locus, is, i, j, *sons, accepted = 0;
   int ntau, iModelNew = -1, s = stree.nspecies;
   double atheta = data.theta_prior[0], btheta = data.theta_prior[1]; /* invG(a,b) */
   double atau = data.tau_prior[0], btau = data.tau_prior[1];
   double lnacceptance = log(PrSplit / (1 - PrSplit)), y, tau, tauU;
   double lnpGd, lnLd, *lnpGinew = space, *lnpDinew = space + data.ngene, lnpSpeciesModelnew = 0;
   double thetafactor = 1, thetai, thetajk0[2], aRJ = mcmc.RJfinetune[0], mRJ = mcmc.RJfinetune[1];
   double pbetatau = 2, qbetatau = 8;

   if (debug == 7) printf("\nUpdateSpeciesJoin\nSpecies tree: %s\n", printDelimitationModel());

   /* count nfeasible for joining at source species tree.  A node is feasible
      for joining if each of its sons is either tip or already joined.
   */
   for (i = s, ntau = 0; i < stree.nnode; i++) {
      if (stree.nodes[i].age == 0) continue;
      ntau++;
      sons = stree.nodes[i].sons;
      if ((sons[0] < s || stree.nodes[sons[0]].age == 0)
         && (sons[1] < s || stree.nodes[sons[1]].age == 0))
         feasibles[nfeasible[0] ++] = i;
   }
   if (nfeasible[0] == 0) return(2);
   copyCoalescent(0);

   /* join species is */
   is = feasibles[(int)(nfeasible[0] * rndu())];
   tau = stree.nodes[is].age;
   thetai = stree.nodes[is].theta;  /* not used if stree.NoTheta */
   sons = stree.nodes[is].sons;

   if (debug == 7) printf("Joining species node %d (sons: %d %d)\n", is, sons[0], sons[1]);

   /* get tauU, for coming back through a split move */
   if (is == stree.root)
      tauU = stree.roottau*0.6;
   else
      tauU = stree.nodes[stree.nodes[is].father].age;

   if (debug == 7)
      printf("spnode %d: tauU = %.6f tau = %.6f\n", is, tauU, tau);
   if (tau >= tauU) return(2);  /* join move disallowed */

   lnacceptance += lnPDFBeta(tau, pbetatau, qbetatau, tauU);

   /* theta_j & theta_k */
   if (!stree.NoTheta) {
      for (i = 0; i < 2; i++)  thetajk0[i] = stree.nodes[sons[i]].theta;
      for (i = 0; i < 2; i++) {
         if (thetajk0[i] <= 0) continue;
         if (mcmc.RJalgorithm == 0) {
            y = exp(aRJ*0.5);
            if (thetajk0[i]<thetai / y || thetajk0[i]>thetai*y) /* move disallowed */
               return(2);
            thetafactor /= aRJ * thetajk0[i];
         }
         else {
            thetafactor *= PDFGamma(thetajk0[i], aRJ, aRJ / (mRJ*thetai));
         }
         /* prior on old theta */
         lnacceptance -= lnPDFinvGamma(thetajk0[i], atheta, btheta);
      }
   }

   stree.nodes[is].age = 0;
   stree.nodes[sons[0]].theta = stree.nodes[sons[1]].theta = -1;
   /* prior on species delimitation model and on tau */
   lnpSpeciesModelnew = lnpriorSpeciesModel(ntau);   /* new model has ntau species */
   lnacceptance += lnpSpeciesModelnew - data.lnpSpeciesModel;
   if (is == stree.root)  /* ntau changed from 1 to 0.  */
      lnacceptance -= lnPDFinvGamma(tau, atau, btau);
   else
      lnacceptance -= log((ntau - 1) / stree.nodes[stree.root].age);  /* Eq 2 in YR2010 */

   /* count nfeasible for splitting at target.  A node is feasible if its father
      is already split.  This is calculated after the species to join (is) is known.
   */
   for (i = s; i < stree.nnode; i++) {
      if (stree.nodes[i].age == 0 &&
         (i == stree.root || stree.nodes[stree.nodes[i].father].age > 0))
         nfeasible[1] ++;
   }
   lnacceptance += log((double)nfeasible[0] / nfeasible[1] * thetafactor);

   /* make a copy of ipop in gene trees.  Change node ipop from sons[0] or sons[1] to is.
      apply the rubber band algorithm
   */

   if (mcmc.saveconP) /* change this to try to save */
      for (i = data.maxns; i < data.maxns * 2 - 1; i++)
         com.oldconP[i] = 0;

   for (locus = 0, lnLd = lnpGd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0); /* copy gtree & use alternative space for conP */
      if (debug == 7) {
         printf("\nbefore rubber band gene tree at locus %d\n", locus + 1);
         printGtree(1);
      }
      j = RubberProportional(is, tauU, tau, 0, &lnacceptance, space + data.ngene * 2);
      if (debug == 7) {
         printf("\nafter rubber band gene tree at locus %d\n", locus + 1);
         printGtree(1);
      }

      lnpDinew[locus] = (j == 0 ? data.lnpDi[locus] : lnpD_locus(locus));
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnpGd + lnLd;

   if (stree.speciesdelimitation == 1 && stree.speciestree == 0)
      iModelNew = GetDmodelIndex();

   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      accepted = 1;
      if (stree.speciesdelimitation == 1 && stree.speciestree == 0)
         stree.iModel = iModelNew;
      for (locus = 0; locus < data.ngene; locus++) {
         data.lnpDi[locus] = lnpDinew[locus];
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
         /* re-do the changes to gtrees: age and ipop */
         UseLocus(locus, 0, 0, 0);
         RubberProportional(is, tauU, tau, 0, &y, space + data.ngene * 2);

         if (debug == 7) {
            printf("\njoin success locus %d\n", locus);
            UseLocus(locus, 0, 0, 0);
            printGtree(1);
         }
      }
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      if (stree.NoTheta)  com.np--;
      else                com.np -= 1 + (thetajk0[0] > 0) + (thetajk0[1] > 0);
      data.lnpSpeciesModel = lnpSpeciesModelnew;
   }
   else {                                            /* reject */
      if (debug == 7) printf("\njoin failure\n");
      stree.nodes[is].age = tau;
      copyCoalescent(1);
      if (!stree.NoTheta) for (i = 0; i < 2; i++) stree.nodes[sons[i]].theta = thetajk0[i];
   }

   if (debug == 7) {
      printf("join end    ages: %9.5f%9.5f\n", gnodes[0][3].age, gnodes[0][4].age);
      printf("stree ages: %9.5f%9.5f\n", stree.nodes[3].age, stree.nodes[4].age);
   }
   return(accepted);
}

double mixing(double* lnL, double* lnpG, double finetune, double space[])
{
   /* This multiplies all tau (& theta), and branch lengths in all gene trees by c.
      This move can bring a node age in a gene tree to be older than OldAge,
      which may not be nice.
   */
   int accepted = 0, locus, is, i, k, ntau = stree.nspecies - 1, ntheta;
   double xold, xnew, *lnpDinew = space, *lnpGinew = lnpDinew + data.ngene;
   double c, lnc, lnacceptance = 0, lnLd, lnpGd;
   double atheta = data.theta_prior[0], btheta = data.theta_prior[1];
   double atau = data.tau_prior[0], btau = data.tau_prior[1];

   if (finetune <= 0) error2("steplength = 0 in mixing.  Check finetune");
   copyCoalescent(0);
   ntheta = (stree.NoTheta ? 0 : stree.npop);
   if (stree.speciesdelimitation == 1 || stree.speciestree == 1) {
      for (is = stree.nspecies, ntau = 0; is < 2 * stree.nspecies - 1; is++)
         ntau += (stree.nodes[is].age > 0);
      if (!stree.NoTheta) {
         for (is = 0, ntheta = 0; is < stree.npop; is++) {
            if (stree.nodes[stree.pops[is]].theta > 0)
               ntheta++;
         }
      }
   }
   lnc = finetune*rndSymmetrical();
   c = exp(lnc);
   for (locus = 0, k = 0; locus < data.ngene; locus++)
      k += data.ns[locus] - 1;
   lnacceptance = (ntheta + ntau + k)*lnc;

   if (!stree.NoTheta) {
      for (i = 0; i < stree.npop; i++) {
         is = stree.pops[i];
         if (stree.nodes[is].theta <= 0) continue;
         xold = stree.nodes[is].theta;
         stree.nodes[is].theta = xnew = xold*c;
         lnacceptance += (-atheta - 1)*lnc - btheta*(1 / xnew - 1 / xold);
      }
   }
   for (i = 0; i < stree.nspecies - 1; i++) {
      is = stree.nspecies + i;
      if (stree.nodes[is].age == 0) continue;
      xold = stree.nodes[is].age;
      stree.nodes[is].age = xnew = xold*c;
      if (is == stree.root)  /* root in species tree */
         lnacceptance += (-atau - 1 - ntau + 1)*lnc - btau*(1 / xnew - 1 / xold);
   }

   if (mcmc.saveconP)
      for (i = 0; i < data.maxns * 2 - 1; i++) com.oldconP[i] = 0;
   for (locus = 0, lnLd = lnpGd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & conP */
      for (i = com.ns; i < tree.nnode; i++)  nodes[i].age *= c;
      lnpDinew[locus] = lnpD_locus(locus);
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnLd + lnpGd;
   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      for (locus = 0; locus < data.ngene; locus++) {
         UseLocus(locus, 0, 1, 0);
         for (i = com.ns; i < tree.nnode; i++)
            nodes[i].age *= c;
         data.lnpDi[locus] = lnpDinew[locus];
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
      }
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      accepted = 1;
   }
   else {
      if (!stree.NoTheta) {
         for (i = 0; i < stree.npop; i++)
            if (stree.nodes[stree.pops[i]].theta > 0)
               stree.nodes[stree.pops[i]].theta /= c;
      }
      for (i = stree.nspecies; i < stree.nspecies * 2 - 1; i++) {
         if (stree.nodes[i].age > 0)
            stree.nodes[i].age /= c;
      }
      copyCoalescent(1);
   }
   return(accepted);
}


double UpdateLocusrateHeredity(double* lnL, double* lnpG, double finetune)
{
   /* This updates locus-specific rates and heredity multipliers.
      There is probably no need to update both when both are estimated for the
      same locus, as the posterior for rates and heredity multipliers should be
      quite flat.
   */
   int  accepted = 0, locus, j, locusref = 0;
   double lnacceptance, lnLd, lnpGd, lnpDinew, lnpGinew;
   double h, hnew, r, rnew, rref, rrefnew, lnpDrefnew;

   if (debug == 2) puts("\nUpdateLocusrateHeredity ");
   if (finetune <= 0) error2("steplength = 0 in UpdateLocusrateHeredity.  Check finetune");

   if (mcmc.saveconP) FOR(j, data.maxns * 2 - 1) com.oldconP[j] = 0;

   if (data.est_locusrate == 1) {
      for (j = 1, locusref = 0; j < data.ngene; j++)
         if (data.npatt[j] > data.npatt[locusref]) locusref = j;

      for (locus = 0; locus < data.ngene; locus++) {
         if (locus == locusref) continue;
         r = data.locusrate[locus];
         rref = data.locusrate[locusref];
         rnew = r + finetune*rndSymmetrical();
         rnew = data.locusrate[locus] = reflect(rnew, 0, r + rref);
         rrefnew = data.locusrate[locusref] -= rnew - r;

         lnacceptance = (data.a_locusrate - 1)*log((rnew*rrefnew) / (r*rref));

         UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & conP */
         lnpDinew = lnpD_locus(locus);
         UseLocus(locusref, 1, mcmc.usedata, 0);
         lnpDrefnew = lnpD_locus(locusref);

         lnLd = lnpDinew - data.lnpDi[locus] + lnpDrefnew - data.lnpDi[locusref];
         lnacceptance += lnLd;
         if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
            accepted++;
            if (mcmc.usedata) {
               *lnL += lnLd;
               data.lnpDi[locus] = lnpDinew;
               data.lnpDi[locusref] = lnpDrefnew;
               CopyconPin(locus, 0);
               CopyconPin(locusref, 0);
            }
         }
         else {                                            /* reject */
            data.locusrate[locus] = r;
            data.locusrate[locusref] = rref;
         }
      }
   }

   if (data.est_heredity == 1) { /* this step does not change likelihood */
      for (locus = 0; locus < data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree but not conP */
         copyCoalescent(0);

         h = data.heredity[locus];
         hnew = h + finetune*rndSymmetrical();
         if (hnew < 0) hnew *= -1;
         data.heredity[locus] = hnew;
         lnacceptance = (data.a_heredity - 1)*log(hnew / h) - data.b_heredity*(hnew - h);
         CountCoalescentLocus(locus);
         if (!stree.NoTheta) {
            lnpGinew = lnpG_S(locus);
            lnpGd = lnpGinew - data.lnpGi[locus];
         }
         else
            lnpGd = lnpG_S_NoTheta() - *lnpG;

         lnacceptance += lnpGd;
         if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
            accepted++;
            if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew;
            *lnpG += lnpGd;
         }
         else {
            data.heredity[locus] = h;
            copyCoalescent(1);
         }
      }
   }
   return(accepted / (data.est_locusrate*(data.ngene - 1.0) + data.est_heredity*data.ngene));
}

double UpdateSequenceErrors(double* lnL, double finetune, double space[])
{
   /* This updates the parameters for sequencing errors
   */
   int  accepted = 0, locus, is, i, j;
   double lnacceptance, lnLd, *lnpDinew = space, eold, enew, eDold, eDnew, *a;

   if (finetune <= 0) error2("steplength = 0 in UpdateSequenceErrors.  Check finetune");
   for (i = com.ns; i < com.ns * 2 - 1; i++)
      com.oldconP[i] = 0;

   for (is = 0; is < stree.nspecies; is++) {
      if (data.iseqerr[is] == 0) continue;
      for (i = 0; i < 4; i++) {
         a = data.a_seqerr[is] + i * 4;
         for (j = 0; j < 4; j++) {
            if (i == j) continue;
            eold = data.e_seqerr[is][i * 4 + j];
            eDold = data.e_seqerr[is][i * 4 + i];
            enew = eold + finetune*rndSymmetrical();
            enew = data.e_seqerr[is][i * 4 + j] = reflect(enew, 1e-20, eDold + eold);
            eDnew = data.e_seqerr[is][i * 4 + i] = eDold + eold - enew;

            lnacceptance = (a[j] - 1)*log(enew / eold) + (a[i] - 1)*log(eDnew / eDold);

            for (locus = 0, lnLd = 0; locus < data.ngene; locus++) {
               UseLocus(locus, 1, mcmc.usedata, 0);   /* copy gtree & conP */
               lnpDinew[locus] = lnpD_locus(locus);
               lnLd += lnpDinew[locus] - data.lnpDi[locus];
               lnacceptance += lnpDinew[locus] - data.lnpDi[locus];
            }

            if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
               accepted++;
               for (locus = 0; locus < data.ngene; locus++)
                  data.lnpDi[locus] = lnpDinew[locus];
               if (mcmc.usedata) SwitchconPin();
               *lnL += lnLd;
            }
            else {
               data.e_seqerr[is][i * 4 + j] = eold;
               data.e_seqerr[is][i * 4 + i] = eDold;
            }
         }
      }
   }
   return((double)accepted / (data.nseqerr*12.0));
}

int BranchWeights(double weight[])
{
   int s = stree.nspecies, i;
   double t, bpower = -0.5;  /* weights for focus branch & target branch c. */

   for (i = s; i < stree.nnode; i++) {
      weight[i - s] = 0;
      if (i != stree.root && stree.nodes[i].age > 0) {
         weight[i - s] = stree.nodes[stree.nodes[i].father].age - stree.nodes[i].age;
         weight[i - s] = pow(weight[i - s], bpower);
      }
   }
   for (i = 0, t = sum(weight, s - 1); i < s - 1; i++) weight[i] /= t;
   if (debug) {
      printf("\nnodes   ");  for (i = 0; i < s - 1; i++) printf("%9d", s + i);
      printf("\nages    ");  for (i = 0; i < s - 1; i++) printf("%9.5f", stree.nodes[s + i].age);
      printf("\nweights ");  for (i = 0; i < s - 1; i++) printf("%9.5f", weight[i]);
   }
   return(0);
}


int UpdateSpeciesTreeNNI(double *lnL, double* lnpG, double space[])
{
   /* This uses NNI to update the species tree topology and SPR to update gene trees
      at the same time (see figure # in Yang & Rannala 2014).
      (a) SPR change of the gene trees.  For each locus, I scan the nodes on the gene tree
      to identity Movednodes, and nodes whose ipop need changing (AB2B), collect the
      necessary info (such as ntarget, nsource, and targetChosen, etc) for each Movednode
      without channging the gene tree.  After all affected nodes for the locus are
      identified, the changes are applied.
      (b) NNI change to the species tree.
      (c) calculation of the prior and likelihood to accept/reject the move.
      Other notes:
      com.oldconP[] is set to 0, so that repeated likelihood calculation may be performed.
   */
   int s = stree.nspecies, ysona, xsony, locus, inode, dad, sons[2], i, j, k, accepted = 0;
   int x, y, a, b, c, nib, nAsons, sizeGtrees;
   int nMoved[NGENE], PrunedSons[NS - 1], targets[NS * 2 - 1] = { 0 }, TargetChosen[NS - 1], nsource, ntarget;
   int sources[NS * 2 - 1] = { 0 };
   char AB2B[NS - 1];
   double tMoved[NS - 1], t, tau0, tau1;  /* tau0>tau1 are ages of x & y. */
   double lnacceptance = 0, lnpGd, lnLd, *lnpGinew = space, *lnpDinew = space + data.ngene;
   double lnpSpeciesModelnew = 0;
   char *flag = (char*)space;   /* flag[ns*2-1] for A & B, using the same space for loci. */
   struct SPECIESTREE stree0;


   error2("NNI needs to be tested, especially after invG priors.  1/11/2016\n");

   for (i = s, nib = 0; i < stree.nnode; i++)
      if (i != stree.root && stree.nodes[i].age > 0)  nib++;

   if (nib < 1) error2("can't apply NNI!  Number of internal branches is < 1");

   /* choose the internal branch in the species tree to apply NNI */
   while (1) {
      y = s + 1 + (int)(rndu()*(s - 2));  /* root is s; non-root is s+1, ...,  */
      if (stree.nodes[y].age) break;
   }

   /* NNI to swap b with c, or to move A to C. */
   ysona = (int)(2 * rndu());    /* this is A.   */
   a = stree.nodes[y].sons[ysona];
   b = stree.nodes[y].sons[1 - ysona];
   x = stree.nodes[y].father;
   xsony = (y == stree.nodes[x].sons[0] ? 0 : 1);
   c = stree.nodes[x].sons[1 - xsony];
   tau0 = stree.nodes[x].age;  tau1 = stree.nodes[y].age;

   if (stree.nodes[c].age > tau1) {
      return(2);
   }

   copyCoalescent(0);

   stree0 = stree;
   for (locus = 0, sizeGtrees = 0; locus < data.ngene; locus++)
      sizeGtrees += (data.ns[locus] * 2 - 1) * sizeof(struct TREEN);
   memmove(gnodes_t[0], gnodes[0], sizeGtrees);

   if (mcmc.usedata && mcmc.saveconP)
      for (j = 0; j < data.maxns * 2 - 1; j++)  com.oldconP[j] = 0;

   if (debug == 11) {
      printf("\n\n******** UpdateSpeciesTreeNNI (moving %d %s to %d %s) ", a, stree.nodes[a].name, c, stree.nodes[c].name);
      i = stree.nodes[4].sons[0]; j = stree.nodes[4].sons[1];
      if (s == 3)
         printf("\nSpecies tree for 3 species: (%s, %s) ", stree.nodes[i].name, stree.nodes[j].name);
   }

   /* SPR move to modify gene tree topologies */
   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, 0, 1);  /* copy gtree & use alternative space for conP */
      nMoved[locus] = 0;
      for (i = 0; i < com.ns - 1; i++)    AB2B[i] = (char)0;
      for (i = 0; i < com.ns * 2 - 1; i++)  sources[i] = -1;
      for (i = 0; i < com.ns * 2 - 1; i++)  flag[i] = 0;
      SetSonNodeFlags(a, x, flag);   /* set flags for populations A & B */

      if (debug == 11) {
         printf("\n\n******* Gene tree before NNI, species node ages (tau)\n");
         printf("tau_%d = %9.6f\ntau_%d = %9.6f\n", x, tau0, y, tau1);
         printf("\nlocus %d (ns = %d)\n", locus + 1, data.ns[locus]);
         if (printGtree(1)) error2("gene tree error before NNI");
         printf("flags\n");
         for (i = com.ns; i < tree.nnode; i++)
            printf("node %3d  flag  %d%d\n", i, flag[i] & 1, flag[i] >> 1);
      }
      /* Identify Moved nodes.  This loop does not change the gene tree. */
      for (i = com.ns; i < tree.nnode; i++) {
         /* Moved node i satisfies the following conditions:
           It is in population AB (y).  It has descendents in both A & B.
           One of its son nodes is pure A and has descendents in A only.
           Node in pop C with age between (tau1, tau0) changes its ipop from C to AC.
           Node with both sons having B flags changes its ipop from AB to B.
         */
         if (nodes[i].ipop != y) continue;
         sons[0] = nodes[i].sons[0];
         sons[1] = nodes[i].sons[1];
         /* AB -> B, square nodes have both sons with B descendents.  When AB disappears,
            nodes in AB change ipop to B.  Nodes in AB ancestors are unaffected.
         */
         if (flag[sons[0]] >> 1 && flag[sons[1]] >> 1) {
            AB2B[i - com.ns] = (char)1;  continue;
         }
         /* Identify Moved nodes, which have 1 A son. */
         nAsons = 0;
         if (flag[sons[0]] == 1) { nAsons++;  j = sons[0]; }
         if (flag[sons[1]] == 1) { nAsons++;  j = sons[1]; }
         if (flag[i] != 3 || nAsons != 1) continue;
         PrunedSons[nMoved[locus]++] = j;
         tMoved[nMoved[locus] - 1] = t = nodes[i].age;
         if (debug == 11)
            printf("\nMovednode %2d (sons %2d %2d), time = %9.5f", i, j, nodes[i].sons[0] + nodes[i].sons[1] - j, t);

         /* count feasible branches in source and target. */
         ntarget = 0;
         /* feasible lineages in target: The branch covers t.  Its ipop is descendent of C. */
         for (j = 0; j < tree.nnode; j++) {
            if (stree.pptable[nodes[j].ipop][c] && nodes[j].age <= t && nodes[nodes[j].father].age > t)
               targets[ntarget++] = j;
         }

         /** If no C seqs exist at the locus, there won't be target branch on gene tree */
         if (ntarget == 0) {
            if (locus)  memmove(gnodes[0], gnodes_t[0], sizeGtrees);  /* gene trees at all loci */
            return(2);
         }

         TargetChosen[nMoved[locus] - 1] = targets[(int)(ntarget*rndu())];
         if (debug == 11) {
            printf("\nnode %d: %d feasible lineages in target (chosen %d): ", i, ntarget, TargetChosen[nMoved[locus] - 1]);
            for (j = 0; j < ntarget; j++) printf(" %2d", targets[j]);  FPN(F0);
         }

         /* feasible lineages in source:
         The branch covers t.  Its ipop is either B or descendent of B. */
         nsource = 1; sources[0] = sons[0] + sons[1] - PrunedSons[nMoved[locus] - 1];
         for (j = 0; j < tree.nnode; j++) {
            if (j == sources[0] || j == nodes[sources[0]].father)     continue;
            if (nodes[j].age >= t || nodes[nodes[j].father].age <= t) continue;
            if (stree.pptable[nodes[j].ipop][b] || (nodes[j].ipop == y && flag[j] != 1))
               sources[nsource++] = j;
         }
         if (debug == 11) {
            printf("node %d: %d feasible lineages in source: ", i, nsource);
            for (j = 0; j < nsource; j++) printf(" %2d", sources[j]);  FPN(F0);
         }

         lnacceptance += log((double)ntarget / nsource);
      }  /* for(i), gene tree nodes */

      /* Now all moved nodes for the locus are identified, I apply the
         SPR changes to the gene tree.
         The Moved node has two sons: sons[0]=PrunedSons[i] is all A and is pruned off and moved to target.
         Then apply the NNI change on the species tree.
      */
      for (i = 0; i < nMoved[locus]; i++) {
         /* Pruning: prune off subtree at sons[0]. */
         sons[0] = PrunedSons[i];  t = tMoved[i];
         inode = nodes[sons[0]].father;
         sons[1] = nodes[inode].sons[0] + nodes[inode].sons[1] - sons[0];
         dad = nodes[inode].father;
         for (j = 0; j < 2; j++)
            if (nodes[dad].sons[j] == inode) break;
         nodes[dad].sons[j] = sons[1];  /* pruning off the subtree at sons[0]. */
         nodes[sons[1]].father = dad;   /* sons[1] changes father from inode to dad */

         /* Regrafting: subtree at sons[0] is inserted on branch k */
         for (k = TargetChosen[i]; ; k = dad) {
            dad = nodes[k].father;
            if (nodes[dad].age > t) break;
         }
         for (j = 0; j < 2; j++)
            if (nodes[dad].sons[j] == k) break;
         nodes[dad].sons[j] = inode;  /* check nodes[inode].age */
         if (nodes[inode].sons[0] == sons[0]) nodes[inode].sons[1] = k;
         else                              nodes[inode].sons[0] = k;
         nodes[inode].father = dad;
         nodes[k].father = inode;
         if (debug == 11) {
            printf("\n*** After Movednode %2d (son %2d) is moved to target %d", inode, sons[0], k);
            printGtree(1);
         }
      }

      /* reset ipop for diamond and square nodes: C -> AC, and AB -> B. */
      for (i = com.ns; i < tree.nnode; i++) {
         if (AB2B[i - com.ns])
            nodes[i].ipop = b;
         else if (nodes[i].ipop == c && nodes[i].age > tau1 && nodes[i].age < tau0)
            nodes[i].ipop = y;
      }
      memcpy(gnodes[locus], nodes_t, (com.ns * 2 - 1) * sizeof(struct TREEN));
   }  /* for(locus) */

   /* Change species tree, update pptable[], concatenate species names, etc. */
   stree.nodes[y].sons[1 - ysona] = c;  stree.nodes[c].father = y;
   stree.nodes[x].sons[1 - xsony] = b;  stree.nodes[b].father = x;
   for (i = s; i < 2 * s - 1; i++) stree.nodes[i].name[0] = '\0';
   DownSptreeSetSpnames(stree.root, 1);
   SetupPopPopTable(debug == 11);

   lnpSpeciesModelnew = lnpriorSpeciesModel(nib + 2);  /* new model has ntau+1 species  */
   lnacceptance += lnpSpeciesModelnew - data.lnpSpeciesModel;

   for (locus = 0, lnpGd = lnLd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 1);  /* copy gtree & use alternative space for conP */
      if (debug == 11) {
         printf("\n\n*********** Gene tree after updating stree in NNI, locus %d\n", locus);
         if (printGtree(1))
            error2("Gene tree error after NNI");
      }
      if (nMoved[locus]) {  /* calculate likelihood only if regrafting in genetree. */
         lnpDinew[locus] = lnpD_locus(locus);
         lnLd += lnpDinew[locus] - data.lnpDi[locus];
      }
      else
         lnpDinew[locus] = data.lnpDi[locus];

      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }  /* for(locus) */
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnpGd + lnLd;

   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      accepted = 1;
      for (locus = 0; locus < data.ngene; locus++) {
         data.lnpDi[locus] = lnpDinew[locus];
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
      }
      data.lnpSpeciesModel = lnpSpeciesModelnew;
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      if (debug == 11) printf("Species tree NNI move accepted");
   }
   else { /* reject */
      stree = stree0;
      memmove(gnodes[0], gnodes_t[0], sizeGtrees);
      copyCoalescent(1);
      if (debug == 11) printf("NNI move rejected");
   }
   return(accepted);
}


int UpdateSpeciesTreeSPR(int NNIonly, double *lnL, double* lnpG, double space[])
{
   /* This uses SPR or NNI to update the species tree topology while modifying the gene trees
      to avoid conflict (figure 2 in Yang & Rannala 2014 MBE).
      (a) SPR change of the gene trees.  For each locus, I scan the nodes on the gene tree
      to identity Movednodes, and nodes whose ipop need changing (AB2B), collect the
      necessary info (such as ntarget, nsource, and targetChosen, etc) for each Movednode
      without channging the gene tree.  After all affected nodes for the locus are
      identified, the changes are applied.
      (b) SPR change to the species tree.
      (c) calculation of the prior and likelihood to accept/reject the move.
      Other notes:
      com.oldconP[] is set to 0, so that repeated likelihood calculation may be performed.

      pathA[] & pathC[].  pathA[] includes Y, X, ..., and pathC has C, ....  Ancestor Z is not included.
      Min number on the two paths is 2, for NNI and with Z = X.
   */
   int s = stree.nspecies, ysona, xsony, locus, inode, dad, sons[2], is, i, j, k, ipopA, ipopC;
   int x, y, z, a, b, c, c0, zt, ndspecies, nAsons, sizeGtrees, accepted = 0, receiver;
   int nMoved[NGENE], PrunedSons[NS - 1], npathA, npathC, pathA[NSPECIES], pathC[NSPECIES];
   int nsource, ntarget, targets[NS * 2 - 1], sources[NS * 2 - 1], TargetChosen[NS - 1], *bs = targets;
   char AB2B[NS - 1];
   double ctheta = 1, atheta = data.theta_prior[0], btheta = data.theta_prior[1];  /* theta for y */
   double lnacceptance = 0, lnpGd, lnLd, *lnpGinew = space, *lnpDinew = space + data.ngene;
   double lnpSpeciesModelnew = 0, tMoved[NS - 1], t, tau1;  /* tau1 is age of y. */
   char *flag = (char*)space;   /* flag[ns*2-1] for A & B, using the same space for loci. */
   double weight[NSPECIES * 2 - 2], r;  /* weights for focus branch & target branch c. */
   struct SPECIESTREE stree0;

   if (NNIonly)
      error2("SPR tested, NNI to be implemented/tested .  26/3/2015\n");

   for (j = s, ndspecies = 1; j < s * 2 - 1; j++)
      if (stree.nodes[j].age > 0)  ndspecies++;

   if (debug == 11)
      printf("\n\n######### UpdateSpeciesTreeSPR\nSpecies tree:\n");

   stree0 = stree;
   for (locus = 0, sizeGtrees = 0; locus < data.ngene; locus++)
      sizeGtrees += (data.ns[locus] * 2 - 1) * sizeof(struct TREEN);
   memmove(gnodes_t[0], gnodes[0], sizeGtrees);
   if (mcmc.usedata && mcmc.saveconP)
      for (j = 0; j < data.maxns * 2 - 1; j++)  com.oldconP[j] = 0;
   copyCoalescent(0);

   /* choose focus branch in the species tree to apply SPR, using bl^power as weights */
   if (debug == 11) printSptreeBPP(F0);
   if (debug == 11) printf("\nBranch weights\n");
   BranchWeights(weight);
   for (i = 0, t = 0, r = rndu(); i < s - 2; i++) if (r < (t += weight[i])) break;
   y = s + i;
   lnacceptance -= log(weight[i]);

   /* SPR to move A to C. */
   ysona = (int)(2 * rndu());    /* this is A.   */
   a = stree.nodes[y].sons[ysona];
   b = stree.nodes[y].sons[1 - ysona];
   x = stree.nodes[y].father;
   xsony = (y == stree.nodes[x].sons[0] ? 0 : 1);
   c0 = stree.nodes[x].sons[1 - xsony];
   tau1 = stree.nodes[y].age;
   if (debug == 11) printf("\nMoving %d %s", a, stree.nodes[a].name);

   /* sample potential targets according to weights (# ancestors on path, min = 2) */
   for (is = 0, ntarget = 0; is < stree.nnode; is++) {
      dad = stree.nodes[is].father;
      if (stree.pptable[is][y] || stree.nodes[is].age >= tau1 || stree.nodes[dad].age <= tau1)
         continue;
      for (zt = dad; zt != -1; zt = stree.nodes[zt].father)
         if (stree.pptable[x][zt]) break;
      weight[ntarget] = 1;  /* should this be 2? */
      for (i = y; i != zt; i = stree.nodes[i].father) weight[ntarget]++; /* pathA */
      for (i = is; i != zt; i = stree.nodes[i].father) weight[ntarget]++; /* pathC */
      targets[ntarget++] = is;
   }

   if (debug == 11) {
      printf("\ntargets & weights");
      matIout(F0, targets, 1, ntarget);
      matout(F0, weight, 1, ntarget);
   }
   for (i = 0, t = 0; i < ntarget; i++) t += weight[i] = 1 / weight[i];
   for (i = 0; i < ntarget; i++) weight[i] /= t;
   for (i = 0, t = 0, r = rndu(); i < ntarget - 1; i++) if (r < (t += weight[i])) break;
   c = targets[i];
   lnacceptance -= log(weight[i]);

   /* c = targets[(int)(ntarget*rndu())]; */

   for (z = stree.nodes[c].father; z != -1; z = stree.nodes[z].father)
      if (stree.pptable[x][z]) break;
   for (is = y, npathA = 0; is != z; is = stree.nodes[is].father)
      pathA[npathA++] = is;
   for (is = c, npathC = 0; is != z; is = stree.nodes[is].father)
      pathC[npathC++] = is;

   if (debug == 11) {
      printf("\n\nMoving %d %s to %d %s (chosen from ", a, stree.nodes[a].name, c, stree.nodes[c].name);
      for (i = 0; i < ntarget; i++) printf(" %d %s", targets[i], stree.nodes[targets[i]].name);
      printf(") MCA: %d %s\npathA (%2d): ", z, stree.nodes[z].name, npathA);
      for (i = 0; i < npathA; i++)  printf(" %2d", pathA[i]);
      printf("\npathC (%2d): ", npathC);
      for (i = 0; i < npathC; i++)  printf(" %2d", pathC[i]);
   }

   /* SPR move to modify gene tree topologies */
   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, 0, 1);  /* copy gtree & use alternative space for conP */
      if (debug == 11) {
         printf("\n\n******* Gene tree before SPR, species node ages (tau)\n");
         printf("tau_%d = %9.6f\ntau_%d = %9.6f\n", x, stree.nodes[x].age, y, tau1);
         printf("\nlocus %d (ns = %d)\n", locus + 1, data.ns[locus]);
         if (printGtree(1)) error2("gene tree error before SPR");
      }
      nMoved[locus] = 0;
      for (i = 0; i < com.ns - 1; i++)    AB2B[i] = (char)0;
      for (i = 0; i < com.ns * 2 - 1; i++)  sources[i] = -1;
      for (i = 0; i < com.ns * 2 - 1; i++)  flag[i] = 0;
      /* SetSonNodeFlags(a, z, flag); */  /* set flags for populations A (or ancestor) & B (or ancestor) */
      for (i = 0; i < com.ns; i++) {           /* 1st bit is for A (or ancestor) */
         if (stree.pptable[nodes[i].ipop][a] == 0) continue;
         for (flag[i] = 1, j = nodes[i].father; flag[j] == 0; j = nodes[j].father) {
            if (nodes[j].ipop == z || stree.pptable[z][nodes[j].ipop])
               break;
            flag[j] = 1;
            if (j == tree.root) break;
         }
      }
      for (i = 0; i < com.ns; i++) {  /* 2nd bit is for sibling, counsin, (B, ..., ) */
         if (stree.pptable[nodes[i].ipop][a]) continue;
         for (j = 0; j < npathA; j++) if (stree.pptable[nodes[i].ipop][pathA[j]]) break;
         if (j == npathA) continue;
         for (flag[i] += 2, j = nodes[i].father; flag[j] < 2; j = nodes[j].father) {
            if (nodes[j].ipop == z || stree.pptable[z][nodes[j].ipop])
               break;
            flag[j] += 2;
            if (j == tree.root) break;
         }
      }

      if (debug == 11) {
         printf("Flags\n");
         for (i = com.ns; i < tree.nnode; i++)
            printf("node %3d  flag  %d%d\n", i, flag[i] & 1, flag[i] >> 1);
      }
      /* Identify Moved nodes.  This loop does not change the gene tree. */
      for (inode = com.ns; inode < tree.nnode; inode++) {
         /* Moved node i satisfies the following conditions:
           It is on pathA, i.e., it is in population AB (y) or ancestor.
           One and only one of its son nodes is pure A and has descendents in A only.
           Node in pop C with age between (tau1, tau0) changes its ipop from C to AC.
           Node with both sons having B flags changes its ipop from AB to B.
         */

         for (i = 0; i < npathA; i++) if (nodes[inode].ipop == pathA[i]) break;
         if (i == npathA)
            continue;
         ipopA = pathA[i];
         sons[0] = nodes[inode].sons[0];
         sons[1] = nodes[inode].sons[1];
         /* AB -> B, square nodes have both sons with B descendents. */
         if (ipopA == y && flag[sons[0]] >> 1 && flag[sons[1]] >> 1) {
            AB2B[inode - com.ns] = (char)1;  continue;
         }
         /* Identify Moved nodes, which have 1 A son. */
         nAsons = 0;
         if (flag[sons[0]] == 1) { nAsons++;  j = sons[0]; }
         if (flag[sons[1]] == 1) { nAsons++;  j = sons[1]; }
         if (flag[inode] != 3 || nAsons != 1) continue;    /* both flags are set and 1 pure-A son  */
         PrunedSons[nMoved[locus]++] = j;
         tMoved[nMoved[locus] - 1] = t = nodes[inode].age;
         for (i = npathC - 1; i >= 0; i--) if (t > stree.nodes[pathC[i]].age) break;
         ipopC = pathC[i];

         if (debug == 11)
            printf("\nMovednode %2d (sons %2d %2d), time = %9.5f (from %d %s to %d %s)",
               inode, j, nodes[inode].sons[0] + nodes[inode].sons[1] - j, t, ipopA, stree.nodes[ipopA].name, ipopC, stree.nodes[ipopC].name);

         /* count feasible branches in source and target. */
         /* determine the pop on pathC at time t.
            target branches must be descendent of ipopC and must be on pathC.
         */
         for (j = 0, ntarget = 0; j < tree.nnode; j++) {
            if (nodes[j].age >= t || nodes[nodes[j].father].age <= t) continue;
            if (stree.pptable[nodes[j].ipop][ipopC])
               targets[ntarget++] = j;
         }
         /** If no C seqs exist at the locus, there won't be target branch on gene tree */
         if (ntarget == 0) {
            if (locus) memmove(gnodes[0], gnodes_t[0], sizeGtrees);  /* all gene trees */
            return(2);
         }

         TargetChosen[nMoved[locus] - 1] = targets[(int)(ntarget*rndu())];
         if (debug == 11) {
            printf("\nnode %d: %d feasible lineages in target (chosen %d): ", inode, ntarget, TargetChosen[nMoved[locus] - 1]);
            for (j = 0; j < ntarget; j++) printf(" %2d", targets[j]);  FPN(F0);
         }

         /* feasible lineages in source:
         The branch covers t.  Its ipop is in pathA or is a descendent of pathA. */
         nsource = 1; sources[0] = sons[0] + sons[1] - PrunedSons[nMoved[locus] - 1];
         for (j = 0; j < tree.nnode; j++) {
            if (j == sources[0] || j == nodes[sources[0]].father) continue;
            if (nodes[j].age >= t || nodes[nodes[j].father].age <= t) continue;
            if (stree.pptable[nodes[j].ipop][ipopA] && flag[j] != 1)  /* not all red */
               sources[nsource++] = j;
         }
         if (debug == 11) {
            printf("node %d: %d feasible lineages in source: ", inode, nsource);
            for (j = 0; j < nsource; j++) printf(" %2d", sources[j]);  FPN(F0);
         }

         lnacceptance += log((double)ntarget / nsource);
      }  /* for(inode), gene tree nodes */

      /* All Moved nodes for the locus are identified.  Now apply SPR to the gene tree.
         Moved node has two sons: sons[0]=PrunedSons[i] is all A and is pruned off and moved to target.
         After gene trees for all loci are changed, apply SPR on the species tree.
      */
      for (i = 0; i < nMoved[locus]; i++) {
         /* Pruning: prune off subtree at sons[0]. */
         sons[0] = PrunedSons[i];  t = tMoved[i];
         for (j = npathC - 1; j >= 0; j--) if (t > stree.nodes[pathC[j]].age) break;
         ipopC = pathC[j];
         inode = nodes[sons[0]].father;
         sons[1] = nodes[inode].sons[0] + nodes[inode].sons[1] - sons[0];
         dad = nodes[inode].father;
         for (j = 0; j < 2; j++)
            if (nodes[dad].sons[j] == inode) break;
         nodes[dad].sons[j] = sons[1];  /* pruning off the subtree at sons[0]. */
         nodes[sons[1]].father = dad;   /* sons[1] changes father from inode to dad */

         /* Regrafting: subtree at sons[0] is inserted on branch receiver */
         for (receiver = TargetChosen[i]; ; receiver = dad) {
            dad = nodes[receiver].father;
            if (nodes[dad].age > t) break;
         }
         for (j = 0; j < 2; j++)
            if (nodes[dad].sons[j] == receiver) break;
         nodes[dad].sons[j] = inode;  /* check nodes[inode].age */
         if (nodes[inode].sons[0] == sons[0]) nodes[inode].sons[1] = receiver;
         else                              nodes[inode].sons[0] = receiver;
         nodes[inode].father = dad;
         nodes[receiver].father = inode;
         nodes[inode].ipop = ipopC;

         if (debug == 11) {
            printf("\n*** After Movednode %2d (son %2d) is moved to target %d", inode, sons[0], receiver);
            printGtree(1);
         }
      }  /* for(i=0; i<nMoved[locus]; ) */

      /* reset ipop.
         square nodes: AB -> B.
         diamond: C -> AC.
         circle&triangle: AB -> AC or ancestor.
         */
      for (inode = com.ns; inode < tree.nnode; inode++) {
         t = nodes[inode].age;
         if (AB2B[inode - com.ns])
            nodes[inode].ipop = b;
         else if (nodes[inode].ipop == c && t > tau1)                       /* diamond: C -> AC */
            nodes[inode].ipop = y;
         else if (flag[inode] == 1 && t > tau1 && t < stree.nodes[z].age) {  /* circle&triangle: AB -> AC or ancestor */
            for (i = npathC - 1; i >= 1; i--) if (t > stree.nodes[pathC[i]].age) break;
            if (i >= 1)                                                   /* species on pathC */
               nodes[inode].ipop = pathC[i];
            else                                                       /* species AC */
               nodes[inode].ipop = y;
         }
      }
      memcpy(gnodes[locus], nodes_t, (com.ns * 2 - 1) * sizeof(struct TREEN));
   }  /* for(locus) */

   /* Change species tree, update pptable[], concatenate species names, etc. */
   dad = stree.nodes[c].father;
   stree.nodes[c].father = y;  stree.nodes[y].sons[1 - ysona] = c;
   stree.nodes[y].father = dad;
   stree.nodes[b].father = x;  stree.nodes[x].sons[xsony] = b;
   for (i = 0; i < 2; i++) if (stree.nodes[dad].sons[i] == c) stree.nodes[dad].sons[i] = y;
   for (i = s; i < 2 * s - 1; i++) stree.nodes[i].name[0] = '\0';
   DownSptreeSetSpnames(stree.root, 1);
   SetupPopPopTable(debug == 11);
   if (debug == 11) printSptreeBPP(F0);
   /* probability of choosing focus branch in reverse move  */
   BranchWeights(weight);
   lnacceptance += log(weight[y - s]);

   /* probability of sampling stree target branches in reverse move  */
   for (is = 0, ntarget = 0; is < stree.nnode; is++) {
      dad = stree.nodes[is].father;
      if (stree.pptable[is][y] || stree.nodes[is].age >= tau1 || stree.nodes[dad].age <= tau1)
         continue;
      if (is == b) k = ntarget;
      for (zt = dad; zt != -1; zt = stree.nodes[zt].father)
         if (stree.pptable[y][zt]) break;      /* y is father of A&C after move */
      weight[ntarget] = 1;
      for (i = y; i != zt; i = stree.nodes[i].father)  weight[ntarget]++; /* pathA */
      for (i = is; i != zt; i = stree.nodes[i].father)  weight[ntarget]++; /* pathC */
      ntarget++;
   }
   if (debug == 11) {
      printf("\nsource branches & weights");
      matout(F0, weight, 1, ntarget);
   }
   for (i = 0, t = 0; i < ntarget; i++) t += weight[i] = 1 / weight[i];
   lnacceptance += log(weight[k] / t);

   lnpSpeciesModelnew = lnpriorSpeciesModel(ndspecies);
   lnacceptance += lnpSpeciesModelnew - data.lnpSpeciesModel;

   for (locus = 0, lnLd = lnpGd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 1);  /* copy gtree & use alternative space for conP */
      if (debug == 11) {
         printf("\n\n*********** Gene tree after updating stree in SPR, locus %d\n", locus);
         if (printGtree(1))  error2("Gene tree error after SPR");
      }
      if (nMoved[locus]) {  /* calculate likelihood only if regrafting in genetree. */
         lnpDinew[locus] = lnpD_locus(locus);
         lnLd += lnpDinew[locus] - data.lnpDi[locus];
      }
      else
         lnpDinew[locus] = data.lnpDi[locus];

      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }  /* for(locus) */
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnpGd + lnLd;
   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      accepted = 1;
      for (locus = 0; locus < data.ngene; locus++) {
         data.lnpDi[locus] = lnpDinew[locus];
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
      }
      data.lnpSpeciesModel = lnpSpeciesModelnew;
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      if (debug == 11) printf("Speciestree SPR move accepted");
   }
   else { /* reject */
      stree = stree0;
      memmove(gnodes[0], gnodes_t[0], sizeGtrees);
      stree.nodes[y].theta /= ctheta;
      copyCoalescent(1);
      if (debug == 11) printf("Speciestree SPR move rejected");
   }
   return(accepted);
}

double lnPDFpow(double x, double b, double lambda)
{
   /* log PDF for power distribution: f(x; b, lambda) = lambda/b *(x/b)^(lambda - 1),  0 < x < b.  */
   return log(lambda / b) + (lambda - 1)*log(x / b);
}
double lnPDFexp(double x, double a, double rate)
{
   /* log PDF for exponential distribution: f(x; a, mean) = rate *exp(-rate*(x - a),  x > a. */
   return log(rate) - (x - a)*rate;
}

/* UpdateSpeciesTreeNodeSlider()
   This uses NodeSlider to update the species tree topology while modifying gene trees
   to avoid conflict.  It prunes off clade A on stree and re-inserts on branch C.
   It consists of the following steps:
   (A) Choose an internal branch on stree using weights and flip a coin to go up or down.
       In the Expand move, the chosen branch is XY.  There will be only one target.
       In the Shrink move, the chosen branch is YB.  Note Y has daughters A & B.
       Move branch YA to a descendent of B, and there will be multiple target branches.
       Generate new age of Y (tau_Y) and slide YA up or down the stree, to find the target
       branch C for re-grafting.
       RWay[] has species on path in reverse move: Y, C-ancestors&root.  It is used to
       identify ipopTarget & to set up ipop within pure-A clade.

   (B) Cycle through all loci.  For each locus, set up flagA, identity all affected nodes
       and their targets in the gene tree.  Then apply the SPR operations for the locus.

       flagA[]: 0 if black (on skeleton); 1: red (pure-A); 2: black-red (pure-A clade origin).
       flagA=2 identifies movednodes.  For those, count ntarget, nsource, and sample target.

       Variables targets[], flagA[ns*2-1] etc. are shared for loci.

       There may be multiple Movednodes, and multiple pruning and regrafting operations.
       Prune off all affected nodes (or A clades), and "lay them on the ground".
       Determine their new ages (rescaled by tau*_Y/tau_Y), and identify and mark
       reattachment points on the genetree skeleton.  Reattach pruned branches back on
       to genetree.  The order at which the affected nodes are pruned off or reattached
       is unimportant.  We do not reattach one pruned branch onto another, but it is
       possible for a pruned branch to be regrafted to the same branch on skeleton,
       in which case only node ages change but not tree topology.  It is possible for
       multiple pruned branches to be inserted onto same branch on skeleton.

   (C) Apply SPR operation on species tree, after all gene trees have been changed.

   (D) Calculate prior and likelihood to accept/reject the move.

   Other notes:
   (1) com.oldconP[] is right now set to 0.  Change this to avoid duplicated calculation
*/

int NodeSlider_ProposeSpeciesTree(int *a, int *c, double *tauYnew, int RWay[], double tauRWay[], double *lnacceptance)
{
   /* This samples branch YA & target C on stree, sets up RWay, and updates lnacceptance.
   */
   int    s = stree.nspecies, ysona, xsony, is, i, k;
   int    x, y, b, ntarget = 1, nsource = 1, targets[NSPECIES * 2 - 1];
   double tau, taufactor, r, t;
   double lambda = log(mcmc.ShrinkRatio) / log(1 - mcmc.ShrinkRatio), xage, bage, weight[NSPECIES - 1];

   if (debug == 12) printSptreeBPP(F0);
   if (debug == 12) printf("\nBranch weights");
   BranchWeights(weight);
   for (is = 0, t = 0, r = rndu(); is < s - 2; is++) if (r < (t += weight[is])) break;
   *lnacceptance -= log(weight[is]);

   if (rndu() > 0.5) {  /* Expand: tau*_Y > tau_X, exponential. ntarget=1, nsource>1 */
      y = s + is;    /* The chosen interior branch is XY. */
      ysona = (int)(2 * rndu());
      *a = stree.nodes[y].sons[ysona];  b = stree.nodes[y].sons[1 - ysona];
      tau = stree.nodes[y].age;
      x = stree.nodes[y].father;
      xsony = (y == stree.nodes[x].sons[0] ? 0 : 1);
      xage = stree.nodes[x].age;
      *tauYnew = xage + rndexp(mcmc.ExpandRatio*xage);
      for (i = x; i != -1; i = stree.nodes[i].father) {
         if (*tauYnew < stree.nodes[i].age) break;
         targets[0] = *c = i;   /* *c is finally set the last time this is visited. */
      }
      *lnacceptance -= lnPDFexp(*tauYnew, xage, 1 / (mcmc.ExpandRatio*xage));      /* tauY PDF */
      *lnacceptance += lnPDFpow(tau, stree.nodes[*c].age, lambda);  /* tauY PDF, reverse */
      /* source, in reverse move, is is descendent of C, not descendent of A or B. */
      for (i = 0; i < stree.nnode; i++) {  /* first source is B. */
         if (stree.pptable[i][*c] && !stree.pptable[i][y]
            && tau > stree.nodes[i].age && tau < stree.nodes[stree.nodes[i].father].age)
            nsource++;
      }
      *lnacceptance -= log(0.5);   /* probability of choosing A */
      if (debug == 12) printf("\nchosen X %d - Y %d with prob %8.5f", x, y, weight[is]);
   }
   else {            /* Shrink: tau*_Y < tau_B.  power, ntarget>=1, nsource=1 */
      b = s + is;    /* The chosen interior branch is YB. */
      bage = stree.nodes[b].age;
      y = stree.nodes[b].father;
      ysona = (stree.nodes[y].sons[0] == b ? 1 : 0);
      x = stree.nodes[y].father;
      *a = stree.nodes[y].sons[ysona];
      tau = stree.nodes[y].age;
      *tauYnew = bage*pow(rndu(), 1 / lambda);
      for (i = 0, ntarget = 0; i < s * 2 - 1; i++) {  /* target i must be descendent of B. */
         if (stree.pptable[i][b] && stree.nodes[i].age< *tauYnew
            && stree.nodes[stree.nodes[i].father].age> *tauYnew)
            targets[ntarget++] = i;
      }
      *c = targets[(int)(ntarget*rndu())];
      xage = stree.nodes[stree.nodes[*c].father].age;           /* tau_X for reverse move */
      *lnacceptance -= lnPDFpow(*tauYnew, bage, lambda);     /* tauY PDF */
      *lnacceptance += lnPDFexp(tau, xage, 1 / (mcmc.ExpandRatio*xage));   /* tauY PDF, reverse */
      *lnacceptance += log(0.5);               /* probability of choosing A in reverse move */
      if (debug == 12) printf("\nchosen Y %d - B %d with prob %8.5f", y, b, weight[is]);
   }
   *lnacceptance += log((double)ntarget / nsource);
   taufactor = *tauYnew / tau;

   if (debug == 12) {
      if (x != -1) printf("\ntau_X(node %d) = %8.5f", x, stree.nodes[x].age);
      printf("\ntau_Y(node %d) = %8.5f -> %8.5f (factor %8.5f)\n", y, tau, *tauYnew, taufactor);
      printf("tau_A(node %d) = %8.5f tau_B(node %d) = %8.5f", *a, stree.nodes[*a].age, b, stree.nodes[b].age);
   }

   /* set up RWay[] & tauRWay[] */
   for (i = 0; i < s; i++) RWay[i] = -1;
   RWay[0] = y;
   tauRWay[0] = *tauYnew;
   for (i = stree.nodes[*c].father, k = 1; i != -1; i = stree.nodes[i].father) {
      if (i == y) continue;
      RWay[k] = i;
      tauRWay[k++] = stree.nodes[i].age;
   }

   if (debug == 12) {
      printf("\nMoving branch %d %s-%d to branch %d (tau*_Y = %8.5f fact=%8.5f)", y, stree.nodes[y].name, *a, *c, *tauYnew, taufactor);
      printf("\n  RWay");   for (i = 0; RWay[i] != -1; i++) printf("%9d", RWay[i]);
      printf("\n  age ");   for (i = 0; RWay[i] != -1; i++) printf("%9.5f", tauRWay[i]);
      printf("\ntarget %d %s chosen from:", *c, stree.nodes[*c].name);
      for (i = 0; i < ntarget; i++) printf(" %d %s", targets[i], stree.nodes[targets[i]].name);
   }
   return(0);
}

int NodeSlider_AcceptSpeciesTree(int a, int c, double taunew, int *nscale, double *lnacceptance)
{
   /* This applies SPR-NodeSlider on stree, scales clade A, & sample branch in reverse move.
      nscale is changed if there are nodes inside A clade.
   */
   int    s = stree.nspecies, x, y, b, ysona, xsony = -1, dad, i, k;
   double tau, taufactor, weight[NSPECIES];  /* weights for stree branches. */

   if (debug == 12)
      printf("\n*** Change stree.  Gene trees for all loci have been changed.\n");

   /* recover known info from a & c */
   y = stree.nodes[a].father;   x = stree.nodes[y].father;
   ysona = (stree.nodes[y].sons[0] == a ? 0 : 1);
   b = stree.nodes[y].sons[1 - ysona];
   if (y != stree.root)
      xsony = (stree.nodes[x].sons[0] == y ? 0 : 1);
   tau = stree.nodes[y].age;
   taufactor = taunew / tau;

   /* apply SPR on stree */
   dad = stree.nodes[c].father;
   stree.nodes[c].father = y;  stree.nodes[y].sons[1 - ysona] = c;
   stree.nodes[y].father = dad;
   stree.nodes[b].father = x;
   if (x != -1) stree.nodes[x].sons[xsony] = b;
   stree.nodes[y].age = taunew;
   if (c != stree.root) {
      for (i = 0; i < 2; i++)
         if (stree.nodes[dad].sons[i] == c)  stree.nodes[dad].sons[i] = y;
   }
   if (y == stree.root)   stree.root = b;
   else if (c == stree.root)   stree.root = y;
   for (i = s; i < 2 * s - 1; i++) { /* Scale nodes inside clade A on species tree by taufactor */
      if (stree.pptable[i][a] && stree.nodes[i].age) {
         stree.nodes[i].age *= taufactor;
         (*nscale)++;    /* this counts a, but not y. */
      }
   }
   for (i = s; i < 2 * s - 1; i++) stree.nodes[i].name[0] = '\0';
   DownSptreeSetSpnames(stree.root, 1);
   SetupPopPopTable(debug == 12);
   if (debug == 12) {
      printf("\nProposed species tree\n");
      printSptreeBPP(F0);
   }

   /* probability of choosing branch YC (k) on stree in reverse move  */
   if (debug == 12) printf("\nBranch weights in reverse");
   BranchWeights(weight);
   k = (taufactor > 1 ? c - s : y - s);
   *lnacceptance += log(weight[k]);   /* probability of choosing YC */
   if (debug == 12)
      printf("\nChosen in reverse %s (node %d) with prob %9.5f", (taufactor > 1 ? "YB" : "XY"), s + k, weight[k]);
   return(0);
}


int NodeSlider_PaintPureA(int inode, char flagA[])
{
   /* flagA[] = 0 if not pure-A (on skeleton); 1: inside pure-A clade; 2: pure-A clade origin.
      daughters ==> inode
            1 1 ==> 1
            0 1 ==> 2
            1 2 ==> 2
            0 2 ==> 0
            0 0 ==> 0
            2 2 ==> 0
   */
   int i, ison, *sons = nodes[inode].sons;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].nson)
         NodeSlider_PaintPureA(ison, flagA);
   }
   if (flagA[sons[0]] == 1 && flagA[sons[1]] == 1)     flagA[inode] = 1;
   else if (abs(flagA[sons[0]] - flagA[sons[1]]) == 1)  flagA[inode] = 2;
   else                                                 flagA[inode] = 0;
   return(0);
}

int NodeSlider_ScaleAClade(int inode, int RWay[], double tauRWay[], double taufactor)
{
   /* This multiplies all node ages inside gene-tree clade inode by taufactor.
      This also sets the ipop for nodes inside the clade (inode).  The ipop is not changed
      if node is younger than tau*_Y (because all genetree nodes are scaled by tau*_Y/tau_Y).
      If node is older than tau*_Y, the population should be on RWay.
   */
   int i, nnodesScaled = 1;  /* age of inode rescaled */

   nodes[inode].age *= taufactor;
   if (nodes[inode].age > tauRWay[0]) {  /* tau_Y */
      for (i = 1; RWay[i] != -1; i++) if (nodes[inode].age < tauRWay[i]) break;
      nodes[inode].ipop = RWay[i - 1];
   }

   for (i = 0; i < nodes[inode].nson; i++) {
      if (nodes[nodes[inode].sons[i]].nson)
         nnodesScaled += NodeSlider_ScaleAClade(nodes[inode].sons[i], RWay, tauRWay, taufactor);
   }
   return nnodesScaled;
}

int UpdateSpeciesTreeNodeSlider(double *lnL, double* lnpG, double finetune, double space[])
{
   int    s = stree.nspecies, ndspecies, ysona, xsony, locus, inode, dad, sons[2], i, j, k, im;
   int    ipopA, ipopC;  /* ipop for moved node at source & target */
   int    x, y, a, b, c, sizeGtrees, accepted = 0, nscale = 0;
   int    ntarget = 0, nsource, source, targets[NS * 2 - 1], target, *bs = targets, isonA, newdad;
   int    Aroot[NS - 1], TargetChosen[NS - 1], ipopTarget[NS - 1], RWay[NSPECIES];
   double tauRWay[NSPECIES], tMoved[NS - 1];
   double lnacceptance = 0, lnpGd, lnLd, *lnpGinew = space, *lnpDinew = space + data.ngene;
   double lnpSpeciesModelnew = 0, t, tnew, tau, taunew, taufactor;
   struct SPECIESTREE stree0 = stree;  /* copy of stree */
   int    *root0 = (int*)(space + 2 * data.ngene), *nMoved = root0 + data.ngene;
   char   *flagA = (char*)(nMoved + data.ngene), *flagGchange = flagA + data.maxns * 2 - 1;
   double atau = data.tau_prior[0], btau = data.tau_prior[1];
   double tau0 = stree.nodes[stree.root].age, tau0new;

   for (j = s, ndspecies = 1; j < s * 2 - 1; j++)
      if (stree.nodes[j].age > 0)  ndspecies++;
   if (ndspecies <= 2) return(2);
   copyCoalescent(0);

   if (debug == 12) {
      printf("\n\n######### UpdateSpeciesTreeNodeSlider\n");
      printf("\n**** (A) Identify branch YA & target branch C on species tree\n");
   }
   NodeSlider_ProposeSpeciesTree(&a, &c, &taunew, RWay, tauRWay, &lnacceptance);
   y = stree.nodes[a].father;
   x = stree.nodes[y].father;
   ysona = (stree.nodes[y].sons[0] == a ? 0 : 1);
   b = stree.nodes[y].sons[1 - ysona];
   if (y != stree.root)
      xsony = (stree.nodes[x].sons[0] == y ? 0 : 1);
   tau = stree.nodes[y].age;
   taufactor = taunew / tau;

   for (locus = 0, sizeGtrees = 0; locus < data.ngene; locus++)
      sizeGtrees += (data.ns[locus] * 2 - 1) * sizeof(struct TREEN);
   memmove(gnodes_t[0], gnodes[0], sizeGtrees);
   for (locus = 0; locus < data.ngene; locus++) root0[locus] = data.root[locus];
   memset(flagGchange, 0, data.ngene * sizeof(char));

   if (mcmc.usedata && mcmc.saveconP)  /* perhaps change this to save computation? */
      for (j = 0; j < data.maxns * 2 - 1; j++)  com.oldconP[j] = 0;

   /* (B) Scan each gene tree to identify Movednodes & apply SPR.  Cycle through loci. */
   if (debug == 12)
      printf("\n\n**** (B) Identify Movednodes & apply SPR operations\n");
   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, 0, 1);  /* copy gtree & use alternative space for conP */
      if (debug == 12) {
         printf("\n**** (B1) locus %d (ns = %d), before SPR\n", locus, data.ns[locus]);
         if (printGtree(1)) error2("gene tree error before SPR-NodeSlider");
      }
      nMoved[locus] = 0;
      for (i = 0; i < com.ns; i++) flagA[i] = stree.pptable[nodes[i].ipop][a];
      NodeSlider_PaintPureA(tree.root, flagA);
      if (debug == 12)
         for (j = 0; j < tree.nnode; j++) printf("\nnode %3d  flag  %d", j, flagA[j]);

      /* (B1) Identify Moved nodes.  This loop does not change the gene tree. */
      for (inode = com.ns; inode < tree.nnode; inode++) {
         if (flagA[inode] != 2) continue;  /*  flagA=2 marks moved nodes */
         ipopA = nodes[inode].ipop;
         sons[0] = nodes[inode].sons[0];
         sons[1] = nodes[inode].sons[1];
         isonA = (flagA[sons[0]] == 1 ? 0 : 1);
         t = nodes[inode].age;
         Aroot[nMoved[locus]] = nodes[inode].sons[isonA];
         tMoved[nMoved[locus]] = tnew = t*taufactor;
         /* Use RWay to determine ipopC, ipop at reattachment point in the new stree. */
         for (i = 1; RWay[i] != -1; i++)  if (tnew < tauRWay[i]) break;
         ipopTarget[nMoved[locus]++] = ipopC = RWay[i - 1];
         if (debug == 12)
            printf("\nMovednode %2d (sons %2d %2d), time = %8.5f ->%8.5f (from %d %s to pop %d)",
               inode, sons[isonA], sons[1 - isonA], t, tnew, ipopA, stree.nodes[ipopA].name, ipopC);

         /* List target branches on gene tree and sample one for reattachment.
            Feasible target branch must pass ipopC, but there are two complications
            because ipopC is for the new stree while pptable is for the old stree.
            If ipopC = Y in new stree, the branch must pass C on old stree.  (C->AC=Y)
            If ipopC = B in new stree, the branch must pass AB on old stree. (AB->B)
         */
         k = ipopC;
         if (k == y) k = c;
         if (k == b) k = y;
         for (j = 0, ntarget = 0; j < tree.nnode; j++) {
            if (nodes[j].age >= tnew || (j != tree.root && nodes[nodes[j].father].age <= tnew))
               continue;
            if (flagA[j] != 1 && stree.pptable[nodes[j].ipop][k])  /* flagA must be 0 or 2. */
               targets[ntarget++] = j;
         }
         /** If no C seqs exist at the locus, there won't be target branch on gene tree,
             and the move is disabled.  This will never happen when tau_Y increases, and may
             happen when tau_Y decreases in the move.
         */
         if (ntarget == 0) {
            if (taufactor > 1)
               error2("jgl: this should not happen, for taufactor>1.  Let me know.");
            if (locus) memmove(gnodes[0], gnodes_t[0], sizeGtrees);  /* all gene trees */
            memmove(data.root, root0, data.ngene * sizeof(int));
            return(2);
         }

         target = targets[(int)(ntarget*rndu())];
         if (debug == 12) {
            printf("\nnode %d: %d feasible target branches (chosen %d): ", inode, ntarget, target);
            for (j = 0; j < ntarget; j++) printf(" %2d", targets[j]);  FPN(F0);
         }

         /* Reset target, if it has flag 2, by tracing towards tips until flag=0. */
         for (k = 0; ; ) {
            if (flagA[target] == 0) break;
            k++;
            j = nodes[target].sons[0];
            target = (flagA[j] == 1 ? nodes[target].sons[1] : j);
         }
         if (debug == 12 && k)  printf("\nMovednode %2d, target reset to %d", inode, target);
         TargetChosen[nMoved[locus] - 1] = target;

         /* count nsource.  Feasible lineages should satisfy the following conditions:
            The branch covers t.  It is not pureA and its ipop passes ipopA. */
         nsource = 1; source = sons[1 - isonA];
         for (j = 0; j < tree.nnode; j++) {
            if (j == source || j == inode) continue;
            if (nodes[j].age >= t || (j != tree.root && nodes[nodes[j].father].age <= t)) continue;
            if (flagA[j] != 1 && stree.pptable[nodes[j].ipop][ipopA])  /* not pureA*/
               nsource++;
         }
         if (debug == 12)
            printf("node %d: nsource = %d", inode, nsource);
         lnacceptance += log((double)ntarget / nsource);
      }  /* for(inode), gene tree nodes */

      /* (B2) All Movednodes for locus have been identified.  Apply SPR on gene tree.
         For each movednode (inode), prune off sons[0]=Aroot[im] and reattach it at target.
      */
      for (im = 0; im < nMoved[locus]; im++) {  /* im loops through movednodes at the locus */
         sons[0] = Aroot[im];  inode = nodes[sons[0]].father;
         nodes[inode].age = tnew = tMoved[im];
         target = TargetChosen[im];
         nodes[inode].ipop = ipopTarget[im];

         sons[1] = nodes[inode].sons[0] + nodes[inode].sons[1] - sons[0];
         dad = nodes[inode].father;

         if (debug == 12)
            printf("\nSPR move branch %d-%d onto target %d\n", inode, sons[0], target);
         for (; ; ) {   /* possibly many insertions on same branch */
            newdad = nodes[target].father;
            if (target == tree.root || nodes[newdad].age > tnew) break;
            target = newdad;
         }
         /* check whether the genetree topology changes. */
         if (target != sons[1] && target != inode) { /* genetree topology changes */
            if (inode != tree.root) {
               if (nodes[dad].sons[0] == inode) nodes[dad].sons[0] = sons[1];
               else                          nodes[dad].sons[1] = sons[1];
               nodes[sons[1]].father = dad;
            }
            if (target != tree.root) {
               for (j = 0; j < 2; j++)
                  if (nodes[newdad].sons[j] == target) break;
               nodes[newdad].sons[j] = inode;
            }
            nodes[inode].father = newdad;
            nodes[target].father = inode;
            if (inode == tree.root)           /* old root disappears. */
               tree.root = sons[1];
            else if (target == tree.root)     /* creates new root. */
               tree.root = inode;
            if (nodes[inode].sons[0] == sons[0]) nodes[inode].sons[1] = target;
            else                              nodes[inode].sons[0] = target;
            nodes[tree.root].father = -1;  nodes[tree.root].branch = 0;
         }
         nscale++;               /* this is for scaling age of Movednode Aroot. */
         if (nodes[sons[0]].age)   /* re-scale node ages and reset ipop in clade A */
            nscale += NodeSlider_ScaleAClade(sons[0], RWay, tauRWay, taufactor);
         if (debug == 12) {
            printf("\n**** (B2) After Movednode %2d (son %2d) is moved to target %d", inode, sons[0], target);
            printf(" (tau & ipop may be incorrect)\n");
            printGtree(1);
         }
      }  /* for(im), moved nodes */

      /* recale ages and reset ipop if the whole tree is a pure-A clade.  */
      if (flagA[tree.root] == 1) {
         nscale += NodeSlider_ScaleAClade(tree.root, RWay, tauRWay, taufactor);
         flagGchange[locus] = (char)1;
      }

      /* reset ipop.  square nodes: AB -> B; diamond: C -> AC.  */
      for (inode = com.ns; inode < tree.nnode; inode++) {
         if (flagA[inode]) continue;    /* flag has to be 0, black on the skeleton. */
         if (nodes[inode].ipop == y)      /* AB -> B */
            nodes[inode].ipop = b;
         else if (nodes[inode].ipop == c && nodes[inode].age > taunew) /* C -> AC */
            nodes[inode].ipop = y;
      }

      data.root[locus] = tree.root;
      memcpy(gnodes[locus], nodes_t, (com.ns * 2 - 1) * sizeof(struct TREEN));
   }  /* for(locus) */


   /* (C) Change stree, update pptable[], concatenate species names.  All gene trees have been changed. */
   NodeSlider_AcceptSpeciesTree(a, c, taunew, &nscale, &lnacceptance);
   lnacceptance += nscale*log(taufactor);
   lnpSpeciesModelnew = lnpriorSpeciesModel(ndspecies);
   lnacceptance += lnpSpeciesModelnew - data.lnpSpeciesModel;

   tau0new = stree.nodes[stree.root].age;
   if (fabs(tau0new - tau0) > 1e-20)  /* prior on tau's (YR2010: Eq. 2) */
      lnacceptance += (-atau - 1 - (ndspecies - 2))*log(tau0new / tau0) - btau*(1 / tau0new - 1 / tau0);

   for (locus = 0, lnpGd = lnLd = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 1);  /* copy gtree & use alternative space for conP */
      if (debug == 12) {
         printf("\n\n**** Gene tree after updating stree (locus %d)\n", locus);
         if (printGtree(1))
            error2("Gene tree error after updating stree");
      }

      lnpDinew[locus] = data.lnpDi[locus];
      if (nMoved[locus] || flagGchange[locus]) {  /* calculate likelihood only if genetree has changed. */
         lnpDinew[locus] = lnpD_locus(locus);
         lnLd += lnpDinew[locus] - data.lnpDi[locus];
      }
      CountCoalescentLocus(locus);
      if (!stree.NoTheta) {
         lnpGinew[locus] = lnpG_S(locus);
         lnpGd += lnpGinew[locus] - data.lnpGi[locus];
      }
   }  /* for(locus) */
   if (stree.NoTheta)
      lnpGd = lnpG_S_NoTheta() - *lnpG;

   lnacceptance += lnpGd + lnLd;
   if (lnacceptance >= -1e-10 || rndu() < exp(lnacceptance)) {
      accepted = 1;
      for (locus = 0; locus < data.ngene; locus++) {
         data.lnpDi[locus] = lnpDinew[locus];
         if (!stree.NoTheta) data.lnpGi[locus] = lnpGinew[locus];
      }
      data.lnpSpeciesModel = lnpSpeciesModelnew;
      if (mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      *lnpG += lnpGd;
      if (debug == 12) printf("Sptree NodeSlider move accepted. ");
   }
   else { /* reject */
      stree = stree0;
      memmove(gnodes[0], gnodes_t[0], sizeGtrees);
      memmove(data.root, root0, data.ngene * sizeof(int));
      copyCoalescent(1);
      if (debug == 12) printf("Sptree NodeSlider move rejected. ");
   }
   return(accepted);
}


void checkGtree(void)
{
   int locus, inode, j, ipop, father;
   double t, tb[2];

   for (locus = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 0, 0, 0);
      for (inode = com.ns; inode < tree.nnode; inode++) {
         t = nodes[inode].age;
         father = nodes[inode].father;
         ipop = nodes[inode].ipop;
         tb[0] = stree.nodes[ipop].age;
         tb[1] = OLDAGE;
         for (j = 0; j < nodes[inode].nson; j++)
            tb[0] = max2(tb[0], nodes[nodes[inode].sons[j]].age);
         if (ipop != stree.root)
            tb[1] = stree.nodes[stree.nodes[ipop].father].age;
         if (inode != tree.root)
            tb[1] = min2(tb[1], nodes[nodes[inode].father].age);

         if (tb[1] < tb[0] || t<tb[0] || t>tb[1]) {
            printf("\n\n****** checkGtree ******\n");
            printf("locus %d node %d t=%9.6f tb: (%9.6f %9.6f)", locus, inode, t, tb[0], tb[1]);
            printf("\nspecies tree:\n");
            for (j = 0; j < 2 * stree.nspecies - 1; j++, FPN(F0)) {
               printf("species %d %-12s ", j + 1, stree.nodes[j].name);
               printf("age %10.6g ", stree.nodes[j].age);
            }
            printf("\ngene tree:\n");
            printGtree(1);
            printf("\ninconsistent gene tree at locus %d node %d\ttime %.5g (%.5g %.5g)\n",
               locus, inode, t, tb[0], tb[1]);
            exit(-1);
         }

         /*
         else if(t<tb[0])
            nodes[inode].age = tb[0];
         else if(t>tb[1])
            nodes[inode].age = tb[1];
         */
      }
   }
}

int DownSptreeSetTime(int inode)
{
   /* This moves down species tree to specify node ages, starting from root age.
      It sets the theta to -1 if the tip species is collapsed.
   */
   int j, ison;
   double prop = (stree.nspecies > 10 ? 0.9 : 0.5);
   /*
   double prop = pow(1/(stree.nspecies-1.0), 1/(stree.nspecies-2.0));
   */
   if (stree.nspecies <= 1) error2("fewer than 2 species?");
   for (j = 0; j < stree.nodes[inode].nson; j++) {
      ison = stree.nodes[inode].sons[j];
      if (stree.nodes[inode].age == 0) {   /* inode collapsed */
         stree.nodes[ison].theta = -1;
      }
      if (stree.nodes[ison].nson) {        /* ison is an interior node */
         if (stree.nodes[inode].age == 0)
            stree.nodes[ison].age = 0;

         if (stree.nodes[ison].age > 0) {
            stree.nodes[ison].age = stree.nodes[inode].age * (prop + (1 - prop - 0.02)*rndu());
         }
         DownSptreeSetTime(ison);
      }
   }
   return(0);
}

int GetInitials(void)
{
   /* This sets the initial values for starting the MCMC, and returns np, the
      number of parameters in the MCMC, to be collected in collectx().

      If (data.est_locusrate == 1), the rate at locus 1 is printed out.
      See also collectx().
   */
   int np, is, i, j, k, ntheta, ntau = stree.nspecies - 1, s = stree.nspecies;
   double atau = data.tau_prior[0], btau = data.tau_prior[1];
   double atheta = data.theta_prior[0], btheta = data.theta_prior[1], y;

   ntheta = (stree.NoTheta ? 0 : stree.npop);
   /* initial theta's.  Some may be set to -1 later.
      if(stree.NoTheta), theta are used to generate initial gene trees & then removed. */
   for (i = 0, k = 0; i < stree.npop; i++)  /* invgamma prior  */
      stree.nodes[stree.pops[i]].theta = btheta / (atheta - 1)*(0.9 + 0.2*rndu());

   /* initial tau's.  Some theta's are reset to -1. */
   if (s > 1) {
      /* initially age is used to indicate the presence of the tau parameter. */
      for (i = s; i < stree.nnode; i++) stree.nodes[i].age = 1;

      if (stree.analysis == A10) {  /* A10: starting model chosen at random from nModels */
         k = (int)(stree.nModels*rndu());
         /**** to specify the starting tree at the keyboard, uncomment the following.  ***/
         /*
         printf("starting tree (1-%d)? ", stree.nModels);
         scanf("%d", &k);
         k--;
         */
         for (j = s; j < s * 2 - 1; j++) {
            stree.nodes[j].age = stree.DelimitationModels[k*s + j - s] - '0';
         }
         printf("\nStarting species-delimitation model: %s\n", printDelimitationModel());
      }
      else if (stree.analysis == A11) {  /* A11: starting model by random collaps of a node */
         k = (int)(s*rndu());
         if (k < s - 1)
            for (i = s; i < s * 2 - 1; i++) stree.nodes[i].age = !stree.pptable[i][s + k];
      }
      if (stree.nodes[s].age)          /* species tree root age */
         stree.nodes[s].age = btau / (atau - 1)*(0.9 + 0.2*rndu());
      DownSptreeSetTime(s);  /* this initializes tau and sets theta = -1 for daughter nodes */
   }

   if (stree.analysis == A11)
      PriorS_IntegerPartitions(stree.nspecies, (stree.nspecies < 20), stree.PriorSA11);
   if (stree.speciesdelimitation == 1 || stree.speciestree == 1) {
      for (i = 0, ntheta = ntau = 0; i < stree.nnode; i++) {
         if (debug) {
            printf("node %2d %-20s ", i + 1, stree.nodes[i].name);
            printf("theta = %9.6f  tau = %9.6f\n", stree.nodes[i].theta, stree.nodes[i].age);
         }
         if (!stree.NoTheta) ntheta += (stree.nodes[i].theta > 0);
         ntau += (stree.nodes[i].age > 0);
      }
      data.lnpSpeciesModel = lnpriorSpeciesModel(ntau + 1);
   }

   if (data.est_heredity == 1) {
      for (i = 0; i < data.ngene; i++)
         data.heredity[i] = data.a_heredity / data.b_heredity*(0.8 + 0.4*rndu());
   }
   if (data.est_locusrate == 1) {
      for (i = 0; i < data.ngene; i++) {
         data.locusrate[i] = 0.8 + 0.4*rndu();
      }
      y = sum(data.locusrate, data.ngene) / data.ngene;
      for (i = 0; i < data.ngene; i++)
         data.locusrate[i] /= y;
   }

   /* sequence errors */
   if (data.nseqerr) {
      printf("\nInitials for seqerrors");

      for (is = 0; is < s; is++) {
         if (data.iseqerr[is]) {
            for (i = 0; i < 4; i++) {
               for (j = 0, y = 0; j < 4; j++) y += data.a_seqerr[is][i * 4 + j];
               for (j = 0; j < 4; j++)
                  data.e_seqerr[is][i * 4 + j] = data.a_seqerr[is][i * 4 + j] / y * (0.8 + 0.4*rndu());
               for (j = 0, y = 0; j < 4; j++) y += data.e_seqerr[is][i * 4 + j];
               for (j = 0; j < 4; j++) data.e_seqerr[is][i * 4 + j] /= y;
            }
            matout(F0, data.e_seqerr[is], 4, 4);
         }
      }
   }

   np = ntheta + ntau + data.nseqerr * 16
      + (mcmc.printlocusrate == 1) * data.ngene
      + (mcmc.printheredity == 1) * data.ngene;

   if (s == 1 && mcmc.printGenetree) np += data.ngene;   /* t_MRCA is printed for 1 species */

   return(np);
}


int collectx(int mode, FILE* fout, double x[])
{
   /* This collects parameters into x[] for printing and summarizing.  It checks the number of parameters.
      mode = 0: prints out the header line into fout.
      mode = 1: Collect stree.theta, stree.age, etc. into x[].  In other words, take an MCMC sample.
      mode = 2: Copy theta and tau from x[] into stree, for printing.  x[] should have the posterior means.

      if(com.np == -1) com.np is assigned a value here.
   */
   int i, j, k = 0, is, ipop;

   if (mode == 0 && fout == NULL)
      return(0);
   if (mode == 0) fprintf(fout, "Gen");
   if (mode == 0 && stree.speciesdelimitation == 1 && stree.speciestree == 0) fprintf(fout, "\tnp\ttree");
   /* theta's  */
   if (!stree.NoTheta) {
      for (i = 0; i < stree.npop; i++) {
         ipop = stree.pops[i];
         if (stree.nodes[ipop].theta == -1)
            continue;
         if (mode == 0) {
            if (ipop < stree.nspecies || stree.nspecies <= 10)
               fprintf(fout, "\ttheta_%d%s", ipop + 1, stree.nodes[ipop].name);
            else
               fprintf(fout, "\ttheta_%d", ipop + 1);
         }
         if (mode == 2)  stree.nodes[ipop].theta = x[k];
         else         x[k] = stree.nodes[ipop].theta;
         k++;
      }
   }
   /* tau's  */
   for (i = stree.nspecies; i < stree.nnode; i++) {
      if (stree.nodes[i].age == 0)  continue;
      if (mode == 0) {
         if (stree.nspecies <= 10) fprintf(fout, "\ttau_%d%s", i + 1, stree.nodes[i].name);
         else                       fprintf(fout, "\ttau_%d", i + 1);
      }
      if (mode == 2)  stree.nodes[i].age = x[k];
      else         x[k] = stree.nodes[i].age;
      k++;
   }
   if (mcmc.printlocusrate) {
      for (i = 0; i < data.ngene; i++) {
         if (mode == 0)       fprintf(fout, "\trate_L%d", i + 1);
         else if (mode == 1)  x[k] = data.locusrate[i];
         k++;
      }
   }
   if (mcmc.printheredity) {
      for (i = 0; i < data.ngene; i++) {
         if (mode == 0)       fprintf(fout, "\theredity_L%d", i + 1);
         else if (mode == 1)  x[k] = data.heredity[i];
         k++;
      }
   }
   if (data.nseqerr) {
      for (is = 0; is < stree.nspecies; is++)
         if (data.iseqerr[is]) {
            for (i = 0; i < 4; i++)  for (j = 0; j < 4; j++) {
               if (mode == 0)        fprintf(fout, "\tseqerr_[%s]%c%c", stree.nodes[is].name, BASEs[i], BASEs[j]);
               else if (mode == 1)   x[k] = data.e_seqerr[is][i * 4 + j];
               k++;
            }
         }
   }
   if (stree.nspecies == 1 && mcmc.printGenetree) /* t_MRCA */
      for (i = 0; i < data.ngene; i++) {
         if (mode == 0)       fprintf(fout, "\tt_MRCA_L%d", i + 1);
         else if (mode == 1)  x[k] = gnodes[i][data.root[i]].age;
         k++;
      }

   if (k != com.np) {
      if (com.np != -1)
         printf("np is changed from %d to %d in collectx().", com.np, k);
      com.np = k;
   }
   if (mode == 0 && mcmc.usedata) fprintf(fout, "\tlnL");
   if (mode == 0) fprintf(fout, "\n");
   return(0);
}


void copySptree(void)
{
   /* This copies stree into nodes = nodes_t, for printing or editing
   */
   int i, j;

   nodes = nodes_t;
   com.ns = stree.nspecies;   tree.root = stree.root;
   tree.nnode = stree.nnode;  tree.nbranch = stree.nbranch;
   for (i = 0; i < stree.nnode; i++) {
      if (i < com.ns) strcpy(com.spname[i], stree.nodes[i].name);
      nodes[i].father = stree.nodes[i].father;
      nodes[i].nson = stree.nodes[i].nson;
      for (j = 0; j < nodes[i].nson; j++)
         nodes[i].sons[j] = stree.nodes[i].sons[j];
      nodes[i].age = stree.nodes[i].age;

      nodes[i].label = stree.nodes[i].theta;
   }
   for (i = 0; i < stree.nnode; i++) {
      if (i == stree.root) nodes[i].branch = 0;
      else               nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
   }
}

void printSptreeBPP(FILE* fout)
{
   copySptree();
   /* OutTreeN(fout, 1, PrBranch); */
   OutTreeN(fout, 1, PrBranch | PrLabel);
   /* OutTreeN(fout, 1, 0); */
}

char *printDelimitationModel(void)
{
   static char cmodel[NSPECIES] = "";
   int i;

   for (i = 0; i < stree.nspecies - 1; i++)
      cmodel[i] = (stree.nodes[stree.nspecies + i].age > 0 ? '1' : '0');
   return(cmodel);
}

int NumberDelimitationModels(int inode)
{
   /* This counts the species models (Yang & Rannala 2010, p.9265 left column, near bottom)
      and also the total number of labeled histories on the guide tree.
   */
   int j, ison, Z = 1;

   for (j = 0; j < stree.nodes[inode].nson; j++) {
      ison = stree.nodes[inode].sons[j];
      if (ison >= stree.nspecies)
         Z *= NumberDelimitationModels(ison);
   }
   return(Z + 1);
}


int EnumerateDelimitationModels(void)
{
   /* This enumerates species models in stree.DelimitationModels[].  It returns nModels.
      If(stree.speciesdelimitation && !stree.speciestree), the models are sorted alphabetically
      (eg, 000, 100, 110, 111).
      If(stree.speciesdelimitation && stree.speciestree), the models are not sorted but the order does not matter.
      The algorithm uses one loop for( ; ; ) to represent s1=s-1 loops.
      Loop indices are ic[], loop limits are nc[], and current loop is Icurrent.
      Index for an older node (outer loop) determines the limits of younger nodes (inner loops).
      For example, if ic[] = 0, then nc[] for all descendent nodes are 1.
   */
   int nModels = 0, s = stree.nspecies, s1 = s - 1;
   int j, Icurrent, offset, ic[NSPECIES - 1] = { 0 }, nc[NSPECIES - 1] = { 2 };
   char debug = 0;

   if (stree.analysis != A10)
      error2("EnumerateDelimitationModels(): should not be here?");
   Icurrent = 0;
   for (; ; ) {
      /* s1-1 is index for last node in pre-order tree traversal, so it is time to print. */
      if (Icurrent == s1 - 1) {
         for (; ic[Icurrent] < nc[Icurrent]; ic[Icurrent]++) {
            offset = nModels++ * s;
            for (j = 0; j < s1; j++)
               stree.DelimitationModels[offset + j] = (char)('0' + ic[j]);
            if (debug)
               printf("\tdelimitation model %2d: %s\n", nModels, stree.DelimitationModels + offset);
         }
      }

      if (Icurrent < 0) break;
      if (ic[Icurrent] >= nc[Icurrent]) {
         if (--Icurrent >= 0)
            ic[Icurrent]++;
      }
      else {
         /*  Set nc[] = 1 or 2 for descendents of Icurrent. */
         for (j = Icurrent + 1; j < s1; j++)
            if (stree.pptable[s + j][s + Icurrent])
               nc[j] = 1 + ic[Icurrent];
         if (Icurrent < s1 - 1)
            ic[++Icurrent] = 0;
      }

      if (debug) {
         printf("ic: "); for (j = 0; j < s1; j++) printf("%d", ic[j]);
         printf("  nc: "); for (j = 0; j < s1; j++) printf("%d", nc[j]);
         printf("  Icurrent = %d\n", Icurrent);
      }
   }

   if (debug) {
      printf("\n%d models: ", nModels);
      for (j = 0; j < nModels; j++)  printf(" %s", stree.DelimitationModels + j*s);
      FPN(F0);
   }

   return(nModels);
}

double CountLHsTree2(void)
{
   /* This counts the number of labeled histories for the current rooted species tree, using
      stree.nodes[].age.  An internal node is collapsed if its age is 0.
      This is modified from CountLHsTree().
   */
   int s = stree.nspecies, i, j, k, change, *sons, LR[NSPECIES - 1][2];
   double nLH, y = 0;
   char debug = 0;

   if (debug)
      printf("\n\n*******CountLHsTree2******\n");

   for (i = 0; i < s - 1; i++)
      LR[i][0] = LR[i][1] = -1;
   for (k = 0; k < s; k++) {
      if (debug) {
         printf("\n*****Round %d\n", k + 1);
         for (i = 0; i < s - 1; i++)
            printf("\nnode %2d %d (%2d %2d): %2d %2d ", i + s, stree.nodes[i + s].age > 0, stree.nodes[i + s].sons[0], stree.nodes[i + s].sons[1], LR[i][0], LR[i][1]);
         FPN(F0);
      }
      for (i = 2 * s - 1 - 1, change = 0; i >= 0; i--) {
         if (stree.nodes[i].age == 0) continue;
         sons = stree.nodes[i].sons;
         for (j = 0; j < 2; j++) {
            if (LR[i - s][j] != -1) continue;
            if (sons[j] < s || stree.nodes[sons[j]].age == 0) { /* sons[j] is tip or collapsed. */
               LR[i - s][j] = 0;
               change = 1;
            }
            else if (LR[sons[j] - s][0] != -1 && LR[sons[j] - s][1] != -1) {
               LR[i - s][j] = LR[sons[j] - s][0] + LR[sons[j] - s][1] + 1;
               change = 1;
            }
         }
      }
      if (!change) break;
   }
   if (debug && k > 2) printf("%d rounds used in CountLHsTree2()\n", k);
   for (i = s, nLH = 1; i < s * 2 - 1; i++) {
      if (stree.nodes[i].age == 0) continue;
      if (debug)
         printf("\nnode %2d %d (%2d %2d): %2d %2d ", i, stree.nodes[i].age > 0, stree.nodes[i].sons[0], stree.nodes[i].sons[1], LR[i - s][0], LR[i - s][1]);

      if (LR[i - s][0] == -1 || LR[i - s][1] == -1)
         error2("LR = -1");
      if (LR[i - s][0] && LR[i - s][1]) {
         nLH *= Binomial((double)(LR[i - s][0] + LR[i - s][1]), LR[i - s][0], &y);
         if (y) error2("y!=0 not expected");
      }
   }
   if (debug)
      printf("\nnLH = %6.0f\n", nLH);
   if (!(nLH >= 1 && nLH < 1e300)) error2("too many labeled histories.");
   return(nLH);
}


int GetDmodelIndex(void)
{
   /* This gets the species delimitation model index, by constructing the model string in cmodel and
      searching in stree.DelimitationModels.
   */
   char cmodel[NSPECIES] = "", *model;
   int s = stree.nspecies, i;

   for (i = 0; i < stree.nspecies - 1; i++)
      cmodel[i] = (stree.nodes[s + i].age > 0 ? '1' : '0');

   model = (char*)bsearch(cmodel, stree.DelimitationModels, stree.nModels, sizeof(char)*s, (int(*)(const void*, const void*))strcmp);

   return (int)(model - stree.DelimitationModels) / (sizeof(char)*s);
}



int ProcessGtrees(FILE* fout)
{
#if (defined(__unix__) || defined(__APPLE__))
   int Unix = 1;
#else
   int Unix = 0;
#endif
   int locus;
   char filename[128], line[1024];
   FILE **fGtree, *finconsensus;

   fGtree = (FILE**)malloc(data.ngene * sizeof(FILE*));
   if (fGtree == NULL) error2("fGtree mem allocation error.");
   for (locus = 0; locus < data.ngene; locus++) {
      sprintf(line, "gtree/L%03d.tre", locus + 1);
      fGtree[locus] = (FILE*)fopen(line, "r");
      if (fGtree[locus] == NULL) {
         printf("gtree file open error for locus %d\n", locus + 1);
         exit(-1);
      }
   }

   fprintf(fout, "\n\nConsensus gene trees with posterior clade probabilities\n");

   printf("\nSummarizing gene trees, %d trees per locus\n", mcmc.nsample);
   for (locus = 0; locus < data.ngene; locus++)
      fclose(fGtree[locus]);
   for (locus = 0; locus < data.ngene; locus++) {
      if (Unix) {
         system("rm outfile outtree");
         sprintf(filename, "gtree/L%03d.tre\0", locus + 1);
      }
      else {
         system("del outfile outtree");
         sprintf(filename, "gtree\\L%03d.tre\0", locus + 1);
      }
      finconsensus = gfopen("in.consense", "w");
      fprintf(finconsensus, "%s\nC\nC\nR\nY\n", filename);
      fclose(finconsensus);

      printf("Locus %3d %s\n", locus + 1, filename);
      system("consense  < in.consense > log");
      fprintf(fout, "Locus%03d: ", locus + 1);
      appendfile(fout, "outtree");
   }

   free(fGtree);

   return(0);
}


int InitializeDelimitationModel(int prinT, double priorPmodel[])
{
   /* This allocates space for stree.pmodel and stree.SpeciedModels and lists all models for speciesdelimitation=1
   */
   int s = stree.nspecies, j, k;

   if (stree.analysis != A10) error2("A10: should not be here?");
   /* count and enumerate species models */
   stree.nModels = k = NumberDelimitationModels(stree.root);
   if (prinT) printf("Number of species-delimitation models = %2d\n", stree.nModels);
   if ((stree.pmodel = (double*)malloc(stree.nModels * sizeof(double))) == NULL) error2("oom for pmodel");
   if ((stree.DelimitationModels = (char*)malloc(k*s * sizeof(char))) == NULL) error2("oom");
   memset(stree.DelimitationModels, 0, k*s * sizeof(char));
   EnumerateDelimitationModels();

   /* prior probs for delimitation models */
   for (k = 0; k < stree.nModels; k++) {
      for (j = s; j < stree.nnode; j++)
         stree.nodes[j].age = stree.DelimitationModels[k*s + j - s] - '0';
      stree.pmodel[k] = exp(lnpriorSpeciesModel(-1)) * CountLHsTree2();  /* prior model prob */
   }
   abyx(1 / sum(stree.pmodel, stree.nModels), stree.pmodel, stree.nModels);
   if (stree.nModels <= 20) { /* list the species models and their prior probs */
      for (k = 0; k < stree.nModels; k++)
         if (prinT) printf("\tdelimitation model %3d: %s  prior %8.5f\n", k + 1, stree.DelimitationModels + k*s, stree.pmodel[k]);
   }

   if (prinT) printf("\n[Note: Ancestral nodes in order: ");
   for (k = s; k < s * 2 - 1; k++)
      if (prinT) printf(" %2d %s", k + 1, stree.nodes[k].name);
   if (prinT) printf("]\n");

   return(0);
}


void SummarizeA10_SpeciesDelimitation(FILE* fout, char mcmcf[])
{
   /* This reads the mcmc sample of species delimitation models (A10) from mcmcf and
      generate the posterior distribution of the models.
   */
   FILE *fmcmc = gfopen(mcmcf, "r");
   int s = stree.nspecies, nr, np, i, j, k, lline = 1000000;
   char *line, cmodel[NSPECIES] = { 0 };
   double PPnodes[NSPECIES - 1] = { 0 }, *postpmodel;

   printf("\n\nSummarizing the species-delimitation sample in file %s\n", mcmcf);
   fprintf(fout, "\n\nSummarizing the species-delimitation sample in file %s\n", mcmcf);

   InitializeDelimitationModel(0, stree.pmodel);
   postpmodel = (double*)malloc(stree.nModels * sizeof(double));
   line = (char*)malloc(lline * sizeof(char));
   if (postpmodel == NULL || line == NULL) error2("oom postpmodel or line");
   zero(postpmodel, stree.nModels);

   fgets(line, lline, fmcmc);  /* pop header line */
   for (nr = 0; ; nr++) {
      k = fscanf(fmcmc, "%d%d%s", &i, &np, cmodel);
      fgets(line, lline, fmcmc);  /* pop rest of line */
      if (k != 3) break;

      for (i = 0; i < stree.nspecies - 1; i++)
         stree.nodes[s + i].age = cmodel[i] - '0';
      k = GetDmodelIndex();
      postpmodel[k] ++;
   }
   for (i = 0; i < stree.nModels; i++) postpmodel[i] /= nr;
   for (j = 0; j < stree.nModels; j++) {
      for (k = 0; k < s - 1; k++)
         if (stree.DelimitationModels[j*s + k] == '1')
            PPnodes[k] += postpmodel[j];
   }
   printf("\nNumber of species-delimitation models = %2d\n", stree.nModels);
   printf("      model    prior  posterior\n");
   for (j = 0; j < stree.nModels; j++)
      printf(" %4d  %s %9.5f %9.5f\n", j + 1, stree.DelimitationModels + j*s, stree.pmodel[j], postpmodel[j]);
   printf("\n[Note: Ancestral nodes in order: ");
   for (k = s; k < s * 2 - 1; k++)
      printf(" %2d %s", k + 1, stree.nodes[k].name);
   printf("]\n");

   fprintf(fout, "\nmodel         prior   posterior\n");
   for (j = 0; j < stree.nModels; j++)
      fprintf(fout, " %4d  %s %9.5f %9.5f\n", j + 1, stree.DelimitationModels + j*s, stree.pmodel[j], postpmodel[j]);
   fprintf(fout, "\n[Note: Ancestral nodes in order: ");
   for (k = s; k < s * 2 - 1; k++)
      fprintf(fout, " %2d %s", k + 1, stree.nodes[k].name);
   fprintf(fout, "]\n");

   copySptree();
   for (k = 0; k < s * 2 - 1; k++) {
      nodes[k].label = (k < s ? 0 : PPnodes[k - s]);
      nodes[k].annotation = NULL;
   }
   printf("\nGuide tree with posterior probability for presence of nodes\n");
   OutTreeN(F0, 1, PrLabel);   FPN(F0);
   fprintf(fout, "\nGuide tree with posterior probability for presence of nodes\n");
   OutTreeN(fout, 1, PrLabel); FPN(fout);

   free(stree.pmodel);
   free(stree.DelimitationModels);
}

void Tree2PartitionDescentTreeBPP(int inode, char split[])
{
   int i, ison;

   for (i = 0; i < nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      if (nodes[ison].label != -1)
         split[(int)nodes[ison].label] = '1';
      else
         Tree2PartitionDescentTreeBPP(ison, split);
   }
}

void Tree2PartitionBPP(char splits[], int ndspecies)
{
   /* This is modified from Tree2Partition.  It generates branch partitions in splits.
      This is called from SummarizeA11_SpeciesTreeDelimitation, and an internal node may represent
      a delimited species, in which case its label is set to species index.  Otherwise label=-1.
      stree.nspecies or s is the number of populations on the guide tree, while ndspecies is
      the number of delimited species.
   */
   int s = stree.nspecies, lsplit = s + 1, nsplit = ndspecies - 2, i, j, k = 0;
   char *split;

   if (ndspecies <= 2) return;
   memset(splits, 0, nsplit*lsplit * sizeof(char));
   for (i = s; i < tree.nnode; i++) {
      if (i == tree.root) continue;
      if (nodes[i].label != -1)  /* i is a collapsed node (delimited species)  */
         continue;
      split = splits + k*lsplit;
      for (j = 0; j < ndspecies; j++) split[j] = '0';
      Tree2PartitionDescentTreeBPP(i, split);
      k++;
   }
}

int SpeciesTreeDelimitationModelRepresentation(struct SMODEL *model, int ndspecies)
{
   /* This converts the species tree, in nodes[], into a model representation in model.
      ndspecies is the number of species in the delimitation model.
      If a node (tip or internal node) represents a delimited species, its label nodes[].label
      is set to the species index.  Otherwise it is -1.  See also Tree2Partition.
      This destroys nodes[].label.
   */
   int s = stree.nspecies, s1 = s + 1, s21 = s * 2 - 1, i, j, k, ndsp = ndspecies, anc, status = 0;
   char *splits = model->splits, visited[NSPECIES * 2 - 1] = { '\0' }, *pptable = NULL;
   int debug = 0;

   if (ndspecies<1 || ndspecies>s)
      error2("#species out of range");
   else if (ndspecies != 1 || ndspecies != s) {
      if ((pptable = (char*)malloc(s21*s21 * sizeof(char))) == NULL)
         error2("oom in Representation");
   }
   if (debug) printf("\nndspecies %d\n", ndspecies);
   memset(model, 0, sizeof(struct SMODEL));

   for (i = 0; i < s * 2 - 1; i++)  nodes[i].label = -1;
   /* species names */
   if (ndspecies == 1) {   /* one species, the name is 123...s */
      for (i = 0; i < s; i++) model->spnames[0][i] = i + 1;
   }
   else if (ndspecies == s) {  /* s species, the names are 1, 2, 3, ..., s */
      for (i = 0; i < s; i++)
      {
         model->spnames[i][0] = i + 1;  nodes[i].label = i;
      }
   }
   else {
      memset(pptable, 0, s21*s21 * sizeof(char));
      for (i = 0; i < s * 2 - 1; i++)
         for (j = nodes[i].father; j != -1; j = nodes[j].father)
            pptable[i*s21 + j] = (char)1;
      if (debug) {
         printf("\npptable\n");
         for (i = 0; i < s21; i++, FPN(F0))
            for (j = 0; j < s21; j++) printf(" %2d", (int)pptable[i*s21 + j]);
      }

      for (i = 0, ndsp = 0; i < s; i++) {
         if (visited[i]) continue;
         if (nodes[i].branch > 0) {  /* tip pop i is a distinct species */
            model->spnames[ndsp][0] = i + 1;  nodes[i].label = ndsp;
         }
         else {
            /* tip pop i is part of a new species (anc).  It consists of k pops, and is
               represented by first tip pop. */
            for (anc = nodes[i].father; anc != tree.root; anc = nodes[anc].father)
               if (nodes[anc].branch > 0) break;
            for (j = i, k = 0; j < s; j++)
               if (pptable[j*s21 + anc]) {
                  model->spnames[ndsp][k++] = j + 1;
                  visited[j] = (char)1;
               }
            nodes[anc].label = ndsp;
            for (j = s; j < s * 2 - 1; j++)
               if (pptable[j*s21 + anc]) nodes[j].label = ndsp;
         }
         ndsp++;
      }
      if (ndsp > ndspecies)
         error2("# delimited species look strange.");
      else if (ndsp < ndspecies - 2)
         puts("# delimited species look strange.");
      else if (ndsp < ndspecies)
         status = 1;
   }

   /* partitions */
   if (ndsp > 2) {
      Tree2PartitionBPP(splits, ndsp);
      if (debug) {
         printf("\n%2d splits before sorting : ", ndsp - 2);
         for (i = 0; i < ndsp - 2; i++) printf(" %s", splits + i*s1);  printf("\n");
      }
      qsort(splits, ndsp - 2, s1, (int(*)(const void *, const void *))strcmp);
      if (debug) {
         printf("%2d splits after sorting  : ", ndsp - 2);
         for (i = 0; i < ndsp - 2; i++) printf(" %s", splits + i*s1);  printf("\n");
      }
   }

   model->nspecies = ndsp;
   if (debug)  PrintSmodel(F0, model, 1);
   if (ndspecies > 1 && ndspecies < s) free(pptable);
   return(status);
}

int PrintSmodel(FILE *fout, struct SMODEL *model, int printPhylogeny)
{
   int s1 = stree.nspecies + 1, ns = model->nspecies, is, i, k, maxspname = LSPNAME * 2, lspname;
   char *p, *space, *spnames0[NS];

   if ((space = (char*)malloc((ns*maxspname) * sizeof(char))) == NULL)
      error2("oom in PrintSmodel");
   memset(space, 0, ns*maxspname * sizeof(char));
   for (i = 0; i < ns; i++) spnames0[i] = com.spname[i];
   for (i = 0; i < ns; i++) com.spname[i] = space + i*maxspname;

   fprintf(fout, "%2d (", ns);
   for (i = 0; i < ns; i++) {
      lspname = 0;
      p = model->spnames[i];

      while (*p) {      /* concatenate the spnames to make up a new spname */
         is = *p - 1;  /* pop is */
         k = (int)strlen(stree.nodes[is].name);
         if (lspname + k + 1 > maxspname) error2("name too long");
         strcpy(com.spname[i] + lspname, stree.nodes[is].name);
         lspname += k;
         p++;
      }
      com.spname[i][lspname++] = 0;
      fprintf(fout, "%s", com.spname[i]);
      if (i < ns - 1) fprintf(fout, " ");
   }
   fprintf(fout, ") ");
   if (printPhylogeny) {
      for (i = 0; i < s1 - ns; i++) fprintf(fout, " ");
      Partition2Tree(model->splits, s1, ns, ns - 2, NULL);
      fprintf(fout, " ");
      OutTreeN(fout, 1, 0);

      if (stree.nspecies < 20) {
         for (i = 0; i < (s1 - ns) * 4; i++) fprintf(fout, " ");
         for (i = 0; i < ns - 2; i++)  fprintf(fout, " %s", model->splits + i*s1);
      }
   }

   for (i = 0; i < ns; i++) com.spname[i] = spnames0[i];
   free(space);
   return(0);
}

int compareSmodel(const struct SMODEL* model1, const struct SMODEL* model2)
{
   int i, j, s1 = stree.nspecies + 1;

   j = model1->nspecies - model2->nspecies;
   if (j) return(j);
   for (i = 0; i < model1->nspecies; i++) {
      j = strcmp(model1->spnames[i], model2->spnames[i]);
      if (j) return(j);
   }
   for (i = 0; i < model1->nspecies - 2; i++) {
      j = strcmp(model1->splits + i*s1, model2->splits + i*s1);
      if (j) return(j);
   }
   return(0);
}

int compareDelimit(const struct SDELIMIT* delimit1, const struct SDELIMIT* delimit2)
{
   int i, j;

   j = delimit1->nspecies - delimit2->nspecies;
   if (j) return(j);
   for (i = 0; i < delimit1->nspecies; i++) {
      j = strcmp(delimit1->spnames[i], delimit2->spnames[i]);
      if (j) return(j);
   }
   return(0);
}


void SummarizeA11_SpeciesTreeDelimitation(FILE* fout, char mcmcf[])
{
   /* This reads the MCMC sample of species delimitation-tree models from mcmcf and
      generate the posterior distribution of the models and of the species.
      Populations are ordered as in control file bpp.ctl.
      Visited models are collected in the list models.
      Delimited species are collected in the list dspecies.
   */
   int i, j, k = 0, im, ntree, sizemodel, found;
   int s = stree.nspecies, s1 = s + 1, *index, *indexspace, lline = 1024;
   char line[1024], *dspecies = NULL, *p;
   int maxnmodel = s, nmodel = 0, maxndspecies = s, ndspecies = 0, maxndelimit = s, ndelimit = 0, nchangends = 0;
   double *countmodel = NULL, *countdspecies = NULL, *countdelimit = NULL, y, cdf, Pspecies[NSPECIES] = { 0 };
   struct SMODEL *models = NULL, model;
   struct SDELIMIT *delimits = NULL, delimit;
   FILE *ft = gfopen(mcmcf, "r");
   int debug = 0;

   printf("\nSummarizing the species-tree sample in file %s\n", mcmcf);
   fprintf(fout, "\nSummarizing the species-tree sample in file %s\n", mcmcf);
   copySptree();  /* this sets com.ns, com.spname[], etc. */
   if (s < 2) error2("need >=2 pops for this to work.");
   if ((nodes = (struct TREEN*)malloc((s * 2 - 1) * sizeof(struct TREEN))) == NULL) error2("oom");
   for (i = 0; i < s * 2 - 1; i++) nodes[i].annotation = NULL;
   sizemodel = sizeof(struct SMODEL);


   /* (A) posterior of models */
   for (ntree = 0; ; ntree++) {
      if (nmodel + s >= maxnmodel) {
         maxnmodel = (int)(maxnmodel*(ntree < 1000 ? 2 : 1.2));
         models = (struct SMODEL*)realloc(models, maxnmodel*sizemodel);
         countmodel = (double*)realloc(countmodel, maxnmodel * sizeof(double));
         if (models == NULL || countmodel == NULL) error2("oom models || countmodel");
         memset(models + nmodel, 0, (maxnmodel - nmodel)*sizemodel);
         memset(countmodel + nmodel, 0, (maxnmodel - nmodel) * sizeof(double));
      }

      if (ReadTreeN(ft, &i, 0, 0)) break;
      fscanf(ft, "%d", &k);  /* k is number of delimited species */
      for (j = 0; j < s; j++)  Pspecies[j] = (Pspecies[j] * ntree + ((j == k - 1))) / (ntree + 1);
      fgets(line, lline, ft);
      if (debug || (ntree + 1) % 5000 == 0) {
         printf("\rread tree %5d  ", ntree + 1);
         if (s < 15) OutTreeN(F0, 1, 0);
      }
      nchangends += SpeciesTreeDelimitationModelRepresentation(&model, k);
      j = binarysearch(&model, models, nmodel, sizemodel, (int(*)(const void *, const void *))compareSmodel, &found);

      if (found)
         countmodel[j]++;
      else {
         if (j < nmodel) {
            memmove(models + j + 1, models + j, (nmodel - j)*sizemodel);
            memmove(countmodel + j + 1, countmodel + j, (nmodel - j) * sizeof(double));
         }
         memmove(models + j, &model, sizemodel);
         nmodel++;
         countmodel[j] = 1;
      }
      if (debug) {
         printf("\n%4d models:\n", nmodel);
         for (i = 0; i < nmodel; i++) {
            PrintSmodel(F0, models + i, 1);
            printf(" (%.0f)\n", countmodel[i]);
         }
      }
   }  /* for(ntree) */
   if (nchangends) printf("\n%2d trees have # species reduced due to near-0 blengths.\n", nchangends);

   /* index+indexspace is used to sort models, dspecies, etc.  Size here is estimate. */
   if ((index = (int*)malloc(max2(nmodel, 2 * s) * 2 * sizeof(int))) == NULL)
      error2("oom index");
   indexspace = index + nmodel;
   memset(index, 0, nmodel * sizeof(int));
   indexing(countmodel, nmodel, index, 1, indexspace);
   printf("\n(A) List of best models (count postP #species SpeciesTree)\n");
   for (k = 0, cdf = 0; k < nmodel; k++) {
      j = index[k];  y = countmodel[j];
      printf("%6.0f %8.6f %8.6f  ", y, y / ntree, (cdf += y / ntree));
      PrintSmodel(F0, models + j, 1);
      FPN(F0);
      if (cdf > 0.99 && y / ntree < 0.001) break;
      if (nmodel > 99 && k >= 50) {
         printf("\nMany models.  See output file for other models.\n");  break;
      }
   }
   fprintf(fout, "\n(A) List of best models (count postP #species SpeciesTree)\n\n");
   for (k = 0, cdf = 0; k < nmodel; k++) {
      j = index[k];  y = countmodel[j];
      fprintf(fout, "%6.0f %8.6f %8.6f  ", y, y / ntree, (cdf += y / ntree));
      PrintSmodel(fout, models + j, 1);
      FPN(fout);
      if (cdf > 0.999 && y / ntree < 0.0001) break;
   }

   /* (B) posterior of delimitations */
   for (im = 0; im < nmodel; im++) {
      if (ndelimit + s >= maxndelimit) {
         maxndelimit = (int)(maxndelimit*(im < 1000 ? 2 : 1.2));
         delimits = (struct SDELIMIT*)realloc(delimits, maxndelimit * sizeof(struct SDELIMIT));
         countdelimit = (double*)realloc(countdelimit, maxndelimit * sizeof(double));
         if (delimits == NULL || countdelimit == NULL) error2("oom delimit || countdelimit");
         memset(delimits + ndelimit, 0, (maxndelimit - ndelimit) * sizeof(struct SDELIMIT));
         memset(countdelimit + ndelimit, 0, (maxndelimit - ndelimit) * sizeof(double));
      }
      delimit.nspecies = models[im].nspecies;
      memcpy(delimit.spnames[0], models[im].spnames, s*(NSPECIES + 1) * sizeof(char));
      j = binarysearch(&delimit, delimits, ndelimit, sizeof(struct SDELIMIT), (int(*)(const void *, const void *))compareDelimit, &found);
      if (found)
         countdelimit[j] += countmodel[im];
      else {
         if (j < ndelimit) {
            memmove(delimits + j + 1, delimits + j, (ndelimit - j) * sizeof(struct SDELIMIT));
            memmove(countdelimit + j + 1, countdelimit + j, (ndelimit - j) * sizeof(double));
         }
         memmove(delimits + j, &delimit, sizeof(struct SDELIMIT));
         ndelimit++;
         countdelimit[j] = countmodel[im];
      }
   }
   indexspace = index + ndelimit;
   memset(index, 0, ndelimit * sizeof(int));
   indexing(countdelimit, ndelimit, index, 1, indexspace);
   printf("\n(B) %2d species delimitations & their posterior probabilities\n", ndelimit);
   for (k = 0; k < ndelimit; k++) {
      j = index[k];
      printf("%7.0f %9.6f  ", countdelimit[j], countdelimit[j] / ntree);
      PrintSmodel(F0, (struct SMODEL *)(delimits + j), 0);   /* is this casting safe */
      FPN(F0);
      if (ndelimit > 99 && k >= 50) {
         printf("\nMany delimitations.  See output file for others.\n");  break;
      }

   }
   fprintf(fout, "\n(B) %2d species delimitations & their posterior probabilities\n\n", ndelimit);
   for (k = 0; k < ndelimit; k++) {
      j = index[k];
      fprintf(fout, "%7.0f %9.6f  ", countdelimit[j], countdelimit[j] / ntree);
      PrintSmodel(fout, (struct SMODEL *)(delimits + j), 0);   /* is this casting safe */
      FPN(fout);
   }


   /* (C) posterior of delimited species */
   for (im = 0; im < nmodel; im++) {
      if (ndspecies + s >= maxndspecies) {
         maxndspecies = (int)(maxndspecies*(im < 1000 ? 2 : 1.2));
         dspecies = (char*)realloc(dspecies, maxndspecies*s1 * sizeof(char));
         countdspecies = (double*)realloc(countdspecies, maxndspecies * sizeof(double));
         if (dspecies == NULL || countdspecies == NULL) error2("oom dspecies || countdspecies");
         memset(dspecies + ndspecies*s1, 0, (maxndspecies - ndspecies)*s1 * sizeof(char));
         memset(countdspecies + ndspecies, 0, (maxndspecies - ndspecies) * sizeof(double));
      }
      for (i = 0; i < models[im].nspecies; i++) {
         j = binarysearch(models[im].spnames[i], dspecies, ndspecies, s1, (int(*)(const void *, const void *))strcmp, &found);
         if (found)
            countdspecies[j] += countmodel[im];
         else {
            if (j < ndspecies) {
               memmove(dspecies + (j + 1)*s1, dspecies + j*s1, (ndspecies - j)*s1);
               memmove(countdspecies + j + 1, countdspecies + j, (ndspecies - j) * sizeof(double));
            }
            memmove(dspecies + j*s1, models[im].spnames[i], s1);
            ndspecies++;
            countdspecies[j] = countmodel[im];
         }
      }
   }
   indexspace = index + ndspecies;
   memset(index, 0, ndspecies * sizeof(int));
   indexing(countdspecies, ndspecies, index, 1, indexspace);
   printf("\n(C) %2d delimited species & their posterior probabilities\n", ndspecies);
   for (k = 0; k < ndspecies; k++) {
      j = index[k];
      printf("%7.0f %9.6f  ", countdspecies[j], countdspecies[j] / ntree);
      p = dspecies + j*s1;
      while (*p) printf("%s", stree.nodes[*p++ - 1].name);
      FPN(F0);
      if (ndspecies > 99 && k >= 99) {
         printf("\nMany delimited species.  See output file for others.\n");  break;
      }
   }
   fprintf(fout, "\n(C) %2d delimited species & their posterior probabilities\n\n", ndspecies);
   for (k = 0; k < ndspecies; k++) {
      j = index[k];
      fprintf(fout, "%7.0f %9.6f  ", countdspecies[j], countdspecies[j] / ntree);
      p = dspecies + j*s1;
      while (*p) fprintf(fout, "%s", stree.nodes[*p++ - 1].name);
      FPN(fout);
   }


   printf("\n(D) Posterior probability for # of species\n");
   if (stree.analysis == A11 && stree.SpeciesModelPrior >= 2)
      for (j = 0; j < s; j++) stree.PriorSA11[j] = 1.0 / s;
   for (j = 0; j < s; j++) if (s <= 60 || Pspecies[j] > .0001)
      printf("P[%2d] = %8.6f  prior[%2d] = %8.6f\n", j + 1, Pspecies[j], j + 1, stree.PriorSA11[j]);
   fprintf(fout, "\n(D) Posterior probability for # of species\n\n");
   for (j = 0; j < s; j++) if (s <= 60 || Pspecies[j] > .0001)
      fprintf(fout, "P[%2d] = %8.6f  prior[%2d] = %8.6f\n", j + 1, Pspecies[j], j + 1, stree.PriorSA11[j]);

   free(nodes);  free(models);  free(dspecies);  free(delimits);
   free(countmodel);  free(countdspecies);  free(countdelimit);
   free(index);
   fclose(ft);
   return;
}


int SummarizeA00_DescriptiveStatisticsSimpleBPP(FILE *fout, char infile[], int SkipColumns)
{
   FILE *fin = gfopen(infile, "r");
   int  n, p, i, j;
   char *fmt = " %9.6f", *fmt1 = " %9.1f", timestr[32];
   double *D, *x, *mean, *median, *minx, *maxx, *x005, *x995, *x025, *x975, *xHPD025, *xHPD975, *var;
   double *Tint, tmp[2], *y, rho1;
   char *line;
   static int lline = 1000000, ifields[MAXNFIELDS], HasHeader = 1;
   static char varstr[MAXNFIELDS][96] = { "" };
   int s = stree.nspecies, nskip = SkipColumns + (stree.NoTheta ? 0 : stree.npop);
   FILE *fFigTree;
   char FigTreef[96] = "FigTree.tre";

   printf("\nSummarizing.");
   fprintf(fout, "\n\nSummary of MCMC results\n");
   if ((line = (char*)malloc(lline * sizeof(char))) == NULL) error2("oom ds");
   scanfile(fin, &n, &p, &HasHeader, line, ifields);
   printf("\n%d records, %d variables\n", n, p);

   if (n < 1) exit(0);

   D = (double*)malloc(p*n * sizeof(double));
   mean = (double*)malloc((p * 13 + n) * sizeof(double));
   if (D == NULL || mean == NULL) error2("oom DescriptiveStatistics.");
   memset(D, 0, p*n * sizeof(double));
   memset(mean, 0, (p * 12 + n) * sizeof(double));
   median = mean + p; minx = median + p; maxx = minx + p;
   x005 = maxx + p; x995 = x005 + p; x025 = x995 + p; x975 = x025 + p; xHPD025 = x975 + p; xHPD975 = xHPD025 + p;
   var = xHPD975 + p;   Tint = var + p;  y = Tint + p;

   if (HasHeader)
      for (i = 0; i < p; i++) sscanf(line + ifields[i], "%s", varstr[i]);
   for (i = 0; i < n; i++)
      for (j = 0; j < p; j++)
         fscanf(fin, "%lf", &D[j*n + i]);
   fclose(fin);

   printf("Collecting mean, median, min, max, percentiles, etc.\n");
   for (j = SkipColumns, x = D + j*n; j < p; j++, x += n) {
      memmove(y, x, n * sizeof(double));
      Tint[j] = 1 / Eff_IntegratedCorrelationTime(y, n, &mean[j], &var[j], &rho1); /* y destroyed */
      qsort(x, (size_t)n, sizeof(double), comparedouble);
      minx[j] = x[0];  maxx[j] = x[n - 1];
      median[j] = (n % 2 == 0 ? (x[n / 2 - 1] + x[n / 2]) / 2 : x[n / 2]);
      x005[j] = x[(int)(n*.005)];    x995[j] = x[(int)(n*.995)];
      x025[j] = x[(int)(n*.025)];    x975[j] = x[(int)(n*.975)];

      HPDinterval(x, n, tmp, 0.05);
      xHPD025[j] = tmp[0];
      xHPD975[j] = tmp[1];
      if ((j + 1) % 2 == 0 || j == p - 1)
         printf("\r\t\t\t%6d/%6d done  %s", j + 1, p, printtime(timestr));
   }

   fprintf(fout, "\n\n       ");
   for (j = SkipColumns; j < p; j++) fprintf(fout, "   %s", varstr[j]);
   fprintf(fout, "\nmean    ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, mean[j]);
   fprintf(fout, "\nmedian  ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, median[j]);
   fprintf(fout, "\nS.D.    ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, sqrt(var[j]));
   fprintf(fout, "\nmin     ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, minx[j]);
   fprintf(fout, "\nmax     ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, maxx[j]);
   fprintf(fout, "\n2.5%%    "); for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, x025[j]);
   fprintf(fout, "\n97.5%%   "); for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, x975[j]);
   fprintf(fout, "\n2.5%%HPD "); for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, xHPD025[j]);
   fprintf(fout, "\n97.5%%HPD"); for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, xHPD975[j]);
   fprintf(fout, "\nESS*    ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt1, n / Tint[j]);
   fprintf(fout, "\nEff*    ");  for (j = SkipColumns; j < p; j++) fprintf(fout, fmt, 1 / Tint[j]);
   FPN(fout);
   fflush(fout);

   /* print results in FigTree.tre */
   com.np = GetInitials();               /* this sets up stree. */
   if (stree.NoTheta)
      for (j = 0; j < stree.nnode; j++) stree.nodes[j].theta = -1;
   collectx(2, NULL, mean + SkipColumns);  /* this moves means of theta and tau into stree */
   nodes = nodes_t;
   copySptree();                         /* copy stree into nodes[] */
   for (j = s; j < stree.nnode; j++) {       /* copy 95%HPD CI into nodes[].annotation for printing */
      nodes[j].annotation = (char*)malloc(64 * sizeof(char));
      if (nodes[j].label > 0)
         sprintf(nodes[j].annotation, "[&height_95%%_HPD={%.6g, %.6g}, theta=%.6g]\0", xHPD025[nskip + j - s], xHPD975[nskip + j - s], nodes[j].label);
      else
         sprintf(nodes[j].annotation, "[&height_95%%_HPD={%.6g, %.6g}]\0", xHPD025[nskip + j - s], xHPD975[nskip + j - s]);
   }
   if (stree.nspecies > 1) {
      if ((fFigTree = (FILE*)fopen(FigTreef, "w")) == NULL)
         error2("FigTree.tre file creation error");
      fprintf(fFigTree, "#NEXUS\nBEGIN TREES;\n\n\tUTREE 1 = ");
      OutTreeN(fFigTree, 1, PrBranch | PrNodeStr);
      fprintf(fFigTree, "\n\nEND;\n");
      fprintf(fFigTree, "\n\n[Species tree with tau as branch lengths and theta as labels, for FigTree.\n");
      fprintf(fFigTree, "In FigTree, choose 95%%HPD for Node Bars and label for Node Labels]\n");
      fclose(fFigTree);
      printf("\nFigTree tree is in FigTree.tre");
   }
   for (j = s; j < stree.nnode; j++)
      free(nodes[j].annotation);
   free(D); free(mean); free(line);
   return(0);
}


int MCMC(FILE* fout)
{
   FILE *fmcmc = NULL, **fGtree;
   char line[16000];
   int nsteps = 5 + (data.est_locusrate == 1 || data.est_heredity == 1) + (data.nseqerr > 0);
   int locus, j, k, ir, Bmodel = 0, lline = 16000, s = stree.nspecies;
   double *x = NULL, *mx = NULL, lnL, lnpG, mlnL = 0, PrSplit = 0.5;
   double PjumpRJ = 0, PjumpSlider = 0, Pjump[7] = { 0 }, nround = 0, nroundRJ = 0, nroundSPR = 0, Dpi = 1;
   double mrootage = 0, mroottheta = 0, Pspecies[NSPECIES] = { 0 };
   int ndspecies = 0, ndspeciesbest = 0;

   if (mcmc.print > 0)
      fmcmc = gfopen(com.mcmcf, "w");

   if (stree.analysis == A10) {
      InitializeDelimitationModel(1, stree.pmodel);  /* pmodel has prior */
      zero(stree.pmodel, stree.nModels);
   }
   mcmc.moveinnode = 1;   /* moves internal nodes in the gene tree as well */
   mcmc.saveconP = 1;
   if (!mcmc.usedata) mcmc.saveconP = 0;

   printf("\nMCMC settings: %d burnin, sampling every %d, %d samples\n", mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   fprintf(fout, "\n\nMCMC settings: %d burnin, sampling every %d, %d samples\n", mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   if (mcmc.usedata) puts("Approximating posterior, using sequence data");
   else             puts("Approximating prior, not using sequence data");

   printf("(Settings: cleandata=%d print=%d saveconP=%d moveinnode=%d)\n",
      com.cleandata, mcmc.print, mcmc.saveconP, mcmc.moveinnode);

   if (stree.analysis == A00)  printf("\nStarting MCMC...\n");
   else if (stree.analysis == A10)  printf("\nStarting rjMCMC...\n");
   else if (stree.analysis == A01)          /* screen output incorrect if NNI_ONLY */
      printf("\nStarting MCMC-SPR(%.0f%%)/Slider(%.0f%%)...\n", (1 - mcmc.pSlider) * 100, mcmc.pSlider * 100);
   else if (stree.analysis == A11)            /* screen output incorrect if NNI_ONLY */
      printf("\nStarting rjMCMC+SPR(%.0f%%)/Slider(%.0f%%)...\n", (1 - mcmc.pSlider) * 100, mcmc.pSlider * 100);
   if (stree.speciestree == 1)
      printf("NodeSlider: p = %.4f ExpandRatio = %.3f ShrinkRatio = %.3f\n", mcmc.pSlider, mcmc.ExpandRatio, mcmc.ShrinkRatio);

   if (stree.speciesdelimitation == 1) {
      printf("PrSplit = %.6f\nrj algorithm %d: new theta from ", PrSplit, mcmc.RJalgorithm);
      if (mcmc.RJalgorithm == 0)       printf("sliding window with c = %.2f\n", mcmc.RJfinetune[0]);
      else if (mcmc.RJalgorithm == 1)  printf("G(a=%.2f, m=%.2f)\n", mcmc.RJfinetune[0], mcmc.RJfinetune[1]);
   }

   com.np = GetInitials();

   stree.iModel = -1;
   if (stree.analysis == A10)
      stree.iModel = GetDmodelIndex();
   if (stree.speciesdelimitation == 1)
      GetRootTau();

   if (!stree.speciestree) {
      k = max2(4, s * 3 + 3);
      if (mcmc.printlocusrate == 1)  k += data.ngene;  /* print locus rates */
      if (mcmc.printheredity == 1)   k += data.ngene;  /* print heredity scalars */
      if (mcmc.printGenetree == 1 && s == 1)   k += data.ngene;  /* print tMRCA for loci */
      if (data.nseqerr) k += 16;
      x = (double*)malloc(k * 2 * sizeof(double));  /* 2 for x & mx */
      mx = x + k;
      if (x == NULL) error2("oom");
      zero(x, k * 2);
   }

   if (mcmc.printGenetree && s > 1) {
      fGtree = (FILE**)malloc(data.ngene * sizeof(FILE*));
      if (fGtree == NULL) error2("fGtree mem allocation error.");
      system("mkdir gtree");
      for (locus = 0; locus < data.ngene; locus++) {
         sprintf(line, "gtree/L%03d.tre", locus + 1);
         fGtree[locus] = (FILE*)fopen(line, "w");
         if (fGtree[locus] == NULL) {
            printf("gtree file open error for locus %d\n", locus + 1);
            exit(-1);
         }
      }
   }

   if (com.checkpoint == 2) {
      ReadMCMCstate(com.checkpointf);
      mcmc.burnin = 0;
      printf("\nInitial parameters, np = %d.\nSptree & Genetrees read from %s.\n", com.np, com.checkpointf);
   }
   else {
      printf("\nInitial parameters, np = %d.\nGenetrees generated from the MSC density.\n", com.np);
      /* initialize Gtrees and likelihood, for each locus */
      for (locus = 0; locus < data.ngene; locus++) {
         GetRandomGtree(locus);
         data.root[locus] = tree.root;
      }
      if (stree.NoTheta)  /* theta values were used to generate Gtrees & are now deleted. */
         for (j = 0; j < s * 2 - 1; j++) stree.nodes[j].theta = -1;
   }
   for (locus = 0, lnpG = 0; locus < data.ngene; locus++) {
      UseLocus(locus, 0, mcmc.usedata, 0);
      CountCoalescentLocus(locus);
      if (!stree.NoTheta) lnpG += (data.lnpGi[locus] = lnpG_S(locus));
   }
   if (stree.NoTheta)
      lnpG = lnpG_S_NoTheta();

   if (stree.analysis == A01) {  /* A01: speciestree.  This is the starting tree. */
      printSptreeBPP(fmcmc);
      fprintf(fmcmc, "\n");
      fprintf(fout, "\nStarting species tree:\n");
      printSptreeBPP(fout);
      fprintf(fout, "\n");
   }
   if (stree.speciestree == 0)  collectx(0, fmcmc, x);
   if (stree.speciesdelimitation == 1 || stree.speciestree == 1)
      checkGtree();

   if (stree.speciestree == 0) {
      k = com.np;
      if ((mcmc.printGenetree == 1 && s == 1)) /* t_MRCA is printed in this case */
         k -= data.ngene;
      if (mcmc.printlocusrate == 1) k -= data.ngene;
      if (mcmc.printheredity == 1)  k -= data.ngene;
      for (j = 0; j < k; j++) printf("%9.5f", x[j]);
   }
   for (j = 0; j < data.maxns * 2 - 1; j++) com.oldconP[j] = 0;
   lnL = lnpData(data.lnpDi);
   printf("\nlnpG0 = %9.6f  lnL0 = %6.6f\n", lnpG, lnL);
   fprintf(fout, "\nlnpG0 = %9.6f  lnL0 =%9.6f\n", lnpG, lnL);
   fflush(fout);

   for (ir = -mcmc.burnin; ir < mcmc.sampfreq*mcmc.nsample; ir++) {
      if (debug) printf("\n\n***** MCMC ir = %6d", ir + 1);
      if (ir == 0 || (mcmc.resetFinetune && nround >= 100 && mcmc.burnin >= 200
         && ir < 0 && ir % (mcmc.burnin / 4) == 0)) {
         /* reset finetune parameters.  Do this four times. */
         if (mcmc.resetFinetune && mcmc.burnin >= 200) {
            ResetFinetuneSteps(fout, Pjump, mcmc.finetune, nsteps);
         }
         nround = nroundRJ = nroundSPR = 0;  PjumpRJ = PjumpSlider = 0;
         mlnL = 0;
         zero(Pjump, nsteps);
         zero(Pspecies, s);
         if (!stree.speciestree)
            zero(mx, com.np);
         if (stree.speciesdelimitation == 1 && stree.speciestree == 0)
            zero(stree.pmodel, stree.nModels);
         mrootage = 0;  mroottheta = 0;
      }
      nround++;
      if (stree.speciesdelimitation == 1) {
         if (rndu() < PrSplit)
            j = UpdateSpeciesSplit(&lnL, &lnpG, com.space, PrSplit);
         else
            j = UpdateSpeciesJoin(&lnL, &lnpG, com.space, PrSplit);

         if (j != 2) {
            nroundRJ++;
            PjumpRJ += j;
         }
      }
      for (j = s, ndspecies = 1; j < s * 2 - 1; j++) if (stree.nodes[j].age > 0)  ndspecies++;
      if (ndspecies > 2 && stree.speciestree == 1) {
#if(defined NNI_ONLY)         /*  NNI only */
         j = UpdateSpeciesTreeNNI(&lnL, &lnpG, com.space);
#else
         if (rndu() > mcmc.pSlider) j = UpdateSpeciesTreeSPR(0, &lnL, &lnpG, com.space);  /* 1 for NNIonly */
         else                       j = UpdateSpeciesTreeNodeSlider(&lnL, &lnpG, 5 * data.tau_prior[0] / data.tau_prior[1], com.space);
#endif
         if (j == 0 || j == 1)  nroundSPR++;
         if (j == 1)          PjumpSlider++;
      }
      Pjump[0] = (Pjump[0] * (nround - 1) + UpdateGB_InternalNode(&lnL, &lnpG, mcmc.finetune[0])) / nround;
      Pjump[1] = (Pjump[1] * (nround - 1) + UpdateGB_SPR(&lnL, &lnpG, mcmc.finetune[1])) / nround;
      if (!stree.NoTheta) Pjump[2] = (Pjump[2] * (nround - 1) + UpdateTheta(&lnpG, mcmc.finetune[2], com.space)) / nround;
      if (s > 1 && stree.nodes[stree.root].age > 0)
         Pjump[3] = (Pjump[3] * (nround - 1) + UpdateTau(&lnL, &lnpG, mcmc.finetune[3], com.space)) / nround;

      Pjump[4] = (Pjump[4] * (nround - 1) + mixing(&lnL, &lnpG, mcmc.finetune[4], com.space)) / nround;

      if (data.est_locusrate == 1 || data.est_heredity == 1)
         Pjump[5] = (Pjump[5] * (nround - 1) + UpdateLocusrateHeredity(&lnL, &lnpG, mcmc.finetune[5])) / nround;
      if (data.nseqerr)
         Pjump[6] = (Pjump[6] * (nround - 1) + UpdateSequenceErrors(&lnL, mcmc.finetune[6], com.space)) / nround;

      if (stree.speciestree == 0) {
         collectx(1, NULL, x);
         for (j = 0; j < (stree.speciesdelimitation == 1 && stree.speciestree == 0 ? 1 : com.np); j++)
            mx[j] = (mx[j] * (nround - 1) + x[j]) / nround;
      }
      mrootage = (mrootage*(nround - 1) + stree.nodes[stree.root].age) / nround;
      mroottheta = (mroottheta*(nround - 1) + stree.nodes[stree.root].theta) / nround;
      mlnL = (mlnL*(nround - 1) + lnL / BFbeta) / nround;

      if (stree.speciesdelimitation == 1) {
         if (stree.speciestree == 0)  stree.pmodel[stree.iModel] ++;
         Pspecies[ndspecies - 1] = (Pspecies[ndspecies - 1] * (nround - 1) + 1) / nround;
         for (j = 0; j < s; j++)
            if (j != ndspecies - 1) Pspecies[j] = (Pspecies[j] * (nround - 1) + 0) / nround;
         for (j = 1, ndspeciesbest = 0; j < s; j++)
            if (Pspecies[j] > Pspecies[ndspeciesbest]) ndspeciesbest = j;
      }

      if (mcmc.print && ir >= 0 && (ir + 1) % mcmc.sampfreq == 0) {     /* mcmc sampling */
         if (stree.speciestree == 0) fprintf(fmcmc, "%d", ir + 1);

         if (stree.analysis == A00) {       /* A00 */
            for (j = 0; j < com.np; j++)    fprintf(fmcmc, "\t%.5g", x[j]);
            if (mcmc.usedata) fprintf(fmcmc, "\t%.3f", lnL / BFbeta);
            fprintf(fmcmc, "\n");
         }
         else if (stree.analysis == A10) {  /* A10: delimitation */
            fprintf(fmcmc, "\t%d\t%s", com.np, printDelimitationModel());
            for (j = 0; j < com.np; j++)    fprintf(fmcmc, "\t%.5g", x[j]);
            if (mcmc.usedata) fprintf(fmcmc, "\t%.3f", lnL / BFbeta);
            fprintf(fmcmc, "\n");
         }
         else if (stree.analysis == A01) {  /* A01: speciestree */
            printSptreeBPP(fmcmc);
            fprintf(fmcmc, "\n");
         }
         else if (stree.analysis == A11) {  /* A11: delimitation & speciestree */
            printSptreeBPP(fmcmc);
            fprintf(fmcmc, " %d \n", ndspecies);
         }

         /* print the gene tree */
         if (mcmc.printGenetree && s > 1)
            for (locus = 0; locus < data.ngene; locus++) {
               UseLocus(locus, 1, 0, 1);
               OutTreeN(fGtree[locus], 1, 1);
               fprintf(fGtree[locus], "\n");
               /* fprintf(fGtree[locus], " [%9.5f]\n", nodes[tree.root].age); */
            }
      }   /* if() for sampling */
      k = mcmc.sampfreq*mcmc.nsample;
      if (noisy && (k <= 500 || (ir + 1) % (k / 200) == 0)) {  /* printing on screen */
         printf("\r%3.0f%%", (ir + 1.499) / (mcmc.nsample*mcmc.sampfreq)*100.);
         for (j = 0; j < nsteps; j++) printf(" %4.2f", Pjump[j]);
         printf(" ");

         fprintf(fout, "\n%3.0f%%", (ir + 1.499) / (mcmc.nsample*mcmc.sampfreq)*100.);
         for (j = 0; j < nsteps; j++) fprintf(fout, " %5.3f", Pjump[j]);
         fprintf(fout, " ");

         if (stree.analysis == A00) {          /* A00: print first few theta's, first 3 tau's */
            if (!stree.NoTheta) {
               for (j = 0; j < min2(stree.npop, 3); j++)  printf(" %6.4f", mx[j]);
               printf(" ");
            }
            k = (stree.NoTheta ? 0 : stree.npop);
            for (j = k; j < min2(com.np, k + 3); j++) printf(" %6.4f", mx[j]);

            if (!stree.NoTheta) {
               for (j = 0; j < min2(stree.npop, 3); j++)  fprintf(fout, " %6.4f", mx[j]);
               fprintf(fout, " ");
            }
            k = (stree.NoTheta ? 0 : stree.npop);
            for (j = k; j < min2(com.np, k + 3); j++) fprintf(fout, " %6.4f", mx[j]);
         }
         else if (stree.analysis == A10) {     /* A10 */
            for (j = 1, Bmodel = 0; j < stree.nModels; j++)
               if (stree.pmodel[j] > stree.pmodel[Bmodel]) Bmodel = j;
            printf(" %2d %6.4f %s", com.np, (nroundRJ ? PjumpRJ / nroundRJ : 0), printDelimitationModel());
            printf(" P[%d]=%6.4f  %6.4f %6.4f", Bmodel + 1, stree.pmodel[Bmodel] / nround, mroottheta, mrootage);

            fprintf(fout, " %2d %6.4f %s", com.np, (nroundRJ ? PjumpRJ / nroundRJ : 0), printDelimitationModel());
            fprintf(fout, " P[%d]=%6.4f  %6.4f %6.4f", Bmodel + 1, stree.pmodel[Bmodel] / nround, mroottheta, mrootage);
         }
         else if (stree.analysis == A01) {     /* A01 */
            printf(" %6.4f  %6.4f %6.4f", (nroundSPR ? PjumpSlider / nroundSPR : 0), mroottheta, mrootage);
            fprintf(fout, " %6.4f  %6.4f %6.4f", (nroundSPR ? PjumpSlider / nroundSPR : 0), mroottheta, mrootage);
         }
         else if (stree.analysis == A11) {     /* A11 */
            printf(" %2d %2d %6.4f %6.4f P(%d)=%6.4f  %6.4f %6.4f", ndspecies, com.np,
               (nroundRJ ? PjumpRJ / nroundRJ : 0), (nroundSPR ? PjumpSlider / nroundSPR : 0),
               ndspeciesbest + 1, Pspecies[ndspeciesbest], mroottheta, mrootage);
            fprintf(fout, " %2d %2d %6.4f %6.4f P(%d)=%6.4f  %6.4f %6.4f", ndspecies, com.np,
               (nroundRJ ? PjumpRJ / nroundRJ : 0), (nroundSPR ? PjumpSlider / nroundSPR : 0),
               ndspeciesbest + 1, Pspecies[ndspeciesbest], mroottheta, mrootage);
         }

         printf(" %7.2f", lnpG);
         fprintf(fout, " %7.2f", lnpG);
         if (mcmc.usedata) {
            printf(" %8.4f", mlnL);
            fprintf(fout, " %8.4f", mlnL);
         }
         if (mcmc.sampfreq*mcmc.nsample >= 50 && (ir + 1) % (mcmc.sampfreq*mcmc.nsample / 20) == 0) {
            testlnL = 1;
            if (fabs(lnL - lnpData(data.lnpDi)) > 0.01) {
               printf("lnL incorrect: %12.6f != %12.6f (ir=%5d)", lnL, lnpData(data.lnpDi), ir + 1);
               exit(0);
            }
            testlnL = 0;
            if (mcmc.print) fflush(fmcmc);
            printf(" %s\n", printtime(timestr));
            fprintf(fout, " %s\n", printtime(timestr));
            fflush(fout);

            if (ir >= 0 && com.checkpoint == 1) { /* save MCMC state */
               SaveMCMCstate(com.checkpointf, lnpG, lnL);
            }

         }
      }
   }  /* for(ir) */

   if (mcmc.print) fclose(fmcmc);
   if (mcmc.printGenetree && s > 1) {
      for (locus = 0; locus < data.ngene; locus++)
         fclose(fGtree[locus]);
      free(fGtree);
   }

   if (BFbeta != 1) printf("\nBFbeta = %8.6f  E_b(lnf(X)) = %9.4f\n", BFbeta, mlnL);
   printf("\nTime used: %s", printtime(timestr));

   return(0);
}

#endif
