// This program samples nucleoid configurations according to a 
// coarse-grained model, version 2
// Michael Feig 2012-2013, Asli Yildirim 2016

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#include "vector.h"
#include "pdb.h"

#define SMALL_NUM 0.000001

#define MAX_THREADS 1024
int num_threads=4; //MJ was 8

int check=0;

int mbps=816394;                      // number of basepairs
int ndomain=400;                       // number of domains: 50-400

double rfac=0.879;
double cgbeadsize=5.0; double rise=0.340;
int cgbpbead=(int)(cgbeadsize/rise+0.5);
int maxcgbp=int(mbps/cgbpbead)*1.1;

double bpssingle=1.0/rise;      // bps per nm
double bpsscoil=2.0/rise/0.879;   // bps per nm

double setringlen=20;       // nm
double rdistk=1.0;    // force constant for deviation from set distance
double setchainlen=50;    // nm
double bdistk=0.05;    // force constant for deviation from set distance

double halfchainlen;

int cgbeadperchain;
int cgbeadperring;
int beadforhalfchain;

double branchradius=6.5;   // nm
double ringradius=1;       // nm

double initkT=0.59*500;  // initial temperature
double tfactor=0.99999;      // reduction of temperature at each step, set to 1 for constant T

double restforcefac=0.001;   // multiplier for restraint forces

double rgyr0=0;            // target radius of gyration
double rgyrk=0.00;        // force constant for radius of gyration restraint

double maxr=2000;          // maximum size (spherical radius)
double maxrk=0.0001;       // force constant for maximum size restraint

double minangle0=90.0;     // minimum angle for connected cylinders
double maxangle0=180.0;    // minimum angle for connected cylinders
double minangle1=0.0;      // minimum angle for connected cylinders
double maxangle1=90.0;     // minimum angle for connected cylinders
double anglek=0.05;         // force constant for angle restraint

double contactmin=2.0;     // minimum contact distance
double contactk=1.0;       // force constant for contact restraint

int    savefreq=25000;      // how often should we save the PDB?
int    eneoutfreq=2500;      // how often should we save ebnergy

double cutofflen;

int    regenlookupfreq=100; // how often should we regenerate entire lookup table

double deltadisp;

double beadradius=5.0;

int basemapoffset=100000;
int natom=0;

class Bond;

class Bead {
  public:
    Vector pos;
    int type;  // 1: ring, 2: branch, 3: long branch
    int index;

    double radius;
   
    double rgyrener;
    double rgyr;

    Bead(Vector p, int t, int i, double r) :
      pos(p), type(t), index(i), radius(r), rgyr(0.0), rgyrener(0.0) {}

    Bead& operator=(const Bead& from) {
      pos=from.pos;
      type=from.type;
      index=from.index;
      radius=from.radius;
      rgyr=from.rgyr;
      rgyrener=from.rgyrener;
      return (*this); 
    }
 
    void dumpInfo() {
      fprintf(stderr,"bead %d is of type %d\n",index,type);
    }

    void calcRgyrEner() {
      double bnorm=pos.norm();
      rgyr=bnorm*bnorm;
      if (bnorm>maxr) {
        double dmax=(bnorm-maxr);
        rgyrener=dmax*dmax*maxrk;
      } else {
        rgyrener=0.0;
      }
    } 
};

class Connection {
  public:
    int inx;
    int jnx;
}; 

class Bond {
  public:
    Bead *bead1;
    Bead *bead2;

    int index;

    int domain;

    int nlink;
    Bond *next[3];
    Bond *prev;

    double radius;
    double bpsnm;
    double lennm;
    double fconst;

    Vector dvec;
    double d;

    int done;

    Bond() : bead1(0), bead2(0), radius(0), bpsnm(0), lennm(0), nlink(0),fconst(0.0),index(0) {}

    Bond(Bead *b1, Bead *b2, double r, double bnm, double lnm, double fc, int inx) :
      bead1(b1), bead2(b2), radius(r), bpsnm(bnm), lennm(lnm), nlink(0), fconst(fc), index(inx) {}

    ~Bond() {}

    Bond& operator=(const Bond& from) {
      bead1=from.bead1;
      bead2=from.bead2;

      radius=from.radius;
      bpsnm=from.bpsnm;
      lennm=from.lennm;
      fconst=from.fconst;

      dvec=from.dvec;
      d=from.d;

      return (*this);
    } 

    void dumpInfo() {
      setDiff();
      fprintf(stderr,"bond between %d and %d : length = %f\n",bead1->index,bead2->index,d);
      fprintf(stderr,"  bead1 %d pos: %lf %lf %lf\n",bead1->index,bead1->pos.x(),bead1->pos.y(),bead1->pos.z());
      fprintf(stderr,"  bead2 %d pos: %lf %lf %lf\n",bead2->index,bead2->pos.x(),bead2->pos.y(),bead2->pos.z());
    }

    void setDiff() {
      dvec=bead2->pos-bead1->pos;
      d=dvec.norm();
      dvec/=d;
    }

    double length() {
      return d;
    }

    double selfEnergy() {
      setDiff();
      return (d-lennm)*(d-lennm)*fconst;
    }

    double angleEnergy(Vector d1, Vector d2, double minangle, double maxangle) {
      double d1n=d1.norm();
      double d2n=d2.norm();
      if (d1n<SMALL_NUM || d2n<SMALL_NUM) return 0.0;

      double cosa=d1*d2/d1n/d2n; 
      if (abs(cosa)  >1) {
      fprintf(stderr,"cosa is out of 1!\n");
      }  
      double a=acos(cosa)*180.0/M_PI;
      double da=0.0;
      if (a<minangle && a>=0.0) {                
        da=(a-minangle);
      } else if (a>-minangle && a<0.0) {
        da=(a+minangle);
      } else if (a>maxangle && a>=0.0) {
        da=(a-maxangle);
      } else if (a<-maxangle && a<0.0) {
        da=(a+maxangle);
      } else {
        da=0.0;
      }
      return da*da*anglek;
    }

    double minimumDistance(Bond *bx) {

      if (bead1->index == bx->bead1->index ||
          bead1->index == bx->bead2->index ||
          bead2->index == bx->bead1->index ||
          bead2->index == bx->bead2->index) {
        return 0.0;
      }

      Vector u=bead2->pos-bead1->pos;
      Vector v=bx->bead2->pos-bx->bead1->pos;
      double un=u.norm();
      double vn=v.norm();
      
      Vector c1=bead1->pos+bead2->pos;
      Vector c2=bx->bead1->pos+bx->bead2->pos;
      Vector cd=c1-c2;
      if (cd.norm()-un/2.0-vn/2.0>cutofflen) {
         return cutofflen+SMALL_NUM;
      } 

      Vector w=bead1->pos-bx->bead1->pos;

      double a=un*un; //u*u;
      double b=u*v;
      double c=vn*vn; //v*v;
      double d=u*w;
      double e=v*w;
      double D=a*c-b*b;
      double sc,sN,sD=D;
      double tc,tN,tD=D;
      
      if (abs(D) < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on norm S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
      } else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        } else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
      }
    
      if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        if (-d < 0.0)         // recompute sc for this edge
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
      } else if (tN > tD) {       // tc > 1 => the t=1 edge is visible
        tN = tD;
        if ((-d + b) < 0.0)   // recompute sc for this edge
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d + b);
            sD = a;
        }
      }
    
      sc = ((abs(sN) < SMALL_NUM || abs(sD)<SMALL_NUM) ? 0.0 : sN / sD);
      tc = ((abs(tN) < SMALL_NUM || abs(tD)<SMALL_NUM) ? 0.0 : tN / tD);

      Vector su=u*sc;
      Vector sv=v*tc;
      Vector duv=su-sv;
      Vector dP=w+duv;

      double mind=dP.norm();
      if (mind>0) {
         double cosau=(u*dP)/un/mind;
         double cosav=(v*dP)/vn/mind;
         double csu=1.0-cosau*cosau;
         double csv=1.0-cosav*cosav;

         double sinau=(csu>SMALL_NUM)?sqrt(csu):0.0; 
         double sinav=(csv>SMALL_NUM)?sqrt(csv):0.0; 
	 
	 mind-=radius*sinau;
         mind-=radius*sinav;
      }
      
      return mind;
    }     

    double cylinderContactEnergy(Bond *bx) {
      double mind=minimumDistance(bx);
      mind-=contactmin; 
      if (mind<0) {
        return mind*mind*contactk;
      } else {
        return 0.0;
      }  
    }

    double interBondEnergy(Bond *b) {
      if (b->index<=index) return 0.0;
      if (bead1==b->bead1) {
        Vector d1=bead2->pos-bead1->pos;
        Vector d2=b->bead2->pos-b->bead1->pos;
        return angleEnergy(d1,d2,minangle0,maxangle0); 
      } else if (bead1==b->bead2) {
        Vector d1=bead2->pos-bead1->pos;
        Vector d2=b->bead2->pos-b->bead1->pos;
        return angleEnergy(d1,d2,minangle1,maxangle1); 
      } else if (bead2==b->bead1) {
        Vector d1=bead2->pos-bead1->pos;
        Vector d2=b->bead2->pos-b->bead1->pos;
        return angleEnergy(d1,d2,minangle1,maxangle1); 
      } else if (bead2==b->bead2) {
        Vector d1=bead2->pos-bead1->pos;
        Vector d2=b->bead2->pos-b->bead1->pos;
	fprintf(stderr,"should not get here!\n");
        return angleEnergy(d1,d2,0,0); 
      } else {
        return cylinderContactEnergy(b);
      }
    }   
    
    double bps() {
      return length()*bpsnm;
    }

    void remove(Bond *b) {
      int tinx=-1;
      for (int i=0; i<nlink; i++) {
        if (next[i]==b) 
          tinx=i;
      }
      if (tinx<0) {
        fprintf(stderr,"cannot find bond?!\n");
        exit(1);
      }
      if (tinx<--nlink) {
        for (int i=tinx; i<nlink; i++) 
          next[i]=next[i+1];
      }
    }  
    
    void add(Bond *b) {
      if (nlink>=3) {
        fprintf(stderr,"cannot add more bonds\n");
        exit(1);
      }
      next[nlink++]=b; 
      b->domain=domain;
      b->prev=this;
    }
};

class Restraint {
 public:
   long inx;
   long jnx;
   double force;
   double dist;
   int type; // 1: harmonic, 2: onesided
 
   Restraint() : inx(-1), jnx(-1), force(1.0), dist(0.0), type(1) {}
   Restraint(long i, long j, double f=1.0, double d=0.0, int t=1) : inx(i), jnx(j), force(f), dist(d), type(t) {}
}; 

void readRestraints(char *fname, Restraint *&rlist, int &nrest, int maxrest, int type) {
  
  FILE *fptr;
  fptr=fopen(fname,"r");
  if (fptr==0) {
    fprintf(stderr,"Cannot open restraint file %s\n",fname);
    exit(1);
  }

  char line[1024];
  while (!feof(fptr)) {
    if (fscanf(fptr,"%ld%ld%lf%lf\n",
          &rlist[nrest].inx,&rlist[nrest].jnx,&rlist[nrest].force,&rlist[nrest].dist)) {
       rlist[nrest].force*=restforcefac;
       rlist[nrest].type=type;
       if (++nrest>=maxrest) {
         fprintf(stderr,"exceeding maximum number of restraints (%d)\n",maxrest);
         exit(1);
       }
     } 
  }
  fclose(fptr);
}

double getGeometry(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil,       
                   double *rgyr_ring, double *rgyr_scoil, double *rgyr_out, 
                   double *averagering_out, double *averagescoil_out) {  
 
  int i;

  double rgyr_nring=0.0;
  for (i=0; i<nring; i++) {
      double bnorm=rbead[i]->pos.norm();
      rgyr_nring+=bnorm*bnorm;
  }

  rgyr_nring/=nring;
  rgyr_nring=sqrt(rgyr_nring);

  double rgyr_nscoil=0.0;
  for (i=0; i<nscoil; i++) {
      double bnorm=sbead[i]->pos.norm();
      rgyr_nscoil+=bnorm*bnorm;
  }

  rgyr_nscoil/=nscoil;
  rgyr_nscoil=sqrt(rgyr_nscoil);

  double rgyr=0.0;
  for (i=0; i<nring; i++) {
      double bnorm=rbead[i]->pos.norm();
      rgyr+=bnorm*bnorm;
  }

  for (i=0; i<nscoil; i++) {
      double bnorm=sbead[i]->pos.norm();
      rgyr+=bnorm*bnorm;
  }

  rgyr/=(nring+nscoil);
  rgyr=sqrt(rgyr);
 
  double length=0.0;
  double averagering=0.0;
  double averagescoil=0.0;
  for (i=0;i<nring;i++) {
   length+=rbond[i]->length();
    averagering = length/nring;
  }
  for (i=0;i<nscoil;i++) {
    length+=sbond[i]->length();
    averagescoil = length/nscoil;
  }

  *rgyr_ring = rgyr_nring;
  *rgyr_scoil = rgyr_nscoil;
  *rgyr_out = rgyr;
  *averagering_out = averagering;
  *averagescoil_out = averagescoil;
}                                                                                                     

void checkRestraints(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, int *bpmap, Restraint *restlist, int nrest) {
  double rener=0.0;
  int i;
  for (i=0; i<nrest; i++) {
    int ri=int(restlist[i].inx/cgbpbead+0.5);
    int rj=int(restlist[i].jnx/cgbpbead+0.5);
    ri=(ri+basemapoffset)%natom;
    rj=(rj+basemapoffset)%natom;
    if (ri<maxcgbp && rj<maxcgbp) {
      int bpi=bpmap[ri];  
      int bpj=bpmap[rj];
      if (bpi>=0 && bpj>=0) {
        Bead *tb1=0;
        if (bpi<nring) {
          tb1=rbead[bpi];
        } else {
          tb1=sbead[bpi-nring];
        } 
 
        Bead *tb2=0; 
        if (bpj<nring) {
          tb2=rbead[bpj];
        } else {
          tb2=sbead[bpj-nring];
        } 
        Vector d=tb1->pos-tb2->pos;
        double dd=d.norm();
        dd-=restlist[i].dist;
        if (dd<0 || restlist[i].type==1) {
          rener+=restlist[i].force*(dd*dd);
        }
      }
    }
  }
}

class thread_data {
  public:
  int tid;
  Bond **rbond;
  Bead **rbead;
  int nring;
  Bond **sbond;
  Bead **sbead;
  int nscoil;
  int *bondlookupinx;
  int *nbondlookup;
  Bond **bondlookup;
  Bead *sb;
  int *bpmap;
  Restraint *restlist;
  int nrest;
  double cutoff;

  double energy;
  double bondenergy;
  double ibondenergy;
  double rmaxenergy;
  double rener;
  double rgyr;

  thread_data(int t, Bond **rbnd, Bead **rbd, int nrng, Bond **sbnd, Bead **sbd, int nsc, int *blupinx, 
              int *nbl, Bond **bl, Bead *s, int *bpm, Restraint *rlist, int nr, double cut=0.0) : 
    tid(t), rbond(rbnd), rbead(rbd), nring(nrng), sbond(sbnd), sbead(sbd), nscoil(nsc), bondlookupinx(blupinx),nbondlookup(nbl),bondlookup(bl),
    sb(s), bpmap(bpm), restlist(rlist), nrest(nr), cutoff(cut), energy(0.0), bondenergy(0.0), ibondenergy(0.0), rmaxenergy(0.0), rener(0.0), rgyr(0.0) {}
};

void *gCutEner(void *tdata) {
   thread_data *td=(thread_data *)tdata;

   int i;
   for (i=td->tid; i<td->nring; i+=num_threads) {
     double tval=td->rbond[i]->selfEnergy();
     td->energy+=tval;
     td->bondenergy+=tval;
     for (int j=0; j<td->nbondlookup[i]; j++) {
       int inx=j+td->bondlookupinx[i];
       double tval=td->rbond[i]->interBondEnergy(td->bondlookup[inx]);
       td->energy+=tval;
       td->ibondenergy+=tval; 
     }
   }

   for (i=td->tid; i<td->nscoil; i+=num_threads) {
     double tval=td->sbond[i]->selfEnergy();
     td->energy+=tval;
     td->bondenergy+=tval;
     for (int j=0; j<td->nbondlookup[i+td->nring]; j++) {
       int inx=j+td->bondlookupinx[i+td->nring];
       double tval=td->sbond[i]->interBondEnergy(td->bondlookup[inx]);
       td->energy+=tval;
       td->ibondenergy+=tval;
     } 
   }

   for (i=td->tid; i<td->nring; i+=num_threads) {
     td->rgyr+=td->rbead[i]->rgyr;
     td->energy+=td->rbead[i]->rgyrener;
     td->rmaxenergy+=td->rbead[i]->rgyrener;
   }

   for (i=td->tid; i<td->nscoil; i+=num_threads) {
     td->rgyr+=td->sbead[i]->rgyr;
     td->energy+=td->sbead[i]->rgyrener;
     td->rmaxenergy+=td->sbead[i]->rgyrener;
   }

   for (i=td->tid; i<td->nrest; i+=num_threads) {
     int ri=int(td->restlist[i].inx/cgbpbead+0.5);
     int rj=int(td->restlist[i].jnx/cgbpbead+0.5);
     ri=(ri+basemapoffset)%natom;
     rj=(rj+basemapoffset)%natom;
     if (ri<maxcgbp && rj<maxcgbp) {
       int bpi=td->bpmap[ri];  
       int bpj=td->bpmap[rj];
       if (bpi>=0 && bpj>=0) {
         Bead *tb1=0;
         if (bpi<td->nring) {
           tb1=td->rbead[bpi];
         } else {
           tb1=td->sbead[bpi-td->nring];
         } 
 
        Bead *tb2=0; 
        if (bpj<td->nring) {
          tb2=td->rbead[bpj];
        } else {
          tb2=td->sbead[bpj-td->nring];
        } 
        Vector d=tb1->pos-tb2->pos;
        double dd=d.norm();
        dd-=td->restlist[i].dist;
        if (dd<0 || td->restlist[i].type==1) {
          td->rener+=td->restlist[i].force*(dd*dd);
        }
      }
    }
   }
   td->rener*=0.5;
   td->energy+=td->rener;
  
   pthread_exit(tdata);
}

double tgetCutoffEnergy(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, 
                        int *bondlookupinx, int *nbondlookup, Bond **bondlookup, int *bpmap, Restraint *restlist, int nrest, int check=0) {
  pthread_t threads[MAX_THREADS];
  thread_data *ta[MAX_THREADS];

  for (int t=0; t<num_threads; t++) {
    ta[t]=new thread_data(t,rbond,rbead,nring,sbond,sbead,nscoil,bondlookupinx,nbondlookup,bondlookup,0,bpmap,restlist,nrest);

    int rc=pthread_create(&threads[t],NULL,gCutEner, (void *)ta[t]);
    if (rc) {
      fprintf(stderr,"Could not create thread. Error code: %d\n",rc);
      exit(1);
    }
  }
  for (int t=0; t<num_threads; t++) {
    pthread_join(threads[t],NULL); 
  }

  double energy=0.0;
  double bondenergy=0.0;
  double ibondenergy=0.0;
  double rmaxenergy=0.0;
  double rener=0.0;
  double rgyr=0.0;
  for (int t=0; t<num_threads; t++) {
    energy+=ta[t]->energy;
    bondenergy+=ta[t]->bondenergy;
    ibondenergy+=ta[t]->ibondenergy;
    rmaxenergy+=ta[t]->rmaxenergy;
    rener+=ta[t]->rener;
    rgyr+=ta[t]->rgyr;
    delete ta[t];
  }

  rgyr/=(nring+nscoil);
  rgyr=sqrt(rgyr);
  double drgyr=(rgyr-rgyr0);
  double rgyrenergy=drgyr*drgyr*rgyrk/2.0; 
  energy+=drgyr*drgyr*rgyrk/2.0;

  if (check) {
    fprintf(stderr,"cutoff energy: %lf, bond: %lf, ibond: %lf, rgyr: %lf, rmax: %lf, rest: %lf\n",energy,bondenergy,ibondenergy,rgyrenergy,rmaxenergy,rener);
  }
  return energy;
}

double getCutoffEnergy(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, 
                       int *bondlookupinx, int *nbondlookup, Bond **bondlookup, int *bpmap, Restraint *restlist, int nrest, int check=0) {
  double energy=0.0;
  double bondenergy=0.0;
  double ibondenergy=0.0;
  double rgyrenergy=0.0;
  double rmaxenergy=0.0;

  int i;
  for (i=0; i<nring; i++) {
    double tval=rbond[i]->selfEnergy();
    energy+=tval;
    bondenergy+=tval;

    for (int j=0; j<nbondlookup[i]; j++) {
      int inx=j+bondlookupinx[i];
      double tval=rbond[i]->interBondEnergy(bondlookup[inx]);
      energy+=tval;
      ibondenergy+=tval;
    } 
  }

  for (i=0; i<nscoil; i++) {
    double tval=sbond[i]->selfEnergy();
    energy+=tval;
    bondenergy+=tval;
    for (int j=0; j<nbondlookup[i+nring]; j++) {
      int inx=j+bondlookupinx[i+nring];
      double tval=sbond[i]->interBondEnergy(bondlookup[inx]);
      energy+=tval;
      ibondenergy+=tval;
    } 
  }

  double rgyr=0.0;
  for (i=0; i<nring; i++) {
    rgyr+=rbead[i]->rgyr;
    energy+=rbead[i]->rgyrener;
    rmaxenergy+=rbead[i]->rgyrener;
  }

  for (i=0; i<nscoil; i++) {
    rgyr+=sbead[i]->rgyr;
    energy+=sbead[i]->rgyrener;
    rmaxenergy+=sbead[i]->rgyrener;
  }

  rgyr/=(nring+nscoil);
  rgyr=sqrt(rgyr);
  double drgyr=(rgyr-rgyr0);
  energy+=drgyr*drgyr*rgyrk/2.0;
  rgyrenergy+=drgyr*drgyr*rgyrk/2.0;
 
  double rener=0.0;
  for (i=0; i<nrest; i++) {
    int ri=int(restlist[i].inx/cgbpbead+0.5);
    int rj=int(restlist[i].jnx/cgbpbead+0.5);
    ri=(ri+basemapoffset)%natom;
    rj=(rj+basemapoffset)%natom;
    if (ri<maxcgbp && rj<maxcgbp) {
      int bpi=bpmap[ri];  
      int bpj=bpmap[rj];
      if (bpi>=0 && bpj>=0) {
        Bead *tb1=0;
        if (bpi<nring) {
          tb1=rbead[bpi];
        } else {
          tb1=sbead[bpi-nring];
        } 
 
        Bead *tb2=0; 
        if (bpj<nring) {
          tb2=rbead[bpj];
        } else {
          tb2=sbead[bpj-nring];
        } 
        Vector d=tb1->pos-tb2->pos;
        double dd=d.norm();
        dd-=restlist[i].dist;
        if (dd<0 || restlist[i].type==1) {
          rener+=restlist[i].force*(dd*dd);
        }
      }
    }
  }
  rener*=0.5;
  energy+=rener;

  if (check) {
    fprintf(stderr," cutoff energy: %lf, bond: %lf, ibond: %lf, rgyr: %lf, rmax: %lf, rest: %lf\n",energy,bondenergy,ibondenergy,rgyrenergy,rmaxenergy,rener);
  }
  return energy;
}

void *gSinCutEner(void *tdata) {
   thread_data *td=(thread_data *)tdata;

   int i;
   for (i=td->tid; i<td->nring; i+=num_threads) {
     if (td->rbond[i]->bead1==td->sb || td->rbond[i]->bead2==td->sb) {
       double tval=td->rbond[i]->selfEnergy();
       td->energy+=tval;
       for (int j=0; j<td->nbondlookup[i]; j++) {
         int inx=j+td->bondlookupinx[i];
 	 double tval=td->rbond[i]->interBondEnergy(td->bondlookup[inx]);
	 td->energy+=tval;
       }
     } else {
       for (int j=0; j<td->nbondlookup[i]; j++) {
         int inx=j+td->bondlookupinx[i];
         if (td->bondlookup[inx]->bead1==td->sb || td->bondlookup[inx]->bead2==td->sb) {
           double tval=td->rbond[i]->interBondEnergy(td->bondlookup[inx]);
	   td->energy+=tval;
         }
       } 
     }
   }

   for (i=td->tid; i<td->nscoil; i+=num_threads) {
    if (td->sbond[i]->bead1==td->sb || td->sbond[i]->bead2==td->sb) {
      double tval=td->sbond[i]->selfEnergy();
      td->energy+=tval;
      for (int j=0; j<td->nbondlookup[i+td->nring]; j++) {
        int inx=j+td->bondlookupinx[i+td->nring];
  	double tval=td->sbond[i]->interBondEnergy(td->bondlookup[inx]);
	td->energy+=tval;
      } 
    } else {
      for (int j=0; j<td->nbondlookup[i+td->nring]; j++) {
        int inx=j+td->bondlookupinx[i+td->nring];
        if (td->bondlookup[inx]->bead1==td->sb || td->bondlookup[inx]->bead2==td->sb) {
  	  double tval=td->sbond[i]->interBondEnergy(td->bondlookup[inx]);
  	  td->energy+=tval;
        }
      }
    } 
   }

   for (i=td->tid; i<td->nring; i+=num_threads) {
     td->rgyr+=td->rbead[i]->rgyr;
     td->energy+=td->rbead[i]->rgyrener;
   }

   for (i=td->tid; i<td->nscoil; i+=num_threads) {
     td->rgyr+=td->sbead[i]->rgyr;
     td->energy+=td->sbead[i]->rgyrener;
   }

 
   double rener=0.0;
   for (i=td->tid; i<td->nrest; i+=num_threads) {
     int ri=int(td->restlist[i].inx/cgbpbead+0.5);
     int rj=int(td->restlist[i].jnx/cgbpbead+0.5);
     ri=(ri+basemapoffset)%natom;
     rj=(rj+basemapoffset)%natom;
     if (ri<maxcgbp && rj<maxcgbp) {
       int bpi=td->bpmap[ri];  
       int bpj=td->bpmap[rj];
       if (bpi>=0 && bpj>=0) {
         Bead *tb1=0;
         if (bpi<td->nring) {
           tb1=td->rbead[bpi];
         } else {
           tb1=td->sbead[bpi-td->nring];
         } 
 
        Bead *tb2=0; 
        if (bpj<td->nring) {
          tb2=td->rbead[bpj];
        } else {
          tb2=td->sbead[bpj-td->nring];
        } 
        if (td->sb==tb1 || td->sb==tb2 ) {
          Vector d=tb1->pos-tb2->pos;
          double dd=d.norm();
          dd-=td->restlist[i].dist;
          if (dd<0 || td->restlist[i].type==1) {
            rener+=td->restlist[i].force*(dd*dd);
          }
        }
      }
    }
   }
   rener*=0.5;
   td->energy+=rener;
  
   pthread_exit(tdata);
}

double tgetSingleCutoffEnergy(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, 
                              int *bondlookupinx, int *nbondlookup, Bond **bondlookup, Bead *sb, int *bpmap, Restraint *restlist, int nrest) {
  pthread_t threads[MAX_THREADS];
  thread_data *ta[MAX_THREADS];

  for (int t=0; t<num_threads; t++) {
    ta[t]=new thread_data(t,rbond,rbead,nring,sbond,sbead,nscoil,bondlookupinx,nbondlookup,bondlookup,sb,bpmap,restlist,nrest);

    int rc=pthread_create(&threads[t],NULL,gSinCutEner, (void *)ta[t]);
    if (rc) {
      fprintf(stderr,"Could not create thread. Error code: %d\n",rc);
      exit(1);
    }
  }
  for (int t=0; t<num_threads; t++) {
    pthread_join(threads[t],NULL); 
  }

  double energy=0.0;
  double rgyr=0.0;
  for (int t=0; t<num_threads; t++) {
    energy+=ta[t]->energy;
    rgyr+=ta[t]->rgyr;
    delete ta[t];
  }

  rgyr/=(nring+nscoil);
  rgyr=sqrt(rgyr);
  double drgyr=(rgyr-rgyr0);
  energy+=drgyr*drgyr*rgyrk/2.0;

  return energy;
}

void updateLookup(Bond **rbond, int nring, Bond **sbond, int nscoil, int *bondlookupinx, int *nbondlookup, Bond **bondlookup, double cutoff, Bond *b) {
  int i;
  int inx=bondlookupinx[b->index];
  int nbl=nbondlookup[b->index];
  for (i=0; i<nbl; i++, inx++) {
    Bond *bt=bondlookup[inx];
    int jnx=bondlookupinx[bt->index];
    int jnbl=nbondlookup[bt->index];
    for (int j=0; j<jnbl; j++) {
      if (bondlookup[jnx+j]==b) {
        jnbl--;
        bondlookup[jnx+j]=bondlookup[jnbl+jnx];
      }
    }
    nbondlookup[bt->index]=jnbl; 
  }
  nbondlookup[b->index]=0;

  for (i=0; i<nring; i++) {
    if (i!=b->index && rbond[i]->minimumDistance(b)<cutoff) {
      bondlookup[bondlookupinx[i]+nbondlookup[i]++]=b;
      bondlookup[bondlookupinx[b->index]+nbondlookup[b->index]++]=rbond[i];
    } 
  }

  for (i=0; i<nscoil; i++) {
    if (i+nring!=b->index && sbond[i]->minimumDistance(b)<cutoff) {
      bondlookup[bondlookupinx[i+nring]+nbondlookup[i+nring]++]=b;
      bondlookup[bondlookupinx[b->index]+nbondlookup[b->index]++]=sbond[i];
    } 
  }
}

void *genLook(void *tdata) {
  thread_data *td=(thread_data *)tdata;
 
  int i;
  for (i=td->tid; i<td->nring; i+=num_threads) {
    for (int j=i+1; j<td->nring; j++) {
      if (td->rbond[i]->minimumDistance(td->rbond[j])<td->cutoff) {
        td->bondlookup[td->bondlookupinx[i]+td->nbondlookup[i]++]=td->rbond[j];
      }
    }
    for (int j=0; j<td->nscoil; j++) {
      if (td->rbond[i]->minimumDistance(td->sbond[j])<td->cutoff) {
        td->bondlookup[td->bondlookupinx[i]+td->nbondlookup[i]++]=td->sbond[j];
      }
    } 
  }

  for (i=td->tid; i<td->nscoil; i+=num_threads) {
    for (int j=i+1; j<td->nscoil; j++) {
      if (td->sbond[i]->minimumDistance(td->sbond[j])<td->cutoff) {
        td->bondlookup[td->bondlookupinx[i+td->nring]+td->nbondlookup[i+td->nring]++]=td->sbond[j];
      }
    }
  }
}

void tgenerateLookup(Bond **rbond, int nring, Bond **sbond, int nscoil, int *bondlookupinx, int *nbondlookup, Bond **bondlookup, double cutoff) {
  pthread_t threads[MAX_THREADS];
  thread_data *ta[MAX_THREADS];

  int i;
  for (i=0; i<nscoil+nring; i++) {
    bondlookupinx[i]=i*nscoil;
    nbondlookup[i]=0;
  }

  for (int t=0; t<num_threads; t++) {
    ta[t]=new thread_data(t,rbond,0,nring,sbond,0,nscoil,bondlookupinx,nbondlookup,bondlookup,0,0,0,0,cutoff);

    int rc=pthread_create(&threads[t],NULL,genLook,(void *)ta[t]);
    if (rc) {
      fprintf(stderr,"Could not create thread. Error code: %d\n",rc);
      exit(1);
    }
  }
  for (int t=0; t<num_threads; t++) {
    pthread_join(threads[t],NULL); 
  }

  for (i=0; i<nscoil+nring; i++) {
    int nbl=nbondlookup[i];
    for (int j=0; j<nbl; j++) {
      Bond *bt=bondlookup[bondlookupinx[i]+j];
      bondlookup[bondlookupinx[bt->index]+nbondlookup[bt->index]++]=bt;
    }
  }
}

void dumpLookupTable(Bond **rbond, int nring, Bond **sbond, int nscoil, int *bondlookupinx, int *nbondlookup, Bond **bondlookup) { 
  for (int i=0; i<nring+nscoil; i++) {
    int inx=bondlookupinx[i];
    for (int j=0; j<nbondlookup[i]; j++) {
      fprintf(stderr,"lookup %d %d %d\n",i,j,inx+j);
    } 
  }
}

void dumpBondTable(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil) {
  for (int i=0; i<nring; i++) {
    fprintf(stderr,"rbond %d inx %d beads %d %d prev %d links %d to ",i,rbond[i]->index,rbond[i]->bead1->index,rbond[i]->bead2->index,rbond[i]->prev->index,rbond[i]->nlink);
    for (int j=0; j<rbond[i]->nlink; j++) {
      fprintf(stderr,"%d (%d:%d) ",rbond[i]->next[j]->index,rbond[i]->next[j]->bead1->index,rbond[i]->next[j]->bead2->index);
    }
    fprintf(stderr,"\n");
  }
  for (int i=0; i<nscoil; i++) {
    fprintf(stderr,"sbond %d inx %d beads %d %d prev %d links %d to ",i,sbond[i]->index,sbond[i]->bead1->index,sbond[i]->bead2->index,sbond[i]->prev->index,sbond[i]->nlink);
    for (int j=0; j<sbond[i]->nlink; j++) {
      fprintf(stderr,"%d (%d:%d) ",sbond[i]->next[j]->index,sbond[i]->next[j]->bead1->index,sbond[i]->next[j]->bead2->index);
    }
    fprintf(stderr,"\n");
  }
}

void writePDB(FILE *fptr,int n, double x, double y, double z, const char *atomname, const char *resname, char c, const char *seg) {
  fprintf(fptr,
     "ATOM %6d  %-3s%c%-4s%c%5d   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
     n,atomname,' ',resname,c,n,x,y,z,0.0,0.0,seg);
}

void dumpPDB(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, char *filename) {
  char pdbfilename[512];
  sprintf(pdbfilename,"%s.pdb",filename);
  FILE *fptr=fopen(pdbfilename,"w");

  for (int i=0; i<nring; i++) {
    Bead *b=rbead[i];
    writePDB(fptr,b->index+1,b->pos.x()/10.0,b->pos.y()/10.0,b->pos.z()/10.0,"C","CDN",'A',"RING");
  }
  for (int i=0; i<nscoil; i++) {
    Bead *b=sbead[i];
    writePDB(fptr,b->index+1,b->pos.x()/10.0,b->pos.y()/10.0,b->pos.z()/10.0,"P","PDN",'B',"COIL");
  }

  fprintf(fptr,"TER\n");
  for (int i=0; i<nring; i++) {
    fprintf(fptr,"CONECT%5d%5d\n",rbond[i]->bead1->index+1,rbond[i]->bead2->index+1);
  }
  for (int i=0; i<nscoil; i++) {
    fprintf(fptr,"CONECT%5d%5d\n",sbond[i]->bead1->index+1,sbond[i]->bead2->index+1);
  }
  fprintf(fptr,"END\n");

  fclose(fptr);
}

void dumpCOOR(Bead **rbead, int nring, Bead **sbead, int nscoil, char *filename) {
  char pdbfilename[512];
  sprintf(pdbfilename,"%s.coord",filename);
  FILE *fptr=fopen(pdbfilename,"w");

  for (int i=0; i<nring; i++) {
    Bead *b=rbead[i];
    fprintf(fptr,"%d %20.12lf %20.12lf %20.12lf\n",
            b->index+1,b->pos.x(),b->pos.y(),b->pos.z());
  }
  for (int i=0; i<nscoil; i++) {
    Bead *b=sbead[i];
    fprintf(fptr,"%d %20.12lf %20.12lf %20.12lf\n",
            b->index+1,b->pos.x(),b->pos.y(),b->pos.z());
  }

  fclose(fptr);
}

void dumpRestraints(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, int *bpmap, Restraint *restlist, int nrest, char *filename) {
  if (nrest>0) {
   char pdbfilename[512];
   sprintf(pdbfilename,"%s.rest",filename);
   FILE *fptr=fopen(pdbfilename,"w");

   fprintf(fptr,"MMFP\n");
   fprintf(fptr,"geo maxgeo %d sphere rcm force 0.01 select all end\n",nrest*5);
   for (int i=0; i<nrest; i++) {
    int ri=int(restlist[i].inx/cgbpbead+0.5);
    int rj=int(restlist[i].jnx/cgbpbead+0.5);
    ri=(ri+basemapoffset)%natom;
    rj=(rj+basemapoffset)%natom;
    if (ri<maxcgbp && rj<maxcgbp) {
      int bpi=bpmap[ri];  
      int bpj=bpmap[rj];
      if (bpi>=0 && bpj>=0) {
        Bead *tb1=0;
        if (bpi<nring) {
          tb1=rbead[bpi];
        } else {
          tb1=sbead[bpi-nring];
        } 
 
        Bead *tb2=0; 
        if (bpj<nring) {
          tb2=rbead[bpj];
        } else {
          tb2=sbead[bpj-nring];
        } 

        if (restlist[i].type == 1) {
          if (tb1->index+1<nring && tb2->index+1<nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d C %d C\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1>=nring && tb2->index+1>=nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d P %d P\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1>=nring && tb2->index+1<nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d P %d C\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1<nring && tb2->index+1>=nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d C %d P\n",tb1->index+1,tb2->index+1);
          }
        } else {
          if (tb1->index+1<nring && tb2->index+1<nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf NEGATIVE -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d C %d C\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1>=nring && tb2->index+1>=nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf NEGATIVE -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d P %d P\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1>=nring && tb2->index+1<nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf NEGATIVE -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d P %d C\n",tb1->index+1,tb2->index+1);
          } else if (tb1->index+1<nring && tb2->index+1>=nring) {
            fprintf(fptr,"RESDistance KVAL %lf RVAL %lf NEGATIVE -\n",restlist[i].force*100.0,restlist[i].dist/10.0);
            fprintf(fptr,"  1.0 %d C %d P\n",tb1->index+1,tb2->index+1);
          }
        }
        }
      }
    }
   fprintf(fptr,"END\n");
   fprintf(fptr,"rgyr force %lf ref %lf select all end\n",rgyrk*100.0,rgyr0);

   fclose(fptr);
  }
}
 
void readPDB(char *pdbfilename, PDBEntry *pdb, int& natom, Connection *connect, int& nconnect) {
  FILE *fptr;

  int lastnum=-1;
  natom=0;

  fptr=fopen(pdbfilename,"r");
  if (fptr==0) {
    fprintf(stderr,"Cannot open PDB file %s\n",pdbfilename);
    exit(1);
  }

  while (!feof(fptr)) {
    if (pdb[natom].read(fptr)>0 && pdb[natom].residueNumber()!=lastnum) {
        lastnum=pdb[natom].residueNumber();
        natom++;
    } 
  }
  fclose(fptr);

  nconnect=0;
  char line[1024];
 
  fptr=fopen(pdbfilename,"r");
  while (!feof(fptr)) {
    if (fgets(line,1024,fptr)) {
      if (!strncmp(line,"CONECT",6)) {
         int inx=atoi(PDBEntry::substr(line,6,5));
         int jnx=atoi(PDBEntry::substr(line,11,5));
         connect[nconnect].inx=inx;
         connect[nconnect].jnx=jnx;
         nconnect++;
      } 
    }
  }
  fclose(fptr);
}

void readCOOR(Bead **rbead, int nring, Bead **sbead, int nscoil, char *filename) {
  char pdbfilename[512];
  sprintf(pdbfilename,"%s.coord",filename);
  FILE *fptr=fopen(pdbfilename,"r");

  int inx;
  double x,y,z;
  for (int i=0; i<nring; i++) {
    Bead *b=rbead[i];
    if (fscanf(fptr,"%d%lf%lf%lf",&inx,&x,&y,&z)) {
      b->pos=Vector(x,y,z);
      b->calcRgyrEner();
    } else {
      fprintf(stderr,"data missing\n");
    }
  }
  for (int i=0; i<nscoil; i++) {
    Bead *b=sbead[i];
    if (fscanf(fptr,"%d%lf%lf%lf",&inx,&x,&y,&z)) {
      b->pos=Vector(x,y,z);
      b->calcRgyrEner();
    } else {
      fprintf(stderr,"data missing\n");
    }
  }

  fclose(fptr);
}

void dumpPSF(Bond **rbond, Bead **rbead, int nring, Bond **sbond, Bead **sbead, int nscoil, char *filename) {
  char psffilename[512];
  sprintf(psffilename,"%s.psf",filename);
  FILE *fptr=fopen(psffilename,"w");

  fprintf(fptr,"PSF EXT CMAP CHEQ\n\n");
  fprintf(fptr,"%10d !NTITLE\n",2);
  fprintf(fptr,"* TITLE\n* automatically generated by mcnuc2\n\n");
  fprintf(fptr,"%10d !NATOM\n",nring+nscoil);

  int i;
  for (i=0; i<nring; i++) {
    Bead *b=rbead[i]; 
    fprintf(fptr,"%10d %-8s %-8d %-8s %-8s %4d %10.5f %13.3f        %4d%10.5f     %12.6e\n",
       b->index+1,"RING",b->index+1,"CDN","C",1,0.0,100.0,0,0.0,-0.0030114);
  } 

  for (i=0; i<nscoil; i++) {
    Bead *b=sbead[i]; 
    fprintf(fptr,"%10d %-8s %-8d %-8s %-8s %4d %10.5f %13.3f        %4d%10.5f     %12.6e\n",
       b->index+1,"COIL",b->index+1,"PDN","P",2,0.0,100.0,0,0.0,-0.0030114);
  } 
  
  fprintf(fptr,"\n%10d !NBOND: bonds\n",nring+nscoil);
  int n=0;
  for (i=0; i<nring; i++) {
    fprintf(fptr,"%10d%10d",rbond[i]->bead1->index+1,rbond[i]->bead2->index+1);
    if (++n%4==0) fprintf(fptr,"\n");
  }
  for (i=0; i<nscoil; i++) {
    fprintf(fptr,"%10d%10d",sbond[i]->bead1->index+1,sbond[i]->bead2->index+1);
    if (++n%4==0) fprintf(fptr,"\n");
  }
  if (n%4!=0) fprintf(fptr,"\n");

  fprintf(fptr,"\n%10d !NTHETA: angles\n\n",0);
  fprintf(fptr,"\n%10d !NPHI: dihedrals\n\n",0);
  fprintf(fptr,"\n%10d !NIMPHI: impropers\n\n",0);
  fprintf(fptr,"\n%10d !NDON: donors\n\n",0);
  fprintf(fptr,"\n%10d !NACC: acceptors\n\n",0);
  fprintf(fptr,"\n%10d !NNB\n\n",0);
  
  n=0;
  for (i=0; i<nring+nscoil; i++) {
    fprintf(fptr,"%10d",0);
    if (++n%8==0) fprintf(fptr,"\n");
  }
  if (n%8!=0) fprintf(fptr,"\n");

  fprintf(fptr,"\n%10d%10d !NGRP NST2\n",nring+nscoil,0);

  for (i=0; i<nring+nscoil; i++) {
    fprintf(fptr,"%10d%10d%10d",i,0,0);
    if (i%3==0) fprintf(fptr,"\n");
  }
  if ((i-1)%3!=0) fprintf(fptr,"\n");

  fprintf(fptr,"\n%10d !MOLNT\n",2);

  n=0;
  for (i=0; i<nring; i++) {
    fprintf(fptr,"%10d",1);
    if (++n%8==0) fprintf(fptr,"\n");
  }
  for (i=0; i<nscoil; i++) {
    fprintf(fptr,"%10d",2);
    if (++n%8==0) fprintf(fptr,"\n");
  }
  if (n%8!=0) fprintf(fptr,"\n");

  fprintf(fptr,"\n%10d%10d !NUMLP NUMLPH\n\n",0,0);
  fprintf(fptr,"\n%10d !NCRTERM: cross-terms\n\n",0);
  fclose(fptr);
}

int enterBond(Bond *b, int &natom, double &bps, int *bpmap) {
  if (b->done) {
    return 0;
  }
  b->done=1;
  
  int nbeads;
  if (b->bead2->type==1) {
    nbeads=cgbeadperring;
    bps+=cgbeadperring*15;
  } else if (b->bead2->type==2) {
    nbeads=cgbeadperchain/2;
    bps+=cgbeadperchain*15;    
  } else if  (b->bead2->type==3) {
    nbeads=(cgbeadperchain+beadforhalfchain+1)/2;
    bps+=(cgbeadperchain+beadforhalfchain)*15;
  }

  int binx1=b->bead1->index;
  int binx2=b->bead2->index;

  int i;
  for (i=0; i<nbeads/2; i++) {
    bpmap[natom++]=binx1;
  } 
  for (;i<nbeads; i++) {
    bpmap[natom++]=binx2;
  }

  if (b->nlink==1) {
    enterBond(b->next[0],natom,bps,bpmap);
  } else if (b->nlink==2) {
    if (b->next[0]->index>b->next[1]->index) {
      enterBond(b->next[0],natom,bps,bpmap);
      enterBond(b->next[1],natom,bps,bpmap);
    } else {
      enterBond(b->next[1],natom,bps,bpmap);
      enterBond(b->next[0],natom,bps,bpmap);
    }
  } 
  if (b->bead2->type == 2) {
    for (i=0; i<nbeads/2; i++) { 
      bpmap[natom++]=b->bead2->index;
    }
    for (;i<nbeads; i++) {
      bpmap[natom++]=b->bead1->index;
    }
  }
  if (b->bead2->type == 3) {
    for (i=0; i<nbeads/2; i++) {
      bpmap[natom++]=b->bead2->index;
    }
    for (;i<nbeads-1; i++) {
      bpmap[natom++]=b->bead1->index;
    }
  }
  return 0;
}

int mapBasePairs(int *bpmap, Bond **ringbond, int ndomain, 
                             Bond **chainbond, int nchain, Bead **ringbead, Bead **chainbead) {

  for (int i=0; i<maxcgbp; bpmap[i++]=-1);
  for (int i=0; i<ndomain; i++) 
    ringbond[i]->done=0;
  for (int i=0; i<nchain; i++)
    chainbond[i]->done=0;
  
  natom=0;
  double bps=0.0;
  enterBond(ringbond[0],natom,bps,bpmap);
  fprintf(stderr,"natom: %d bps %.2f\n",natom,bps);

  return 0;
}  

int main(int argc, char **argv) {
  int maxtrials=10000;
  char fname[512];
  char rfnameharm[512];
  char rfnamelower[512];

  srandom(time(NULL));

  maxtrials=atoi(argv[1]);
  strcpy(fname,argv[2]);

  Restraint *restlist; 
  int nrest=0;
 
  if (argc>3) {
    strcpy(rfnameharm,argv[3]);
    strcpy(rfnamelower,argv[4]);

    int maxrest=1000000;
    restlist=new Restraint[maxrest];

    readRestraints(rfnameharm,restlist,nrest,maxrest,1); 
    fprintf(stderr,"read %d restraints from %s\n",nrest,rfnameharm); 
    readRestraints(rfnamelower,restlist,nrest,maxrest,2); 
    fprintf(stderr,"read %d restraints from %s\n",nrest,rfnamelower); 
  }

  if (argc>5) {
    ndomain=atoi(argv[5]);
  }

  if (argc>6) {
    setchainlen=atoi(argv[6]);
  }

  int totalcgbeads=(int)(mbps/(int)(cgbeadsize*bpssingle+0.5)); // totalcgbeads: 269529
  int totalringcgbeads=ndomain*setringlen/cgbeadsize; // 200*20/5
  int totalchaincgbeads=totalcgbeads-totalringcgbeads; //269529-200*20/5

  cgbeadperring=totalringcgbeads/ndomain; // 20/5
  cgbeadperchain=(int)((setchainlen*bpsscoil/(int)(cgbeadsize*bpssingle+0.5))+0.5); //40*6.69/15
  if (cgbeadperchain%2!=0) {
   cgbeadperchain++;
  }

  cutofflen=setchainlen*2.0;
  deltadisp=setchainlen/2.0;
 
  Bead **ringbead;
  Bead **chainbead;

  int nchain=0;

  Bond **ringbond;
  Bond **chainbond;

  char oldpdbfilename[512];
  sprintf(oldpdbfilename,"%s.pdb",fname);
 
  PDBEntry *dnapdb;
  dnapdb = new PDBEntry[100000];
  int natompdb;

  Connection *connectpdb;
  connectpdb=new Connection[100000];
  int nconnect;

  readPDB(oldpdbfilename, dnapdb, natompdb, connectpdb, nconnect);
  
  fprintf(stderr,"read %d atoms from %s\n", natompdb,oldpdbfilename);

  nchain=natompdb-ndomain;
  beadforhalfchain=totalcgbeads-(nchain*cgbeadperchain+ndomain*cgbeadperring);
  ringbead=new Bead*[ndomain];
  chainbead=new Bead*[nchain];

  for (int i=0; i<ndomain; i++) {
    ringbead[i]=new Bead(dnapdb[i].coordinates()*10.0,1,i,beadradius);
    ringbead[i]->calcRgyrEner();
  }

  for (int i=0; i<nchain; i++) {
    if (i == nchain-1) {
      chainbead[i]=new Bead(dnapdb[i+ndomain].coordinates()*10.0,3,i+ndomain,beadradius);
    } else {
      chainbead[i]=new Bead(dnapdb[i+ndomain].coordinates()*10.0,2,i+ndomain,beadradius);
    }
    chainbead[i]->calcRgyrEner();
  }

  readCOOR(ringbead,ndomain,chainbead,nchain,fname);
  
  fprintf(stderr,"read coordinates\n");

  ringbond=new Bond*[ndomain];
  chainbond=new Bond*[nchain];

  for (int i=0; i<ndomain; i++) {
    Bead *next=(i<ndomain-1)?ringbead[i+1]:ringbead[0];
    ringbond[i]=new Bond(ringbead[i],next,ringradius,bpssingle,setringlen,rdistk,i);
    ringbond[i]->domain=i;
  }
 
  for (int i=0; i<ndomain; i++) {
    Bond *nextring=(i<ndomain-1)?ringbond[i+1]:ringbond[0];
    ringbond[i]->add(nextring);
  }

  int nchainbond=0; 
  for (int i=0; i<nconnect; i++) {
    if (connectpdb[i].jnx>ndomain) {
      int inx=connectpdb[i].inx;
      int jnx=connectpdb[i].jnx;

      Bead *prev=(inx<=ndomain)?ringbead[inx-1]:chainbead[inx-1-ndomain];
      Bead *next=chainbead[jnx-1-ndomain];
      if (next->type==2) {
        chainbond[nchainbond]=new Bond(prev,next,branchradius,bpsscoil,setchainlen,bdistk,nchainbond+ndomain);
      } else {
        double halfchainlen=(double)(beadforhalfchain*setchainlen/cgbeadperchain);
        chainbond[nchainbond]=new Bond(prev,next,branchradius,bpsscoil,setchainlen+halfchainlen,bdistk,nchainbond+ndomain);
      }
      nchainbond++;
    }
  }

  for (int i=0; i<ndomain; i++) {
    for (int j=0; j<nchainbond; j++) {     
      if (ringbond[i]->bead2->index == chainbond[j]->bead1->index) {
        ringbond[i]->add(chainbond[j]);
      }
    }
  }
  for (int i=0; i<nchainbond; i++) {
    for (int j=0; j<nchainbond; j++) {     
      if (chainbond[i]->bead2->index == chainbond[j]->bead1->index) {
        chainbond[i]->add(chainbond[j]);
      }
    }
  }

  delete connectpdb;
  delete dnapdb;
  fprintf(stderr,"setup complete\n");
  
  // now let's do Monte Carlo sampling

  int *bpmap=new int[maxcgbp];
  if (nrest>0) mapBasePairs(bpmap,ringbond,ndomain,chainbond,nchain,ringbead,chainbead);

  Bond **bondlookup=new Bond*[(nchain+ndomain)*(nchain+ndomain)];
  int *bondlookupinx=new int[nchain+ndomain+1];
  int *nbondlookup=new int[nchain+ndomain+1];

  tgenerateLookup(ringbond,ndomain,chainbond,nchain,bondlookupinx,nbondlookup,bondlookup,cutofflen);
  double lastenergy=getCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bpmap,restlist,nrest,1);
  fprintf(stderr,"initial energy: %10.5lf\n",lastenergy);

  double kT=initkT;
  int mctrial=1;

  int naccept=0;
  int step;
  while (mctrial++<=maxtrials) {
    Bead *bb=0;
    Bond *bd=0;
    Bond *bdt=0;
    Bond *bdp=0;

    double rx=0.0;
    double ry=0.0;
    double rz=0.0;

    double deltaE=0.0;

    int doffset=0;
    
    step=mctrial%10;

    if (step>0 && step<7) { // 1 2 3 4 5 6 - move branch point
      int ibead=random()%nchain;
      bb=chainbead[ibead];

      double delta=deltadisp;
      rx=(double((random()%9999999))/9999999.0-0.5)*delta;
      ry=(double((random()%9999999))/9999999.0-0.5)*delta;
      rz=(double((random()%9999999))/9999999.0-0.5)*delta;

      double n1=tgetSingleCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bb,bpmap,restlist,nrest);
      bb->pos+=Vector(rx,ry,rz);
      bb->calcRgyrEner();
      double n2=tgetSingleCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bb,bpmap,restlist,nrest);
      
      deltaE=n2-n1;
    } else if (step>6) {  //  7 8 9 - move ring point 
      int ibead=random()%ndomain;
      bb=ringbead[ibead];

      double delta=deltadisp/5.0;
      rx=(double((random()%9999999))/9999999.0-0.5)*delta;
      ry=(double((random()%9999999))/9999999.0-0.5)*delta;
      rz=(double((random()%9999999))/9999999.0-0.5)*delta;

      double n1=tgetSingleCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bb,bpmap,restlist,nrest);
      bb->pos+=Vector(rx,ry,rz);
      bb->calcRgyrEner();
      double n2=tgetSingleCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bb,bpmap,restlist,nrest);
      deltaE=n2-n1;     
    } else { // 0 - offset change
      doffset=random()%200+100;
      doffset=0;
      basemapoffset+=doffset;
      double newenergy=tgetCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bpmap,restlist,nrest);
      deltaE=newenergy-lastenergy;
    }

    double rnum=double((random()%9999999))/9999999.0;
    if (deltaE<=0 || exp(-(deltaE)/kT)>=rnum) {
      // accept
      lastenergy+=deltaE;

//MJ too much printing
      //fprintf(stderr,"%d accepted: %10.5lf temp: %10.5lf\n",mctrial,lastenergy,kT);
      //if (doffset!=0) {
      // fprintf(stderr," new basemapoffset: %d\n",basemapoffset);
      //}
      if (++naccept%regenlookupfreq==0) {
         tgenerateLookup(ringbond,ndomain,chainbond,nchain,bondlookupinx,nbondlookup,bondlookup,cutofflen);
         double nenergy=tgetCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bpmap,restlist,nrest,0); //last argument changed to 0 for less printing
      } 
    } else {
      // reject
      //fprintf(stderr,"rejected %d with %10.5lf\n",mctrial,lastenergy+deltaE);
      
      if (bb!=0) {
        bb->pos-=Vector(rx,ry,rz);
        bb->calcRgyrEner();
      } else if (doffset!=0) {
        basemapoffset-=doffset;
      }
    }    

    if (mctrial%savefreq==0) {
      fprintf(stderr,"saving conformation\n");
      dumpPDB(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,fname);
      dumpCOOR(ringbead,ndomain,chainbead,nchain,fname);
      dumpPSF(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,fname);
      dumpRestraints(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bpmap,restlist,nrest,fname);
    }


    if (mctrial%eneoutfreq==0) {
      //fprint(stderr, "step: %d, Temp: %lf\n", mctrial, kT);
      fprintf(stderr, "step: %d, Temp: %lf\n", mctrial, kT);
      tgetCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bpmap,restlist,nrest,1);
      double rgyr_nring,rgyr_nscoil,rgyr,averagering,averagescoil;
      getGeometry(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,&rgyr_nring,&rgyr_nscoil,&rgyr,&averagering,&averagescoil);  
      fprintf(stderr,"rgyr_nring: %lf rgyr_nscoil: %lf rgyr: %lf averagering: %lf averagescoil: %lf\n",rgyr_nring,rgyr_nscoil,rgyr,averagering,averagescoil); 
    } 
    kT*=tfactor;
  }
 
  // all finished

  if (nrest>0) mapBasePairs(bpmap,ringbond,ndomain,chainbond,nchain,ringbead,chainbead);
  double totcut=getCutoffEnergy(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bondlookupinx,nbondlookup,bondlookup,bpmap,restlist,nrest,1);
  dumpPDB(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,fname);
  dumpCOOR(ringbead,ndomain,chainbead,nchain,fname);
  dumpPSF(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,fname);
  dumpRestraints(ringbond,ringbead,ndomain,chainbond,chainbead,nchain,bpmap,restlist,nrest,fname);

  for (int i=0; i<ndomain; i++) {
    delete ringbead[i];
    delete ringbond[i];
  }
  delete ringbead;
  delete ringbond;

  for (int i=0; i<nchain; i++) {
    delete chainbead[i];
    delete chainbond[i];
  }
  delete chainbead;
  delete chainbond;

  delete bondlookup;
  delete bondlookupinx;
  delete bpmap;

  if (nrest>0) {
     delete restlist;
  }
}

