/*
 * Author:  Oye  Date:  11/26/90
 * $Id: xcedcio.c,v 1.3 2003/01/22 16:54:15 oye Exp $
 *
 *
 * Updated cvt_toLatLon with algorithm that will work better
 * at higher latitudes.  DFF, April 1, 2009.
 */

/* c------------------------------------------------------------------------ */


# include <stdio.h>
# include <time.h>
# include <math.h>
# include <string.h>
# ifdef NETCDF
# include <netcdf.h>
# endif
# include <stdlib.h>
# include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

# define RADIANS(x) ((x)*0.017453292)
# define MAX_FILES 25
# define SOMETHING_LEN 56
# define CED_BIG_ENDIAN 0
# define YES 1
# define NO 0

# include <endian.h>
# if __BYTE_ORDER == __LITTLE_ENDIAN
# define LITTLENDIAN
# endif

struct ced_file_head {
    char descriptor_name[4];	/* should contain "CED1" */
    int byte_ordering;
    int size_file;
    int reserved_word_4;
    int file_offset[MAX_FILES];
    char say_something[MAX_FILES*SOMETHING_LEN];
    int reserved_word_55;
    int reserved_word_56;
    int reserved_word_57;
    int reserved_word_58;
    int reserved_word_59;
    int reserved_word_60;
    /* this is 1540 bytes long on the Sun */
};


struct date_st
{
	long	ds_yymmdd;	/* Day portion	*/
	long	ds_hhmmss;	/* time portion */
};

typedef struct _ZebTime
{
	long	zt_Sec;		/* Seconds since 1/1/70		*/
	long	zt_MicroSec;	/* Microseconds added to zt_Sec */
} ZebTime;

typedef struct date_st UItime;	/* Different from UI "date" so we can 
				   change it. */
/*
 * Field info struct.
 */
typedef struct _FldInfo
{
	char	fi_Name[8];	/* Name of the field	*/
	short	fi_Scale;	/* Scale factor		*/
} FldInfo;

typedef struct _CoordInfo
{
	float	ci_MinVal;
	float	ci_MaxVal;
	float	ci_Spacing;
	int	ci_NStep;
} CoordInfo;

/*
 * The mudras header looks like this:
 */
typedef struct _mudras_header
{
	char	mh_Id[20];		/* 1 Ident information		*/
	short	mh_Fill1[22];		/* 11 fill */
	
	short	mh_OrLatDeg;		/* 33 Latitude part 1		*/
	short	mh_OrLatMin;		/* 34 Latitude part 2		*/
	short	mh_OrLatSec;		/* 35 Latitude part 3 (scaled)	*/
	short	mh_OrLonDeg;		/* 36 Longitude part 1		*/
	short	mh_OrLonMin;		/* 36 Longitude part 2		*/
	short	mh_OrLonSec;		/* 38 Longitude part 3 (scaled)	*/
	
	short	mh_Fill39[23];
	char	mh_Machine[2];		/* 62 2 char creating machine	*/
	short	mh_Fill63[2];
	short	mh_BlkSize;		/* 65 blocking	(3200)		*/
	short	mh_Fill66;
	short	mh_BadVal;		/* Bad value flag		*/
	short	mh_Scale;		/* 68 Scale factor		*/
	
	short	mh_Fill69[27];

	short	mh_BlkPerRec;		/* 96 Blocks per record		*/
	short	mh_Fill97[9];

	short	mh_NLevel;		/* 106 Number of levels (nZ?)	*/
	short	mh_Fill107[9];

	short	mh_Year;		/* 116 Year			*/
	short	mh_Month;		/* 117 Month			*/
	short	mh_Day;		/* 118 Day			*/
	short	mh_Hour;		/* 119 Hour			*/

	short	mh_Minute;		/* 120 Minute			*/
	short	mh_Second;		/* 121 Second			*/
	short	mh_Year1;		/* 122 Year			*/
	short	mh_Month1;		/* 123 Month			*/
	short	mh_Day1;		/* 124 Day			*/
	short	mh_Hour1;		/* 125 Hour			*/
	short	mh_Minute1;		/* 126 Minute			*/
	short	mh_Second1;		/* 127 Second			*/

	short	mh_Fill128[32];

	short	mh_XMin;		/* 160 Minimum X position	*/
	short	mh_XMax;		/* 161 Max X position		*/
	short	mh_nX;			/* 162 X resolution		*/
	short	mh_XSpacing;		/* 163 X spacing		*/
	short	mh_Fill164;
	short	mh_YMin;		/* 165 Min y			*/
	short	mh_YMax;		/* 166 Max   Y			*/
	short   mh_nY;			/* 167 Y resolution		*/
	short	mh_YSpacing;		/* 168 Y spacing		*/
	short	mh_Fill169;
	
	short	mh_ZMin;		/* 170 Min z			*/
	short	mh_ZMax;		/* 171 Max Z			*/
	short	mh_nZ;			/* 172 Z resolution		*/
	short	mh_ZSpacing;		/* 173 Z spacing		*/
	short	mh_Fill174;
	short	mh_NField;		/* 175 number of fields		*/
	FldInfo mh_Fields[25];		/* 176-300 Field information	*/
	
	short	mh_WPerPlane;		/* 301 Words per plane		*/
	short	mh_NLand1;		/* 302 number of landmarks (1)	*/
	short	mh_NLand2;		/* 303 number of landmarks (2)	*/
	
	short	mh_Fill[207];		/* Fill out the struct		*/
} MudHdr;

# define VOL_HDR_SIZE 510
# define LVL_HRR_SIZE 10


MudHdr *mh, *mhc;
char * mudbuf = NULL;
char * mudcbuf = NULL;


ZebTime zt;

struct ced_fixes {
    int fix_alt;
    int fix_lat;
    int fix_lon;
    int secs_to_UTC;
    float alt;			/* meters */
    float lat;			/* degrees */
    float lon;			/* degrees */
};


/*
 * Field info.
 */
# define MAXFLD 25

int Nfield = 0;				/* How many?		*/
char SrcFlds[MAXFLD][40];		/* Mudras field names	*/
char DstFlds[MAXFLD][40];		/* Netcdf field names	*/
int VFields[MAXFLD];			/* Netcdf var ID's	*/
int MFields[MAXFLD];			/* Mudras field numbers	*/
int VBTime;				/* Base time		*/
int VTOff;				/* Time offset		*/
int VLat, VLon, VAlt;			/* Origin position	*/
int VXs, VYs, VZs;			/* Grid spacing		*/

int Nfile;
int DTime;			/* Time dimension		*/

# ifndef PI
# define PI 3.141592654
# endif

# ifdef LITTLENDIAN
# endif

#define MAX_PATH_LEN 1024
#define PATH_DELIM "/"

/*
 * Radius of the earth, in km
 */
# define R_EARTH  6378.137

/*
 * Origin latitude and longitude (radians)
 */
static double Origin_lat = -99.0, Origin_lon = -99.0;
static double r_earth=R_EARTH;

# define BADVAL -32768

struct ced_file_info {
    int num_files;
    int offsets[MAXFLD];
    int vol_size[MAXFLD];
};

static struct ced_file_info *cfi=NULL;
static struct ced_file_head *cfhead=NULL;
static FILE *cedstrm = NULL;
static FILE *scrstrm = NULL;
static int rlens;
static char scratch_name[MAX_PATH_LEN];
static char scratch_dir[MAX_PATH_LEN];

static int makedir(const char *path);
static int ta_makedir(const char *path);
static int ta_makedir_recurse(const char *path);
static int ta_makedir_for_file(const char *file_path);

static void cvt_Origin (double lat, double lon);
static void cvt_ToLatLon (double x, double y, double *lat, double *lon);
static void cvt_ToXY (double lat, double lon, double *x, double *y);

static void setScratchDir();

/********************************************************/

/* set scratch dir */

static void setScratchDir()
{

  if (getenv("SCRATCH") == NULL) {
    printf("\n\n");
    printf("**** WARNING: SCRATCH var not set.\n");
    printf("***  Writing scratch files to /tmp\n");
    sprintf(scratch_dir, "/tmp");
  } else {
    sprintf(scratch_dir, getenv("SCRATCH"));
    if (ta_makedir_recurse(scratch_dir)) {
      perror(scratch_dir);
      fprintf(stderr, "Cannot make output dir: %s\n", scratch_dir);
      exit(1);
    }
  }

}

/* c------------------------------------------------------------------------ */

int reonetcdf_(int *prefix,
               int *np)
{
    /* c...mark */
    int ii, nn, nfields=0;	/* just use the fields and names as is */
    int vol_num=1;
    char file_name_prefix[MAX_PATH_LEN], str[MAX_PATH_LEN], *a, *b, *c;
    char **src_ids, **dst_ids;
    float rr_earth=6300.;
    struct ced_fixes cf;

# ifdef NETCDF

    a = str;
    dst_ids = src_ids = &a;
    cedrwd_();			/* position to absorb cedric data */

    strncpy(str, (char *)prefix, *np);
    str[*np] = '\0';
    b = a +strlen(a) -1;
    if(!strrchr(a, '/')) {
	/* directory name contains no slashes...punt!
	 */
	printf("Bad directory name: %s\n", a);
	return;
    }
    else if(*b != '/') {
	strcat(a, "/");
    }
    if (ta_makedir_recurse(a)) {
      perror(a);
      fprintf(stderr, "Cannot make output dir: %s\n", a);
      exit(1);
    }

    ced_netcdf(cedstrm, nfields, src_ids, dst_ids, a
	       , &cf, vol_num);
# endif

}
/* c------------------------------------------------------------------------ */

# ifdef NETCDF
int ced_netcdf(FILE *cedstrm, 
               int nfields, 
               char **src_ids, 
               char **dst_ids, 
               char *file_name_prefix,
               struct ced_fixes *cf,
               int vol_num)

{
    int ii, nn, fld, nflds, bv, dx, dy, dz, nb, rs, dims[4], level, nxy;
    int mark;
    float header_scale, alt, rcp_scale, f_lat, f_lon;
    static short *grid = NULL, *xgrid = NULL, *ygrid = NULL;
    static float *f_grid = NULL;
    double d, lat0, lat, lon0, lon;
    CoordInfo xi, yi, zi;
    char fname[256], str[256], *a, *b;
    char ced_names[MAXFLD][12], cdf_names[MAXFLD][12];
    char **sptr=src_ids, **dptr=dst_ids, *end_string();
    long start[4], count[4], file_offset, btime, zero=0;
    short lh[LVL_HRR_SIZE], *srs, *stop, *dst;
    float bad_val = BADVAL, *fgrd;
    struct tm t;
    struct ced_file_head cfh_struct, *cfh;
    size_t st;
    char tz[20];

    if( !mudbuf ) {
       mudbuf = (char *)malloc( VOL_HDR_SIZE * sizeof(short) );
# ifdef LITTLENDIAN
       mudcbuf = (char *)malloc( VOL_HDR_SIZE * sizeof(short) );
# endif
    }

    /* read in the cedric file header */
    cfh = &cfh_struct;
    rs = fseek(cedstrm, 0L, 0L);
    rs = fread((char *)cfh, (long)sizeof(struct ced_file_head), 1L, cedstrm);
# ifdef LITTLENDIAN
    st = swap32_( cfh->file_offset[vol_num-1] );
# else
    st = ( cfh->file_offset[vol_num-1] );
# endif
    rs = fseek(cedstrm, st, 0L);

    /*
     * read in the volume header
     */
    
# ifdef LITTLENDIAN
    mhc = (MudHdr *)mudcbuf;	/* for keeping character alignment */
    mh =  (MudHdr *)mudbuf;	/* swapped 16-bit data */
# else
    mh = mhc = (MudHdr *)mudbuf; /* big endian alignment */
# endif

    rs = fread((char *)mhc, (size_t)(VOL_HDR_SIZE*sizeof(short)), 1L, cedstrm);

# ifdef LITTLENDIAN
    swack_short( mhc, mh, (int)VOL_HDR_SIZE ); /* swap the bytes */
# endif

    /* get to Zeb time
     */
    t.tm_year = mh->mh_Year > 1900 ? mh->mh_Year-1900 : mh->mh_Year;
    t.tm_mon = mh->mh_Month - 1;
    t.tm_mday = mh->mh_Day;
    t.tm_hour = mh->mh_Hour;
    t.tm_min = mh->mh_Minute;
    t.tm_sec = mh->mh_Second;

# ifdef SYSV
# define KNOWN_OST
    strcpy (tz, "TZ=GMT");
    putenv (tz);
    timezone = 0;
    daylight = 0;
    t.tm_wday = t.tm_yday = 0;
    t.tm_isdst = -1;
    zt.zt_Sec = (long) mktime (&t);
# endif

#ifdef SVR4
# define KNOWN_OST
    strcpy (tz, "TZ=GMT");
    putenv (tz);
    timezone = 0;
    altzone = 0;
    daylight = 0;
    t.tm_wday = t.tm_yday = 0;
    t.tm_isdst = -1;
    zt.zt_Sec = (long) mktime (&t);
# endif

# ifdef linux
# define KNOWN_OST
    strcpy (tz, "TZ=GMT");
    putenv (tz);
    timezone = 0;
    daylight = 0;
    t.tm_wday = t.tm_yday = 0;
    t.tm_isdst = -1;
    zt.zt_Sec = (long) mktime (&t);
# endif

# ifndef KNOWN_OST
    t.tm_zone = (char *) 0;
    t.tm_wday = t.tm_isdst = t.tm_yday = 0;
    
    zt.zt_Sec = timegm (&t);
# endif

    zt.zt_MicroSec = 0;

    header_scale = mh->mh_Scale;
    
    /* get the origin */
    lat0 = mh->mh_OrLatDeg + mh->mh_OrLatMin/60.0 +
	  ((float)mh->mh_OrLatSec/(header_scale*3600.0));
    lon0 = mh->mh_OrLonDeg + mh->mh_OrLonMin/60.0 +
	  ((float)mh->mh_OrLonSec/(header_scale*3600.0));
    alt = mh->mh_ZMin*.001;
    alt = mh->mh_Fill39[0]*.001;

    cvt_Origin(lat0, lon0);

    /* get coordinate info */
    xi.ci_MinVal  = mh->mh_XMin/header_scale;
    xi.ci_MaxVal  = mh->mh_XMax/header_scale;
    xi.ci_Spacing = mh->mh_XSpacing/1000.0;
    xi.ci_NStep   = mh->mh_nX;
    yi.ci_MinVal  = mh->mh_YMin/header_scale;
    yi.ci_MaxVal  = mh->mh_YMax/header_scale;
    yi.ci_Spacing = mh->mh_YSpacing/1000.0;
    yi.ci_NStep   = mh->mh_nY;
    zi.ci_MinVal  = mh->mh_ZMin/(header_scale);
    zi.ci_MaxVal  = mh->mh_ZMax/(header_scale);
    zi.ci_Spacing = mh->mh_ZSpacing/1000.0;
    zi.ci_NStep   = mh->mh_nZ;

    cvt_ToLatLon ((double)xi.ci_MinVal, (double)yi.ci_MinVal, &lat, &lon);
    f_lat = lat; f_lon = lon;	/* back to float */
/*
 * Field info.
 */
    nflds = mh->mh_NField;
    nxy = xi.ci_NStep * yi.ci_NStep;
    bv = mh->mh_BadVal;
    for(fld=0; fld < nflds; fld++) {
	if(nfields) {
	    cdf_names[fld][0] = '\0';
	    /* relate the fields */
	    sptr = src_ids;
	    dptr = dst_ids;
	    for(ii=0; ii < nfields; ii++,sptr++,dptr++) {
		if(strstr(end_string(mhc->mh_Fields[fld],8,str), *sptr)) {
		    strcpy(cdf_names[fld], *dptr);
		}
	    }
	}
	else {
	    end_string(mhc->mh_Fields[fld].fi_Name, 8, cdf_names[fld]);
	}
    }

    /*
     * Set up the netCDF file
     */
    sprintf(fname, "%sddop.%02d%02d%02d.%02d%02d%02d.cdf"
	    , file_name_prefix
	    , t.tm_year, mh->mh_Month, mh->mh_Day
	    , mh->mh_Hour, mh->mh_Minute, mh->mh_Second);
    
    Nfile = nccreate (fname, NC_CLOBBER);
    fprintf(stderr, "Writing NetCDF file: %s\n", fname);

    /*
     * Make some dimensions.
     */
    DTime = ncdimdef (Nfile, "time", NC_UNLIMITED);
    dx = ncdimdef (Nfile, "x", xi.ci_NStep);
    dy = ncdimdef (Nfile, "y", yi.ci_NStep);
    dz = ncdimdef (Nfile, "z", zi.ci_NStep);

    /*
     * Make the variables.
     */
    dims[0] = DTime;
    dims[1] = dz;
    dims[2] = dy;
    dims[3] = dx;
    
    /*
     * Times
     */
    VBTime = ncvardef (Nfile, "base_time", NC_LONG, 0, 0);
    VTOff = ncvardef (Nfile, "time_offset", NC_FLOAT, 1, dims);
    /*
     * And positions.
     */
    VLat = ncvardef (Nfile, "lat", NC_FLOAT, 0, 0);
    VLon = ncvardef (Nfile, "lon", NC_FLOAT, 0, 0);
    VAlt = ncvardef (Nfile, "alt", NC_FLOAT, 0, 0);
    VXs = ncvardef (Nfile, "x_spacing", NC_FLOAT, 0, 0);
    VYs = ncvardef (Nfile, "y_spacing", NC_FLOAT, 0, 0);
    VZs = ncvardef (Nfile, "z_spacing", NC_FLOAT, 0, 0);
    /*
     * and fields
     */
    nn = nfields > 0 ? nfields : nflds;
    for (ii = 0; ii < nn; ii++)
	  {
	      VFields[ii] = ncvardef (Nfile, cdf_names[ii], NC_FLOAT, 4
				      , dims);
	      (void) ncattput (Nfile, VFields[ii], "missing_value",
			       NC_FLOAT, 1, &bad_val);
	  }
    ncendef (Nfile);

    ncvarput1 (Nfile, VLat, 0, &f_lat);
    ncvarput1 (Nfile, VLon, 0, &f_lon);
    ncvarput1 (Nfile, VAlt, 0, &alt);
    ncvarput1 (Nfile, VXs, 0, &xi.ci_Spacing);
    ncvarput1 (Nfile, VYs, 0, &yi.ci_Spacing);
    ncvarput1 (Nfile, VZs, 0, &zi.ci_Spacing);
    ncvarput1 (Nfile, VBTime, 0, &zt.zt_Sec);
    ncvarput1 (Nfile, VTOff, &zero, &zero);

    /* c...mark */
    if( !xgrid ) {
       grid = xgrid = (short *)malloc(nxy * sizeof(short));
    }

# ifdef LITTLENDIAN
    if( !ygrid ) {
       ygrid = (short *)malloc(nxy * sizeof(short));
    }
# endif
    if( !f_grid )
      { f_grid = (float *)malloc(nxy * sizeof(float)); }

    for(level=0; level < zi.ci_NStep; level++){
	/*
	 * read level header but ignore it
	 */
	fread((char *) lh, (size_t)(LVL_HRR_SIZE * sizeof(short))
	      , 1L, cedstrm);
	
	for(fld=0; fld < nflds; fld++) {
	    rcp_scale = 1.0/mh->mh_Fields[fld].fi_Scale;
# ifdef LITTLENDIAN
	    grid = xgrid;
# endif
	    /*
	     * read in the next field
	     */
	    mark = fread(grid, (long)(nxy * sizeof(short)), 1L, cedstrm);

# ifdef LITTLENDIAN
	    swack_short((char *)xgrid, (char *)ygrid, nxy );
	    grid = ygrid;
# endif

	    if(strlen(cdf_names[fld])) {
		for(srs=grid,stop=grid+nxy,fgrd=f_grid; srs < stop;
		    srs++, fgrd++) {
		    *fgrd = *srs == bv ? bad_val : *srs * rcp_scale;
		}
		/*
		 * Dump it into the netcdf file.
		 */
		start[0] = start[2] = start[3] = 0;
		start[1] = level;
		count[0] = 1;		/* time		*/
		count[1] = 1;		/* level	*/
		count[2] = yi.ci_NStep;
		count[3] = xi.ci_NStep;
		ncvarput (Nfile, VFields[fld], start, count, f_grid);
	    }
	}
    }
    ncclose (Nfile);
    free(grid);
    free(f_grid);
}
# endif
/* c------------------------------------------------------------------------ */

char *end_string(char *srs, int n, char *dst)

{
    /* take "n" characters of srs
     * remove trailing blanks
     * and terminate it
     */
    char *a=dst, *b=dst+n;

    if(n < 0)
	  return(0);

    for(; a < b;)
	  *a++ = *srs++;

    for(; b > dst && *(b-1) == ' '; b--);

    *b = '\0';

    return(dst);
}
/* c------------------------------------------------------------------------ */

int cedopn_(char *name,  int *n, int *isize)
{

    setScratchDir();

    char real_name[MAX_PATH_LEN];
    int ii, jj, kk, rslt;
    size_t sof, offs;
    char *a=name;

    name[*n] = '\0';
    if(cfhead == NULL) {
	cfhead = (struct ced_file_head *)malloc(sizeof(struct ced_file_head));
	memset((char *)cfhead, 0, sizeof(struct ced_file_head));
	
	cfi = (struct ced_file_info *)malloc(sizeof(struct ced_file_info));
	memset((char *)cfi, 0, sizeof(struct ced_file_info));

	if( !strchr(name, '/')) {	/* construct the name */
	    a = real_name;
            sprintf(a, "%s/", scratch_dir);
	    strcat( a, name );
	    strcat( a, ".ced" );
	}
    
	printf(" Real cedric file name is %s\n", a );
        ta_makedir_for_file(a);
	cedstrm = fopen( a, "w+" );
        if (cedstrm == NULL) {
          perror(a);
          fprintf(stderr, "ERROR - cannot open cedric file: %s\n", a);
          fprintf(stderr, "  Make sure directory exists\n");
          exit(1);
        }
	
	strncpy(cfhead->descriptor_name, "CED1", 4 );
	cfhead->byte_ordering = CED_BIG_ENDIAN;
	
	a = "This file was created by REORDER";
	strcpy(cfhead->say_something, a );
    }
# ifdef LITTLENDIAN
    sof = swap32_( cfhead->size_file );
# else
    sof = cfhead->size_file;
# endif
    /*
     * find the first zero offset
     */
    for(jj=0; jj < MAX_FILES && cfhead->file_offset[jj]; jj++);

    if(jj == MAX_FILES) {
	printf("Number of volumes has exceeded %d!\n", MAXFLD);
	exit(1);
    }
    if( jj ) {
# ifdef LITTLENDIAN
       offs = swap32_( cfhead->file_offset[jj -1] );
# else
       offs = cfhead->file_offset[jj -1];
# endif
    }
    sof += sizeof( struct ced_file_head );
    offs = sof;
    sof += *isize;

# ifdef LITTLENDIAN
    cfhead->file_offset[jj] = swap32_( offs );
    cfhead->size_file = swap32_( sof );
# else
    cfhead->file_offset[jj] = offs;
    cfhead->size_file = sof;
# endif
    ii = fseek(cedstrm, 0L, 0L);
    kk = sizeof( struct ced_file_head );
    ii = cedwrt_((char *)cfhead, &kk);
    ii = fseek(cedstrm, offs, 0L);
    return(ii);
}
/* c------------------------------------------------------------------------ */

int opnscr_(int *dummy1, int *lenrec, int *dummy2 )
{

    char pid[128];
    const char *scratchDir = NULL;

    sprintf(pid, "%d", getpid());

    setScratchDir();

    sprintf(scratch_name, "%s/reo-", scratch_dir);
    strcat(scratch_name, pid );
    scrstrm = fopen( scratch_name, "w+" );
    if (scrstrm == NULL) {
      perror(scratch_name);
      fprintf(stderr, "ERROR - cannot open scratch file: %s\n", scratch_name);
      fprintf(stderr, "  Make sure directory exists\n");
      exit(1);
    }
    rlens = *lenrec*4;
}
/* c------------------------------------------------------------------------ */

int cedcls_()
{
    fclose(cedstrm);
}
/* c------------------------------------------------------------------------ */

int cedloc_()
{
    /* return the offset to the current postion on the disk
     */
    int i=fseek(cedstrm, 0L, 1L);
    return(i);
}
/* c------------------------------------------------------------------------ */

int clsscr_()
{
    char command[99];

    if (scrstrm) {
      fclose(scrstrm);
    }
    sprintf( command, "rm " ); 
    strcat( command, scratch_name );
    system( command );
}
/* c------------------------------------------------------------------------ */

int cedrd_( int *buf, int *isize )
{
    /*
     * position to write logical record lrec of isize 16-bit words
     */
    int i=fread( (char *)buf, *isize, 1, cedstrm );
    return( i );
}
/* c------------------------------------------------------------------------ */

int cedwrt_( int *buf, int *isize )
{
    /*
     * write isize bytes
     */
    int i=fwrite( (char *)buf, *isize, 1, cedstrm );
    return( i );
}
/* c------------------------------------------------------------------------ */

int cedrwd_()
{
    /*
     * Position to the start of the last volume
     *
     * rewind the cedric file and read the file header
     * but not the volume header.
     */
    int i, j, k;
    size_t offs;

    if(!cfhead)
	  return(-1);
    
    for(i=0; i < MAX_FILES-1 && cfhead->file_offset[i]; i++);
# ifdef LITTLENDIAN
    offs = swap32_( cfhead->file_offset[i-1] );
# else
    offs = ( cfhead->file_offset[i-1] );
# endif
    
    j = fseek(cedstrm, offs, 0L);
    return( j );
}
/* c------------------------------------------------------------------------ */

int branwt_( int *lrec, int *buf, int *isize )

{
    /*
     * position to write logical record lrec of isize 16-bit words
     */
  if (scrstrm == NULL) {
    return 0;
  }
    int i = fseek( scrstrm, (long)((*lrec -1)*rlens), (int)0 );
    i = fwrite( (char *)buf, rlens, 1, scrstrm );
    return( i );
}
/* c------------------------------------------------------------------------ */

int branrd_( int *lrec, int *buf, int *isize )
{
    int i;
  if (scrstrm == NULL) {
    return 0;
  }
    /*
     * position to write logical record lrec of isize 16-bit words
     */
    i = ((*lrec -1)*rlens);
    i = fseek( scrstrm, (long)i, (int)0 );
    i = fread( (char *)buf, rlens, 1, scrstrm );
    return( i );
}
/* c------------------------------------------------------------------------ */

long TC_FccToSys (UItime *fcc)
/*
 * Convert an FCC time into a system time.
 */
{
	struct tm t;
	char tz[20];


	t.tm_year = fcc->ds_yymmdd/10000;
	t.tm_mon = (fcc->ds_yymmdd/100) % 100 - 1;
	t.tm_mday = fcc->ds_yymmdd % 100;
	t.tm_hour = fcc->ds_hhmmss/10000;
	t.tm_min = (fcc->ds_hhmmss/100) % 100;
	t.tm_sec = fcc->ds_hhmmss % 100;

# ifdef SYSV
# define KNOWN_OST
	strcpy (tz, "TZ=GMT");
	putenv (tz);
	timezone = 0;
	daylight = 0;
	t.tm_wday = t.tm_yday = 0;
	t.tm_isdst = -1;

	return( (long) mktime (&t) );
# endif
	
#ifdef SVR4
# define KNOWN_OST
	strcpy (tz, "TZ=GMT");
	putenv (tz);
	timezone = 0;
	altzone = 0;
	daylight = 0;
	t.tm_wday = t.tm_yday = 0;
	t.tm_isdst = -1;

	return( (long) mktime (&t) );
# endif
	
# ifdef linux
# define KNOWN_OST
	strcpy (tz, "TZ=GMT");
	putenv (tz);
	timezone = 0;
	daylight = 0;
	t.tm_wday = t.tm_yday = 0;
	t.tm_isdst = -1;

	return( (long) mktime (&t) );
# endif

# ifndef KNOWN_OST
	t.tm_zone = (char *) 0;
	t.tm_wday = t.tm_isdst = t.tm_yday = 0;
	
	return (timegm (&t));
#endif

}

/* c------------------------------------------------------------------------ */

long TC_ZtToSys (ZebTime *zt)
/*
 * Convert a zeb format time into a basic system format representation.
 */
{
	return (zt->zt_Sec);
}

/* c----------------------------------------------------------------------- */

static void cvt_Origin (double lat, double lon)

/*
 * Use lat,lon (deg) as the reference location for 
 * latitude,longitude <-> x,y conversions
 *
 * Return TRUE if we set the origin successfully, otherwise FALSE.
 */
{
/*
 * Store the values in radians
 */
	Origin_lat = lat * PI / 180.0;
	Origin_lon = lon * PI / 180.0;


	printf("cvt_Origin:\n");
	printf("\torigin lat: %8.4f\n",lat);
	printf("\torigin lon: %8.4f\n",lon);

}
/* c------------------------------------------------------------------------ */

static void cvt_ToLatLon (double x, double y, double *lat, double *lon)

/*
 * Convert x and y (km) to lat and lon (deg)
 */
{
    double range;
    double angleRadians;
    double alpha;
    double latRadians;
    double deltaPhi;
    double deltaLat;
    double deltaLon;
    double lonRadians;

    range = sqrt((x * x) + (y * y));
    angleRadians = atan2(y,x);

    alpha = range/r_earth;
    
    latRadians = Origin_lat + (alpha * cos(angleRadians));
    
    deltaPhi = log(tan((latRadians/2) + M_PI/4)/tan((Origin_lat/2) + M_PI/4));
    deltaLat = latRadians - Origin_lat;

    double q;
    if(deltaPhi == 0.0){
	q = cos(latRadians);
    }else{
	q = deltaLat / deltaPhi;
    }

    deltaLon = alpha * sin(angleRadians)/q;
    
    lonRadians = Origin_lon + deltaLon;
    
    *lat = latRadians * 180.0 / M_PI;
    *lon = lonRadians * 180.0 / M_PI;

    *lon = fmod(*lon,360.0);

    if(*lon < -180.0){
	*lon += 360.0;
    }

    if(*lon > 180.0){
	*lon -= 360.0;
    }

    printf("cvt_ToLatLon:\n");
    printf("\tx  : %8.4f\n",x);
    printf("\ty  : %8.4f\n",y);
    printf("\tlat: %8.4f\n",*lat);
    printf("\tlon: %8.4f\n",*lon);

}
/* c----------------------------------------------------------------------- */

static void cvt_ToXY (double lat, double lon, double *x, double *y)

/* 
 * Convert lat and lon (deg) to x and y (km) using azimuthal 
 * orthographic projection
 */
{
    double sin(), cos();
	double	del_lat, del_lon;
/*
 * Convert the lat,lon to x,y
 */
	lat *= PI / 180.0;
	lon *= PI / 180.0;

	del_lat = lat - Origin_lat;
	del_lon = lon - Origin_lon;

	*x = r_earth * cos ((double)lat) * sin (del_lon);
	*y = r_earth * sin (del_lat);
}
/* c------------------------------------------------------------------------ */

/* c----------------------------------------------------------------------- */

/* c------------------------------------------------------------------------ */

/* c----------------------------------------------------------------------- */

/*************************************************
 * Directory creation routine
 *
 * Nancy Rehak, Rap, NCAR, Boulder, CO, 80303, USA
 *
 * Updated by Mike Dixon, Feb 1999.
 *
 * March 1996
 */

// copy string safely

static char *STRncopy(char *s1, const char *s2, int maxs1)
{
  if (!s1 || !s2)
    return NULL;
  
  if (maxs1 > 0)
  {
    strncpy(s1, s2, (size_t) (maxs1-1));
    s1[maxs1-1] = '\0';
  }
  return(s1);
}

/********************************************************
 * makedir()
 *
 * Utility routine to create a directory.  If the directory
 * already exists, does nothing.
 *
 * Returns -1 on error, 0 otherwise.
 */

static int makedir(const char *path)
{
  return (ta_makedir(path));
}

/********************************************************
 * ta_makedir()
 *
 * Utility routine to create a directory.  If the directory
 * already exists, does nothing.
 *
 * Returns -1 on error, 0 otherwise.
 */

static int ta_makedir(const char *path)
{
  struct stat stat_buf;
  
  /*
   * Status the directory to see if it already exists.
   */

  if (stat(path, &stat_buf) == 0) {
    return(0);
  }
  
  /* this check seems to cause problems on Linux because errno is
   * not set correctly after stat - dixon
   *
   * if (errno != ENOENT) {
   *   return(-1); 
   * }
   */
  
  /*
   * Directory doesn't exist, create it.
   */

  if (mkdir(path,
	    S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
    /*
     * failed
     * check if dir has been made bu some other process
     * in the mean time, in which case return success
     */
    if (stat(path, &stat_buf) == 0) {
      return(0);
    }
    return(-1);
  }
  
  return(0);
}


/********************************************************
 * ta_makedir_recurse()
 *
 * Utility routine to create a directory recursively.
 * If the directory already exists, does nothing.
 * Otherwise it recurses through the path, making all
 * needed directories.
 *
 * Returns -1 on error, 0 otherwise.
 */

static int ta_makedir_recurse(const char *path)
{

  char up_dir[MAX_PATH_LEN];
  char *last_delim;
  struct stat dir_stat;
  int delim = PATH_DELIM[0];
  
  /*
   * Status the directory to see if it already exists.
   * '/' dir will always exist, so this stops the recursion
   * automatically.
   */
  
  if (stat(path, &dir_stat) == 0) {
    return(0);
  }
  
  /*
   * create up dir - one up the directory tree -
   * by searching for the previous delim and removing it
   * from the string.
   * If no delim, try to make the directory non-recursively.
   */
  
  STRncopy(up_dir, path, MAX_PATH_LEN);
  last_delim = strrchr(up_dir, delim);
  if (last_delim == NULL) {
    return (ta_makedir(up_dir));
  }
  *last_delim = '\0';
  
  /*
   * make the up dir
   */
  
  if (ta_makedir_recurse(up_dir)) {
    return (-1);
  }

  /*
   * make this dir
   */

  if (ta_makedir(path)) {
    return(-1);
  } else {
    return(0);
  }

}

/********************************************************
 * ta_makedir_for_file()
 *
 * Utility routine to create a directory recursively,
 * given a file path. The directory name is determined
 * by stripping the file name off the end of the file path.
 * If the directory already exists, does nothing.
 * Otherwise it recurses through the path, making all
 * needed directories.
 *
 * Returns -1 on error, 0 otherwise.
 */

static int ta_makedir_for_file(const char *file_path)
{
  
  char *delim, *lastDelim;
  int delimLen = strlen(PATH_DELIM);

  /*
   * get dir path from file path
   */

  delim = strstr(file_path, PATH_DELIM);
  if (delim == NULL) {
    /*
     * no delim, so file goes in current directory
     * no need to make dir
     */
    return 0;
  }

  lastDelim = delim; 
  while (delim != NULL) {
    delim = strstr(lastDelim + delimLen, PATH_DELIM);
    if (delim != NULL) {
      lastDelim = delim;
    }
  }

  /*
   * dir path extends to last delim
   */

  int dirPathLen = lastDelim - file_path;
  char *dir_path = (char *) malloc(dirPathLen + 1);
  strncpy(dir_path, file_path, dirPathLen);
  dir_path[dirPathLen] = '\0';

  if (ta_makedir_recurse(dir_path)) {
    free(dir_path);
    return -1;
  }

  free(dir_path);
  return 0;

}

