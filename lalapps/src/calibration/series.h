typedef enum { Time, Freq } domain;

#define IS_TIME( domain_ ) ( domain_ == Time )
#define IS_FREQ( domain_ ) ( domain_ == Freq )

typedef struct { int sec; int nan; } epoch;

struct series
{
  const char *name;
  epoch       tbeg;
  epoch       tend;
  domain      type;
  float       step;
  const char *unit;
  size_t      size;
  float      *data;
};

double epoch_diff( const epoch *t2, const epoch *t1 );
void epoch_add( epoch *t1, epoch *t0, double dt );
int write_ilwd( const char *fname, const struct series *ser );
struct FrameH *fr_add_proc_data( struct FrameH *frame, const struct series *ser );
