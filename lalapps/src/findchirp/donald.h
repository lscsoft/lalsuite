#define DELM 0.0001
#define H1L1DISTANCE 0.018
#define COINWINDOW 0.011

#define MAXSTR       256

#define RESPONSEC_ENORM  0
#define RESPONSEC_ESUB   1
#define RESPONSEC_EARG   2
#define RESPONSEC_EVAL   3
#define RESPONSEC_EFILE  4
#define RESPONSEC_EINPUT 5
#define RESPONSEC_EMEM   6

#define RESPONSEC_MSGENORM  "Normal exit"
#define RESPONSEC_MSGESUB   "Subroutine failed"
#define RESPONSEC_MSGEARG   "Error parsing arguments"
#define RESPONSEC_MSGEVAL   "Input argument out of valid range"
#define RESPONSEC_MSGEFILE  "Could not open file"
#define RESPONSEC_MSGEINPUT "Error reading file"
#define RESPONSEC_MSGEMEM   "Out of memory"

/* Usage format string. */
#define USAGE "Usage: %s [options]\n"\
        "--help                           Print this help message\n" \
        "--ifofiles ifo1.cfg ifo2.cfg     Configurat\n"


