#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define GTYPE CONCAT2(TYPE,Grid)
#define FUNC(f) CONCAT3(LAL,TYPECODE,f)
#define CFUNCGRID FUNC(CreateGrid)
#define DFUNCGRID FUNC(DestroyGrid)
#define CFUNCARRAY FUNC(CreateArray)
#define DFUNCARRAY FUNC(DestroyArray)

void
CFUNCGRID(LALStatus * stat, GTYPE ** grid, UINT4Vector * dimLength,
          UINT4 dimension)
{
    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input arguments. */
    ASSERT(dimLength, stat, GRIDH_ENUL, GRIDH_MSGENUL);
    ASSERT(dimLength->data, stat, GRIDH_ENUL, GRIDH_MSGENUL);
    ASSERT(grid, stat, GRIDH_ENUL, GRIDH_MSGENUL);
    ASSERT(!*grid, stat, GRIDH_EOUT, GRIDH_MSGEOUT);
    ASSERT(dimension <= dimLength->length, stat, GRIDH_ENUL,
           GRIDH_MSGENUL);

    /* Allocate the grid. */
    if (!(*grid = (GTYPE *) LALMalloc(sizeof(GTYPE)))) {
        ABORT(stat, GRIDH_EMEM, GRIDH_MSGEMEM);
    }
    memset(*grid, 0, sizeof(GTYPE));

    /* Allocate the dimension units array. */
    if (!((*grid)->dimUnits = (LALUnit *)
          LALMalloc(dimension * sizeof(LALUnit)))) {
        LALFree(*grid);
        *grid = NULL;
        ABORT(stat, GRIDH_EMEM, GRIDH_MSGEMEM);
    }
    memset((*grid)->dimUnits, 0, dimension * sizeof(LALUnit));

    /* Allocate the offset and interval vectors. */
    LALDCreateVector(stat->statusPtr, &((*grid)->offset), dimension);
    BEGINFAIL(stat) {
        LALFree((*grid)->dimUnits);
        LALFree(*grid);
        *grid = NULL;
    }
    ENDFAIL(stat);
    LALDCreateVector(stat->statusPtr, &((*grid)->interval), dimension);
    BEGINFAIL(stat) {
        TRY(LALDDestroyVector(stat->statusPtr, &((*grid)->offset)), stat);
        LALFree((*grid)->dimUnits);
        LALFree(*grid);
        *grid = NULL;
    }
    ENDFAIL(stat);

    /* Allocate the data array. */
    CFUNCARRAY(stat->statusPtr, &((*grid)->data), dimLength);
    BEGINFAIL(stat) {
        TRY(LALDDestroyVector(stat->statusPtr, &((*grid)->interval)),
            stat);
        TRY(LALDDestroyVector(stat->statusPtr, &((*grid)->offset)), stat);
        LALFree((*grid)->dimUnits);
        LALFree(*grid);
        *grid = NULL;
    }
    ENDFAIL(stat);

    /* Done. */
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}


void DFUNCGRID(LALStatus * stat, GTYPE ** grid)
{
    INITSTATUS(stat);
    ATTATCHSTATUSPTR(stat);

    /* Check for valid input argument. */
    ASSERT(grid, stat, GRIDH_ENUL, GRIDH_MSGENUL);
    ASSERT(*grid, stat, GRIDH_ENUL, GRIDH_MSGENUL);

    /* Destroy the internal arrays and vectors. */
    TRY(DFUNCARRAY(stat->statusPtr, &((*grid)->data)), stat);
    TRY(LALDDestroyVector(stat->statusPtr, &((*grid)->interval)), stat);
    TRY(LALDDestroyVector(stat->statusPtr, &((*grid)->offset)), stat);
    LALFree((*grid)->dimUnits);

    /* Destroy the structure itself, and exit. */
    LALFree(*grid);
    *grid = NULL;
    DETATCHSTATUSPTR(stat);
    RETURN(stat);
}
