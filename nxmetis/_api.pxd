cdef extern from "src/include/parmetis.h":
    ctypedef Py_ssize_t idx_t
    ctypedef float real_t

    int ParMETIS_V3_PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm)
    int ParMETIS_V3_PartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
	     idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
	     idx_t *edgecut, idx_t *part, MPI_Comm *comm)
    int ParMETIS_V3_PartGeom(
             idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, MPI_Comm *comm)
    int ParMETIS_V3_RefineKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, 
	     idx_t *part, MPI_Comm *comm)
    int ParMETIS_V3_AdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, 
	     idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist, 
	     idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
    int ParMETIS_V3_Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, 
	     idx_t *ncommonnodes, idx_t **xadj, idx_t **adjncy, MPI_Comm *comm)
    int ParMETIS_V3_PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
	     idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm)
    int ParMETIS_V3_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
	     idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm)
    int ParMETIS_V32_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
             idx_t *numflag, idx_t *mtype, idx_t *rtype, idx_t *p_nseps, idx_t *s_nseps,
             real_t *ubfrac, idx_t *seed, idx_t *dbglvl, idx_t *order, 
             idx_t *sizes, MPI_Comm *comm)
    int ParMETIS_SerialNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
             idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm)


cdef int PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm) nogil:
    return ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag,
        numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm)


cdef int PartGeomKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ndims, real_t *xyz, 
	     idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
	     idx_t *edgecut, idx_t *part, MPI_Comm *comm) nogil:
    return ParMETIS_V3_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag,
        numflag, ndims, xyz, ncon, nparts, tpwgts, ubvec, options, edgecut,
        part, comm)


cdef int PartGeom(idx_t *vtxdist, idx_t *ndims, real_t *xyz, idx_t *part, MPI_Comm *comm) nogil:
    return ParMETIS_V3_PartGeom(vtxdist, ndims, xyz, part, comm)


cdef int RefineKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, 
	     idx_t *part, MPI_Comm *comm) nogil:
    return ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag,
        numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm)


cdef int AdaptiveRepart(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
	     idx_t *vsize, idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, 
	     idx_t *nparts, real_t *tpwgts, real_t *ubvec, real_t *ipc2redist, 
	     idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm) nogil:
    return ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt,
        wgtflag, numflag, ncon, nparts, tpwgts, ubvec, ipc2redist, options,
        edgecut, part, comm)


cdef int Mesh2Dual(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag, 
	     idx_t *ncommonnodes, idx_t **xadj, idx_t **adjncy, MPI_Comm *comm) nogil:
    return ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, numflag, ncommonnodes,
    xadj, adjncy, comm)


cdef int PartMeshKway(idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt, 
	     idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts, 
	     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, 
	     MPI_Comm *comm) nogil:
    return ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, elmwgt, wgtflag,
        numflag, ncon, ncommonnodes, nparts, tpwgts, ubvec, options, edgecut,
        part, comm)


cdef int NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
	     idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm) nogil:
    return ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, numflag, options, order,
        sizes, comm)


cdef int V32_NodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
             idx_t *numflag, idx_t *mtype, idx_t *rtype, idx_t *p_nseps, idx_t *s_nseps,
             real_t *ubfrac, idx_t *seed, idx_t *dbglvl, idx_t *order, 
             idx_t *sizes, MPI_Comm *comm) nogil:
    return ParMETIS_V32_NodeND(vtxdist, xadj, adjncy, vwgt, numflag, mtype,
        rtype, p_nseps, s_nseps, ubfrac, seed, dbglvl, order, sizes, comm)


cdef int SerialNodeND(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *numflag, 
             idx_t *options, idx_t *order, idx_t *sizes, MPI_Comm *comm) nogil:
    return ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order,
        sizes, comm)


    enum:
        METIS_VER_MAJOR
        METIS_VER_MINOR
        METIS_VER_SUBMINOR
        METIS_NOPTIONS

    enum rstatus_et:
        METIS_OK
        METIS_ERROR_INPUT
        METIS_ERROR_MEMORY
        METIS_ERROR

    enum moptype_et:
        METIS_OP_PMETIS
        METIS_OP_KMETIS
        METIS_OP_OMETIS

    enum moptions_et:
        METIS_OPTION_PTYPE
        METIS_OPTION_OBJTYPE
        METIS_OPTION_CTYPE
        METIS_OPTION_IPTYPE
        METIS_OPTION_RTYPE
        METIS_OPTION_DBGLVL
        METIS_OPTION_NITER
        METIS_OPTION_NCUTS
        METIS_OPTION_SEED
        METIS_OPTION_NO2HOP
        METIS_OPTION_MINCONN
        METIS_OPTION_CONTIG
        METIS_OPTION_COMPRESS
        METIS_OPTION_CCORDER
        METIS_OPTION_PFACTOR
        METIS_OPTION_NSEPS
        METIS_OPTION_UFACTOR
        METIS_OPTION_NUMBERING

    enum mptype_et:
        METIS_PTYPE_RB
        METIS_PTYPE_KWAY

    enum mgtype_et:
        METIS_GTYPE_DUAL
        METIS_GTYPE_NODAL

    enum mctype_et:
        METIS_CTYPE_RM
        METIS_CTYPE_SHEM

    enum miptype_et:
        METIS_IPTYPE_GROW
        METIS_IPTYPE_RANDOM
        METIS_IPTYPE_EDGE
        METIS_IPTYPE_NODE

    enum mrtype_et:
        METIS_RTYPE_FM
        METIS_RTYPE_GREEDY
        METIS_RTYPE_SEP2SIDED
        METIS_RTYPE_SEP1SIDED

    enum mdbglvl_et:
        METIS_DBG_INFO
        METIS_DBG_TIME
        METIS_DBG_COARSEN
        METIS_DBG_REFINE
        METIS_DBG_IPART
        METIS_DBG_MOVEINFO
        METIS_DBG_SEPINFO
        METIS_DBG_CONNINFO
        METIS_DBG_CONTIGINFO
        METIS_DBG_MEMORY

    enum mobjtype_et:
        METIS_OBJTYPE_CUT
        METIS_OBJTYPE_VOL
