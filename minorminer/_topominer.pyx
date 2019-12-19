# distutils: language = c++
# cython: language_level=2
"""
Topological `minorminer` or `topominer` is a topology-aware minor embedding tool.
`topominer` adds a global placement step to `minorminer` by parsing 2D location
information of the problem graph and assigning initial_chains to them.

This is an implementation of the work of Pinilla and Wilton in [2].

Detailed description of `minorminer` can be found in `_minorminer.pyx`

[1] https://arxiv.org/abs/1406.2741
[2] https://link.springer.com/chapter/10.1007/978-3-030-20656-7_7
"""
include "_minorminer_h.pxi"
import os

def topo_embedding(S, T, source_layout, target_layout, **params):
    """
    topo_embedding(S, T, source_layout, target_layout, **params)

    Args::

        S: an iterable of label pairs representing the edges in the source graph, or a NetworkX Graph

        T: an iterable of label pairs representing the edges in the target graph, or a NetworkX Graph

        source_layout: Dictionary of 2D locations of nodes in source
            graph. Values should be len-2 numpy.ndarray x,y = [-1.0, 1.0]. If
            not given topo_embedding creates an artificial layout.

        target_layout: Dictionary of 2D locations of nodes in target
            graph. Values should be len-2 numpy.ndarray x,y = [-1.0, 1.0]. If
            not given topo_embedding creates an artificial layout.

        **params (optional): see 'minorminer'

    Returns::

        When return_overlap = False (the default), returns a dict that maps labels in S to lists of labels in T.
            If the heuristic fails to find an embedding, an empty dictionary is returned

        When return_overlap = True, returns a tuple consisting of a dict that maps labels in S to lists of
            labels in T and a bool indicating whether or not a valid embedding was found

        When interrupted by Ctrl-C, returns the best embedding found so far

        Note that failure to return an embedding does not prove that no embedding exists

    """
    cdef _input_parser _in
    try:
        _in = _input_parser(S, T, source_layout, target_layout, params)
    except EmptySourceGraphError:
        return {}

    cdef chainmap candidates
    findCandidates(_in.Sg, _in.Tg, _in.Sloc, _in.Tloc, candidates)
    _in.opts.initial_chains =  candidates

    cdef vector[int] chain
    cdef vector[vector[int]] chains
    cdef int success = findEmbedding(_in.Sg, _in.Tg, _in.opts, chains)

    cdef int nc = chains.size()

    rchain = {}
    if chains.size():
        for v in range(nc-_in.pincount):
            chain = chains[v]
            rchain[_in.SL.label(v)] = [_in.TL.label(z) for z in chain]

    if _in.opts.return_overlap:
        return rchain, success
    else:
        return rchain

class EmptySourceGraphError(RuntimeError):
    pass

cdef class _input_parser:
    cdef input_graph Sg, Tg
    cdef labeldict SL, TL
    cdef locmap Sloc, Tloc
    cdef optional_parameters opts
    cdef int pincount
    def __init__(self, S, T, source_layout, target_layout, params):
        cdef uint64_t *seed
        cdef object z

        self.opts.localInteractionPtr.reset(new LocalInteractionPython())

        names = {"max_no_improvement", "random_seed", "timeout", "max_beta",
                 "tries", "inner_rounds", "chainlength_patience", "max_fill",
                 "threads", "return_overlap", "skip_initialization", "verbose"}

        for name in params:
            if name not in names:
                raise ValueError("%s is not a valid parameter for find_embedding"%name)

        z = params.get("max_no_improvement")
        if z is not None:
            self.opts.max_no_improvement = int(z)

        z = params.get("skip_initialization")
        if z is not None:
            self.opts.skip_initialization = int(z)

        z = params.get("chainlength_patience")
        if z is not None:
            self.opts.chainlength_patience = int(z)

        z = params.get("random_seed")
        if z is not None:
            self.opts.seed( long(z) )
        else:
            seed_obj = os.urandom(sizeof(uint64_t))
            seed = <uint64_t *>(<void *>(<uint8_t *>(seed_obj)))
            self.opts.seed(seed[0])

        z = params.get("tries")
        if z is not None:
            self.opts.tries = int(z)

        z = params.get("verbose")
        if z is not None:
            self.opts.verbose = int(z)

        z = params.get("inner_rounds")
        if z is not None:
            self.opts.inner_rounds = int(z)

        z = params.get("timeout")
        if z is not None:
            self.opts.timeout = float(z)

        z = params.get("max_beta")
        if z is not None:
            self.opts.max_beta = float(z)

        z = params.get("return_overlap")
        if z is not None:
            self.opts.return_overlap = int(z)

        z = params.get("max_fill")
        if z is not None:
            self.opts.max_fill = int(z)

        z = params.get("threads")
        if z is not None:
            self.opts.threads = int(z)

        self.SL = _read_graph(self.Sg, S)
        if not self.SL:
            raise EmptySourceGraphError

        self.TL = _read_graph(self.Tg, T)
        if not self.TL:
            raise ValueError("Cannot embed a non-empty source graph into an empty target graph.")

        _get_locmap(source_layout, self.Sloc, self.SL)
        _get_locmap(target_layout, self.Tloc, self.TL)

cdef int _get_chainmap(C, chainmap &CMap, SL, TL, parameter) except -1:
    cdef vector[int] chain
    CMap.clear();
    try:
        for a in C:
            chain.clear()
            if C[a]:
                for x in C[a]:
                    if x in TL:
                        chain.push_back(<int> TL[x])
                    else:
                        raise RuntimeError, "%s uses target node labels that weren't referred to by any edges"%parameter
                if a in SL:
                    CMap.insert(pair[int,vector[int]](SL[a], chain))
                else:
                    raise RuntimeError, "%s uses source node labels that weren't referred to by any edges"%parameter

    except (TypeError, ValueError):
        try:
            nc = next(C)
        except:
            nc = None
        if nc is None:
            raise ValueError("initial_chains and fixed_chains must be mappings (dict-like) from ints to iterables of ints; C has type %s and I can't iterate over it"%type(C))
        else:
            raise ValueError("initial_chains and fixed_chains must be mappings (dict-like) from ints to iterables of ints; C has type %s and next(C) has type %s"%(type(C), type(nc)))

cdef int _get_locmap(loc, locmap &LMap, L):
    LMap.clear();
    try:
        for x in L:
            LMap.insert(pair[int,intpair](L[x],loc[x]))
    except KeyError as e:
        raise RuntimeError("All nodes in graph should have locations.")

cdef _read_graph(input_graph &g, E):
    cdef labeldict L = labeldict()
    if hasattr(E, 'edges'):
        E = E.edges()
    for a,b in E:
        g.push_back(L[a],L[b])
    return L

__all__ = ["topo_embedding"]
