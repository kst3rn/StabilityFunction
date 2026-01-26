r"""
The component graph of a semistable curve
=========================================

Given a connected, semistable projective curve `X` over an algebraically closed
field `k`, we associate to `X` its *component graph* `G` as follows:

- the vertices of `G` correspond to the irreducible components of `X`;
- the edges of `G` correspond to the *nodes* of `X`, i.e. the singular points,
  which are, by assumption, ordinary double points. Such an edge connects the
  two vertices corresponding to the components intersecting in the node.

A component may intersect itself, and two components may intersect in several
nodes. Therefore `G` is an (undirected) graph with multiple edges and loops
allowed.

Vertex genera
-------------

In addition to the underlying graph we record, for each vertex `v` of `G`, a
nonnegative integer `g_v`, namely the geometric genus of the component
corresponding to `v` (i.e. the genus of its normalization).  The arithmetic
genus of `X` can then be recovered from the familiar formula

.. MATH::

    g(X) = \sum_v g_v + |E(G)| - |V(G)| + 1

for connected graphs.

Purely combinatorial design
---------------------------

This module is deliberately *purely combinatorial*.  An instance of
:class:`ComponentGraph` stores only:

- an abstract multigraph with loops,
- the genus attached to each vertex.

In particular, vertices are *not* modeled as curve objects.  This keeps the
class lightweight and reusable: the same data structure can represent component
graphs coming from many different sources (stable models, semistable models,
special fibres of arithmetic surfaces, etc.), independent of how the underlying
curve is presented.

One-tails and the core
----------------------

For applications to stable curves of small genus it is convenient to single out
*one-tails*, i.e. components of arithmetic genus one attached to the rest of the
curve by a single node.  Combinatorially, one-tails correspond to leaf vertices
of arithmetic genus one; there are two important cases:

- an *elliptic tail*: a leaf vertex of geometric genus one with no loop;
- a *pigtail*: a leaf vertex of geometric genus zero with one loop.

Removing all one-tails leaves the *core* of the curve; in the genus three
applications this typically reduces the size of the relevant graph drastically.

Isomorphism and canonical signatures
------------------------------------

Two semistable curves can have isomorphic component graphs even if their
components are numbered differently.  For this reason :class:`ComponentGraph`
provides a method :meth:`canonical_signature` which assigns to a component graph
a canonical invariant (a tuple of integers) built from the combinatorial data of
the core together with the distribution of one-tails.  Two component graphs are
isomorphic (as graphs with vertex genera and one-tail structure) if and only if
their canonical signatures agree.  This gives a simple and reliable way to test
isomorphism and to classify component graphs in small genus.


EXAMPLES::

    We construct a few small component graphs and compare their canonical
    signatures. Two component graphs are isomorphic (as weighted multigraphs
    with tails) if and only if their canonical signatures agree.

    We start with two graphs which are isomorphic.  The first one has a core
    consisting of two components of genera 1 and 0 meeting in three nodes
    (type ``1---0``), and the second one is obtained from the first by creating
    the vertices in a different order::

        sage: from component_graphs import ComponentGraph
        sage: G1 = ComponentGraph()
        sage: v1 = G1.add_vertex(1)      # core vertex of genus 1
        sage: w1 = G1.add_vertex(0)      # core vertex of genus 0
        sage: G1.add_edge(v1, w1); G1.add_edge(v1, w1); G1.add_edge(v1, w1)
        sage: G2 = ComponentGraph()
        sage: w2 = G2.add_vertex(0)      # same graph, but vertices created in opposite order
        sage: v2 = G2.add_vertex(1)
        sage: G2.add_edge(v2, w2); G2.add_edge(v2, w2); G2.add_edge(v2, w2)
        sage: G1.canonical_signature() == G2.canonical_signature()
        True

    Next we modify the first graph by attaching an elliptic tail to the genus-0
    component. The resulting graph is no longer isomorphic to the previous one::

        sage: H = ComponentGraph()
        sage: v = H.add_vertex(1)
        sage: w = H.add_vertex(0)
        sage: H.add_edge(v, w); H.add_edge(v, w); H.add_edge(v, w)
        sage: _ = H.add_elliptic_tail(w)
        sage: H.canonical_signature() == G1.canonical_signature()
        False

    Finally, we distinguish an elliptic tail from a pigtail.  Starting from a
    single core component of genus 1, we attach two one-tails, once as two
    elliptic tails and once as two pigtails.  These graphs are not isomorphic::

        sage: A = ComponentGraph()
        sage: cA = A.add_vertex(1)
        sage: _ = A.add_elliptic_tail(cA)
        sage: _ = A.add_elliptic_tail(cA)
        sage: B = ComponentGraph()
        sage: cB = B.add_vertex(1)
        sage: _ = B.add_pigtail(cB)
        sage: _ = B.add_pigtail(cB)
        sage: A.canonical_signature() == B.canonical_signature()
        False


"""

from sage.all import SageObject, Graph
from itertools import permutations


class ComponentGraph(SageObject):
    r""" The component graph of a semistable curve over a field.

    OUTPUT:

    An empty component graph.

    To this component graph one can add new vertices and edges. 

    """

    def __init__(self):
        # upon creation, the component graph is empty 
        self._graph = Graph(0, loops=True, multiedges=True)
        self._genus = {}
        # omitted for now:
        # self._tail_kind = {}    
        # for a vertex v, this should be "symbolic", "elliptic", "pigtail" or None

    def __repr__(self):
        return f"abstract component graph of a semistable projective curve"
    
    def genus_of_vertex(self, v):
        r"""Return the geometric genus attached to a vertex.

        INPUT:

        - ``v`` -- a vertex of the component graph

        OUTPUT:

        the nonnegative integer giving the geometric genus of the component
        corresponding to ``v``.
        """
        return self._genus[v]
    
    def genus(self):
        r""" Return the genus of this component graph.
        
        """
        if self.is_connected():
            return (sum(self.genus_of_vertex(v) for v in self.vertices()) 
                + self.number_of_edges() - self.number_of_vertices() + 1)
        else:
            raise NotImplementedError("the curve must be connected")
        
    def vertices(self):
        r""" Return the list of vertices of this component graph.

        """
        return self.graph().vertices()

    def number_of_vertices(self):
        r""" Return the number of vertices of this component graph.
        
        """
        return self.graph().num_verts()
           
    def edges(self):
        r""" Return the list of edges.
        
        """
        return self.graph().edges(labels=False)
    
    def number_of_edges(self):
        r""" Return the number of edges of this component graph.
        
        """
        return self.graph().num_edges()
    
    def graph(self):
        r""" Return the abstract graph of this component graph.
        
        """
        return self._graph
    
    def is_connected(self):
        r""" Return whether this component graph is connected.
        
        """
        return self.graph().is_connected()
    
    def add_vertex(self, genus):
        r""" Add a vertex to this component graph.
        
        INPUT:

        - ``genus`` -- a nonnegative integer

        OUTPUT:

        the vertex that has been created (a nonnegative integer)
        
        """
        v = self.graph().add_vertex()
        self._genus[v] = genus
        return v

    def add_edge(self, u, v):
        r"""Add an edge between two vertices.

        INPUT:

        - ``u, v`` -- vertices of the component graph

        This adds an (undirected) edge between ``u`` and ``v``.  Multiple edges
        and loops (when ``u == v``) are allowed.
        """
        self.graph().add_edge(u, v)

    def add_elliptic_tail(self, v):
        r""" Add an elliptic tail to this component graph.
        
        INPUT:

        - ``v`` -- a vertex of the graph
        
        OUTPUT:
        
        A new vertex of genus one is added to the graph, and connected to `v`
        by a unique edge. This is an *elliptic tail*.

        The newly created vertex is returned.
        """
        t = self.add_vertex(1)
        self.add_edge(v, t)
        return t

    def add_one_tail(self, v, tail_type):
        r""" Add a one-tail to this component.

        INPUT:

        - ``v`` -- a vertex of the graph
        - ``tail_type`` -- either "e" or "m"
        
        OUTPUT:

        A one-tail of type ``tail_type`` is attached to the vertex `v`,
        i.e. an elliptic tail if ``tail_type`` is "e", and a pigtail
        if ``tail_type`` is "m".

        The newly created vertex is returned.                       
        """
        if tail_type == "e":
            return self.add_elliptic_tail(v)
        elif tail_type == "m":
            return self.add_pigtail(v)
        else:
            ValueError()

    def add_pigtail(self, v):
        r""" Add a pigtail to this component graph.
        
        INPUT:

        - ``v`` -- a vertex of the graph
        
        OUTPUT:

        A new vertex of genus zero is added to the graph, connected to `v`
        by a unique edge, and a loop on the new vertex is added.
        This is a *pigtail*. 

        The newly created vertex is returned.       
        """
        t = self.add_vertex(0)
        self.add_edge(v, t)
        self.add_edge(t, t)   # loop
        return t
    
    def n_loops(self, v):
        r""" Return the number of loops on this vertex.
        
        INPUT:

        - ``v`` -- a vertex of the graph

        OUTPUT:

        the number of loops on `v`, i.e. edges of the form `(v,v)`.

        """
        return sum(1 for (a,b) in self.graph().edges(labels=False) if a == v and b == v)

    def nonloop_neighbors(self, v):
        r"""Return the neighbors of ``v`` excluding loops.

        INPUT:

        - ``v`` -- a vertex of the component graph

        OUTPUT:

        a list of vertices adjacent to ``v`` by a non-loop edge.
        """
        return [u for u in self.graph().neighbors(v) if u != v]

    def nonloop_degree(self, v):
        r"""Return the number of non-loop edges incident to ``v``.

        INPUT:

        - ``v`` -- a vertex of the component graph

        OUTPUT:

        the number of edges connecting ``v`` to *distinct* vertices.
        """
        return len(self.nonloop_neighbors(v))

    def edge_multiplicity(self, u, v, G=None):
        r"""Return the multiplicity of edges between two vertices.

        INPUT:

        - ``u, v`` -- vertices of the component graph
        - ``G`` -- an optional graph (default: the core graph)

        OUTPUT:

        the number of edges between ``u`` and ``v`` in ``G``.  If ``u == v``,
        this is the number of loops on ``u``.
        """
        H = self.core_graph() if G is None else G
        if u == v:
            return sum(1 for (a,b) in H.edges(labels=False) if a==u and b==u)
        return sum(1 for (a,b) in H.edges(labels=False)
                   if (a==u and b==v) or (a==v and b==u))
   
    def is_one_tail(self, v):
        r""" Return whether `v` corresponds to a one-tail.
        
        A vertex `v` of a component graph is a *one-tail* if it has arithmetic
        genus one, and is connected to the rest of the graph by a single edge. 

        There are two types of one-tails:

        - elliptic tails: the component is smooth and has geometric genus one,
        - pigtail: the component has geometric genus zero, and has one loop.

        """
        return self.is_elliptic_tail(v) or self.is_pigtail(v)
    
    def is_elliptic_tail(self, v):
        r""" Return whether `v` corresponds to an elliptic tail.
        
        A vertex `v` of a component graph is a *one-tail* if it has arithmetic
        genus one, and is connected to the rest of the graph by a single edge. 

        There are two types of one-tails:

        - elliptic tails: the component is smooth and has geometric genus one,
        - pigtail: the component has geometric genus zero, and has one loop.

        """
        return self.genus_of_vertex(v)==1 and self.n_loops(v)==0 and self.nonloop_degree(v)==1

    def is_pigtail(self, v):
        r""" Return whether `v` corresponds to a pigtail.
        
        A vertex `v` of a component graph is a *one-tail* if it has arithmetic
        genus one, and is connected to the rest of the graph by a single edge. 

        There are two types of one-tails:

        - elliptic tails: the component is smooth and has geometric genus one,
        - pigtail: the component has geometric genus zero, and has one loop.

        """
        return self.genus_of_vertex(v)==0 and self.n_loops(v)==1 and self.nonloop_degree(v)==1
    
    def n_elliptic_tails(self, v):
        r""" Return the number of elliptic tails attached to this vertex.
        
        INPUT:

        - `v` -- a vertex of the underlying graph, which is not a one-tail

        OUTPUT:

        the number of elliptic tails attached to the component corresponding to `v`.

        """
        assert not self.is_one_tail(v), "v must not be a one-tail"
        return sum(1 for t in self.graph()[v] if self.is_elliptic_tail(t))
    
    def n_pigtails(self, v):
        r""" Return the number of pigtails attached to this vertex.
        
        INPUT:

        - `v` -- a vertex of the underlying graph

        OUTPUT:

        the number of pigtails attached to the component corresponding to `v`.

        """
        assert not self.is_one_tail(v), "v must not be a one-tail"
        return sum(1 for t in self.graph()[v] if self.is_pigtail(t))

    def vertex_label(self, v):
        r""" Return the label of a vertex.
        
        INPUT:

        - ``v`` -- a vertex of the underlying graph

        OUTPUT:

        the *label* of `v`, which is a three-tuple `(g, e, p)` where

        - `g` is the geometric genus of the component corresponding to `v`,
        - `e` is the number of elliptic tails attached, and  
        - `p` is the number of pigtails attached.
        
        """
        return (self.genus_of_vertex(v), self.n_elliptic_tails(v), self.n_pigtails(v))
    
    def core_vertices(self):
        r""" Return the list of vertices belonging to the core.
        
        """
        return sorted(v for v in self.vertices() if not self.is_one_tail(v))

    def core_graph(self):
        r""" Return the subgraph corresponding to the core.
        
        """
        return self.graph().subgraph(self.core_vertices())
    
    def signature(self, sigma):
        r"""Return the signature of this component graph, with respect to a permutation.

        INPUT:

        - ``sigma`` -- a permutation of ``range(n)``, where ``n`` is the number
          of core vertices. (So ``sigma`` is a tuple/list of length ``n``.)

        OUTPUT:

        A tuple ``(L, A)`` where

        - ``L`` is the tuple of vertex labels of the core vertices in the order
          given by ``sigma``; each label is ``(g,e,p)``.
        - ``A`` is the tuple encoding the upper-triangular part of the adjacency
          matrix of the core graph in that same order, including the diagonal.
          The diagonal entries are the numbers of loops.

        This is designed so that ``canonical_signature()`` is the lexicographic
        minimum of these values over all permutations.
        """
        core = self.core_vertices()
        n = len(core)

        # normalize sigma to a tuple
        sigma = tuple(sigma)

        if len(sigma) != n or set(sigma) != set(range(n)):
            raise ValueError("sigma must be a permutation of range(n), where n=#core vertices")

        # core vertices in permuted order
        Vp = [core[i] for i in sigma]

        # L: ordered vertex labels (g,e,p)
        L = tuple(self.vertex_label(v) for v in Vp)

        # A: upper triangular adjacency, including diagonal (loops)
        H = self.core_graph()
        A = []
        for i in range(n):
            for j in range(i, n):
                A.append(self.edge_multiplicity(Vp[i], Vp[j], G=H))
        A = tuple(A)

        return (L, A)

    def canonical_signature(self):
        r"""Return the canonical signature of this component graph.

        OUTPUT:

        The canonical signature, defined as the lexicographically minimal value
        of ``self.signature(sigma)`` over all permutations ``sigma`` of the core
        vertices.

        Since in the genus 3 applications the number of core vertices is small,
        this brute-force canonicalization is practical and robust.
        """
        core = self.core_vertices()
        n = len(core)

        best = None
        for sigma in permutations(range(n)):
            s = self.signature(sigma)
            if best is None or s < best:
                best = s
        return best

    def as_text(self):
        r"""
        Return a human-readable description of this component graph.

        The output describes the combinatorial structure of the graph:
        vertices with their genera, loops, and attached one-tails, and the
        adjacency structure of the core graph.

        OUTPUT:

        A multiline string.

        EXAMPLES::

            sage: from component_graphs import ComponentGraph
            sage: G = ComponentGraph()
            sage: v = G.add_vertex(1)
            sage: w = G.add_vertex(0)
            sage: G.add_edge(v, w); G.add_edge(v, w); G.add_edge(v, w)
            sage: _ = G.add_elliptic_tail(w)
            sage: print(G.as_text())
            Component graph of a semistable curve
            ------------------------------------
            Connected: True
            Total genus: 3

            Vertices:
              v0: genus=1, loops=0, elliptic_tails=0, pigtails=0  [core]
              v1: genus=0, loops=0, elliptic_tails=1, pigtails=0  [core]
              v2: genus=1, loops=0                               [elliptic tail]

            Core adjacency (with multiplicities):
              v0 -- v1 : 3
        """
        lines = []
        lines.append("Component graph of a semistable curve")
        lines.append("-" * 36)
        lines.append(f"Connected: {self.is_connected()}")
        if self.is_connected():
            lines.append(f"Total genus: {self.genus()}")
        else:
            lines.append("Total genus: undefined (graph not connected)")
        lines.append("")
        lines.append("Vertices:")

        for v in self.vertices():
            g = self.genus_of_vertex(v)
            loops = self.n_loops(v)

            if self.is_one_tail(v):
                if self.is_elliptic_tail(v):
                    kind = "elliptic tail"
                elif self.is_pigtail(v):
                    kind = "pigtail"
                else:
                    kind = "one-tail"
                lines.append(
                    f"  v{v}: genus={g}, loops={loops}                               [{kind}]"
                )
            else:
                e = self.n_elliptic_tails(v)
                p = self.n_pigtails(v)
                lines.append(
                    f"  v{v}: genus={g}, loops={loops}, elliptic_tails={e}, pigtails={p}  [core]"
                )

        core = self.core_vertices()
        if len(core) > 0:
            lines.append("")
            lines.append("Core adjacency (with multiplicities):")
            H = self.core_graph()
            for i, u in enumerate(core):
                for v in core[i:]:
                    m = self.edge_multiplicity(u, v, G=H)
                    if m == 0:
                        continue
                    if u == v:
                        lines.append(f"  v{u} -- v{u} : {m} (loops)")
                    else:
                        lines.append(f"  v{u} -- v{v} : {m}")
        else:
            lines.append("")
            lines.append("Core adjacency: (empty)")

        return "\n".join(lines)
