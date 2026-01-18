r"""
The component graph of a semistable curve
=========================================

Let `X` be a projective plane curve over a perfect field `k`. We say that `X` is 
*semistable* if it is reduced and has at most ordinary double points (nodes) as 
singularities.

To such a semistable curve we can associate a *component graph* `G`, as follows. 
We first assume that all irreducible components are absolutely irreducible and all 
singular points are *split* ordinary double points (this can be achieved with
a finite extension of the base field `k`). Then the vertices of `G` correspond to
irreducible components of `X` and the edges correspond to the nodes. 

We also consider the following generalization, which is relevant for our paper
`Semistable reduction of plane quartics`. Assume that `X` is a projective plane
curve over a field `k` which is reduced and has at most *nodes* and *cusps*
as singularities. Then we can attach to `X` a component graph as before, where
we attach for each cusp formally a *one-tail* as a leaf of the component graph.  

In this module we realize the computation of the component graph of a semistable
curve `X` (possibly with cusps).


EXAMPLES:

We compute the component graph of a triangle:

    sage: k = GF(2)
    sage: R.<x,y,z> = k[]
    sage: X = Curve(x*y*z)
    sage: G = component_graph(X); G
    component graph of a projective curve over Finite Field of size 2

    sage: G.genus()
    1

    sage: print(G.as_text())
    Component graph over Finite Field of size 2
    Vertices: 3, Edges: 3
    v0: g=0, component: z
    v1: g=0, component: y
    v2: g=0, component: x
    Adjacency:
    v0: v1, v2
    v1: v0, v2
    v2: v0, v1

The example from § 5.1 of our paper. This is a reduced and integral quartic over `\mathbb{F}_2`
with one loop and two cusps (defined over `\mathbb{F}_4$).

    sage: X = Curve(x^4 + x*y^2*z + x^2*z^2 + z^4)
    sage: G = component_graph(X)
    sage: G.genus()
    3

    sage: print(G.as_text())
    Component graph over Finite Field in z2 of size 2^2
    Vertices: 3, Edges: 3
    v0: g=0, component: Projective Plane Curve over Finite Field in z2 of size 2^2 defined by x^4 + x*y^2*z + x^2*z^2 + z^4
    v1: g=1, component: symbolic tail
    v2: g=1, component: symbolic tail
    Adjacency:
    v0: v0, v1, v2
    v1: v0
    v2: v0


"""

from semistable_model.curves.plane_curves import ProjectivePlaneCurve
from sage.all import SageObject, Graph, Curve, Set, lcm, GF
from itertools import product


def component_graph(X, allow_cusps=True):
    r""" Return the component graph of a projective plane curve.
    
    INPUT:

    - `X` -- a projective plane curve over a field
    - ``allow_cusps`` -- a boolean

    OUTPUT:

    the component graph of `X`.

    The curve `X` 
    - must be defined over a finite field, 
    - must be reduced, 
    - if ``allow_cusps`` is False, `X` must be semistable, and
    - if it is ``True``, nodes and cusps are also allowed as singularities

    
    EXAMPLES:

        sage: from semistable_model.curves.component_graphs import component_graph
        sage: k = GF(2)
        sage: R.<x,y,z> = k[]
        sage: X = Curve(x*y*z)
        sage: G = component_graph(X); G
        component graph of a projective curve over Finite Field of size 2

        sage: G.vertices()
        {0: Projective Plane Curve over Finite Field of size 2 defined by z,
        1: Projective Plane Curve over Finite Field of size 2 defined by y,
        2: Projective Plane Curve over Finite Field of size 2 defined by x}

        sage: G.genus()
        1

        sage: k = GF(3)
        sage: R.<x,y,z> = k[]
        sage: X = Curve(y^2*z+x^3+x^2*z)
        sage: G = component_graph(X); G
        component graph of a projective curve over Finite Field in z2 of size 3^2

    Note that the base field was extended to have a split node at (0:0:1).

        sage: G.vertices()
        {0: Projective Plane Curve over Finite Field in z2 of size 3^2 defined by x^3 + x^2*z + y^2*z}
    
        sage: G.edges()
        [[(0,
        Place (y + 2*z2 + 2) on Projective Plane Curve over Finite Field in z2 of size 3^2 defined by x^3 + x^2*z + y^2*z),
        (0,
        Place (y + z2 + 1) on Projective Plane Curve over Finite Field in z2 of size 3^2 defined by x^3 + x^2*z + y^2*z)]]

        sage: G.genus()
        1

        sage: X = Curve((x^2+y^2-z^2)*((x-z)^2+y^2+z^2))
        sage: G = component_graph(X); G
        component graph of a projective curve over Finite Field in z2 of size 3^2

        sage: G.genus()
        3

    If the curve is not semistable, but has cusps, then each cusp is represented by
    a one-tail:

        sage: X = Curve(x^3 - y^2*z)
        sage: G = component_graph(X); G
        component graph of a projective curve over Finite Field of size 3

        sage: print(G.as_text())
        Component graph over Finite Field of size 3
        Vertices: 2, Edges: 1
          v0: g=0, component: Projective Plane Curve over Finite Field of size 3 defined by x^3 - y^2*z
          v1: g=1, component: symbolic tail
        Adjacency:
          v0: v1
          v1: v0

    """
    # we want X to be an object of the native Curve class
    if isinstance(X, ProjectivePlaneCurve):
        X = Curve(X.defining_polynomial())
    assert X.defining_polynomial().is_squarefree(), "X must be reduced"
    k = X.base_ring()
    assert k.is_field() and k.is_finite(), "X has to be defined over a finite field"
    k1  = splitting_field(X)
    if not k1 == k:
        X = Curve(X.base_extend(k1))
        k = k1

    components = [Curve(Y) for Y in X.irreducible_components()]
    branches = singular_branches(X, components)
    singular_points = Set(b.rational_point() for b in branches)
    nodes = []
    cusps = []
    for P in singular_points:
        # test whether P is a node; there is a problem here to be solved!
        if is_node(X, P):
            nodes.append(P)
        elif allow_cusps and is_cusp(X, P):
            cusps.append(P)
        else:
            raise ValueError("X is not semistable")

    G = ComponentGraph(X.base_ring())
    for Y in components:
        G.add_vertex(CurveComponent(Y))
    for P in nodes:
        # there must be exactly two branches passing through P 
        e = [b for b in branches if b.rational_point() == P]
        assert len(e) == 2
        G.add_node(e)
    if allow_cusps:
        for P in cusps:
            b1 = [b for b in branches if b.rational_point() == P][0]
            tail_component = TailComponent()
            G.add_tail(b1, tail_component)
    return G


def splitting_field(X):
    r""" Return the splitting field of a plane curve.
    
    INPUT:

    - `X` -- a plane projective curve 

    OUTPUT:

    a finite extension of `k` over which 

    - all irreducible components
    - all singular points, and
    - all branches at singular points 
    
    are rational.

    
    EXAMPLES:

    sage: from semistable_model.curves.component_graphs import splitting_field
    sage: k = GF(7)
    sage: R.<x,y,z> = k[]
    sage: F = y^2*z + x^3 + x^2*z
    sage: X = Curve(F)
    sage: splitting_field(X)
    Finite Field in z2 of size 7^2

    """
    S = singular_branches(X)
    p = X.base_ring().characteristic()
    n = lcm(s.residue_field().degree() for s in S)
    return GF(p**n)

    
def singular_branches(X, components=None):
    r""" Return the list of all singular branches of a projective plane curve.

    INPUT:

    - ``X`` -- a projective plane curve over a field
    - ``components`` -- the list of irreducible components of `X`, ``None``

    OUTPUT:

    the list of all branches of `X` passing through the singular points.


    EXAMPLES:

    sage: from semistable_model.curves.component_graphs import singular_branches
    sage: k = GF(7)
    sage: R.<x,y,z> = k[]
    sage: F = x^2 + y^2
    sage: X = Curve(F)
    sage: singular_branches(X)
    [Place (1/z) on Projective Plane Curve over Finite Field of size 7 defined by x^2 + y^2]

    """
    if components is None:
        components = [Curve(Y) for Y in X.irreducible_components()]
    S = X.singular_subscheme().reduce()
    I = S.defining_ideal()
    R = I.ring()
    ret = []
    for Y in components:
        RY = Y.coordinate_ring()
        coordinates = []
        for x in R.gens():
            xb = RY(x)
            if not xb.is_zero() and not xb.is_unit():
                coordinates.append(x)
        places = Set()
        for f, x in product(I.gens(), coordinates):
            g = Y.function(f/x**f.degree())
            if not g.is_zero():
                places = places.union(Set(g.zeros()))
        ret += [CurveBranch(Y, v) for v in places 
                if all(f in Y.place_to_closed_point(v).prime_ideal() for f in I.gens())] 
    return ret

def is_node(X, P):
    r""" Return whether `P` is an ordinary double point of `X`.
    
    INPUT:

    - `X` -- a projective plane curve over a field
    - `P` -- a rational point on `X`

    OUTPUT:

    whether `P` is an ordinary double point on `X`, i.e.
    whether the tangent cone is a reduced conic (i.e.
    two distinct lines, possibley after a quadratic extension
    of th ground field).
    
    EXAMPLES:

        sage: k = GF(2)
        sage: R.<x,y,z> = k[]
        sage: X = Curve(x*y*z)
        sage: P = X.point([1,0,0])
        sage: is_node(X, P)
        True

    Note that Sage's native method :meth:`is_ordinary_singularity` gives a false result:

        sage: X.is_ordinary_singularity(P)
        False

        sage: P = X.point([1,1,0])
        sage: is_node(X, P)
        False

    A nonsplit ordinary double point counts as a node, too:

        sage: X = Curve(x^2 + x*y + y^2)
        sage: P = X.point([0,0,1])
        sage: is_node(X, P)
        True

    A cusp is not an ordinary double point:

        sage: k = GF(3)
        sage: R.<x,y,z> = k[]
        sage: X = Curve(x^3 + y^2*z)
        sage: P = X.point([0,0,1])
        sage: is_node(X, P)
        False


    """
    assert P in X, "P must be a point on X"
    x, y, z = X.ambient_space().coordinate_ring().gens()
    F = X.defining_polynomial()(x-P[0], y-P[1], z-P[2]).homogeneous_components()
    return not 1 in F.keys() and 2 in F.keys() and F[2].is_squarefree()


def is_cusp(X, P):
    r""" Return whether `P` is a cusp of `X`.
    
    INPUT:

    - ``X`` -- a projective plane curve over a field
    - ``P`` -- a rational point on `X`

    OUTPUT:

    whether `P` is a cusp, i.e. a singularity of type `A_2`.


    EXAMPLES:

        sage: k = GF(2)
        sage: R.<x,y,z> = k[]
        sage: X = Curve(y^2*z + x^3)
        sage: P = X.point([0,0,1])
        sage: is_cusp(X, P)
        True

    """ 
    assert P in X, "P must be a point on X"
    x, y, z = X.ambient_space().coordinate_ring().gens()
    F = X.defining_polynomial()(x-P[0], y-P[1], z-P[2]).homogeneous_components()
    if 1 in F.keys() or not 2 in F.keys():
        # P is not a double point:
        return False
    F2 = F[2]
    if F2.is_squarefree():
        # P is an ordinary double point
        return False
    # we must test whether P is really a cusp:
    if not 3 in F.keys():
        return False
    L = F2.factor()[0][0]
    return not F[3] in L.parent().ideal(L)


# -----------------------------------------------------------------------------

#                      classes


class CurveBranch(SageObject):
    r""" Return a curve branch.

    INPUT:

    - ``Y`` -- an integral curve over a field
    - ``v`` -- a place of the function of `Y`
    
    OUTPUT:

    an object representing the pair `(Y, v)`.

    """
    def __init__(self, Y, v):
        self._curve = Y
        self._place = v
        self._point = Y.place_to_closed_point(v)

    def __repr__(self):
        return f"{self.place()} on {self.curve()}"

    def curve(self):
        return self._curve
    
    def place(self):
        return self._place
    
    def point(self):
        r""" Return the point on the curve corresponding to this branch.
        
        Note that what is returned is a closed point of the curve. You cannot
        use it directly to compare it with a point on another curve. 
        """
        return self._point
    
    def rational_point(self):
        r""" Return the corresponding rational point of the projective plane.
        
        If the point is not rational, an error is raised.
        """
        assert self.point().degree() == 1, "the center of this point is not rational"
        return self.curve().ambient_space().point(self.point().rational_point())
    
    def residue_field(self):
        return self.place().residue_field()[0]


# -----------------------------------------------------------------------------


class Component(SageObject):
    r""" Return a component object which can be used as vertex of a component graph.
    
    This is an abstract base class. Subclasses must, on initilization, define at
    least the following secret attributes:

    - self._curve : an integral curve, or ``None``
    - self._genus : the (geometric) genus of the component
    - self._label : a string giving a short description of the component
    - self._kind  : a string revealing the subclass (e.g. "curve" or "tail")
    
    """
    def curve(self):
        return self._curve

    def genus(self):
        r""" Return the genus of this component.
        
        """
        return self._genus

    def label(self):
        r""" Return the label of this component.
        
        This is a string which gives a short description of the component.

        """
        return self._label

    def kind(self):
        r""" Return the kind of this ccomponent.
        
        The *kind* must be a short string reavealing the subclass this object
        belongs to, e.g. "curve" or "tail"

        """
        return self._kind


class CurveComponent(Component):
    r""" Return the component corresponding to to an integral projective curve over a field.
    
    INPUT:

    - ``curve`` -- an integral projective curve over a field

    At the moment, only plane curves are allowed.
    
    """
    def __init__(self, curve):
        # we have to make sure that ``curve`` lives in the right class
        curve = Curve(curve)
        assert curve.is_irreducible(), "the curve must be irreducible"
        self._curve = curve
        self._genus = curve.geometric_genus()
        # this only works for plane curves:
        self._label = str(curve.defining_polynomial())
        self._kind = "curve"

    def __repr__(self):
        return str(self.curve())

    def curve(self):
        return self._curve

    def coordinate_ring(self):
        return self.curve().coordinate_ring()


class TailComponent(Component):
    def __init__(self, tail_kind="symbolic"):
        self._curve = None
        self._tail_kind = tail_kind   # "elliptic","pig","symbolic"
        self._genus = 1
        self._label = f"{tail_kind} tail"
        self._kind = "tail"

    def __repr__(self):
        return self.label()


# -----------------------------------------------------------------------------


class ComponentGraph(SageObject):
    r""" The component graph of a semistable curve over a field.

    INPUT:

    - ``base_field`` -- a finite field

    OUTPUT:

    An empty component graph over `k`.

    To this component graph one can add new vertices and edges. 
    
    
    """

    def __init__(self, base_field):
        assert base_field.is_field() and base_field.is_finite(), "the base field must be a finite field"
        self._base_field = base_field
        # upon creation, the component graph is empty 
        self._vertices = {}
        self._edges = []
        self._graph = Graph(0, loops=True, multiedges=True)

    def __repr__(self):
        return f"component graph of a projective curve over {self.base_field()}"
    
    def base_field(self):
        r""" Return the base field of this component graph.
        
        """
        return self._base_field
    
    def genus(self):
        r""" Return the genus of this component graph.
        
        """
        if self.is_connected():
            return (sum(Y.genus() for Y in self.components()) 
                + self.number_of_edges() - self.number_of_vertices() + 1)
        else:
            raise NotImplementedError("the curve must be connected")
        
    def vertices(self):
        r""" Return the dictionary of vertices of this component graph.

        The keys of this dictionary are the names of the vertices of the
        underlying bare graph. The correspondig value, the *vertex*, is 
        an integral curve on the ambient projective plane.

        """
        return self._vertices

    def number_of_vertices(self):
        r""" Return the number of vertices of this component graph.
        
        """
        return len(self.vertices())
    
    def components(self):
        r""" Return the list of irreducible components.
        
        """
        return list(self.vertices().values())
        
    def edges(self):
        r""" Return the list of edges.
        
        An *edge* is a pair of distinct branches with the same center,
        an ordinary double point of the curve.

        """
        return self._edges
    
    def number_of_edges(self):
        r""" Return the number of edges of this component graph.
        
        """
        return len(self.edges())
    
    def graph(self):
        r""" Return the abstract graph of this component graph.
        
        """
        return self._graph
    
    def is_connected(self):
        r""" Return whether this component graph is connected.
        
        """
        return self.graph().is_connected()
        
    def add_vertex(self, Y):
        r""" Add a vertex to this component graph.

        INPUT:

        - ``Y`` -- a component object


        OUTPUT:

        The method adds a new vertex to the underlying graph,
        corresponding to the component `Y`. The newly created index
        for this vertex is returned.
        
        """
        index = self.graph().add_vertex()
        self._vertices[index] = Y
        return index
    
    def add_node(self, node_branches):
        r""" Add an edge to this component graph corresponding to a node of the curve.

        INPUT:

        ``node_branches`` -- a pair of distinct branches with the same center, corresponding
                             to a node of a plane curve

        OUTPUT:

        This function adds an edge to this component graph which corresponds to
        the node represented by ``node_branches``. Nothing is returned. 

        Note that both branches need to be defined on the same ambient space
        (for the moment: a projective plane over a field) to make sense of the
        condition that both have the same center. 

        """
        b1 = node_branches[0]
        b2 = node_branches[1]
        assert b1.rational_point() == b2.rational_point()
        a1 = [a for a, Y in self.vertices().items() if Y.curve() == b1.curve()][0]
        a2 = [a for a, Y in self.vertices().items() if Y.curve() == b2.curve()][0]
        self._edges.append([(a1, b1), (a2, b2)])
        self.graph().add_edge(a1, a2)

    def add_tail(self, cusp_branch, tail_component):
        r""" Add a tail to this component graph corresponding to a cusp of the curve.
        
        INPUT:

        - ``cusp_branch`` -- a curve branch
        - ``tail_component`` -- a tail component

        
        OUTPUT:

        This function adds a tail to the component graph, i.e. a leaf of the
        underlying graph, 

        """
        a1 = [a for a, Y in self.vertices().items() if Y.curve() == cusp_branch.curve()][0]
        a2 = self.add_vertex(tail_component)
        self._edges.append([(a1, cusp_branch), (a2, tail_component)])
        self.graph().add_edge(a1, a2)
        
    def as_text(self, show_edges=False, show_adjacency=True):
        r"""Return a text-only representation of the component graph.

        The output is meant to be stable and readable in docstrings.

        INPUT:

        - ``show_edges`` -- bool (default: False); include an edge list
        - ``show_adjacency`` -- bool (default: True); include adjacency lists

        OUTPUT:

        A string describing this component graph.

        """
        lines = []
        lines.append(f"Component graph over {self.base_field()}")
        lines.append(f"Vertices: {self.number_of_vertices()}, Edges: {self.number_of_edges()}")

        # vertices (deterministic order)
        for v in sorted(self.vertices().keys()):
            Y = self.vertices()[v]
            g = Y.genus()
            # short component descriptor
            try:
                comp = Y.defining_polynomial()
            except Exception:
                comp = Y
            lines.append(f"  v{v}: g={g}, component: {comp}")

        G = self.graph()

        if show_adjacency:
            lines.append("Adjacency:")
            for v in sorted(self.vertices().keys()):
                nbrs = list(G.neighbors(v))
                nbrs.sort()
                lines.append(f"  v{v}: " + (", ".join(f"v{w}" for w in nbrs) if nbrs else "-"))

        if show_edges:
            # Graph.edges(labels=False) gives (u,v) or (u,v,label) depending on Sage version;
            # we normalize to pairs.
            raw = G.edges(labels=False, sort=True)
            pairs = []
            for e in raw:
                if len(e) >= 2:
                    u, v = e[0], e[1]
                    pairs.append((u, v))
            edge_str = ", ".join(f"(v{u},v{v})" for (u, v) in pairs) if pairs else "-"
            lines.append("Edges:")
            lines.append(f"  {edge_str}")

        return "\n".join(lines)



