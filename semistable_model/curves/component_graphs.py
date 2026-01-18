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

In this module we realize the computation of the component graph of a semistable
curve `X`.


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

"""

from semistable_model.curves.plane_curves import ProjectivePlaneCurve
from sage.all import SageObject, Graph, Curve, Set, lcm, GF
from itertools import product


def component_graph(X):
    r""" Return the component graph of a projective plane curve.
    
    INPUT:

    - `X` -- a projective plane curve over a field

    OUTPUT:

    the component graph of `X`.

    The curve `X` 
    - must be defined over a finite field, 
    - must be reduced, and
    - its geometric points have at most two branches.

    
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

    If the curve is not semistable, an error is raised:

        sage: X = Curve(x^3 - y^2*z)
        sage: component_graph(X)
        ---------------------------------------------------------------------------
        ValueError                                Traceback (most recent call last)
        ...
        ValueError: X is not semistable

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
    # cusps = []
    for P in singular_points:
        # test whether P is a node; there is a problem here to be solved!
        if is_node(X, P):
            nodes.append(P)
        # elif is_cusp(X, P):
        #     cusps.append(P)
        else:
            raise ValueError("X is not semistable")

    G = ComponentGraph(X.ambient_space())
    for Y in components:
        G.add_vertex(Y)
    for P in nodes:
        e = [b for b in branches if b.rational_point() == P]
        G.add_edge(e)
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
                # print(f"zeros of g: {g.zeros()}")
                places = places.union(Set(g.zeros()))
        ret += [Branch(Y, v) for v in places 
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


class Branch(SageObject):
    r""" Return a branch of a component graph.

    INPUT:

    - ``Y`` -- an irreducible plane curve
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
        
        If the point is not rational, ann error is raised.
        """
        assert self.point().degree() == 1, "the center of this point is not rational"
        return self.curve().ambient_space().point(self.point().rational_point())
    
    def residue_field(self):
        return self.place().residue_field()[0]


class ComponentGraph(SageObject):
    r""" The component graph of a semistable curve over a field.

    INPUT:

    - ``base_field`` -- a field

    OUTPUT:

    An empty component graph over the given field.

    To this component graph one can add new vertices and edges. 
    
    
    """

    def __init__(self, projective_plane):
        self._projective_plane = projective_plane
        self._base_field = projective_plane.base_ring()
        # upon creation, the component graph is empty 
        self._vertices = {}
        self._edges = []
        self._graph = Graph(0, loops=True, multiedges=True)

    def __repr__(self):
        return f"component graph of a projective curve over {self.base_field()}"
    
    def projective_plane(self):
        return self._projective_plane

    def base_field(self):
        r""" Return the base field of this component graph.
        
        """
        return self._base_field
    
    def genus(self):
        r""" Return the genus of this component graph.
        
        """
        if self.is_connected():
            return (sum(Y.geometric_genus() for Y in self.components()) 
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
    
    def branches(self):
        r""" Return the list of branches of this component graph.
        
        """
        return self._branches
    
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

        - ``Y`` -- an integral curve on the ambient projective plane


        The method adds a new vertex to the underlying graph,
        corresponding to the component `Y`
        
        """
        assert Y.ambient_space() == self.projective_plane()
        assert Y.is_irreducible()
        G = self._graph
        index = G.add_vertex()
        self._vertices[index] = Y
    
    def add_edge(self, e):
        r""" Add an edge to this component graph.

        INPUT:

        ``e`` -- a pair of distinct branches with the same center 

        
        """
        G = self.graph()
        b1 = e[0]
        b2 = e[1]
        assert b1.rational_point() == b2.rational_point()
        a1 = [a for a, Y in self.vertices().items() if Y == b1.curve()][0]
        a2 = [a for a, Y in self.vertices().items() if Y == b2.curve()][0]
        self._edges.append([(a1, b1), (a2, b2)])
        # print(f"vertices = {G.vertices()}")
        # print(f"a1 = {a1}, a2 = {a2}")
        G.add_edge(a1, a2)

    def as_text(self, show_edges=False, show_adjacency=True):
        r"""Return a text-only representation of the component graph.

        The output is meant to be stable and readable in docstrings.

        INPUT:

        - ``show_edges`` -- bool (default: False); include an edge list
        - ``show_adjacency`` -- bool (default: True); include adjacency lists

        OUTPUT:

        A string.

        """
        lines = []
        lines.append(f"Component graph over {self.base_field()}")
        lines.append(f"Vertices: {self.number_of_vertices()}, Edges: {self.number_of_edges()}")

        # vertices (deterministic order)
        for v in sorted(self.vertices().keys()):
            Y = self.vertices()[v]
            try:
                g = Y.geometric_genus()
            except Exception:
                g = "?"
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

    

class Component(SageObject):
    r""" Return a vertex of a component graph.

    INPUT:

    - ``curve`` -- an absolutely reduced and irreducible projective curve over a field
    
    Here ``curve`` has to be (for the moment) an instance of 
    `class`:semistable_model.curves.plane_curves.ProjectivePlaneCurves, and satisfy 
    the following additional conditions:
    
    - its base field must be a finite field
    - it has to be absolutely reduced and irreducible

    OUTPUT:

    An object representing the smooth projective model of ``curve``, as a component
    of a semistable curve. 

    To register this component as a vertex of a component graph, it has to be
    explicitly added to the graph, via `meth`:ComponentGraph.add_component.

    We do not check or even demand that the curve is itself semistable, as we only
    consider its smooth projective model.  

    """
    def __init__(self, curve):
        assert isinstance(curve, ProjectivePlaneCurve), \
            "the curve must be an instance of *ProjectivePlaneCurve*"
        assert curve.is_reduced(), "the curve must be reduced"
        assert curve.is_irreducible(), "the curve must be irreducible"
        self._curve = curve
        self._base_field = curve.base_field()
        self._genus = curve.geometric_genus()

    def __repr__(self):
        if self.curve().is_smooth():
            return self.curve()
        else:
            return f"smooth projective model of {self.curve()}"
        
    def base_field(self):
        return self._base_field
    
    def curve(self):
        return self._curve
    
    def genus(self):
        return self._genus
        


class Node(SageObject):
    r""" Return a node of a component graph. 
    
    INPUT:

    - ``branches`` -- a pair of branches

    """