# genus3_reduction_types.py
r"""
Reduction types for stable curves of genus 3
============================================

This module classifies the component graph of a *stable* curve of genus 3.

INPUT:
  - a purely combinatorial ``ComponentGraph`` (from ``component_graphs.py``)

OUTPUT:
  - a string identifying the reduction type (one of the 42 types in Prop. 1.13)

Strategy
--------

We build a catalogue mapping canonical signatures to type strings:

    catalog[ G_template.canonical_signature() ] = "type_name"

Classification of an input graph ``G`` is then a dictionary lookup:

    catalog[ G.canonical_signature() ]

The catalogue is computed once and cached.

EXAMPLES::

    We classify component graphs of stable genus-3 curves by reduction type.

    First, an irreducible core of geometric genus 1 with two elliptic tails
    corresponds to the type ``"1ee"``::

        sage: from component_graphs import ComponentGraph
        sage: from genus3_reduction_types import classify_genus3_type
        sage: G = ComponentGraph()
        sage: v = G.add_vertex(1)          # core vertex
        sage: _ = G.add_elliptic_tail(v)
        sage: _ = G.add_elliptic_tail(v)
        sage: classify_genus3_type(G)
        '1ee'

    Next, a core of type ``"0---0n"``: two rational components meeting in three
    nodes, with a self-node on one component::

        sage: H = ComponentGraph()
        sage: a = H.add_vertex(0)
        sage: b = H.add_vertex(0)
        sage: H.add_edge(a, b); H.add_edge(a, b); H.add_edge(a, b)  # three nodes
        sage: H.add_edge(a, a)                                      # one self-node
        sage: classify_genus3_type(H)
        '0---0n'

"""

from functools import lru_cache

from semistable_model.curves.component_graphs import ComponentGraph


# --------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------

def classify_genus3_type(G):
    r"""
    Return the reduction type (string) of a stable genus-3 component graph.

    INPUT:
      - ``G`` -- a ``ComponentGraph`` representing a stable curve of genus 3

    OUTPUT:
      A string among the 42 types (e.g. ``"1ee"``, ``"BRAID"``, ...).

    Raises:
      KeyError if the canonical signature is not in the catalogue.
    """
    sig = G.canonical_signature()
    return genus3_catalog()[sig]


@lru_cache(maxsize=1)
def genus3_catalog():
    r"""
    Return the catalogue mapping canonical signatures to type strings.

    The catalogue is computed once and cached.
    """
    catalog = {}

    for typ in GENUS3_TYPES_42:
        T = template_graph_of_type(typ)
        sig = T.canonical_signature()
        if sig in catalog and catalog[sig] != typ:
            raise ValueError(
                "Signature collision: %s and %s have same signature."
                % (typ, catalog[sig])
            )
        catalog[sig] = typ

    if len(catalog) != 42:
        missing = set(GENUS3_TYPES_42) - set(catalog.values())
        raise ValueError(
            "Catalogue size is %s (expected 42). Missing: %s"
            % (len(catalog), sorted(missing))
        )
    return catalog


# --------------------------------------------------------------------
# List of the 42 type strings
# --------------------------------------------------------------------

_IRR = [
    "3", "2n", "1nn", "0nnn",
    "2e", "2m", "1ne", "1nm", "0nne", "0nnm",
    "1ee", "1me", "1mm", "0nee", "0nme", "0nmm",
    "0eee", "0mee", "0mme", "0mmm",
]

_INSEP_RED = ["1---0", "0---0n", "0----0", "CAVE", "BRAID", "0---0e", "0---0m"]

_EQ = [
    "1=1", "1=0n", "0n=0n",
    "1=0e", "1=0m",
    "0n=0e", "0n=0m",
    "0e=0e", "0m=0e", "0m=0m",
]

_Z = ["Z=1", "Z=0n", "Z=Z", "Z=0e", "Z=0m"]

GENUS3_TYPES_42 = _IRR + _INSEP_RED + _EQ + _Z


# --------------------------------------------------------------------
# Template graph construction
# --------------------------------------------------------------------

def template_graph_of_type(typ):
    r"""
    Construct a template ``ComponentGraph`` for the given type string.
    """
    if _looks_like_irreducible_label(typ):
        return _template_irreducible(typ)

    if typ in ("1---0", "0---0n", "0----0", "0---0e", "0---0m"):
        return _template_two_vertex_inseparable(typ)

    if typ == "CAVE":
        return _template_cave()

    if typ == "BRAID":
        return _template_braid()

    if "=" in typ and typ.startswith(("1", "0n", "0e", "0m")):
        return _template_equals(typ)

    if typ.startswith("Z="):
        return _template_Z_equals(typ)

    raise ValueError("Unknown genus-3 type label: %s" % typ)


def _add_one_tail(G, v, kind):
    """
    Attach a one-tail of kind 'e' or 'm' to core vertex v.

    Uses add_one_tail if present, otherwise falls back to
    add_elliptic_tail / add_pigtail.
    """
    if hasattr(G, "add_one_tail"):
        return G.add_one_tail(v, kind)
    if kind == "e":
        return G.add_elliptic_tail(v)
    if kind == "m":
        return G.add_pigtail(v)
    raise ValueError("tail kind must be 'e' or 'm'")


def _looks_like_irreducible_label(typ):
    """
    Return True if typ is one of the core-irreducible labels encoded by:
      digit + optional 'n'* + optional tail letters in {'e','m'}*
    e.g. '3', '2n', '0nnm', '1ee', ...
    """
    if typ in ("CAVE", "BRAID"):
        return False
    if typ.startswith("Z=") or "=" in typ or "---" in typ:
        return False
    if not typ or not typ[0].isdigit():
        return False
    for ch in typ[1:]:
        if ch not in ("n", "e", "m"):
            return False
    return True


def _template_irreducible(typ):
    """
    Build template for irreducible core type:
      'g' + 'n'* + tail letters ('e'/'m')*
    """
    g_core = int(typ[0])
    rest = typ[1:]

    n_loops = 0
    while rest.startswith("n"):
        n_loops += 1
        rest = rest[1:]

    tails = list(rest)  # each char is 'e' or 'm'

    G = ComponentGraph()
    v = G.add_vertex(g_core)
    for _ in range(n_loops):
        G.add_edge(v, v)

    for t in tails:
        _add_one_tail(G, v, t)

    return G


def _template_two_vertex_inseparable(typ):
    """
    Two-vertex 2-inseparable reducible core types:
      1---0, 0---0n, 0----0, 0---0e, 0---0m
    """
    G = ComponentGraph()

    if typ == "1---0":
        a = G.add_vertex(1)
        b = G.add_vertex(0)
        for _ in range(3):
            G.add_edge(a, b)
        return G

    if typ == "0---0n":
        a = G.add_vertex(0)
        b = G.add_vertex(0)
        for _ in range(3):
            G.add_edge(a, b)
        G.add_edge(a, a)
        return G

    if typ == "0----0":
        a = G.add_vertex(0)
        b = G.add_vertex(0)
        for _ in range(4):
            G.add_edge(a, b)
        return G

    if typ in ("0---0e", "0---0m"):
        a = G.add_vertex(0)
        b = G.add_vertex(0)
        for _ in range(3):
            G.add_edge(a, b)
        kind = "e" if typ.endswith("e") else "m"
        _add_one_tail(G, a, kind)
        return G

    raise ValueError("Not a two-vertex inseparable type: %s" % typ)


def _template_cave():
    """
    CAVE: 3 genus-0 core vertices, edge multiplicities (2,2,1) on the triangle.
    """
    G = ComponentGraph()
    a = G.add_vertex(0)
    b = G.add_vertex(0)
    c = G.add_vertex(0)

    G.add_edge(a, b); G.add_edge(a, b)
    G.add_edge(b, c); G.add_edge(b, c)
    G.add_edge(c, a)

    return G


def _template_braid():
    """
    BRAID: core is the simple K4 on 4 genus-0 vertices.
    """
    G = ComponentGraph()
    v = [G.add_vertex(0) for _ in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            G.add_edge(v[i], v[j])
    return G


def _vertex_block_from_token(G, token):
    """
    Build a single-vertex block of type token in {1,0n,0e,0m}.
    Return the core vertex.
    """
    if token == "1":
        return G.add_vertex(1)
    if token == "0n":
        v = G.add_vertex(0)
        G.add_edge(v, v)
        return v
    if token == "0e":
        v = G.add_vertex(0)
        _add_one_tail(G, v, "e")
        return v
    if token == "0m":
        v = G.add_vertex(0)
        _add_one_tail(G, v, "m")
        return v
    raise ValueError("Unknown token: %s" % token)


def _template_equals(typ):
    """
    Two-core-vertex 2-separable types: 'A=B' with A,B in {1,0n,0e,0m}.
    Core has two vertices joined by exactly two edges.
    """
    left, right = typ.split("=")
    G = ComponentGraph()
    a = _vertex_block_from_token(G, left)
    b = _vertex_block_from_token(G, right)
    G.add_edge(a, b)
    G.add_edge(a, b)
    return G


def _template_Z_equals(typ):
    """
    'Z=...' types.

    Z is a two-vertex genus-0 block with two parallel edges. Connect the Z block
    to the RHS block by two parallel edges (the separating pair).

    Allowed RHS: 1, 0n, Z, 0e, 0m
    """
    rhs = typ.split("=", 1)[1]
    G = ComponentGraph()

    # Z block
    z1 = G.add_vertex(0)
    z2 = G.add_vertex(0)
    G.add_edge(z1, z2)
    G.add_edge(z1, z2)

    # RHS block
    if rhs == "Z":
        w1 = G.add_vertex(0)
        w2 = G.add_vertex(0)
        G.add_edge(w1, w2)
        G.add_edge(w1, w2)
        anchor = w1
    else:
        anchor = _vertex_block_from_token(G, rhs)

    # separating pair: two edges between z1 and anchor
    G.add_edge(z1, anchor)
    G.add_edge(z1, anchor)

    return G


def sanity_check_templates():
    """
    Developer helper: build all templates and ensure they classify correctly.
    """
    cat = genus3_catalog()
    for typ in GENUS3_TYPES_42:
        T = template_graph_of_type(typ)
        sig = T.canonical_signature()
        if cat[sig] != typ:
            raise AssertionError("Template mismatch for %s" % typ)


if __name__ == "__main__":
    sanity_check_templates()
