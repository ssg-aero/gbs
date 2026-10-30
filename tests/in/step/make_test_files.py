"""Generates the small STEP test files of tests/in/step (AP214 BREP subset).

Usage: python make_test_files.py tests/in/step
"""
import math, sys, os

OUT = sys.argv[1]

def r(x):
    x = float(x)
    if x == int(x) and abs(x) < 1e15:
        return f"{int(x)}."
    s = f"{x:.15G}"
    if "E" in s:
        m, e = s.split("E")
        if "." not in m:
            m += "."
        return f"{m}E{e}"
    return s

class Step:
    def __init__(self, name, units="mm", angle="rad"):
        self.l = []
        self.name = name
        u = []
        if units == "mm":
            u.append(self.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))"))
            self.len = u[0]
        else:  # inch
            mm = self.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
            m = self.add(f"LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),{mm})")
            d = self.add("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)")
            u.append(self.add(f"(CONVERSION_BASED_UNIT('INCH',{m})LENGTH_UNIT()NAMED_UNIT({d}))"))
            self.len = u[0]
        rad = self.add("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        if angle == "rad":
            u.append(rad)
        else:
            m = self.add(f"PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.0174532925199433),{rad})")
            d = self.add("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)")
            u.append(self.add(f"(CONVERSION_BASED_UNIT('DEGREE',{m})NAMED_UNIT({d})PLANE_ANGLE_UNIT())"))
        u.append(self.add("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())"))
        unc = self.add(f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),{self.len},'distance_accuracy_value','confusion accuracy')")
        self.ctx = self.add(f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT(({unc}))GLOBAL_UNIT_ASSIGNED_CONTEXT(({','.join(u)}))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))")

    def add(self, s):
        self.l.append(s)
        return f"#{len(self.l)}"

    # geometry
    def pt(self, p): return self.add(f"CARTESIAN_POINT('',({','.join(r(x) for x in p)}))")
    def dir(self, d): return self.add(f"DIRECTION('',({','.join(r(x) for x in d)}))")
    def ax2(self, o, z=(0, 0, 1), x=(1, 0, 0)): return self.add(f"AXIS2_PLACEMENT_3D('',{self.pt(o)},{self.dir(z)},{self.dir(x)})")
    def line(self, p, d, mag=1.):
        return self.add(f"LINE('',{self.pt(p)},{self.add(f'VECTOR({chr(39)}{chr(39)},{self.dir(d)},{r(mag)})')})")
    def circle(self, o, rad, z=(0, 0, 1), x=(1, 0, 0)): return self.add(f"CIRCLE('',{self.ax2(o, z, x)},{r(rad)})")
    # topology
    def vertex(self, p): return self.add(f"VERTEX_POINT('',{self.pt(p)})")
    def edge(self, v1, v2, crv, same=True): return self.add(f"EDGE_CURVE('',{v1},{v2},{crv},{'.T.' if same else '.F.'})")
    def oe(self, e, fwd=True): return self.add(f"ORIENTED_EDGE('',*,*,{e},{'.T.' if fwd else '.F.'})")
    def loop(self, coedges): return self.add(f"EDGE_LOOP('',({','.join(self.oe(e, f) for e, f in coedges)}))")
    def bound(self, coedges, outer=True, orient=True):
        t = "FACE_OUTER_BOUND" if outer else "FACE_BOUND"
        return self.add(f"{t}('',{self.loop(coedges)},{'.T.' if orient else '.F.'})")
    def face(self, bounds, srf, same=True, name=""):
        return self.add(f"ADVANCED_FACE('{name}',({','.join(bounds)}),{srf},{'.T.' if same else '.F.'})")

    def product(self, rep):
        app = self.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
        self.add(f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{app})")
        pc = self.add(f"PRODUCT_CONTEXT('',{app},'mechanical')")
        p = self.add(f"PRODUCT('{self.name}','{self.name}','',({pc}))")
        pdf = self.add(f"PRODUCT_DEFINITION_FORMATION('','',{p})")
        pdc = self.add(f"PRODUCT_DEFINITION_CONTEXT('part definition',{app},'design')")
        pd = self.add(f"PRODUCT_DEFINITION('design','',{pdf},{pdc})")
        pds = self.add(f"PRODUCT_DEFINITION_SHAPE('','',{pd})")
        self.add(f"SHAPE_DEFINITION_REPRESENTATION({pds},{rep})")

    def solid(self, faces, name=""):
        sh = self.add(f"CLOSED_SHELL('',({','.join(faces)}))")
        so = self.add(f"MANIFOLD_SOLID_BREP('{name}',{sh})")
        rep = self.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('',({so},{self.ax2((0, 0, 0))}),{self.ctx})")
        self.product(rep)

    def surface_model(self, faces):
        sh = self.add(f"OPEN_SHELL('',({','.join(faces)}))")
        sm = self.add(f"SHELL_BASED_SURFACE_MODEL('',({sh}))")
        rep = self.add(f"MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',({sm},{self.ax2((0, 0, 0))}),{self.ctx})")
        self.product(rep)

    def write(self, fname, comment):
        with open(os.path.join(OUT, fname), "w") as f:
            f.write("ISO-10303-21;\nHEADER;\n")
            f.write(f"/* {comment} */\n")
            f.write("FILE_DESCRIPTION(('gbs STEP reader test'),'2;1');\n")
            f.write(f"FILE_NAME('{fname}','2026-10-10T00:00:00',(''),(''),'','','');\n")
            f.write("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n")
            for i, s in enumerate(self.l):
                f.write(f"#{i + 1}={s};\n")
            f.write("ENDSEC;\nEND-ISO-10303-21;\n")


def box(fname="box.stp", broken=False):
    """Box 2 x 3 x 4: half of the planes have an inward axis (same_sense .F.), some lines run backwards."""
    s = Step("BOX")
    L = (2., 3., 4.)
    c = lambda i: tuple(L[k] * ((i >> k) & 1) for k in range(3))
    V = [s.vertex(c(i)) for i in range(8)]
    E = {}
    for i in range(8):
        for b in range(3):
            if not i & (1 << b):
                j = i | (1 << b)
                d = [0., 0., 0.]
                d[b] = 1.
                if (i + b) % 2:  # line from j to i, used backwards
                    E[(i, j)] = s.edge(V[i], V[j], s.line(c(j), [-x for x in d]), False)
                else:
                    E[(i, j)] = s.edge(V[i], V[j], s.line(c(i), d), True)
    co = lambda p, q: (E[(p, q)], True) if p < q else (E[(q, p)], False)
    faces = []
    for axis in range(3):
        for side in range(2):
            a1, a2 = (axis + 1) % 3, (axis + 2) % 3
            idx = lambda u, v: (side << axis) | (u << a1) | (v << a2)
            cs = [idx(0, 0), idx(1, 0), idx(1, 1), idx(0, 1)]  # counter-clockwise around +axis
            if side == 0:
                cs.reverse()  # counter-clockwise seen from outside (-axis)
            n = [0., 0., 0.]
            n[axis] = 1.
            x = [0., 0., 0.]
            x[a1] = 1.
            inward = (axis + side) % 2 == 0  # plane axis against the outward normal
            outward = n if side == 1 else [-v for v in n]
            z = [-v for v in outward] if inward else outward
            srf = s.add(f"PLANE('',{s.ax2(c(cs[0]), z, x)})")
            if broken and axis == 2 and side == 1:
                srf = s.add(f"SURFACE_REPLICA('',{srf},$)")
            name = "TOP" if (axis == 2 and side == 1) else ""
            faces.append(s.face([s.bound([co(cs[k], cs[(k + 1) % 4]) for k in range(4)])], srf, not inward, name))
    s.solid(faces, "")
    s.write(fname, "box 2 x 3 x 4 mm" + (", top face on an unsupported surface" if broken else ""))


def cylinder():
    s = Step("CYLINDER")
    R, h = 5., 10.
    vb, vt = s.vertex((R, 0, 0)), s.vertex((R, 0, h))
    bottom = s.edge(vb, vb, s.circle((0, 0, 0), R))
    top = s.edge(vt, vt, s.circle((0, 0, h), R))
    seam = s.edge(vb, vt, s.add(f"SEAM_CURVE('',{s.line((R, 0, 0), (0, 0, 1))},(),.CURVE_3D.)"))
    lat = s.face([s.bound([(bottom, True), (seam, True), (top, False), (seam, False)])], s.add(f"CYLINDRICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True, "LATERAL")
    fb = s.face([s.bound([(bottom, False)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), False)
    ft = s.face([s.bound([(top, True)])], s.add(f"PLANE('',{s.ax2((0, 0, h))})"), True)
    s.solid([lat, fb, ft])
    s.write("cylinder.stp", "cylinder R 5 h 10 mm, seam edge used twice")


def sphere():
    s = Step("SPHERE")
    R = 3.
    S, N = s.vertex((0, 0, -R)), s.vertex((0, 0, R))
    meridian = s.circle((0, 0, 0), R, (0, -1, 0), (1, 0, 0))  # plane xz, from S through +x to N
    seam = s.edge(S, N, s.add(f"SEAM_CURVE('',{meridian},(),.CURVE_3D.)"))
    f = s.face([s.bound([(seam, True), (seam, False)])], s.add(f"SPHERICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    s.solid([f])
    s.write("sphere.stp", "sphere R 3 mm bounded by its seam only")


def cone():
    s = Step("CONE", angle="deg")
    R, a = 2., 30.
    H = R / math.tan(math.radians(a))
    A, B = s.vertex((0, 0, -H)), s.vertex((R, 0, 0))
    base = s.edge(B, B, s.circle((0, 0, 0), R))
    seam = s.edge(A, B, s.line((0, 0, -H), (R, 0, H), math.hypot(R, H) / math.hypot(R, H)))
    lat = s.face([s.bound([(seam, True), (base, False), (seam, False)])], s.add(f"CONICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)},{r(a)})"), True)
    top = s.face([s.bound([(base, True)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), True)
    s.solid([lat, top])
    s.write("cone.stp", "cone R 2 mm, semi-angle 30 degrees (angles in degrees), apex without edge")


def plate():
    s = Step("PLATE", units="inch")
    W, Hh, cx, cy, rr = 4., 2., 1., 1., 0.5
    pts = [(0, 0, 0), (W, 0, 0), (W, Hh, 0), (0, Hh, 0)]
    V = [s.vertex(p) for p in pts]
    E = []
    for k in range(4):
        p, q = pts[k], pts[(k + 1) % 4]
        d = [q[i] - p[i] for i in range(3)]
        n = math.sqrt(sum(x * x for x in d))
        E.append(s.edge(V[k], V[(k + 1) % 4], s.line(p, [x / n for x in d])))
    vh = s.vertex((cx + rr, cy, 0))
    hole = s.edge(vh, vh, s.circle((cx, cy, 0), rr))  # counter-clockwise
    f = s.face([s.bound([(e, True) for e in E]), s.bound([(hole, True)], outer=False, orient=False)],
               s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), True, "PLATE")
    s.surface_model([f])
    s.write("plate_inch.stp", "plate 4 x 2 in with a hole of radius 0.5 in, lengths in inches, open shell")


def bspline_patch():
    """Quarter of cylinder R 1, height 2, as a rational B-spline surface bounded by two rational arcs and two lines."""
    s = Step("PATCH")
    w = math.sqrt(0.5)
    g = [[(1, 0, 0), (1, 0, 2)], [(1, 1, 0), (1, 1, 2)], [(0, 1, 0), (0, 1, 2)]]  # [u][v]
    P = [[s.pt(p) for p in row] for row in g]
    grid = ",".join("(" + ",".join(row) + ")" for row in P)
    srf = s.add(f"(BOUNDED_SURFACE()B_SPLINE_SURFACE(2,1,({grid}),.UNSPECIFIED.,.F.,.F.,.F.)"
                f"B_SPLINE_SURFACE_WITH_KNOTS((3,3),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()"
                f"RATIONAL_B_SPLINE_SURFACE(((1.,1.),({r(w)},{r(w)}),(1.,1.)))REPRESENTATION_ITEM('')SURFACE())")
    def arc(z):
        cps = ",".join(s.pt((x, y, z)) for x, y, _ in [(1, 0, 0), (1, 1, 0), (0, 1, 0)])
        return s.add(f"(BOUNDED_CURVE()B_SPLINE_CURVE(2,({cps}),.CIRCULAR_ARC.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)"
                     f"CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,{r(w)},1.))REPRESENTATION_ITEM(''))")
    v00, v10, v11, v01 = s.vertex((1, 0, 0)), s.vertex((0, 1, 0)), s.vertex((0, 1, 2)), s.vertex((1, 0, 2))
    e0 = s.edge(v00, v10, arc(0))
    e1 = s.edge(v10, v11, s.line((0, 1, 0), (0, 0, 1)))
    e2 = s.edge(v01, v11, arc(2))
    e3 = s.edge(v00, v01, s.line((1, 0, 0), (0, 0, 1)))
    f = s.face([s.bound([(e0, True), (e1, True), (e2, False), (e3, False)])], srf, True)
    s.surface_model([f])
    s.write("bspline_patch.stp", "quarter cylinder as a rational B-spline surface, open shell")


box()
box("box_unsupported_face.stp", broken=True)
cylinder()
sphere()
cone()
plate()
bspline_patch()


def band_cylinder():
    """Cylinder whose lateral face is bounded by its two circles only (no seam edge), vertices at different angles."""
    s = Step("BAND")
    R, h, a = 5., 10., 1.
    vb, vt = s.vertex((R, 0, 0)), s.vertex((R * math.cos(a), R * math.sin(a), h))
    bottom = s.edge(vb, vb, s.circle((0, 0, 0), R))
    top = s.edge(vt, vt, s.circle((0, 0, h), R, (0, 0, 1), (math.cos(a), math.sin(a), 0)))
    lat = s.face([s.bound([(bottom, True)]), s.bound([(top, False)], outer=False)],
                 s.add(f"CYLINDRICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    fb = s.face([s.bound([(bottom, False)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), False)
    ft = s.face([s.bound([(top, True)])], s.add(f"PLANE('',{s.ax2((0, 0, h))})"), True)
    s.solid([lat, fb, ft])
    s.write("band_cylinder.stp", "cylinder R 5 h 10 mm, lateral face without seam edge, circle vertices at 0 and 1 rad")


def sphere_cap():
    """Cap of a sphere R 3 above z = 1.5, bounded by one circle, closed by a disk."""
    s = Step("CAP")
    R, z = 3., 1.5
    rc = math.sqrt(R * R - z * z)
    v = s.vertex((rc, 0, z))
    c = s.edge(v, v, s.circle((0, 0, z), rc))
    cap = s.face([s.bound([(c, True)])], s.add(f"SPHERICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    disk = s.face([s.bound([(c, False)])], s.add(f"PLANE('',{s.ax2((0, 0, z))})"), False)
    s.solid([cap, disk])
    s.write("sphere_cap.stp", "spherical cap R 3 mm above z = 1.5, without seam nor pole")


def cone_noseam():
    s = Step("CONE2")
    R, a = 2., math.pi / 6
    H = R / math.tan(a)
    B = s.vertex((R, 0, 0))
    base = s.edge(B, B, s.circle((0, 0, 0), R))
    lat = s.face([s.bound([(base, False)])], s.add(f"CONICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)},{r(a)})"), True)
    top = s.face([s.bound([(base, True)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), True)
    s.solid([lat, top])
    s.write("cone_noseam.stp", "cone R 2 mm, semi-angle pi/6, lateral face bounded by its base circle only")


def band_with_hole():
    """Cylinder band without seam edge, with a rectangular hole across the angle 0."""
    s = Step("BAND_HOLE")
    R, h, a, z1, z2 = 5., 10., 0.3, 3., 6.
    vb, vt = s.vertex((R, 0, 0)), s.vertex((R, 0, h))
    bottom = s.edge(vb, vb, s.circle((0, 0, 0), R))
    top = s.edge(vt, vt, s.circle((0, 0, h), R))
    P = lambda t, z: (R * math.cos(t), R * math.sin(t), z)
    v1, v2, v3, v4 = s.vertex(P(-a, z1)), s.vertex(P(a, z1)), s.vertex(P(a, z2)), s.vertex(P(-a, z2))
    arc1 = s.edge(v1, v2, s.circle((0, 0, z1), R))  # across the angle 0
    arc2 = s.edge(v4, v3, s.circle((0, 0, z2), R))
    l1 = s.edge(v1, v4, s.line(P(-a, z1), (0, 0, 1)))
    l2 = s.edge(v2, v3, s.line(P(a, z1), (0, 0, 1)))
    hole = [(arc1, False), (l1, True), (arc2, True), (l2, False)]  # clockwise seen from outside
    f = s.face([s.bound([(bottom, True)]), s.bound([(top, False)], outer=False), s.bound(hole, outer=False)],
               s.add(f"CYLINDRICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    s.surface_model([f])
    s.write("band_with_hole.stp", "cylinder band R 5 h 10 mm without seam edge, rectangular hole from -0.3 to 0.3 rad")


def half_spheres():
    s = Step("HALF_SPHERES")
    R = 3.
    S, N = s.vertex((0, 0, -R)), s.vertex((0, 0, R))
    srf = s.add(f"SPHERICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})")
    # meridians in the plane y = 0, from S to N: x > 0 (the seam of the sphere) and x < 0
    e0 = s.edge(S, N, s.circle((0, 0, 0), R, (0, -1, 0), (1, 0, 0)))
    epi = s.edge(S, N, s.circle((0, 0, 0), R, (0, 1, 0), (-1, 0, 0)))
    # y > 0: counter-clockwise seen from outside. y < 0: the face normal points inwards (same_sense .F.),
    # its loop is counter-clockwise around it, the shell takes the face reversed
    fa = s.face([s.bound([(epi, True), (e0, False)])], srf, True)
    fb = s.face([s.bound([(epi, True), (e0, False)])], srf, False)
    s.solid([fa, s.add(f"ORIENTED_FACE('',*,{fb},.F.)")])
    s.write("half_spheres.stp", "sphere R 3 mm cut by the plane y = 0 through its seam; the half y < 0 has same_sense .F.")


def sphere_wedge(w=4.):
    """Wedge of a sphere R 3 between the meridians at 0 and w > pi rad: a lune and two half disks."""
    s = Step("WEDGE")
    R = 3.
    S, N = s.vertex((0, 0, -R)), s.vertex((0, 0, R))
    d0, dw = (1, 0, 0), (math.cos(w), math.sin(w), 0)
    meridian = lambda d: s.edge(S, N, s.circle((0, 0, 0), R, (d[1], -d[0], 0), d))
    m0, mw = meridian(d0), meridian(dw)
    axis = s.edge(S, N, s.line((0, 0, -R), (0, 0, 1)))
    lune = s.face([s.bound([(mw, True), (m0, False)])], s.add(f"SPHERICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    p0 = s.face([s.bound([(m0, True), (axis, False)])], s.add(f"PLANE('',{s.ax2((0, 0, 0), (0, -1, 0), d0)})"), True)
    pw = s.face([s.bound([(mw, False), (axis, True)])], s.add(f"PLANE('',{s.ax2((0, 0, 0), (-dw[1], dw[0], 0), dw)})"), True)
    s.solid([lune, p0, pw])
    s.write("sphere_wedge.stp", f"wedge of a sphere R 3 mm between the meridians at 0 and {r(w)} rad (wider than pi)")


def sphere_cut(h, fname):
    """Sphere R 3 (axis z) cut by the plane x = h: the spherical face, bounded by one circle, holds both poles."""
    s = Step("SPHERE_CUT")
    R = 3.
    rc = math.sqrt(R * R - h * h)
    v = s.vertex((h, rc, 0))
    c = s.edge(v, v, s.circle((h, 0, 0), rc, (1, 0, 0), (0, 1, 0)))
    big = s.face([s.bound([(c, False)])], s.add(f"SPHERICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    disk = s.face([s.bound([(c, True)])], s.add(f"PLANE('',{s.ax2((h, 0, 0), (1, 0, 0), (0, 1, 0))})"), True)
    s.solid([big, disk])
    s.write(fname, f"sphere R 3 mm (axis z) cut by the plane x = {r(h)}, the part x < {r(h)}: its spherical face holds both poles")


def closed_bspline(across, fname):
    """Disk of z = 0 bounded by a closed cubic B-spline (area 10.4) used by two edges."""
    s = Step("CLOSED_BSPLINE")
    poles = [(2, 0, 0), (2, 2, 0), (-2, 2, 0), (-2, -2, 0), (2, -2, 0), (2, 0, 0)]
    crv = s.add(f"B_SPLINE_CURVE_WITH_KNOTS('',3,({','.join(s.pt(p) for p in poles)}),.UNSPECIFIED.,.T.,.F.,(4,1,1,4),(0.,1.,2.,3.),.UNSPECIFIED.)")
    # points of the curve at t = 0.7, 1.5 and 2.2 (de Boor, computed by the reader's B-spline)
    def at(t):
        k = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
        def basis(i, p, t):
            if p == 0:
                return 1. if k[i] <= t < k[i + 1] else 0.
            a = 0. if k[i + p] == k[i] else (t - k[i]) / (k[i + p] - k[i]) * basis(i, p - 1, t)
            b = 0. if k[i + p + 1] == k[i + 1] else (k[i + p + 1] - t) / (k[i + p + 1] - k[i + 1]) * basis(i + 1, p - 1, t)
            return a + b
        return tuple(sum(basis(i, 3, t) * poles[i][j] for i in range(6)) for j in range(3))
    # across: vertices at t = 2.2 and 0.7, the second edge runs through the closure point at t = 0 = 3;
    # otherwise at the closure point and t = 1.5, the second edge ends there
    va, vb = (s.vertex(at(2.2)), s.vertex(at(0.7))) if across else (s.vertex(poles[0]), s.vertex(at(1.5)))
    e1, e2 = (s.edge(vb, va, crv), s.edge(va, vb, crv))
    f = s.face([s.bound([(e1, True), (e2, True)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), True)
    s.surface_model([f])
    what = "an edge running through its closure point" if across else "an edge ending at its closure point"
    s.write(fname, f"disk bounded by a closed cubic B-spline curve used by two edges, {what}")


def torus_band():
    """Rounded disk: torus R 5 r 1 band from z = 0 (radius 6) to z = 1 (radius 5) without seam edge, two disks."""
    s = Step("TORUS_BAND")
    a = 1.
    v0, v1 = s.vertex((6, 0, 0)), s.vertex((5 * math.cos(a), 5 * math.sin(a), 1))
    c0 = s.edge(v0, v0, s.circle((0, 0, 0), 6.))
    c1 = s.edge(v1, v1, s.circle((0, 0, 1), 5., (0, 0, 1), (math.cos(a), math.sin(a), 0)))
    band = s.face([s.bound([(c0, True)]), s.bound([(c1, False)], outer=False)], s.add(f"TOROIDAL_SURFACE('',{s.ax2((0, 0, 0))},5.,1.)"), True)
    bottom = s.face([s.bound([(c0, False)])], s.add(f"PLANE('',{s.ax2((0, 0, 0), (0, 0, -1))})"), True)
    top = s.face([s.bound([(c1, True)])], s.add(f"PLANE('',{s.ax2((0, 0, 1))})"), True)
    s.solid([band, bottom, top])
    s.write("torus_band.stp", "torus R 5 r 1 mm band (a fillet) between two circles without seam edge, closed by two disks")


def torus_elbow():
    """Quarter of a torus R 5 r 1 (a pipe elbow) bounded by two circles turning in v, without seam edge, two disks."""
    s = Step("TORUS_ELBOW")
    v0, v1 = s.vertex((6, 0, 0)), s.vertex((0, 6, 0))
    c0 = s.edge(v0, v0, s.circle((5, 0, 0), 1., (0, -1, 0), (1, 0, 0)))
    c1 = s.edge(v1, v1, s.circle((0, 5, 0), 1., (1, 0, 0), (0, 1, 0)))
    elbow = s.face([s.bound([(c1, True)]), s.bound([(c0, False)], outer=False)], s.add(f"TOROIDAL_SURFACE('',{s.ax2((0, 0, 0))},5.,1.)"), True)
    d0 = s.face([s.bound([(c0, True)])], s.add(f"PLANE('',{s.ax2((5, 0, 0), (0, -1, 0), (1, 0, 0))})"), True)
    d1 = s.face([s.bound([(c1, False)])], s.add(f"PLANE('',{s.ax2((0, 5, 0), (-1, 0, 0), (0, 1, 0))})"), True)
    s.solid([elbow, d0, d1])
    s.write("torus_elbow.stp", "quarter of a torus R 5 r 1 mm (a pipe elbow) between two circles without seam edge, closed by two disks")


band_cylinder()
sphere_cap()
cone_noseam()
band_with_hole()
half_spheres()
sphere_wedge()
sphere_cut(2.4, "sphere_cut_off_axis.stp")
sphere_cut(0., "sphere_cut_through_poles.stp")
closed_bspline(False, "closed_bspline_end.stp")
closed_bspline(True, "closed_bspline_across.stp")
torus_band()
torus_elbow()


def box_solid_items(s, L=(2., 3., 4.)):
    """Faces of a box (all plane axes outward) in the Step s; returns the closed shell's faces."""
    c = lambda i: tuple(L[k] * ((i >> k) & 1) for k in range(3))
    V = [s.vertex(c(i)) for i in range(8)]
    E = {}
    for i in range(8):
        for b in range(3):
            if not i & (1 << b):
                j = i | (1 << b)
                d = [0., 0., 0.]
                d[b] = 1.
                E[(i, j)] = s.edge(V[i], V[j], s.line(c(i), d), True)
    co = lambda p, q: (E[(p, q)], True) if p < q else (E[(q, p)], False)
    faces = []
    for axis in range(3):
        for side in range(2):
            a1, a2 = (axis + 1) % 3, (axis + 2) % 3
            idx = lambda u, v: (side << axis) | (u << a1) | (v << a2)
            cs = [idx(0, 0), idx(1, 0), idx(1, 1), idx(0, 1)]
            if side == 0:
                cs.reverse()
            n = [0., 0., 0.]
            n[axis] = 1. if side else -1.
            x = [0., 0., 0.]
            x[a1] = 1.
            srf = s.add(f"PLANE('',{s.ax2(c(cs[0]), n, x)})")
            faces.append(s.face([s.bound([co(cs[k], cs[(k + 1) % 4]) for k in range(4)])], srf, True))
    return faces


def cylinder_solid_items(s, R=1., h=5.):
    vb, vt = s.vertex((R, 0, 0)), s.vertex((R, 0, h))
    bottom = s.edge(vb, vb, s.circle((0, 0, 0), R))
    top = s.edge(vt, vt, s.circle((0, 0, h), R))
    seam = s.edge(vb, vt, s.line((R, 0, 0), (0, 0, 1)))
    lat = s.face([s.bound([(bottom, True), (seam, True), (top, False), (seam, False)])], s.add(f"CYLINDRICAL_SURFACE('',{s.ax2((0, 0, 0))},{r(R)})"), True)
    fb = s.face([s.bound([(bottom, False)])], s.add(f"PLANE('',{s.ax2((0, 0, 0))})"), False)
    ft = s.face([s.bound([(top, True)])], s.add(f"PLANE('',{s.ax2((0, 0, h))})"), True)
    return [lat, fb, ft]


class Assembly(Step):
    """Products, shape representations and placed occurrences, as OCCT writes assemblies."""

    def __init__(self, units="mm"):
        super().__init__("", units)
        self.app = self.add("APPLICATION_CONTEXT('core data for automotive mechanical design processes')")
        self.add(f"APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,{self.app})")
        self.pc = self.add(f"PRODUCT_CONTEXT('',{self.app},'mechanical')")
        self.pdc = self.add(f"PRODUCT_DEFINITION_CONTEXT('part definition',{self.app},'design')")

    def product(self, name, items=None):
        """A product; `items` (list of solids) makes it a part with a BREP representation.
        Returns (product_definition, shape_representation)."""
        p = self.add(f"PRODUCT('{name}','{name}','',({self.pc}))")
        pdf = self.add(f"PRODUCT_DEFINITION_FORMATION('','',{p})")
        pd = self.add(f"PRODUCT_DEFINITION('design','',{pdf},{self.pdc})")
        pds = self.add(f"PRODUCT_DEFINITION_SHAPE('','',{pd})")
        origin = self.ax2((0, 0, 0))
        sr = self.add(f"SHAPE_REPRESENTATION('',({origin}),{self.ctx})")
        self.add(f"SHAPE_DEFINITION_REPRESENTATION({pds},{sr})")
        if items:
            brep = self.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('',({','.join(items)},{self.ax2((0, 0, 0))}),{self.ctx})")
            self.add(f"SHAPE_REPRESENTATION_RELATIONSHIP('','',{sr},{brep})")
        return pd, sr, origin

    def occurrence(self, parent, child, o, z=(0, 0, 1), x=(1, 0, 0), name="instance"):
        """Child product placed in the parent at origin o with axes z, x."""
        nauo = self.add(f"NEXT_ASSEMBLY_USAGE_OCCURRENCE('{name}','{name}','',{parent[0]},{child[0]},$)")
        pds = self.add(f"PRODUCT_DEFINITION_SHAPE('','',{nauo})")
        place = self.ax2(o, z, x)
        idt = self.add(f"ITEM_DEFINED_TRANSFORMATION('','',{child[2]},{place})")
        rel = self.add(f"( REPRESENTATION_RELATIONSHIP('','',{child[1]},{parent[1]}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION({idt}) SHAPE_REPRESENTATION_RELATIONSHIP() )")
        self.add(f"CONTEXT_DEPENDENT_SHAPE_REPRESENTATION({rel},{pds})")


def assembly():
    """ASM: two instances of BOX (one moved, one rotated a quarter turn about z) and SUB, holding one CYL."""
    s = Assembly()
    box = s.product("BOX", [s.add(f"MANIFOLD_SOLID_BREP('',{s.add('CLOSED_SHELL(' + chr(39) + chr(39) + ',(' + ','.join(box_solid_items(s)) + '))')})")])
    cyl = s.product("CYL", [s.add(f"MANIFOLD_SOLID_BREP('',{s.add('CLOSED_SHELL(' + chr(39) + chr(39) + ',(' + ','.join(cylinder_solid_items(s)) + '))')})")])
    asm = s.product("ASM")
    sub = s.product("SUB")
    s.occurrence(asm, box, (10, 0, 0), name="box1")
    s.occurrence(asm, box, (0, 20, 0), (0, 0, 1), (0, 1, 0), name="box2")  # x -> y: quarter turn
    s.occurrence(asm, sub, (0, 0, 10), name="sub1")
    s.occurrence(sub, cyl, (5, 5, 0), name="cyl1")
    s.write("assembly.stp", "assembly ASM: BOX at (10,0,0), BOX turned a quarter about z at (0,20,0), SUB at (0,0,10) holding CYL at (5,5,0)")


def mapped_items():
    """Two instances of a box representation placed by mapped items in one shape representation."""
    s = Step("MAPPED")
    shell = s.add("CLOSED_SHELL('',(" + ",".join(box_solid_items(s)) + "))")
    solid = s.add(f"MANIFOLD_SOLID_BREP('',{shell})")
    brep = s.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('BOX',({solid}),{s.ctx})")
    rmap = s.add(f"REPRESENTATION_MAP({s.ax2((0, 0, 0))},{brep})")
    m1 = s.add(f"MAPPED_ITEM('',{rmap},{s.ax2((-5, 0, 0))})")
    m2 = s.add(f"MAPPED_ITEM('',{rmap},{s.ax2((5, 0, 0), (1, 0, 0), (0, 1, 0))})")  # z -> x
    rep = s.add(f"SHAPE_REPRESENTATION('',({m1},{m2}),{s.ctx})")
    Step.product(s, rep)
    s.write("mapped_items.stp", "two instances of a box placed by mapped items, at (-5,0,0) and turned at (5,0,0)")


assembly()
mapped_items()
