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
