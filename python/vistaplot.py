import pygbs.gbs as gbs
import pyvista as pv
import numpy as np
from pyvista.plotting import colors


def vtk_col_name_to_hex(vtk_col_name: str) -> str:
    """Convert vtk colors names to hexadecimal string

    Args:
        vtk_col_name (str): [description]

    Returns:
        str: [description]
    """
    col = colors.GetColor3ub('Peacock')
    return '#%02x%02x%02x' % (col[0], col[1], col[2])


def mesh_curves(crv_lst: list, npt_min=30, deviation=0.01, npt_max=5000) -> list:
    """build a list of KochanekSpline mesh to be displayed by pyvista

    Args:
        crv_lst (list): [description]
        npt_min (int, optional): [description]. Defaults to 30.
        deviation (float, optional): [description]. Defaults to 0.01.
        npt_max (int, optional): [description]. Defaults to 5000.

    Returns:
        list: [description]
    """
    msh_lst = []
    for crv in crv_lst:
        if isinstance(crv, gbs.Curve2d):
            if isinstance(crv, gbs.BSCurve2d):
                crv = gbs.to_bscurve_3d(crv)
            else:
                crv = gbs.approx(crv, deviation=deviation,
                                   tol=1e-6, p=5, bp=npt_min)
                crv = gbs.to_bscurve_3d(crv)
        crv_pt = gbs.discretize_curve(crv, npt_min, deviation, npt_max)
        msh_lst.append(pv.KochanekSpline(
            np.array(crv_pt), n_points=len(crv_pt)))
    return msh_lst


def add_curves_to_plotter(crv_lst: list, plotter: pv.Plotter, col_lst: list = [], def_col="Tomato") -> None:
    """Add curves to plotter

    Args:
        crv_lst (list): [description]
        plotter (pv.Plotter): [description]
        col_lst (list, optional): [description]. Defaults to [].
        def_col (str, optional): [description]. Defaults to "Tomato".
    """
    col_dft = []
    for i in range(len(col_lst), len(crv_lst)):
        col_dft.append(def_col)

    for msh, color in zip(mesh_curves(crv_lst), col_lst+col_dft):
        plotter.add_mesh(msh, color=color)


def plot_curves(crv_lst: list, col_lst: list = [], def_col="Tomato", jupyter_backend='pythreejs') -> pv.Plotter:
    """Create a pyVista plotter and add curves

    Args:
        crv_lst (list): Curves to plot.
        col_lst (list, optional): List of collors, if missing, filled with default. Defaults to [].
        def_col (str, optional): Default curve colors. Defaults to "Tomato".
        jupyter_backend (str, optional): Jupyter notebook plotting backend to use.  One of the
            following:
            * ``'none'`` : Do not display in the notebook.
            * ``'pythreejs'`` : Show a ``pythreejs`` widget
            * ``'static'`` : Display a static figure.
            * ``'ipygany'`` : Show a ``ipygany`` widget
            * ``'panel'`` : Show a ``panel`` widget.
    Returns:
        pv.Plotter: [description]
    """    
    plotter = pv.Plotter()
    # plotter.add_axes_at_origin()
    add_curves_to_plotter(crv_lst, plotter, col_lst, def_col)
    # plotter.add_axes()
    plotter.show(jupyter_backend=jupyter_backend)
    return plotter


def mesh_surfaces(srf_lst: list) -> list:
    msh_lst = []
    for srf in srf_lst:
        if isinstance(srf, gbs.BSSurface2d):
            srf = gbs.to_bssurface_3d(srf)
        srf_actor = gbs.make_surf3d_actor(srf, nu=100, nv=100)
        mapper = srf_actor.GetMapper()
        poly = mapper.GetInput()
        msh = pv.PolyData(poly)
        msh_lst.append(msh)
    return msh_lst


def add_surfaces_to_plotter(srf_lst: list, plotter: pv.Plotter, col_lst: list = [], def_col='#33a1c9', per=0, use_transparency=False,opacity=1.0) -> None:

    for i in range(len(col_lst), len(srf_lst)):
        col_lst.append(def_col)

    for msh, color in zip(mesh_surfaces(srf_lst), col_lst):
        plotter.add_mesh(msh, color=color, smooth_shading=True, culling=False, use_transparency=use_transparency, opacity=opacity)
        for i in range(1, per):
            msh_per = msh.copy(False)  # avoid mesh duplication
            msh_per.rotate_z((360.*i)/per)
            plotter.add_mesh(msh_per, color=color,
                             smooth_shading=True, culling=False)


def plot_surfaces(srf_lst: list, col_lst: list = [], def_col='#33a1c9', jupyter_backend='pythreejs', per=0, use_transparency=False,opacity=1.0) -> pv.Plotter:
    plotter = pv.Plotter()
    # plotter.add_axes_at_origin()
    add_surfaces_to_plotter(srf_lst, plotter, col_lst, def_col, per, use_transparency, opacity)
    # plotter.add_axes()
    # plotter.show_axes()
    plotter.show(jupyter_backend=jupyter_backend, interactive_update=True)
    return plotter
