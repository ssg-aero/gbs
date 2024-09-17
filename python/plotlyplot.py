from pygbs.gbs import Curve2d, BSCurve2d, BSCurveRational2d, discretize_curve, norm, deviation_based_params
import plotly.graph_objects as go
from typing import Union, Sequence, Dict
from math import log, log10


def pts2d_to_xy(pts, proj_func=None):
    x = []
    y = []
    if proj_func is None:
        for pt in pts:
            x.append(pt[0])
            y.append(pt[1])
    else:
        pt = proj_func(pt)
        for pt in pts:
            x.append(pt[0])
            y.append(pt[1])
    return x, y


def xy_from_pt2d(pts):
    x = []
    y = []
    for pt in pts:
        x.append(pt[0])
        y.append(pt[1])
    return x, y


def xy_from_crv_eval(u, curve: Curve2d, use_log10=False, d=0):
    x = []
    y = []
    for u_ in u:
        pt = curve(u_, d)
        if use_log10:
            x.append(u_)
            cu = norm(pt)
            if cu > 0.:
                y.append(log10(cu))
        else:
            x.append(pt[0])
            y.append(pt[1])
    return x, y


def plot_curve_curvature(crv, width=800, height=600):
    u = deviation_based_params(crv)
    y = []
    for u_ in u:
        pt = crv(u_, 2)
        cu = norm(pt)
        if cu > 0.:
            y.append(log10(cu))

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=u,
            y=y,
            mode='markers+lines',
        )
    )
    fig.update_layout(
        width=width,
        height=height,
        title={
            'text': "Curvature",
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
        },
    )

    fig.update_yaxes(type="log")

    return fig


def add_bs_curve_to_fig(curve: Curve2d, fig: go.Figure, name="", ctrl_pts_on=True, crv_pts_on=False, proj_func=None, dash=None):

    if isinstance(curve, BSCurve2d):
        poles = curve.poles()
        weights = [1.] * len(poles)
    elif isinstance(curve, BSCurveRational2d):
        poles = curve.polesProjected()
        weights = curve.weights()
    else:
       ctrl_pts_on = False 

    pts = discretize_curve(curve)

    x_pts, y_pts = pts2d_to_xy(pts, proj_func)

    if crv_pts_on:
        mode = 'lines+markers'
    else:
        mode = 'lines'

    fig.add_trace(
        go.Scatter(
            x = x_pts,
            y = y_pts,
            mode = mode,
            name = name,
            line = dict(dash=dash)
        )
    )

    if ctrl_pts_on:
        
        x_poles, y_poles = pts2d_to_xy(poles)

        size = [10 * w for w in weights]
        color = [5] * len(poles)

        fig.add_trace(
            go.Scatter(
                x = x_poles,
                y = y_poles,
                mode = 'lines+markers',
                line = dict(
                    color='grey',
                    dash='dash',
                ),
                marker = dict(
                    size=size,
                    color=color,
                ),
                name = f"{name}_poles",
            )
        )


def plot_bs_curve_2d(
    curves: Union[Curve2d, Sequence[Curve2d], Dict[str, Curve2d]],
    width=800,
    height=600,
    ctrl_pts_on=True,
    crv_pts_on=False,
    showlegend=False,
    proj_func=None,
):
    """Plot 2D curve with plotly.

    Parameters:
    -----------
    curves (Curve2d | sequence[Curve2d] | dict[str, Curve2d]):
        GBS curve(s) to plot. If a dictionary is given, keys are used as curve names.
    width (int, optional):
        Width of the figure (default: 800).
    height (int, optional):
        Height of the figure (default: 600).
    ctrl_pts_on (bool, optional):
        If `True` (default), control points are shown.
    crv_pts_on (bool, optional):
        Should curve points be shown? Default: `False`.
    showlegend (bool, optional):
        Should legend be shown? Default: `False`.
    proj_func (callable, optional):
        Projection function. Default: `None`.
    
    Returns:
    --------
    `plotly.graph_objects.Figure` object
    """
    fig = go.Figure()

    options = dict(ctrl_pts_on=ctrl_pts_on, crv_pts_on=crv_pts_on, proj_func=proj_func)

    def add_curve(curve: Curve2d, name="") -> None:
        add_bs_curve_to_fig(curve, fig, name, **options)

    if isinstance(curves, dict):
        for name, curve in curves.items():
            add_curve(curve, name)

    elif isinstance(curves, Sequence):
        for curve in curves:
            add_curve(curve)

    else:
        add_curve(curves)

    fig.update_yaxes(
        scaleanchor="x",
        scaleratio=1.0,
    )
    fig.update_layout(
        width=width,
        height=height,
        showlegend=showlegend,
    )
    return fig
