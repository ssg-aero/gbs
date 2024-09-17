from pygbs.gbs import Curve2d, BSCurve2d, BSCurveRational2d, discretize_curve, norm, deviation_based_params
import plotly.graph_objects as go
from typing import Union
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
        weights = [ 1. for p in poles ]
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

        size  = [10*w for w in weights]
        color = [5 for p in poles]

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


def plot_bs_curve_2d(curves, names=None, width=800, height=600, ctrl_pts_on=True,  crv_pts_on=False, showlegend=False, proj_func=None):

    fig = go.Figure()

    if isinstance(curves, list):
        index = 0
        for curve in curves:
            name = ''
            if names is not None and index < len(names):
                name = names[index]
            add_bs_curve_to_fig(curve, fig, name,ctrl_pts_on=ctrl_pts_on, crv_pts_on=crv_pts_on, proj_func=proj_func)
            index+=1
    else:
        name = ''
        if names:
            name = names
        add_bs_curve_to_fig(curves, fig, name, ctrl_pts_on=ctrl_pts_on, crv_pts_on=crv_pts_on, proj_func=proj_func)

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

