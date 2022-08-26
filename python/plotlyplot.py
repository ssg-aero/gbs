from pygbs.gbs import Curve2d, BSCurve2d, BSCurveRational2d, discretize_curve, norm, deviation_based_params
import plotly.graph_objects as go
from typing import Union
from math import log, log10


def pts2d_to_xy(pts):
    x = []
    y = []
    for pt in pts:
        x.append( pt[0] )
        y.append( pt[1] )
    return x, y


def xy_from_pt2d(pts):
    x =[]
    y =[]
    for pt in pts:
        x.append(pt[0])
        y.append(pt[1])
    return x, y


def xy_from_crv_eval(u, crv: Curve2d, use_log10: bool = False, d:int = 0):
    x =[]
    y =[]
    for u_ in u:
        pt = crv(u_,d)
        if use_log10:
            x.append(u_)
            y.append(log10( norm( pt ) ))
        else:
            x.append(pt[0])
            y.append(pt[1])
    return x, y

def plot_curve_curvature(crv, width :int = 800, height : int = 600):
    u = deviation_based_params(crv)
    y =[]
    for u_ in u:
        pt = crv(u_,2)
        y.append(log10( norm( pt ) ))

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x = u,
            y = y,
            mode = 'markers+lines',
        )
    )

    fig.update_layout(
        width=width,
        height=height,
        title={
        'text': "Curvature",
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'}
    )

    fig.update_yaxes(type="log")

    return fig

def add_bs_curve_to_fig(crv: Union[BSCurve2d, BSCurveRational2d], fig: go.Figure, name = '', ctrl_pts_on: bool = True, crv_pts_on: bool = False):

    if isinstance(crv, BSCurve2d):
        poles = crv.poles()
        weights = [ 1. for p in poles ]
    else:
        poles = crv.polesProjected()
        weights = crv.weights()

    pts = discretize_curve(crv)

    x_poles, y_poles = pts2d_to_xy(poles)
    x_pts, y_pts = pts2d_to_xy(pts)

    if crv_pts_on :
        mode = 'lines+markers'
    else:
        mode = 'lines'

    fig.add_trace(
        go.Scatter(
            x = x_pts,
            y = y_pts,
            mode = mode,
            name=name,
        )
    )

    if ctrl_pts_on:
        fig.add_trace(
            go.Scatter(
                x = x_poles,
                y = y_poles,
                line=dict(
                    color = 'black',
                    dash = 'dash',
                ),
                name=name+'_poles',
            )
        )

        size  = [15*w for w in weights]
        color = [5 for p in poles]

        fig.add_trace(
            go.Scatter(
                x = x_poles,
                y = y_poles,
                mode = 'markers',
                marker = dict(
                    size = size,
                    color = color,
                ),
                name=name+'_poles',
            )
        )

def plot_bs_curve_2d(crv, names =None, width :int = 800, height : int = 600, ctrl_pts_on: bool = True,  crv_pts_on: bool = False, showlegend=False):


    fig = go.Figure()

    if isinstance(crv, list):
        index = 0
        for crv_ in crv:
            name = ''
            if names is not None and index < len(names):
                name = names[index]
            add_bs_curve_to_fig(crv_, fig, name,ctrl_pts_on=ctrl_pts_on, crv_pts_on=crv_pts_on)
            index+=1
    else:
        name = ''
        if names:
            name = names
        add_bs_curve_to_fig(crv, fig, name, ctrl_pts_on=ctrl_pts_on, crv_pts_on=crv_pts_on)

    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1,
    )

    fig.update_layout(
        width=width,
        height=height,
        showlegend=showlegend,
    )

    return fig

