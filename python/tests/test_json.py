import json
import ast
import numpy as np
from pygbs import gbs
from pytest import approx
import pathlib
from pygbs import vistaplot as gbv
import pyvista as pv

convertible_args = ['function', 'curve', 'surface', 'curve2d']
template_args    = ['dim', 'type']

def from_json(data: dict):
    
    entity_type = gbs.entity_type(data["type"])
    dim = data['dim']
    class_obj = getattr(gbs, f'{entity_type.name}{dim}d')

    # Create a new dictionary excluding 'dim' and 'type'
    constructor_args = {key: value for key, value in data.items() if key not in template_args}
    constructor_args = {key: from_json(json.loads(value)) if key in convertible_args else value 
                    for key, value in constructor_args.items()}

    return class_obj(**constructor_args)

def build_tip(side, roundness=1., p = 3):
    u1, u2, v1, v2 = side.bounds()

    le = side.isoU(u1)
    te = side.isoU(u2)
    te.reverse()

    p_begin = le.end()
    p_end = te.begin()
    t_begin = le.end(1)
    t_end = te.begin(1)
    c_begin = le.end(2)
    c_end = te.begin(2)

    u1 = 0.
    u2 = roundness

    pt1 = (u1,p_begin) # curve begin
    pt2 = (u2,p_end) # curve end
    cstr_lst = [
        (u1,t_begin,1),
        (u2,t_end,1),
        (u1,c_begin,2),
        (u2,c_end,2),
    ]
    
    return gbs.interpolate(pt1,pt2,cstr_lst,p), le, te

def unif_bounds(crv_lst):
    for crv in crv_lst[1:]:
        crv.changeBounds(crv_lst[0].bounds())

def test_from_json_file():
    
    path = str(pathlib.Path(__file__).parent.resolve()) + '/in/'
    unix_style_path = str(path).replace('\\', '/')
    file_name = unix_style_path+'prop3.json'
    with open(file_name, 'r') as f:
        data = json.load(f)

    s1 = from_json(json.loads(data['s1']))
    assert isinstance(s1, gbs.BSSurface3d)
    s2 = from_json(json.loads(data['s2']))
    assert isinstance(s2, gbs.BSSurface3d)

    tip1, le1, te1 = build_tip(s1, roundness=0.25)
    tip2, le2, te2 = build_tip(s2, roundness=0.25)

    j1 = gbs.join(gbs.join(le1, tip1), te1)
    j2 = gbs.join(gbs.join(le2, tip2), te2)

    le = from_json(json.loads(data['le']))
    te = from_json(json.loads(data['te']))

    spans =  ast.literal_eval(data['spans'])

    le_lst = [le.isoV(s) for s in spans]
    s1_lst = [s1.isoV(s) for s in spans]
    te_lst = [te.isoV(s) for s in spans]
    s2_lst = [s2.isoV(s) for s in spans]

    u_crv_lst = [le_lst[-1], te_lst[-1]]
    v_crv_lst = [tip1, tip2]
    unif_bounds(u_crv_lst)
    unif_bounds(v_crv_lst)
    T = gbs.gordonbs(u_crv_lst, v_crv_lst)

    T1 = gbs.loftbs([s1_lst[-1], tip1],2)
    T2 = gbs.loftbs([s2_lst[-1], tip2],2)

    # u_crv_lst = le_lst+te_lst[::-1]
    # v_crv_lst = [j1, j2]
    # unif_bounds(u_crv_lst)
    # unif_bounds(v_crv_lst)
    # T2 = gbs.gordonbs(u_crv_lst, v_crv_lst)

    # u_crv_lst = s1_lst+[tip1]
    # te1.reverse()
    # v_crv_lst = [le1, te1]
    # unif_bounds(u_crv_lst)
    # unif_bounds(v_crv_lst)
    # S1 = gbs.gordonbs(u_crv_lst, v_crv_lst)

    plotter = pv.Plotter()

    tip_lst = [
            le1,
            te1,
            tip1,
            le2,
            te2,
            tip2,
            # j1,
            # j2,
        ]

    gbv.add_surfaces_to_plotter([s1, s2, le, te], plotter, per=13)
    gbv.add_surfaces_to_plotter([T, T1, T2], plotter, per=13)
    gbv.add_curves_to_plotter(s1_lst+s2_lst+ le_lst+ te_lst+tip_lst, plotter)
    # gbv.add_curves_to_plotter([j1, j2], plotter, per=13)
    plotter.show()
