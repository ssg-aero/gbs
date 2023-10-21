import numpy as np
from pygbs import gbs
import ipywidgets as widgets


class BsCurveEditor:
    def __init__(self, crv, ranges = None) -> None:
        self.curve = crv
        self.poles = crv.poles()
        self.mults = crv.mults()
        self.knots = crv.knots()
        self.degree= crv.degree()

        try:
            l = gbs.length(crv)
        except:
            l = 1.
        

        self.dim = len(self.poles[0])

        self.dim_labels = ['x', 'y', 'z', 'w']

        rational = isinstance(crv, gbs.BSCurveRational1d) or isinstance(crv, gbs.BSCurveRational2d) or isinstance(crv, gbs.BSCurveRational3d)

        if rational:
            self.ctors = [
                gbs.BSCurveRational1d,
                gbs.BSCurveRational2d,
                gbs.BSCurveRational3d,
            ]
        else:
            self.ctors = [
                gbs.BSCurve1d,
                gbs.BSCurve2d,
                gbs.BSCurve3d,
            ]
        
        if ranges is not None:
            self.min = [ r[0] for r in ranges]
            self.max = [ r[1] for r in ranges]
        else:
            self.min = np.ones(self.dim) *-l
            self.max = np.ones(self.dim) * l

        self._observers = []

    def register_observer(self, observer):
        if observer not in self._observers:
            self._observers.append(observer)

    def remove_observer(self, observer):
        self._observers.remove(observer)

    def notify_observers(self):
        for observer in self._observers:
            observer.update()

    def build_curve(self):
        self.curve = self.ctors[self.dim-1](self.poles, self.knots, self.mults, self.degree)
        self.notify_observers()

    def build_on_pole_change(self, i, j):
        def on_pole_change(change):
            if change['type'] == 'change' and change['name'] == 'value':
                self.poles[i][j] = change['new']
                self.build_curve()

        return on_pole_change
    
    def build_pole_slider(self, i, d):
        slider =  widgets.FloatSlider( 
            value=self.poles[i][d],
            min=self.min[d], max=self.max[d], step=(self.max[d]-self.min[d])/500,
            description=f'{self.dim_labels[d]}{i}:', continuous_update=True,
            orientation='horizontal',
        )
        slider.observe(self.build_on_pole_change(i, d))
        return slider

    def build_pole_sliders(self, i):

        return [self.build_pole_slider(i, d) for d in range(self.dim)]
    
    def build_poles_sliders(self):
        sliders = []
        for i in range(len(self.poles)):
            sliders += self.build_pole_sliders(i)
        return sliders