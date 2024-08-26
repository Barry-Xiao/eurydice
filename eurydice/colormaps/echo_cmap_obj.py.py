
import matplotlib.colors as mpc
from matplotlib import colormaps

import seaborn as sn


class CategoricalColorPalette():

    _markers = ['o', 's', '^', 'v', 'P', 'D', 'X', "<", ">"]
    
    def __init__(self, column_name, order, palette, markers=None, sizes=None, 
                 labels=None, default_marker='.'):
        """
        Puts together color palette and markers
        """
        self.hue = column_name
        self.order = order
        self.palette = self._check_cat_color(order, palette)
    
    def _check_cat_color(self, order, palette):
        """
        Verifies color palette and converts to a dictionary
        """
        if isinstance(palette, str):
            palette =  sn.color_palette(palette, n_colors=len(order))
            colors = {v: to_hex(c) for v, c in zip(*(order, palette))}
        elif isinstance(palette, (list, tuple, np.ndarray)):
            colors = {v: c for v, c in zip(*(order, palette))}
        elif isinstance(palette, dict):
            colors = palette
        else:
            raise ValueError('The color palette must be a dictionary, '
                             'list, or color palette string.')
        return palette

    def _check_markers(self, order, markers):
        if markers is None:
            markers = {o: default_marker for o in order}
        elif markers:
            markers = {o: m for o, m in zip(*(order, self._markers))}
            if len(markers) < len(order):
                raise ValueError('There are 9 default markers, and you have '
                                 'more than 9 categories.\n'
                                 'Please specify specific markers')
        elif isinstance(markers, list):
            markers = {o: m for o, m in zip(*(order, markers))}
        elif isinstance(markers, dict):
            continue
        else:
           raise ValueError('If markers are supplied, they must be a list '
                            'or dictionary') 

    def for_pcoa(self):
        """
        Returns a color palette for a PCoA
        
        Returns
        -------
        dict
        """
        poca_colors = dict(
            hue=self.hue,
            hue_order=self.order,
            palette=self.palette,
        )
        return pcoa_colors

    def for_scatter(self):
        """
        Returns a color palette for a scatter plot
        """
