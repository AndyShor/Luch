from bokeh.plotting import figure, output_file, show
from bokeh.models import PrintfTickFormatter, HoverTool, Range1d, LabelSet, Label, WheelZoomTool, Legend
from bokeh.palettes import Category10_10 as default_palette  # import bokeh palette for
import numpy as np

def default_graph():
    bokehfig = figure(width=800, height=400, sizing_mode='scale_both',
                      tools=['pan', 'box_zoom', 'reset', 'save', 'crosshair'],
                      toolbar_location='right', y_axis_type="linear", x_axis_type="linear",
                      x_axis_label='Coordinate, [m]',
                      y_axis_label='deviation from reference track')  # basic diagram layout

    # diagram typesetting

    bokehfig.title.text = 'Trajectory'  # diagram name
    bokehfig.title.align = 'left'  # setting layout
    bokehfig.title.text_font_size = "14pt"
    bokehfig.xaxis.axis_label_text_font_size = "16pt"
    bokehfig.yaxis.axis_label_text_font_size = "16pt"
    bokehfig.yaxis.major_label_text_font_size = "14pt"
    bokehfig.xaxis.major_label_text_font_size = "14pt"

    # custom cursour for data readout
    custom_hover = HoverTool()
    custom_hover.tooltips = """
        <style>
            .bk-tooltip>div:not(:first-child) {display:none;}
        </style>

        <b>Z [m]: </b> @x <br>
        <b>x,y [m]: </b> @y
    """

    bokehfig.add_tools(custom_hover)

    # setup toolbar
    scroll = WheelZoomTool()
    bokehfig.add_tools(scroll)
    bokehfig.toolbar.active_scroll = scroll

    # setup legend
    legend = Legend(items=[])
    bokehfig.add_layout(legend)
    bokehfig.legend.location = "top_left"
    bokehfig.legend.click_policy = "mute"

    return bokehfig

def color_picker(total_items, current_item, palette):
    """ pick color for charge states"""
    if total_items < len(palette):
        return palette[current_item]
    if total_items >= len(palette):
        return palette[current_item % len(palette)]

def plot_trajectory(particle_list):
    plot=default_graph()
    for i,particle in enumerate(particle_list):
        current_color=color_picker(len(particle_list),i,default_palette)
        plot_x = np.array(particle.trajectory[:, 6]).reshape(-1, )
        plot_y = np.array(particle.trajectory[:, 0]).reshape(-1, )
        plot_y1 = -1 * np.array(particle.trajectory[:, 2]).reshape(-1, )
        plot.line(plot_x, plot_y, line_width=2, muted_alpha=0.2, color=current_color,
                  legend_label='trace '+str(i+1))  # add curve as a line to the figure
        plot.line(plot_x, plot_y1, line_width=2, muted_alpha=0.2, color=current_color,
                  legend_label='trace '+str(i+1))  # add curve as a line to the figure
    return plot