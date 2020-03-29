from bokeh.plotting import figure, output_file, show
from bokeh.models import PrintfTickFormatter, HoverTool, Range1d, LabelSet, Label, WheelZoomTool, Legend

def default_graph():
    bokehfig = figure(width=800, height=400, sizing_mode='scale_both',
                      tools=['pan', 'box_zoom', 'reset', 'save', 'crosshair'],
                      toolbar_location='right', y_axis_type="linear", x_axis_type="linear",
                      x_axis_label='Coordinate, [m]',
                      y_axis_label='deflection from reference trajectory')  # basic diagram layout

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