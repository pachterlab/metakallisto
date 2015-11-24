# Module containing frequently-used plotting functions

import matplotlib
import matplotlib.pylab as pyp
import pylab


def plot_setup_pre( title = None, xaxislabel = None, yaxislabel = None,
					xticks = None, xlabels = None, xrotation = 0,
					yticks = None, ylabels = None, yrotation = 0,
					figure_number = 0):
	"""
	Sets up a figure with all the aspects that will remain constant over
	multiple plots: labels and tick marks
	"""

	# Resizes graphs so the x labels fit on the screen
	pylab.rcParams['figure.figsize'] = 20, 10

	# Creates new figure or switches to existing figure
	pyp.figure(figure_number)

	if title:
		pyp.title(title, fontsize='30')
	if xaxislabel:
		pyp.xlabel(xaxislabel)
	if yaxislabel:
		pyp.ylabel(yaxislabel, fontsize='20')
	if xticks:
		pyp.xticks(xticks, (xlabels if xlabels else xticks),
					rotation = xrotation)
	if yticks:
		pyp.yticks(yticks, (ylabels if ylabels else yticks),
					rotation = yrotation)

	return


def plot_generic(x_values, y_values, plot_func = pyp.plot, hold = True,
					figure_number = None, *args, **kargs):
	"""
	Actual plotting function that passes arguments to the correct plot method.
	"""

	if figure_number:
		pyp.figure(figure_number)

	pyp.hold(hold) # keeps existing plots on the graph if true

	plot_func(x_values, y_values, *args, **kargs)


def plot(x_values, y_values, plot_type = "line", **kargs):
	"""
	Function that gets called by user directly. Determines correct function,
	given the type of plot indicated.

	Good keyword arguments:
		color
		linestyle (dashed, dotted, etc))
		marker (o, x, +, etc)
		label
		linewidth
		width (of bars)
	"""

	lookup_plot_type = {"bar": pyp.bar,
						"scatter": pyp.scatter,
						"histogram": pyp.hist,
						"line": pyp.plot}
	try:
		plot_func = lookup_plot_type[plot_type]
	except KeyError:
		plot_func = pyp.plot

	plot_generic(x_values, y_values, plot_func, **kargs)


def plot_setup_post(figure_number = None, show = True, save_file = None,
					legend = True, legend_location = 0):
	"""
	Handles post-figure setup, including legends, file saving (save_file is
	desired filename), showing the figure, and clearing it.
	"""

	if figure_number:
		pyp.figure(figure_number)

	pyp.subplots_adjust(bottom=.5) # adjustment to give more xlabel space

	# change limits
	# pyp.xlim( xmin = 0, xmax = 10000 )
	# pyp.ylim( ymin = 0, ymax = 10000 )
	# pyp.ylim( (0,10000) ) # equivalent to the line above

	if legend:
		pyp.legend(loc = legend_location)
	if save_file:
		pyp.savefig(save_file)
	if show:
		pyp.show()
	else:
		pyp.clf() # clears figure if not plotted


