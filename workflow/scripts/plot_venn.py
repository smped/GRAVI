import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles

# from numpy import block

parser = argparse.ArgumentParser(
	description="Generate Venn Diagrams for 2 or 3 groups",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-a1", type = int, help = "Area 1")
parser.add_argument("-a2", type = int, help = "Area 2")
parser.add_argument("-a3", type = int, help = "Area 3")
parser.add_argument("-a12", type = int, help = "Intersection between areas 1 & 2")
parser.add_argument("-a13", type = int, help = "Intersection between areas 1 & 3")
parser.add_argument("-a23", type = int, help = "Intersection between areas 2 & 3")
parser.add_argument("-a123", type = int, help = "Intersection between areas 1, 2 & 3")
parser.add_argument("-s1", type = str, default = "X", help = "Label for Set 1")
parser.add_argument("-s2", type = str, default = "Y", help = "Label for Set 2")
parser.add_argument("-s3", type = str, default = "Z", help = "Label for Set 3")
parser.add_argument("-c1", type = str, default = "white", help = "Colour for Set 1")
parser.add_argument("-c2", type = str, default = "white", help = "Colour for Set 2")
parser.add_argument("-c3", type = str, default = "white", help = "Colour for Set 3")
parser.add_argument("-al", "--alpha", type=float, default = 0.3, help = "Alpha for fills")
parser.add_argument("-lc", type = str, default = "#4d4d4d", help = "Line Colour")
parser.add_argument("-lw", type = float, default = 1, help = "Line Width")
parser.add_argument("-ht", type = float, default = 8, help = "Figure Height (in inches)")
parser.add_argument("-w", type = float, default = 10, help = "Figure Width (in inches)")
parser.add_argument('-dpi', type = int, default = 300, help = "Output DPI")
parser.add_argument('-ff', type = str, default='serif', help = 'Font Family')
parser.add_argument('-fs', type = int, default=12, help = "Font Size")
parser.add_argument('--hide-areas', action='store_true', help = 'Hide numbers for all areas')
parser.add_argument('-sf', type = float, default=1.05, help="Scale factor for shifting 2way intersection labels. Set to 1 for default positions")
parser.add_argument("-o", "--output", default = "venn.png", help = 'Output File')

## Make a dictionary of all arguments
args = parser.parse_args()
config = vars(args)

## Set global plotting parameters
plt.rcParams["font.family"] = config['ff']
plt.rcParams["font.size"] = config['fs']

## Determine if only a 2way plot is required
draw2 = (not config['a3']) & (not config['a23']) & (not config['a13']) & (not config['a123'])

## Common values
areas = [config['a1'], config['a2'], config['a12']]
labels = [config['s1'], config['s2']]
colours = [config['c1'], config['c2']]
## Draw 2-way if required
if draw2:
	for index, item in enumerate(areas):
		if item is None:
			areas[index] = 0

	v=venn2(areas, set_labels = labels, set_colors=colours, alpha=config['alpha'])
	for ID in ["10", "01", "11"]:
		if v.get_label_by_id(ID).get_text() == "0":
			v.get_label_by_id(ID).set_text("")
		if config['hide_areas']:
			v.get_label_by_id(ID).set_text("")
	
	venn2_circles(subsets = areas, color=config['lc'], linewidth=config['lw'])

## Draw 3-way if required
if not draw2:
	areas.append(config['a3'])
	areas.append(config['a13'])
	areas.append(config['a23'])
	areas.append(config['a123'])
	for index, item in enumerate(areas):
		if item is None:
			areas[index] = 0
	
	labels.append(config['s3'])
	colours.append(config['c3'])
	v=venn3(areas, set_labels = labels, set_colors=colours, alpha=config['alpha'])
	if config['hide_areas']:
		for ID in ["100", "010", "110", "001", "101", "011", "111"]:
			v.get_label_by_id(ID).set_text("")
	
	if not config['hide_areas']:
		sf = config['sf']
		x0,y0 = v.get_label_by_id('111').get_position()
		for ID in ['110', '101', '011']:
		  pos = v.get_label_by_id(ID)
		  if bool(pos):
			  x1,y1 =pos.get_position()
			  v.get_label_by_id(ID).set_position((x0 + sf*(x1-x0), y0 + sf*(y1-y0)))

	venn3_circles(subsets = areas,color=config['lc'], linewidth=config['lw'])

## Modify dimensions & export
plt.gcf().set_size_inches(config['w'], config['ht'])
plt.gcf().tight_layout()
plt.savefig(config['output'], dpi = config['dpi'])

