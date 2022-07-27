import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2

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
parser.add_argument("-c1", type = str, default = "blue", help = "Colour for Set 1")
parser.add_argument("-c2", type = str, default = "red", help = "Colour for Set 2")
parser.add_argument("-c3", type = str, default = "green", help = "Colour for Set 3")
parser.add_argument("-ht", type = float, default = 8, help = "Figure Height (in inches)")
parser.add_argument("-w", type = float, default = 10, help = "Figure Width (in inches)")
parser.add_argument('-dpi', type = int, default = 300, help = "Output DPI")
parser.add_argument("-o", "--output", default = "venn.png", help = 'Output File')

## Make a dictionary of all arguments
args = parser.parse_args()
config = vars(args)

## Determine if only a 2way plot is required
draw2 = (not config['a3']) & (not config['a23']) & (not config['a13']) & (not config['a123'])

## Common values
areas = [config['a1'], config['a2'], config['a12']]
labels = [config['s1'], config['s2']]
## Draw 2-way if required
if draw2:
	for index, item in enumerate(areas):
		if item is None:
			areas[index] = 0

	v = venn2(areas, set_labels = labels)
	v.get_patch_by_id('10').set_color(config['c1'])
	v.get_patch_by_id('01').set_color(config['c2'])

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
	v = venn3(areas, set_labels = labels)
	v.get_patch_by_id('100').set_color(config['c1'])
	v.get_patch_by_id('010').set_color(config['c2'])
	v.get_patch_by_id('001').set_color(config['c3'])

## Modify dimensions & export
plt.gcf().set_size_inches(config['w'], config['ht'])
plt.savefig(config['output'], dpi = config['dpi'])

