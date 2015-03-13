#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'myourshaw'

import os
from shutil import copyfile
import sys
import argparse
import csv
from collections import namedtuple
import my

#-i /scratch1/tmp/myourshaw/mus/rnaseq/annotated/all_annotated_genes_isc_fpkm_gt_0.cdt --gene_colors /scratch1/tmp/myourshaw/mus/rnaseq/annotated/gene_annotation_colors.txt --experiment_colors  /scratch1/tmp/myourshaw/mus/rnaseq/annotated/experiment_colors.txt


class TreeViewColorsError(Exception): pass


Color = namedtuple('Color', 'fg, bg')


colormap = {
    'aliceblue': '#F0F8FF',
    'antiquewhite': '#FAEBD7',
    'aqua': '#00FFFF',
    'aquamarine': '#7FFFD4',
    'azure': '#F0FFFF',
    'beige': '#F5F5DC',
    'bisque': '#FFE4C4',
    'black': '#000000',
    'blanchedalmond': '#FFEBCD',
    'blue': '#0000FF',
    'blueviolet': '#8A2BE2',
    'brown': '#A52A2A',
    'burlywood': '#DEB887',
    'cadetblue': '#5F9EA0',
    'chartreuse': '#7FFF00',
    'chocolate': '#D2691E',
    'coral': '#FF7F50',
    'cornflowerblue': '#6495ED',
    'cornsilk': '#FFF8DC',
    'crimson': '#DC143C',
    'cyan': '#00FFFF',
    'darkblue': '#00008B',
    'darkcyan': '#008B8B',
    'darkgoldenrod': '#B8860B',
    'darkgray': '#A9A9A9',
    'darkgreen': '#006400',
    'darkkhaki': '#BDB76B',
    'darkmagenta': '#8B008B',
    'darkolivegreen': '#556B2F',
    'darkorange': '#FF8C00',
    'darkorchid': '#9932CC',
    'darkred': '#8B0000',
    'darksalmon': '#E9967A',
    'darkseagreen': '#8FBC8F',
    'darkslateblue': '#483D8B',
    'darkslategray': '#2F4F4F',
    'darkturquoise': '#00CED1',
    'darkviolet': '#9400D3',
    'deeppink': '#FF1493',
    'deepskyblue': '#00BFFF',
    'dimgray': '#696969',
    'dodgerblue': '#1E90FF',
    'firebrick': '#B22222',
    'floralwhite': '#FFFAF0',
    'forestgreen': '#228B22',
    'fuchsia': '#FF00FF',
    'gainsboro': '#DCDCDC',
    'ghostwhite': '#F8F8FF',
    'gold': '#FFD700',
    'goldenrod': '#DAA520',
    'gray': '#808080',
    'green': '#008000',
    'greenyellow': '#ADFF2F',
    'honeydew': '#F0FFF0',
    'hotpink': '#FF69B4',
    'indianred': '#CD5C5C',
    'indigo': '#4B0082',
    'ivory': '#FFFFF0',
    'khaki': '#F0E68C',
    'lavender': '#E6E6FA',
    'lavenderblush': '#FFF0F5',
    'lawngreen': '#7CFC00',
    'lemonchiffon': '#FFFACD',
    'lightblue': '#ADD8E6',
    'lightcoral': '#F08080',
    'lightcyan': '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgray': '#D3D3D3',
    'lightgreen': '#90EE90',
    'lightpink': '#FFB6C1',
    'lightsalmon': '#FFA07A',
    'lightseagreen': '#20B2AA',
    'lightskyblue': '#87CEFA',
    'lightslategray': '#778899',
    'lightsteelblue': '#B0C4DE',
    'lightyellow': '#FFFFE0',
    'lime': '#00FF00',
    'limegreen': '#32CD32',
    'linen': '#FAF0E6',
    'magenta': '#FF00FF',
    'maroon': '#800000',
    'mediumaquamarine': '#66CDAA',
    'mediumblue': '#0000CD',
    'mediumorchid': '#BA55D3',
    'mediumpurple': '#9370DB',
    'mediumseagreen': '#3CB371',
    'mediumslateblue': '#7B68EE',
    'mediumspringgreen': '#00FA9A',
    'mediumturquoise': '#48D1CC',
    'mediumvioletred': '#C71585',
    'midnightblue': '#191970',
    'mintcream': '#F5FFFA',
    'mistyrose': '#FFE4E1',
    'moccasin': '#FFE4B5',
    'navajowhite': '#FFDEAD',
    'navy': '#000080',
    'oldlace': '#FDF5E6',
    'olive': '#808000',
    'olivedrab': '#6B8E23',
    'orange': '#FFA500',
    'orangered': '#FF4500',
    'orchid': '#DA70D6',
    'palegoldenrod': '#EEE8AA',
    'palegreen': '#98FB98',
    'paleturquoise': '#AFEEEE',
    'palevioletred': '#DB7093',
    'papayawhip': '#FFEFD5',
    'peachpuff': '#FFDAB9',
    'peru': '#CD853F',
    'pink': '#FFC0CB',
    'plum': '#DDA0DD',
    'powderblue': '#B0E0E6',
    'purple': '#800080',
    'red': '#FF0000',
    'rosybrown': '#BC8F8F',
    'royalblue': '#4169E1',
    'saddlebrown': '#8B4513',
    'salmon': '#FA8072',
    'sandybrown': '#F4A460',
    'seagreen': '#2E8B57',
    'seashell': '#FFF5EE',
    'sienna': '#A0522D',
    'silver': '#C0C0C0',
    'skyblue': '#87CEEB',
    'slateblue': '#6A5ACD',
    'slategray': '#708090',
    'snow': '#FFFAFA',
    'springgreen': '#00FF7F',
    'steelblue': '#4682B4',
    'tan': '#D2B48C',
    'teal': '#008080',
    'thistle': '#D8BFD8',
    'tomato': '#FF6347',
    'turquoise': '#40E0D0',
    'violet': '#EE82EE',
    'wheat': '#F5DEB3',
    'white': '#FFFFFF',
    'whitesmoke': '#F5F5F5',
    'yellow': '#FFFF00',
    'yellowgreen': '#9ACD32',
}


def hex_color(h):
    test = h[1:] if h.startswith('#') else h
    if len(test) != 6:
        return None
    else:
        try:
            int(test, 16)
        except ValueError:
            return None
        else:
            return '#'+test.upper()


def get_color_hex(c, default=None):
    h = hex_color(c)
    if h:
        return h
    else:
        return colormap.get(c.lower(), default)


def run(input, gene_colors, experiment_colors, names_colors, fgcolor, bgcolor):

    tmp_cdt = input+'.tmp'

    if not(gene_colors or experiment_colors):
        raise TreeViewColorsError('At least one gene or experiment color file required.')

    if names_colors:
        with csv.reader(open(names_colors,'rb'), dialect=csv.excel_tab) as f:
            colormap = {c[0].lower(): hex_color(c[1]) for c in f}

    gene_color_map, experiment_color_map = {}, {}

    if gene_colors:
        with open(gene_colors,'rb') as f:
            gene_color_map = {c[0].upper(): Color(get_color_hex(c[1], fgcolor), get_color_hex(c[2], bgcolor)) for c in csv.reader(f, dialect=csv.excel_tab) if not c[0].startswith('#')}

    if experiment_colors:
        with open(experiment_colors,'rb') as f:
            experiment_color_map = {c[0].upper(): Color(get_color_hex(c[1], fgcolor), get_color_hex(c[2], bgcolor)) for c in csv.reader(f, dialect=csv.excel_tab) if not c[0].startswith('#')}

    if my.file_exists(input+'.orig'):
        copyfile(input+'.orig', input)
    else:
        copyfile(input, input+'.orig')
    with open(input,'rb') as cdt, open(tmp_cdt,'wb') as tmp:
        writer = csv.writer(tmp, dialect=csv.excel_tab)
        row_count = 0
        in_data = False
        for row in csv.reader(cdt, dialect=csv.excel_tab):
            row_count += 1
            if row_count == 1:
                if 'FGCOLOR' in row or 'BGCOLOR' in row:
                    raise TreeViewColorsError('Colors already added to {}'.format(input))
                #last annotation column
                gweight_index = row.index('GWEIGHT')
                if gene_colors:
                    row = row[:gweight_index] + ['FGCOLOR', 'BGCOLOR'] + row[gweight_index:]
                    fgcolor_index, bgcolor_index = row.index('FGCOLOR'), row.index('BGCOLOR')
                row_len = len(row)
                experiment_range = range(row.index('GWEIGHT')+1, row_len)
                experiment_dict = {i: row[i] for i in experiment_range}
            else:
                if gene_colors:
                    row = row[:gweight_index] + [fgcolor, bgcolor] + row[gweight_index:]
                if not in_data:
                    #last annotation row
                    if row[0] == 'EWEIGHT':
                        if experiment_colors:
                            fgcolor_row = ['']*row_len
                            fgcolor_row[0] = 'FGCOLOR'
                            for i in experiment_range:
                                fgcolor_row[i] = experiment_color_map.get(experiment_dict[i].upper(), Color(fgcolor,bgcolor)).fg
                            writer.writerow(fgcolor_row)
                            bgcolor_row = ['']*row_len
                            bgcolor_row[0] = 'BGCOLOR'
                            for i in experiment_range:
                                bgcolor_row[i] = experiment_color_map.get(experiment_dict[i].upper(), Color(fgcolor,bgcolor)).bg
                            writer.writerow(bgcolor_row)
                        in_data = True
                #data rows
                else:
                    row[fgcolor_index] = gene_color_map.get(row[1].upper(), Color(fgcolor,bgcolor)).fg
                    row[bgcolor_index] = gene_color_map.get(row[1].upper(), Color(fgcolor,bgcolor)).bg
            writer.writerow(row)

    os.rename(input, input+'.no_color')
    os.rename(tmp_cdt, input)
    print('Done')


def main():

    #command line arguments
    parser = argparse.ArgumentParser(
        description = 'install dbSNP files table on mysql server',
        epilog = 'vax.dbsnp_installer 1.0β1 ©2011-2014 Michael Yourshaw all rights reserved')
    #required
    parser.add_argument('--input', '-i', type=str, required=True,
        help='Cluster .cdt input file with gene symbols in coumn 2')
    parser.add_argument('--gene_colors', type=str,
        help='Tab-delimited file mapping gene symbols to gene label foreground and background colors (color names or hex RGB color codes)')
    parser.add_argument('--experiment_colors', type=str,
        help='Tab-delimited file mapping experiment ids to experiment label foreground and background colors (color names or hex RGB color codes)')
    parser.add_argument('--names_colors', '-c', type=str,
        help='Tab-delimited file mapping color names to hex codes (default 140 HTML (17) and CSS (123) color names from http://www.w3schools.com/html/html_colornames.asp')
    parser.add_argument('--fgcolor', type=str, default = '#000000',
        help='Default foreground text color (color names or hex RGB color codes)')
    parser.add_argument('--bgcolor', type=str, default = '#FFFFFF',
        help='Default background color (color names or hex RGB color codes)')
    #parse args
    args = parser.parse_args()

    run(input=args.input, gene_colors=args.gene_colors, experiment_colors=args.experiment_colors,
        names_colors=args.names_colors, fgcolor=args.fgcolor, bgcolor=args.bgcolor)

if __name__ == "__main__": sys.exit(main())
