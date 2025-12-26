import plot
import composite
from composite import Composite
from composite import SimpleComposite
import re
import xml.dom.minidom as minidom
import argparse
import math
import parseComposite
import sys
from enum import Enum
import svgFactory
import json
import os
from io import BytesIO

document = minidom.Document()

def main():
    return True

if __name__ == "__main__":
    # Remove 'plotter' from sys.argv
    sys.argv.pop(0)
    # Load subcommands into appropriate arrays
    i = -1
    k = -1
    composite_commands = []
    ref_line_commands = []
    plot_command = ""
    current = ""
    for word in sys.argv:
        if word == "composite":
            i += 1
            composite_commands.append("")
            current = "composite"
        elif word == "reference-line":
            k += 1
            ref_line_commands.append("")
            current = "ref"
        elif word == "plot":
            current = "plot"
        elif current == "composite":
            composite_commands[i] += f" {word}"
        elif current == "ref":
            ref_line_commands[k] += f" {word}"
        elif current == "plot":
            plot_command += f" {word}"
    # Create parser for plot subcommand
    plot_parser = argparse.ArgumentParser()
    plot_parser.add_argument("--smoothing", type=int)
    plot_parser.add_argument("--bp-shift", type=int)
    plot_parser.add_argument("--opacity", type=float)
    plot_parser.add_argument("--title", nargs="+")
    plot_parser.add_argument("--xmin",type=int)
    plot_parser.add_argument("--xmax",type=int)
    plot_parser.add_argument("--xlabel", nargs="+")
    plot_parser.add_argument("--ymin", type=int)
    plot_parser.add_argument("--ymax", type=int)
    plot_parser.add_argument("--ylabel", nargs="+")
    plot_parser.add_argument("--color-trace", action="store_true", default=False)
    plot_parser.add_argument("--combined", action="store_true", default=False)
    plot_parser.add_argument("--hide-legend", action="store_true", default=False)
    plot_parser.add_argument("--no-resize", action="store_true", default=False)
    plot_parser.add_argument("--no-shrink", action="store_true", default=False)
    plot_parser.add_argument("--out")
    plot_parser.add_argument("--export-json")
    plot_parser.add_argument("--import-json")
    plot_parser.add_argument("--import-settings-json")

    # Create plot based on plot subcommand, default values in Plot class will be used if argument is not specified
    plot_args = plot_parser.parse_args(plot_command.split())
    p = plot.Plot(title=" ".join(plot_args.title) if plot_args.title is not None else None, xmin=plot_args.xmin, xmax=plot_args.xmax, ymin=plot_args.ymin, ymax=plot_args.ymax, xlabel=" ".join(plot_args.xlabel) if plot_args.xlabel is not None else None, 
                  ylabel=" ".join(plot_args.ylabel) if plot_args.ylabel is not None else None, opacity=plot_args.opacity, smoothing=plot_args.smoothing, bp_shift=plot_args.bp_shift, combined=plot_args.combined, color_trace=plot_args.color_trace, hide_legend=plot_args.hide_legend)

    # Create arrays for default composite names and colors
    names = range(1, len(composite_commands) + 1)
    colors = ["#FF0000","#FF9100","#D7D700","#07E200","#00B0F0","#0007FF","#A700FF","#FF00D0","#BFBFBF","#000000"]
    # Create parser for composite subcommands
    composite_parser = argparse.ArgumentParser()
    composite_parser.add_argument("files")
    composite_parser.add_argument("--name")
    composite_parser.add_argument("--color")
    composite_parser.add_argument("--secondary-color")
    composite_parser.add_argument("--scale", type=float)
    composite_parser.add_argument("--shift-occupancy", type=float)
    composite_parser.add_argument("--smoothing", type=int)
    composite_parser.add_argument("--opacity", type=float)
    composite_parser.add_argument("--bp-shift", type=int)
    composite_parser.add_argument("--hide-sense", action="store_true", default=False)
    composite_parser.add_argument("--hide-anti", action="store_true", default=False)
    composite_parser.add_argument("--swap-strands", action="store_true", default=False)
    # Parse composite subcommands, use values values in Composite class if not specified 
    i = 0
    for command in composite_commands:
        args = composite_parser.parse_args(command.split())
        composite = Composite(scale=args.scale, color=args.color if args.color is not None else colors[i % len(colors)], secondary_color=args.secondary_color, 
                                         smoothing=args.smoothing, bp_shift=args.bp_shift, hide_sense= args.hide_sense, hide_anti= args.hide_anti, baseline=args.shift_occupancy,
                                         name=args.name if args.name is not None else names[i], opacity=args.opacity,)
        
        composite_files = args.files.split(":")
        for c in composite_files:
            #Check if composite file contains multiple composites
            if sum(1 for line in open(c) if len(line.strip()) != 0) <= 3:
                sc = parseComposite.parse_simple(c)
                composite.load_simple_composite(sc)
            else:
                prefixes = parseComposite.get_prefixes_from_multiple_composites(c)
                cd = parseComposite.parse_multiple_composite(c, prefixes[0])
                composite.load_composite_dict(cd)        
        p.add_composite_group(composite)
        i += 1    
 
    # Import settings and composites from plot, preserving options specified in this call
    if plot_args.import_json:
        p.import_data(plot_args.import_json, plot_args, True)
    elif plot_args.import_settings_json:
        p.import_data(plot_args.import_settings_json, plot_args, False)

    # If --no-shrink is specified, don't change y-axis but resize x-axis
    if plot_args.no_shrink:
        p.autoscale_axes(False)
    # If --no-resize is specified, don't change either axis
    elif not plot_args.no_resize:
        p.autoscale_axes(True)

    p.plot_composites()

    # Create parser for reference-line subcommand
    reference_parser = argparse.ArgumentParser()
    reference_parser.add_argument("axis")
    reference_parser.add_argument("--style")
    reference_parser.add_argument("--color")
    reference_parser.add_argument("--val", type=float)
    reference_parser.add_argument("--opacity",type=float)
    # Add reference lines to plot
    for command in ref_line_commands:
        args = reference_parser.parse_args(command.split())
        p.plot_reference_line(axis=args.axis, val=args.val, style=args.style, color=args.color, opacity=args.opacity)

    # produce svg DOM as before
    svg = svgFactory.generateSVG(p)
    svg_elem = svg.documentElement if hasattr(svg, "documentElement") else svg

    def _parse_viewbox(vb):
        parts = re.split(r'[\s,]+', vb.strip())
        if len(parts) == 4:
            return [float(x) for x in parts]
        return None

    def _parse_length_attr(s, fallback):
        if not s:
            return fallback
        m = re.match(r'([0-9.]+)', s)
        return float(m.group(1)) if m else fallback

    # Obtain inner content bounds (minx, miny, width, height)
    vb_attr = svg_elem.getAttribute('viewBox')
    if vb_attr:
        minx, miny, inner_w, inner_h = _parse_viewbox(vb_attr)
    else:
        minx = 0.0
        miny = 0.0
        inner_w = _parse_length_attr(svg_elem.getAttribute('width'), 800.0)
        inner_h = _parse_length_attr(svg_elem.getAttribute('height'), 600.0)

    # USER TUNABLE: how much bigger you want the figure
    scale_factor = 1.3    # 1.0 = same size, 1.2 = 20% larger
    # USER TUNABLE: margin fraction relative to scaled content size
    margin_frac = 0.1   # 8% margin around scaled content

    # compute scaled content size and margins
    scaled_w = inner_w * scale_factor
    scaled_h = inner_h * scale_factor
    margin_x = scaled_w * margin_frac
    margin_y = scaled_h * margin_frac

    # new canvas dims (viewBox internal units)
    new_w = scaled_w + 2 * margin_x
    new_h = scaled_h + 2 * margin_y

    # compute tx, ty such that after "translate(tx,ty) scale(s)", the old minx,miny
    # maps to the margin: (minx + tx) * s = margin_x  => tx = margin_x/s - minx
    s = scale_factor
    tx = margin_x / s - minx
    ty = margin_y / s - miny

    # create new svg document wrapper
    newdoc = minidom.Document()
    newsvg = newdoc.createElement('svg')
    newsvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
    if svg_elem.getAttribute('xmlns:xlink'):
        newsvg.setAttribute('xmlns:xlink', svg_elem.getAttribute('xmlns:xlink'))

    # set pixel width/height to integer values (optional), keep precise viewBox
    newsvg.setAttribute('width', str(int(new_w)))
    newsvg.setAttribute('height', str(int(new_h)))
    newsvg.setAttribute('viewBox', f"0 0 {new_w} {new_h}")
    newsvg.setAttribute('preserveAspectRatio', 'xMidYMid meet')

    # group: translate then scale (so scaled content's top-left is at margin_x,margin_y)
    g = newdoc.createElement('g')
    g.setAttribute('transform', f"translate({tx},{ty}) scale({s})")

    # import existing children into the scaled group
    for node in list(svg_elem.childNodes):
        g.appendChild(newdoc.importNode(node, deep=True))

    newsvg.appendChild(g)
    newdoc.appendChild(newsvg)

    # write result
    out_path = plot_args.out if getattr(plot_args, "out", None) else "out.svg"
    out_dir = os.path.dirname(out_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        newdoc.writexml(f, addindent='  ', newl='\n')

    # JSON export unchanged
    if plot_args.export_json:
        dump_str = json.dumps(p.export(), indent=2)
        with open(plot_args.export_json, 'w', encoding='utf-8') as f:
            f.write(dump_str)


    