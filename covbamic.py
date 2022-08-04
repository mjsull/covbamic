import argparse, sys, os
import pysam
from collections import defaultdict

def get_depth(pos, samfile):
    basefreq = defaultdict(lambda: 0)
    for pileupcolumn in samfile.pileup("MN908947.3", pos):
        if pileupcolumn.pos != pos:
            continue
        depth = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                basefreq["-"] += 1
                depth += 1
            elif pileupread.is_refskip:
                pass
            else:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                base = base.upper()
                basefreq[base] += 1
                depth += 1
    return(basefreq)



def get_sites(variantA, variantB, variant_file, diff_only=True):
    poslist = []
    with open(variant_file) as f:
        header = f.readline().rstrip().split("\t")
        aa_col, pos_col, a_col, b_col, a_col_pres, b_col_pres = None, None, None, None, None, None
        for num, i in enumerate(header):
            if i == "Nucleotide position":
                pos_col = num
            elif i == "aa_SNP":
                aa_col = num
            elif i == variantA:
                a_col = num
            elif i == variantB:
                b_col = num
            elif i == variantA + "_present":
                a_col_pres = num
            elif i == variantB + "_present":
                b_col_pres = num
        for line in f:
            splitline = line.rstrip().split("\t")
            a_pres, b_pres, pos, a, b, aa = splitline[a_col_pres], splitline[b_col_pres], splitline[pos_col], \
                                            splitline[a_col], splitline[b_col], splitline[aa_col]
            pos = int(pos)
            a = a.upper()
            b = b.upper()
            if a_pres == "1" or b_pres == "1":
                if not diff_only or a != b:
                    poslist.append([aa, pos, a, b])
    return(poslist)

def colorstr(rgb): return "#%02x%02x%02x" % (rgb[0],rgb[1],rgb[2])

class scalableVectorGraphics:

    def __init__(self, height, width):
        self.height = height
        self.width = width
        self.out = '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <svg
       xmlns:dc="http://purl.org/dc/elements/1.1/"
       xmlns:cc="http://creativecommons.org/ns#"
       xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
       xmlns:svg="http://www.w3.org/2000/svg"
       xmlns="http://www.w3.org/2000/svg"
       xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
       xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
       height="%f"
       width="%f"
       id="svg2"
       version="1.1"
       inkscape:version="0.48.4 r9939"
       sodipodi:docname="easyfig">
      <metadata
         id="metadata122">
        <rdf:RDF>
          <cc:Work
             rdf:about="">
            <dc:format>image/svg+xml</dc:format>
            <dc:type
               rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
            <dc:title>Easyfig</dc:title>
          </cc:Work>
        </rdf:RDF>
      </metadata>
      <defs
         id="defs120" />
      <sodipodi:namedview
         pagecolor="#ffffff"
         bordercolor="#666666"
         borderopacity="1"
         objecttolerance="10"
         gridtolerance="10"
         guidetolerance="10"
         inkscape:pageopacity="0"
         inkscape:pageshadow="2"
         inkscape:window-width="640"
         inkscape:window-height="480"
         id="namedview118"
         showgrid="false"
         inkscape:zoom="0.0584"
         inkscape:cx="2500"
         inkscape:cy="75.5"
         inkscape:window-x="55"
         inkscape:window-y="34"
         inkscape:window-maximized="0"
         inkscape:current-layer="svg2" />
      <title
         id="title4">Easyfig</title>
      <g
         style="fill-opacity:1.0; stroke:black; stroke-width:1;"
         id="g6">''' % (self.height, self.width)

    def drawLine(self, x1, y1, x2, y2, th=1, cl=(0, 0, 0), alpha = 1.0):
        self.out += '  <line x1="%f" y1="%f" x2="%f" y2="%f"\n        stroke-width="%f" stroke="%s" stroke-opacity="%f" stroke-linecap="round" />\n' % (x1, y1, x2, y2, th, colorstr(cl), alpha)

    def drawPath(self, xcoords, ycoords, th=1, cl=(0, 0, 0), alpha=0.9):
        self.out += '  <path d="M%f %f' % (xcoords[0], ycoords[0])
        for i in range(1, len(xcoords)):
            self.out += ' L%f %f' % (xcoords[i], ycoords[i])
        self.out += '"\n        stroke-width="%f" stroke="%s" stroke-opacity="%f" stroke-linecap="butt" fill="none" z="-1" />\n' % (th, colorstr(cl), alpha)


    def writesvg(self, filename):
        outfile = open(filename, 'w')
        outfile.write(self.out + ' </g>\n</svg>')
        outfile.close()

    def drawRightArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + wid - ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%f"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f" />\n' % (x, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x1, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%f"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%f,%f %f,%f %f,%f" />\n' % (x, y, x, y+ht, x + wid, y1)

    def drawLeftArrow(self, x, y, wid, ht, fc, oc=(0,0,0), lt=1):
        if lt > ht /2:
            lt = ht/2
        x1 = x + wid
        y1 = y + ht/2
        x2 = x + ht / 2
        ht -= 1
        if wid > ht/2:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%f"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%f,%f %f,%f %f,%f %f,%f %f,%f %f,%f %f,%f" />\n' % (x1, y+ht/4, x2, y+ht/4,
                                                                                                x2, y, x, y1, x2, y+ht,
                                                                                                x2, y+3*ht/4, x1, y+3*ht/4)
        else:
            self.out += '  <polygon fill="%s" stroke="%s" stroke-width="%f"\n' % (colorstr(fc), colorstr(oc), lt)
            self.out += '           points="%f,%f %f,%f %f,%f" />\n' % (x, y1, x1, y+ht, x1, y)


    def drawOutRect(self, x1, y1, wid, hei, fill=(255, 255, 255), outfill=(0, 0, 0), lt=1, alpha=1.0, alpha2=1.0):
        self.out += '  <rect stroke="%s" stroke-width="%f" stroke-opacity="%f"\n' % (colorstr(outfill), lt, alpha)
        self.out += '        fill="%s" fill-opacity="%f"\n' % (colorstr(fill), alpha2)
        self.out += '        x="%f" y="%f" width="%f" height="%f" />\n' % (x1, y1, wid, hei)

    def create_pattern(self, id, fill, pattern, width, line_width):
        fill = colorstr(fill)
        if pattern == 'horizontal':
            self.out += '    <defs>\n'
            self.out += '      <pattern id="%s" width="%d" height="%d" patternUnits="userSpaceOnUse">\n' % (id, width, width)
            self.out += '        <line x1="0" y1="0" x2="%d" y2="0" style="stroke:%s; stroke-width:%d; fill:#FFFFFF" />\n' % (width, fill, line_width)
            self.out += '      </pattern>\n'
            self.out += '   </defs>\n'
        elif pattern == 'forward_diag':
            self.out += '    <defs>\n'
            self.out += '      <pattern id="%s" width="%d" height="%d" patternTransform="rotate(45 0 0)" patternUnits="userSpaceOnUse">\n' % (id, width, width)
            self.out += '        <line x1="0" y1="0" x2="0" y2="%d" style="stroke:%s; stroke-width:%d; fill:#FFFFFF" />\n' % (width, fill, line_width)
            self.out += '      </pattern>\n'
            self.out += '   </defs>\n'
        elif pattern == 'reverse_diag':
            self.out += '    <defs>\n'
            self.out += '      <pattern id="%s" width="%d" height="%d" patternTransform="rotate(135 0 0)" patternUnits="userSpaceOnUse">\n' % (id, width, width)
            self.out += '        <line x1="0" y1="0" x2="0" y2="%d" style="stroke:%s; stroke-width:%d; fill:#FFFFFF" />\n' % (width, fill, line_width)
            self.out += '      </pattern>\n'
            self.out += '   </defs>\n'

    # create a rectangle with a patterned fill. Pattern needs to be created with  create_pattern first
    def drawPatternRect(self, x, y, width, height, id, fill, title="none", webpage="none", lt=1):
        fill = colorstr(fill)
        self.out += '  <rect style="fill:#FFFFFF; stroke: %s; stroke-alignment: inner;"\n' % fill
        self.out += '        x="%d" y="%d" width="%d" height="%d"  />\n' % (x, y, width, height)
        self.out += '  <a href="%s">' % webpage
        self.out += '  <rect style="fill:url(#%s); stroke: %s; stroke-alignment: inner;"\n' % (id, fill)
        self.out += '        x="%d" y="%d" width="%d" height="%d"> <title>%s</title> </rect>\n' % (x, y, width, height, title)
        self.out += '  </a>'


    def writeString(self, thestring, x, y, size, ital=False, bold=False, rotate=0, justify='left', color=(0,0,0)):
        if rotate != 0:
            x, y = y, x
        self.out += '  <text\n'
        self.out += '    style="font-size:%fpx;font-style:normal;font-weight:normal;z-index:10\
;line-height:125%%;letter-spacing:0px;word-spacing:0px;fill:%s;fill-opacity:1;stroke:none;font-family:sans-serif"\n' % (size, colorstr(color))
        if justify == 'right':
            self.out += '    text-anchor="end"\n'
        elif justify == 'middle':
            self.out += '    text-anchor="middle"\n'
        if rotate == 1:
            self.out += '    x="-%f"\n' % x
        else:
            self.out += '    x="%f"\n' % x
        if rotate == -1:
            self.out += '    y="-%f"\n' % y
        else:
            self.out += '    y="%f"\n' % y
        self.out += '    sodipodi:linespacing="125%"'
        if rotate == -1:
            self.out += '\n    transform="matrix(0,1,-1,0,0,0)"'
        if rotate == 1:
            self.out += '\n    transform="matrix(0,-1,1,0,0,0)"'
        self.out += '><tspan\n      sodipodi:role="line"\n'
        if rotate == 1:
            self.out += '      x="-%f"\n' % x
        else:
            self.out += '      x="%f"\n' % x
        if rotate == -1:
            self.out += '      y="-%f"' % y
        else:
            self.out += '      y="%f"' % y
        if ital and bold:
            self.out += '\nstyle="font-style:italic;font-weight:bold"'
        elif ital:
            self.out += '\nstyle="font-style:italic"'
        elif bold:
            self.out += '\nstyle="font-style:normal;font-weight:bold"'
        self.out += '>' + thestring + '</tspan></text>\n'


def draw_output(positions, proportion, depths, names, bases, output_file, varA, varB, refseq, panelC):
    left_buffer = 5
    page_width = max([210, len(positions)*10])
    if panelC:
        page_height = 600
    else:
        page_height = 320
    svg = scalableVectorGraphics(page_height, page_width)
    right_buffer = 40
    width = page_width - left_buffer - right_buffer
    y = 20
    y2 = 50
    height = 30
    line_width = 6
    svg.drawLine(left_buffer, y, left_buffer+width, y, line_width)
    length = 29903
    genes = [
		["ORF1a", 266, 13468],
        ["ORF1b", 13468, 21555],
        ["spike", 21563, 25384],
        ["ORF3a", 25393, 26220],
        ["E", 26245, 26472],
        ["M", 26523, 27191],
        ["ORF6", 27202, 27387],
        ["ORF7a", 27394, 27759],
		["ORF8", 27894, 28259],
        ["N", 28274, 29533],
        ["ORF10", 29558,29674]
    ]
    for i in genes:
        name, start, stop = i
        x = left_buffer + start/length*width
        gene_width = (stop-start)/length*width
        svg.drawRightArrow(x, y-height/2+0.5, gene_width, height, (3, 166, 41))
        if gene_width > 10:
            svg.writeString(name, x+2, y+8/3, 8, color=(255, 255, 255))
    spacer = 5
    column_width = (width - (len(positions)-1)*spacer)/len(positions)
    colors = [(127,163,74),
              (149, 115, 198),
              (153, 86, 107),
              (153, 153, 153)]
    colors2 = [(153,86,107),
               (127,163,74),
               (149,115,198),
               (204,109,65),
               (153, 153, 153)]

    proportion_height = 100

    y3 = proportion_height + y2 + height + 1
    y4 = y3 + 30
    y5 = y4+130
    if len(positions) > 15:
        font_size = 4
    elif len(positions) > 10:
        font_size = 5
    else:
        font_size = 6
    for num, i in enumerate(positions):
        print(num, i)
        x1 = left_buffer + i/length*width
        x2 = left_buffer + num/len(positions) * width + column_width/2
        svg.drawPath([x1, x1, x2, x2], [y-height/2, y+height/2, y2, y2+height/2])
        svg.writeString(str(i), x2-font_size/3, y2+height/2+1, font_size, rotate=-1)# justify="middle")
        prop_y = y2 + height
        for num2, j in enumerate(proportion[num]):
            col_height = j * proportion_height
            svg.drawOutRect(x2-column_width/2, prop_y, column_width, col_height, colors[num2], lt=0)
            prop_y += col_height
        svg.writeString(names[num], x2-font_size/3, y3, font_size, rotate=-1)# justify="middle")
        col_height = min([depths[num]/1000 * 100, 100])
        svg.drawOutRect(x2-column_width/2, y4, column_width, col_height, (52, 116, 235), lt=0)
        svg.writeString(str(depths[num]) + "x", x2-font_size/3, y4+col_height+1, font_size, rotate=-1)#, justify="middle")
        if panelC:
            svg.writeString(refseq[i - 1].upper(), x2, y5-1, font_size, justify="middle")
            prop_y = y5
            for num2, j in enumerate(bases[num]):
                col_height = j * proportion_height
                svg.drawOutRect(x2 - column_width / 2, prop_y, column_width, col_height, colors2[num2], lt=0)
                prop_y += col_height
    if panelC:
        svg.writeString("Ref base", left_buffer+width+2, y5-1, font_size)
        svg.writeString("Legend", page_width - right_buffer, y5 + height, 10)
        svg.drawOutRect(page_width - right_buffer, y5 + height + 6, 10, 10, colors2[0], lt=0)
        svg.writeString("A", page_width - right_buffer + 11, y5 + height + 14, 8)
        svg.drawOutRect(page_width - right_buffer, y5 + height + 18, 10, 10, colors2[1], lt=0)
        svg.writeString("T", page_width - right_buffer + 11, y5 + height + 26, 8)
        svg.drawOutRect(page_width - right_buffer, y5 + height + 30, 10, 10, colors2[2], lt=0)
        svg.writeString("C", page_width - right_buffer + 11, y5 + height + 38, 8)
        svg.drawOutRect(page_width - right_buffer, y5 + height + 42, 10, 10, colors2[3], lt=0)
        svg.writeString("G", page_width - right_buffer + 11, y5 + height + 50, 8)
        svg.drawOutRect(page_width - right_buffer, y5 + height + 54, 10, 10, colors2[4], lt=0)
        svg.writeString("-", page_width - right_buffer + 11, y5 + height + 62, 8)

    svg.drawLine(left_buffer, y4, left_buffer+width, y4)
    svg.writeString("0x", left_buffer+width+2, y4, 6)
    svg.drawLine(left_buffer, y4+10, left_buffer + width, y4+10)
    svg.writeString("100x", left_buffer+width+2, y4+10, 6)
    svg.drawLine(left_buffer, y4+100, left_buffer + width, y4+100)
    svg.writeString("1000x", left_buffer+width+2, y4+100, 6)
    svg.writeString("Depth", page_width-right_buffer, y4 + 50, 10, justify="middle", rotate=-1)
    svg.writeString("Legend", page_width-right_buffer, y2+height, 10)
    svg.drawOutRect(page_width-right_buffer, y2+height+6, 10, 10, colors[0], lt=0)
    svg.writeString(varA, page_width-right_buffer+11, y2+height+14, 8)
    svg.drawOutRect(page_width-right_buffer, y2+height+18, 10, 10, colors[1], lt=0)
    svg.writeString(varB, page_width-right_buffer+11, y2+height+26, 8)
    svg.drawOutRect(page_width-right_buffer, y2+height+30, 10, 10, colors[2], lt=0)
    svg.writeString("both", page_width-right_buffer+11, y2+height+38, 8)
    svg.drawOutRect(page_width-right_buffer, y2+height+42, 10, 10, colors[3], lt=0)
    svg.writeString("other", page_width-right_buffer+11, y2+height+50, 8)





    svg.writesvg(output_file)


def get_minor(sites, alignment, refseq, all_minor_fraction, all_minor_depth):
    site_dict = {}
    for i in sites:
        site_dict[i[1]] = i
    sites = []
    for pileupcolumn in alignment.pileup("MN908947.3"):
        pos = pileupcolumn.reference_pos
        ref_base = refseq[pos].upper()

        basefreq = defaultdict(lambda: 0)
        depth = 0
        if pos in site_dict:
            sites.append(site_dict[pos])
        else:
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    basefreq["-"] += 1
                    depth += 1
                elif pileupread.is_refskip:
                    pass
                else:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    base = base.upper()
                    basefreq[base] += 1
                    depth += 1
            if depth < all_minor_depth:
                continue
            good = 0
            for i in basefreq:
                if basefreq[i] / depth >= all_minor_fraction:
                    good += 1
            if good >= 2:
                sites.append(["ref", pos+1, ref_base, ref_base])
    return(sites)


def __main__(bam_file, variantA, variantB, output_file, all_variants, all_minor, all_minor_fraction, all_minor_cov, panel3):
    dirname = os.path.dirname(__file__)
    variant_file = os.path.join(dirname, 'data', "variants.tsv")
    ref = os.path.join(os.path.dirname(__file__), 'data', "nCoV-2019.reference.fasta")
    refseq = ""
    with open(ref) as f:
        for line in f:
            if not line.startswith(">"):
                refseq += line.rstrip()
    sites = get_sites(variantA, variantB, variant_file, not all_variants)
    alignment = pysam.AlignmentFile(bam_file, "rb")
    if all_minor:
        sites = get_minor(sites, alignment, refseq, all_minor_fraction, all_minor_cov)
    proportion = []
    depths = []
    names = []
    positions = []
    bases = []
    for i in sites:
        aa, pos, varA, varB = i
        basefreq = get_depth(pos-1, alignment)
        varA_count, varB_count, both, other = 0, 0, 0, 0
        total = 0
        for i in basefreq:
            if i == varA == varB:
                both = basefreq[i]
            elif i == varA:
                varA_count = basefreq[i]
            elif i == varB:
                varB_count = basefreq[i]
            else:
                other += basefreq[i]
            total += basefreq[i]
        proportion.append([varA_count/total, varB_count/total, both/total, other/total])
        depths.append(total)
        names.append(aa)
        positions.append(pos)
        bases.append([basefreq["A"]/total, basefreq["T"]/total, basefreq["C"]/total, basefreq["G"]/total, basefreq["-"]/total])
    draw_output(positions, proportion, depths, names, bases, output_file, variantA, variantB, refseq, panel3)



parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="output svg file", required=True)
parser.add_argument("-b", "--bam_file", help="sorted and indexed bam file", required=True)
parser.add_argument("-1", "--variant_1", help="variant 1 (BA.2, BA.4 or BA.5)", required=True)
parser.add_argument("-2", "--variant_2", help="variant 2 (BA.2, BA.4 or BA.5)", required=True)
parser.add_argument("-a", "--all", action="store_true", help="List all sites different from reference "
                                                             "(as opposed to only sites that differ between the two variants selected).")
parser.add_argument("-m", "--all_minor", action="store_true", help="List all sites where the minor allele reaches defined threshold.")
parser.add_argument("-f", "--minor_fraction", type=float, default=0.2, help="Fraction of reads with minor allele to report (when -m set)")
parser.add_argument("-d", "--minor_depth", type=int, default=20, help="minimum depth to report minor allele site (when -m set)")
parser.add_argument("-p3", "--panel3", action="store_true", help="Draw panel 3.")

args = parser.parse_args()

__main__(args.bam_file, args.variant_1, args.variant_2, args.output, args.all, args.all_minor, args.minor_fraction, args.minor_depth, args.panel3)