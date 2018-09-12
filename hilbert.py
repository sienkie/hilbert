from PIL import Image, ImageDraw, ImageFont, ImageOps
import numpy as np
import argparse
from scipy.interpolate import interp1d

depth = 11

# key: hilbert curve n=m-1 orientation (where 'm' is current step)
# value: dict with curve consecutive part as key and value as follows:
# (pixel coords, orientation of this part in next iterration)

hilbert_map = {
    'a': {0: ((0, 0), 'd'), 1: ((0, 1), 'a'), 3: ((1, 0), 'b'), 2: ((1, 1), 'a')},
    'b': {2: ((0, 0), 'b'), 1: ((0, 1), 'b'), 3: ((1, 0), 'a'), 0: ((1, 1), 'c')},
    'c': {2: ((0, 0), 'c'), 3: ((0, 1), 'd'), 1: ((1, 0), 'c'), 0: ((1, 1), 'b')},
    'd': {0: ((0, 0), 'a'), 3: ((0, 1), 'c'), 1: ((1, 0), 'd'), 2: ((1, 1), 'd')}
}


def convertToQuaternary(n, order=depth):
    """ Converts decimal number to quaternary for id_to_xy function.

    Args:
        n (int): Decimal number.
        order (int): Order of Hilbert curve, implicating how many numbers should be returned (see examples).

    Returns:
        List of integers. If quaternary has more digits then order value returns list of quaternary digits.
        Otherwise adds 0 at the beginning of the list until it has length equal to order value.

    Examples:
        >>> print(convertToQuaternary(17, 1))
        [1, 0, 1]
        >>>print(convertToQuaternary(17,5))
        [0, 0, 1, 0, 1]

    """
    digits = []
    while n > 0:
        digits.insert(0, n % 4)
        n = n // 4
    return (order - len(digits)) * [0] + digits if len(digits) <= order else digits


def id_to_xy(point_id, order=depth):
    """ Maps point on line on Hilbert curve of given order.

    Args:
        point_id: ID of point on line (numbered from 0).
        order: Order of Hilbert curve.

    Returns:
        Tuple with coordinates of point (numbered point_id) on Hilbert curve of given order.

    """
    global hilbert_map
    x0, y0 = 0, 0
    xn, yn = 2 ** order - 1, 2 ** order - 1

    quaternary = convertToQuaternary(point_id, order)
    current_square = 'a'

    while x0 != xn and y0 != yn:
        quad_nr = quaternary.pop(0)
        ((quad_x, quad_y), current_square) = hilbert_map[current_square][quad_nr]
        x0 += quad_x * (2 ** (order - 1))
        y0 += quad_y * (2 ** (order - 1))
        xn -= abs(quad_x - 1) * (2 ** (order - 1))
        yn -= abs(quad_y - 1) * (2 ** (order - 1))
        order = order - 1
    return (x0, y0)


def get_spaced_colors(n):
    max_value = 255 ** 3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


def hilbert(chrom_sizes_file, files, output="hilbert.png", colors_list=None, color_map="", invert=False,
            color_legend="", threshold=None, size=2 ** depth, file_names=None):
    """ Creates visualization of genomic data on Hilbert curve.

    Args:
        chrom_sizes_file: Tab delimited file containing chromosomes sizes.
            Chromosome names should be compliant with names used in data files.
        files: List of paths to .bed or .bedgraph files.
        output: Path to output .png file.
        colors_list: List of colors for consecutive files (0-red, 1-green, 2-blue).
        color_map: Optional path to output color map .png file.
        invert: True if pixels should be colored with invert colors (0-cyan, 1-magenta, 2-yellow).
        color_legend: Optional path to output text file with color legend.
        threshold: Optional threshold above which regions are colored. List of floats.
        size: Size of output image corresponding to order of Hilbert Curve. Using default order value (11) results with
            creation of image with resolution 2048 x 2048px (2^11 x 2^11px).
        file_names: Optional list of alternative names for files. Useful when runing from Galaxy which
            converts all files to generic dataset (.dat) files.

    Returns:

    """
    if not colors_list:
        colors_list = range(len(files))
    if not threshold:
        threshold = [0.0 for i in files]
    assert len(files) <= 3
    assert len(threshold) >= len(files)
    assert len(colors_list) >= len(files)
    if not file_names:
        file_names = files

    filetypes = [kind.split(".")[-1] for kind in file_names]

    for kind in filetypes:
        assert kind in ['bed', 'bedgraph', 'bdg']

    with open(chrom_sizes_file, 'r') as f:
        chromlines = f.readlines()
        chrom_sizes = {str(line.split()[0]): int(line.split()[1]) for line in chromlines}
        chrs = [line.split()[0] for line in chromlines]
    whole_genome = sum([v for (k, v) in chrom_sizes.items()])

    data = np.zeros((size, size, 3), dtype=np.uint8)  # default black pixels

    interpolate2HC = interp1d([1, whole_genome],
                              [0, size ** 2 - 1])  # interpolate genome positions to points id in curve

    for i in range(len(files)):

        f = open(files[i], 'r')
        flines = f.readlines()
        while True:
            if flines[0].split()[0] in ["track", "browser", "#"]:
                flines = flines[1:]
            else:
                break
        f.close()

        for line in flines:
            l = line.split()
            minlen = 5 if filetypes[i] == "bed" else 4
            if len(l) >= minlen:
                chr = l[0]
                if chr in chrs:
                    relative_position = sum([chrom_sizes[chrs[c]] for c in range(0, chrs.index(chr))])
                    score = float(l[4]) if filetypes[i] == "bed" else float(l[3])
                    start = int(l[1]) + relative_position
                    stop = int(l[2]) + relative_position
                    if score > threshold[i]:
                        points = [id_to_xy(px) for px in
                                  range(int(interpolate2HC(start)), int(interpolate2HC(stop)) + 1)]
                        for point in points:
                            data[point][colors_list[i]] = 255

        img = Image.fromarray(data, 'RGB')
        if invert:
            img = ImageOps.invert(img)
        img.save(output, format="png")

        if color_legend:
            col = ["red", "green", "blue"]
            if invert:
                col = ["light blue", "pink", "yellow"]
            with open(color_legend, "w") as f:
                for i in range(len(file_names)):
                    f.write(col[colors_list[i]] + '\t' + file_names[i] + '\n')

        if color_map:
            data_colors = np.zeros((size, size, 3), dtype=np.uint8)  # default black pixels
            colors = get_spaced_colors(len(chrs) + 1)[1:]

            for i in range(len(chrs)):
                chr = chrs[i]
                start = sum([chrom_sizes[chrs[c]] for c in range(0, chrs.index(chr))]) + 1
                stop = start + chrom_sizes[chrs[i]] - 1
                points = [id_to_xy(px) for px in
                          range(int(interpolate2HC(start)), int(interpolate2HC(stop)) + 1)]
                for point in points:
                    data_colors[point] = colors[i]

            img_col = Image.fromarray(data_colors, 'RGB')

            padding = 5
            font = ImageFont.load_default()
            width = max([font.getsize(c)[0] for c in chrs])
            height = max([font.getsize(c)[1] for c in chrs])
            img_collist = Image.new('RGB', (width + 2 * padding, size), (0, 0, 0))
            draw = ImageDraw.Draw(img_collist)

            start = padding
            for i in range(len(chrs)):
                draw.text((padding, start), chrs[i], colors[i], font=font)
                draw = ImageDraw.Draw(img_collist)
                start += height

            img = Image.new("RGB", (width + 2 * padding + size, size), (0, 0, 0))

            img.paste(img_col, (0, 0))
            img.paste(img_collist, (size, 0))
            img.save(color_map, format="png")


def run_from_galaxy():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chr_sizes', type=str, required=True, help='tabular file with chromosome sizes')
    parser.add_argument('--files', type=str, required=True, help='input datasets')
    parser.add_argument('--output', type=str, required=True, help='output file')
    parser.add_argument('--output_map', type=str, default="", help='output galaxy dataset for color map')
    parser.add_argument('--output_legend', type=str, default="", help='output galaxy dataset for color legend')
    parser.add_argument('--file_names', type=str, required=False, help='file names for galaxy datasets')
    parser.add_argument('--colors', type=str, required=True, help='colors')
    parser.add_argument('--threshold', type=str, required=True, help='thresholds')
    parser.add_argument('--invert', default=False, help='invert colors')
    parser.add_argument('--size', default=11, help='depth of Hilbert Curve')
    args = parser.parse_args()

    inputs = args.files.split()
    inputs_names = [x.strip() for x in args.file_names.split(';') if x.strip()]
    colors = [int(i) for i in args.colors.split()]  # color list for consecutive files 0-red, 1-green, 2-blue
    thresholds = [float(i) for i in args.threshold.split()]
    size = 2 ** int(args.size)

    hilbert(chrom_sizes_file=args.chr_sizes, files=inputs, output=args.output, colors_list=colors,
            color_map=args.output_map, invert=args.invert, threshold=thresholds,
            color_legend=args.output_legend, size=size, file_names=inputs_names)


if __name__ == "__main__":
    run_from_galaxy()
