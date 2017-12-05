from PIL import Image, ImageDraw, ImageFont, ImageOps
import numpy as np
from scipy.interpolate import interp1d

global depth
depth = 11

# key: hilbert curve n=1 orientation
# value: dict with pixel coords as key and value as follows:
# (curve consecutive part, orientation of this part in next iterration)

hilbert_map_bin = {
    'a': {(0, 0): (0, 'd'), (0, 1): (1, 'a'), (1, 0): (3, 'b'), (1, 1): (2, 'a')},
    'b': {(0, 0): (2, 'b'), (0, 1): (1, 'b'), (1, 0): (3, 'a'), (1, 1): (0, 'c')},
    'c': {(0, 0): (2, 'c'), (0, 1): (3, 'd'), (1, 0): (1, 'c'), (1, 1): (0, 'b')},
    'd': {(0, 0): (0, 'a'), (0, 1): (3, 'c'), (1, 0): (1, 'd'), (1, 1): (2, 'd')},
}


def xy_to_id(x, y, order=depth):
    current_square = 'a'
    position = 0
    for i in range(order - 1, -1, -1):
        position <<= 2
        quad_x = 1 if x & (1 << i) else 0
        quad_y = 1 if y & (1 << i) else 0
        quad_position, current_square = hilbert_map_bin[current_square][(quad_x, quad_y)]
        position |= quad_position
    return position


def convertToQuaternary(n, order=depth):
    digits = []
    while n > 0:
        digits.insert(0, n % 4)
        n = n // 4
    return (order - len(digits)) * [0] + digits if len(digits) <= order else digits


hilbert_map = {
    'a': {0: ((0, 0), 'd'), 1: ((0, 1), 'a'), 3: ((1, 0), 'b'), 2: ((1, 1), 'a')},
    'b': {2: ((0, 0), 'b'), 1: ((0, 1), 'b'), 3: ((1, 0), 'a'), 0: ((1, 1), 'c')},
    'c': {2: ((0, 0), 'c'), 3: ((0, 1), 'd'), 1: ((1, 0), 'c'), 0: ((1, 1), 'b')},
    'd': {0: ((0, 0), 'a'), 3: ((0, 1), 'c'), 1: ((1, 0), 'd'), 2: ((1, 1), 'd')}
}


def id_to_xy(id, order=depth):
    x0, y0 = 0, 0
    xn, yn = 2 ** order - 1, 2 ** order - 1

    quaternary = convertToQuaternary(id, order)
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


def hilbert(chromSizes, files, output_file, create_map=False, revert=False, threshold=0, size=2 ** depth):
    assert len(files) <= 3

    filetypes = [kind.split(".")[-1] for kind in files]
    for kind in filetypes:
        assert kind in ['bed', 'bedgraph', 'bdg']

    with open(chromSizes, 'r') as f:
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
        if flines[0].split()[0] in ["track", "browser"]:  # TODO optional header ??
            flines = flines[1:]
        f.close()

        for line in flines:
            l = line.split()
            minlen = 5 if filetypes[i] == "bed" else 4
            if len(l) >= minlen:
                print(l)  # TODO remove
                chr = l[0]
                if chr in chrs:
                    relative_position = sum([chrom_sizes[chrs[c]] for c in range(0, chrs.index(chr))])
                    score = float(l[4]) if filetypes[i] == "bed" else float(l[3])
                    start = int(l[1]) + relative_position
                    stop = int(l[2]) + relative_position
                    if score > threshold:
                        points = [id_to_xy(px) for px in
                                  range(int(interpolate2HC(start)), int(interpolate2HC(stop)) + 1)]
                        for point in points:
                            data[point][i] = 255

        img = Image.fromarray(data, 'RGB')
        if revert:
            img = ImageOps.invert(img)
        img.save(output_file)

        if create_map:
            data_colors = np.zeros((size, size, 3), dtype=np.uint8)  # default black pixels
            colors = get_spaced_colors(len(chrs) + 1)[1:]

            for i in range(len(chrs)):
                print(i)
                chr = chrs[i]
                start = sum([chrom_sizes[chrs[c]] for c in range(0, chrs.index(chr))]) + 1
                stop = start + chrom_sizes[chrs[i]] - 1
                print(start, stop)  # TODO remove
                print(id_to_xy(int(interpolate2HC(start))), id_to_xy(int(interpolate2HC(stop)) + 1))
                points = [id_to_xy(px) for px in
                          range(int(interpolate2HC(start)), int(interpolate2HC(stop)) + 1)]
                for point in points:
                    data_colors[point] = colors[i]

            img_col = Image.fromarray(data_colors, 'RGB')

            padding = 5
            font_size = 50  # TODO: this is default font size, resize if all chromosomes wont fit
            font = ImageFont.truetype("fonts/THSarabunNew.ttf", font_size)  # TODO: change licence --> GNU font
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
            img.save(output_file.split(".")[0] + "_chrs_territories.png")
