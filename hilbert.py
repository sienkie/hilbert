from PIL import Image
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
    max_value = 16581375  # 255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


def hilbert(chromSizes, files, output_file, create_map=False, threshold=0, size=2 ** depth):
    assert len(files) <= 3

    chrom_sizes = {}
    with open(chromSizes, 'r') as f:
        chrom_sizes = {str(line.split()[0]): int(line.split()[1]) for line in f}
    whole_genome = sum([v for (k, v) in chrom_sizes.items()])
    chrs = list(chrom_sizes.keys())
    chrs.sort(key=lambda item: (len(item), item))
    # TODO other chromosome order?? now: chr1 .. chr9, chrX, chrY, chr10 .. chr21

    data = np.zeros((size, size, 3), dtype=np.uint8)  # default black pixels

    interpolate2HC = interp1d([1, whole_genome],
                              [0, size ** 2 - 1])  # interpolate genome positions to points id in curve

    for i in range(len(files)):

        f = open(files[i], 'r')
        flines = f.readlines()[1:]  # TODO optional header ??
        f.close()

        max_score = max([float(line.split()[4]) for line in flines])
        interpolateColor = lambda x, y: interp1d([threshold, max_score], [x, y])  # interpolate to colors on x,y scale

        for line in flines:
            l = line.split()
            if len(l) > 4:
                print(l)
                chr = l[0]
                if chr in chrs:
                    relative_position = sum([chrom_sizes[chrs[c]] for c in range(0, chrs.index(chr))])
                    score = float(l[4])
                    start = int(l[1]) + relative_position
                    stop = int(l[2]) + relative_position
                    if score > threshold:
                        points = [id_to_xy(px) for px in
                                  range(int(interpolate2HC(start)), int(interpolate2HC(stop)) + 1)]
                        for point in points:
                            data[point][i] = 255

        if create_map:
            pass

    img = Image.fromarray(data, 'RGB')
    img.save(output_file)
    # img.show()
