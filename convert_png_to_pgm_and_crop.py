from PIL import Image
import os
import sys


def resize_dir(inpath, outpath):
    f = []
    for (dirpath, dirnames, filenames) in os.walk(inpath):
        f.extend(filenames)
        break
    
    if not os.path.exists(outpath):
        try:
            os.mkdir(outpath)
        except OSError:
            print("failed to create dir")

    for img in f:
    
        im = Image.open(inpath + img)
        region = im.crop((0, 1, 704, 485))
        img = img[:-3]
        outname = outpath + img + "pgm"
        region.convert('L')
        region.save(outname)
        print(outname)


path_to_files = str(sys.argv[1])
path_to_place_files = str(sys.argv[2])


if not os.path.exists(path_to_place_files):
    try:
        os.makedirs(path_to_place_files)
    except OSError:
        print("Failed to create output dir")
dirs = []
for (dirpath, dirnames, filenames) in os.walk(path_to_files):
    print(dirnames)
    dirs.extend(dirnames)
    break

for name in dirs:
    print(name)
    resize_dir(path_to_files + name + '/', path_to_place_files + name + '/')
