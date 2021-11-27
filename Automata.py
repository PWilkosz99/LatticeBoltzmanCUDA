import numpy

def generateMatrix(ImageMatrix, OldMatrix):
    for x in range(400):
        for y in range(600): # fill with value only part of matrixes
            ImageMatrix[y][x] = [255, 0, 0]
            OldMatrix[y][x] = [0, 255, 0]

    return ImageMatrix, OldMatrix