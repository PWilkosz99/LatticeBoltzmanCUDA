import numpy

def generateMatrix(ImageMatrix, OldMatrix):
    for x in range(600):
        for y in range(300): # fill with value only part of matrixes
            ImageMatrix[x][y] = [255, 0, 0]
            OldMatrix[x][y] = [0, 255, 0]

    return ImageMatrix, OldMatrix