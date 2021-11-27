import numpy

def generateMatrix(ImageMatrix, OldMatrix):
    for x in range(800):
        for y in range(600): # fill with value only part of matrixes
            ImageMatrix[x][y]=1
            OldMatrix[x][y]=1

            return ImageMatrix, OldMatrix